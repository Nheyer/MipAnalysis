#!/usr/bin/env python
# Name: Taylor Real & Nicholas Heyer (tdreal)
# Group Members: NOTCH2NL Group
# note the script is compatible with both python 2.7 and 3.5+

import sys
import pysam
import matplotlib.pyplot as plt
import numpy as np
from operator import add


def add_seq(dict, region, ts):
    if region in dict.keys():
        dict[region].append(ts)
    else:
        dict[region] = [ts]
    return dict


def parse_fasta(fasta_path):
    is_seq = False
    first_seq = True
    seq = ""
    name_region = []
    seq_dict = {}
    name_dict = {}

    fasta = open(fasta_path, 'r')
    # loop through the lines of the fasta
    for line in fasta:
        if not is_seq and line[0] == '>':
            if not first_seq:
                print(name_region)
                seq_dict = add_seq(seq_dict, name_region[1], seq)
                name_dict[seq] = name_region[0]
            name_region = line.rstrip("\n").lstrip('>').split("\t")
            first_seq = False
            is_seq = True
        elif is_seq:
            seq = line.rstrip("\n")
            is_seq = False
        elif line[0].isalpha():
            seq += line.rstrip("\n")
    # add the last sequence in the fasta
    print(name_region)
    seq_dict = add_seq(seq_dict, name_region[1], seq)
    name_dict[seq] = name_region[0]
    fasta.close()
    return seq_dict, name_dict


class ARGS:
    def __init__(self):
        self.make_graphs = True
        self.output_file = sys.stdout
        self.out_path = "/dev/stdout"
        self.input_files = []
        self.ref_path = ""
        self.log_file = sys.stderr
        self.log_path = "/dev/stderr"
        self.bam_path = ""

    def io(self):
        help_str = "Usage :\n"
        help_str += "Python MipAnalysis.py -i <path to bam> \n"
        help_str += "-h\t--help\t\tPrint detailed help information\n"
        if len(sys.argv) < 2:
            sys.stderr.write(help_str)
            sys.exit(-1)
        for i in range(80):
            help_str += "~"
        help_str += "\nFlags:\n"
        help_str += "-i\t--input      \t[path],[path]\tA comma separated list of paths to the sorted\n"
        help_str += "  \t             \t      \tindexed bam file to analyse\n"
        help_str += "-g\t--no-graphics\t[bool]\tIf this flag is included the program sill skip generating graphics\n"
        help_str += "-o\t--output     \t[path]\tPath to output information, default will print to stdout\n"
        help_str += "-b\t--out-bam    \t[path]\tPath to a bam file to put all novel reads in a region checked\n"
        help_str += "-r\t--ref        \t[path]\tPath to a fasta with read name being location (samtools like) and a \n"
        help_str += "  \t             \t      \tcomment of what to call it\n"
        i = 1
        while i < len(sys.argv):
            if sys.argv[i] in ['-i', "--input"]:
                self.input_files = sys.argv[i + 1].split(",")
                i = i + 2
            elif sys.argv[i] in ['-g', "--no-graphics"]:
                self.make_graphs = False
                i = i + 1
            elif sys.argv[i] in ['-o', "--output"]:
                self.output_file = open(sys.argv[i + 1], 'w')
                self.out_path = sys.argv[i + 1]
                i = i + 2
            elif sys.argv[i] in ['-h', "--help"]:
                sys.stderr.write(help_str)
                sys.exit(0)
            elif sys.argv[i] in ['-r', "--ref"]:
                self.ref_path = sys.argv[i + 1]
                i = i + 2
            elif sys.argv[i] in ['-l', "--log"]:
                self.log_path = sys.argv[i + 1]
                self.log_file = open(sys.argv[i + 1], "a")
                i = i + 2
            elif sys.argv[i] in ['-b', "--out-bam"]:
                self.bam_path = sys.argv[i + 1]
                i = i + 1
            else:
                sys.stderr.write("Unrecognised input:\t" + sys.argv[i] + "\n")
                i = i + 1
        return True


# Takes BAM file data using Pysam ans outputs quality reports and graphs.
def runAnalysis(myReader, user_args, name):

    # Performs analysis on the entire run.

    mapQual = []
    readQual = []
    rqData = []
    readSize = []
    baseCoverage = []
    allBX = set([])
    avRS = 0.0
    for read in myReader.fetch("NOTCH2NL-consensus"):   # get all reads ##todo add dynamic contigs
        mapQual.append(read.mapping_quality)  # adds mapping qualities
        rqData.append(read.get_tag("RQ"))
        readSize.append(read.query_alignment_length)  # adds read size
        if read.is_read2 and "N" not in read.get_tag("BX"):  # read 2 contains barcode
            allBX.add(read.get_tag("BX"))
    sys.stderr.write("Finished getting MQ, RQ, RS, and num bx, working on coverage ... \n")
    # get coverage of the entire region
    num_a, num_t, num_c, num_g = myReader.count_coverage("NOTCH2NL-consensus")
    baseCoverage = list(map(add, map(add, num_a, num_t), map(add, num_c, num_g)))
    avRS = sum(readSize) / len(readSize)  # finds average
    avMQ = sum(mapQual) / len(mapQual)  # finds average Mapping quality
    avRQ = sum(rqData) / len(rqData)  # finds average Read Quality
    avCoverage = sum(baseCoverage) / len(baseCoverage)  # finds average
    totBar = len(allBX)
    if user_args.make_graphs:
        plt.hist(mapQual, facecolor='xkcd:lightblue')  # Plot Mapping Quality
        plt.title("Overall Mapping Quality")
        plt.ylabel("Number of Reads")
        plt.xlabel("Mapping Quality of Reads")
        #plt.axis([0, 60, 0, 1500])
        plt.savefig(name + "_MQ.pdf", format='pdf')
        plt.close()
        plt.hist(rqData, facecolor = 'xkcd:lavender')  # Plot Read Quality
        plt.title("Overall Read Quality")
        plt.ylabel("Number of Reads")
        plt.xlabel("Quality Score of Reads")
        #plt.axis([0, 42, 0, 1500])
        plt.savefig(name + "_RQ.pdf", format='pdf')
        plt.close()
        plt.fill_between(range(1, len(baseCoverage)+1), baseCoverage)  # histogram plot
        plt.title("Overall Sequencing Depth")
        plt.ylabel("Number of Reads")
        plt.xlabel("Position")
        #plt.axis([0, 1001, 0, 1501])
        #plt.xticks(np.arange(0, 1001, 100))
        plt.savefig(name + "_BC.pdf", format='pdf')
        plt.close()

    return avMQ, avRQ, avRS, avCoverage, totBar


# generates summery stats for the region of interest
def analysis_of_region(alignmentfile, region, alleles, user_args, fq_dict, name):
    out_align = pysam.AlignmentFile(user_args.bam_path.rstrip(".bam") + region + ".bam", 'wb', template=alignmentfile)
    contig = region.split(":")[0]
    start_stop = region.split(":")[1].split("-")
    # generate a list of sets, one set for each allile to check for ( including ref )
    BXs = []
    all_bx = set()
    BX_freq = []
    num_unq_bx = 0
    for dummy in alleles:
        BXs.append(set())
    # get an array of num a/t/c/g and then combine them to an overall coverage
    num_a, num_t, num_c, num_g = alignmentfile.count_coverage(contig, int(start_stop[0]), int(start_stop[1]))
    baseCoverage = list(map(add, map(add, num_a, num_t), map(add, num_c, num_g)))
    if user_args.make_graphs:
        # Create a histogram
        plt.fill_between(range(len(baseCoverage)), baseCoverage)
        plt.title(region)
        plt.xlabel("Relitive position")
        plt.ylabel("number of reads")
        plt.savefig(name + "_" + region + "_coverage.pdf", format="pdf")
        plt.close()
    average_coverage = (float(sum(baseCoverage))/float(len(baseCoverage)))
    # grab all the reads in the defined region
    for read in alignmentfile.fetch(contig, int(start_stop[0]), int(start_stop[1])):
        if read.is_read2:
            # if it is in read 2, check to see if the sequence is a match ofr any input allele
            found_read = False
            for i_allele in range(len(alleles)):
                if alleles[i_allele] in read.query_sequence:
                    # if it is a match for any of them try and add it to the proper set and a universal set
                    if "N" not in read.get_tag("BX"):
                        BXs[i_allele].add(read.get_tag("BX"))
                        all_bx.add(read.get_tag("BX"))
            # if you find the allele we can go to the next read

            if not found_read:
                out_align.write(read)

    out_align.close()

    # convert the sets we have into frequencies, and append them to a list not we have var ref pairs as % 0 and 1
    # then we will zip them to a dictionary, where the input allele is the key to it's freq and return it
    num_unq_bx = float(len(all_bx))
    for bx_set_index in range(len(BXs)):
        Bx_0 = float(len(BXs[bx_set_index]))
        if num_unq_bx > 0:
            BX_freq.append(Bx_0 / num_unq_bx)
        else:
            print("no BX in region!!")
            BX_freq.append(0.0)
    # we need to add the values now, so doing that here
    for i in range(len(alleles)):
        fq_dict[alleles[i]] = str(BX_freq[i])
    fq_dict[region + "_cov"] = str(average_coverage)
    fq_dict[region + "_bx"] = str(num_unq_bx)
    return fq_dict


# runs all the analysis for a particular bam file
def run_all(run_dict, seq_dict, bam_file, arg):
    output_values = {}
    myReader = pysam.AlignmentFile(bam_file, 'rb')
    MQ, RQ, size, coverage, barcodes = runAnalysis(myReader, arg, bam_file.rstrip(".bam"))
    arg.output_file.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%8d" % (bam_file, MQ, RQ, size, coverage, barcodes))
    sys.stderr.write("QC for file " + bam_file + " Completed! \nRunning Variant Analysis...\n")

    for region in run_dict.keys():
        output_values = analysis_of_region(myReader, region, run_dict[region],
                                           arg, output_values, bam_file.rstrip(".bam"))
        arg.output_file.write(output_values[region + "_cov"] + "\t")
        arg.output_file.write(output_values[region + "_bx"] + "\t")
        sys.stderr.write("Finished analysing region " + region + "in file " + bam_file + "\n")
    for key in seq_dict:
        arg.output_file.write("\t" + output_values[key])
    myReader.close()


def main():
    arguments = ARGS()
    arguments.io()
    run_data, seq_nm = parse_fasta(arguments.ref_path)
    # writhe a header line that is dynamic for each line it needs
    arguments.output_file.write("File/Individual\t"       + "Average Mapping Quality\t" + "Average Read Quality\t" +
                                "Average Size of Reads\t" + "Average Coverage\t"        + "Total Number of Barcodes\t")
    for region in run_data.keys():
        arguments.output_file.write(region + "-number of barcodes \t" + region + "average coverage\t")
    for known_var in seq_nm.keys():
        arguments.output_file.write(seq_nm[known_var] + "\t")
    # finished outputing a header, with all values from the fasta
    for file in arguments.input_files:
        arguments.output_file.write("\n")
        run_all(run_data, seq_nm, file, arguments)

    arguments.output_file.close()
    sys.stderr.write("Variant Analysis Completed!\nAll Graphs Saved as .pdf.\n")


if __name__ == "__main__":
    main() 
