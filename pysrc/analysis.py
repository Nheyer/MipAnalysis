#!/usr/bin/env python3
# Name: Taylor Real & Nicholas Heyer (tdreal)
# Group Members: NOTCH2NL Group
# note the script is compatible with  3.5+

import sys
import pysam
import os
import numpy as np
import matplotlib.pyplot as plt
from operator import add
DEBUG = 1
BITWISE_READ1 = 64
BITWISE_READ2 = 128

class ARGS:
    def __init__(self):
        self.make_graphs = True
        self.trigger_split = False
        self.bwa_ref = ""
        self.bx_len = [-1,-1]
        self.output_file = sys.stdout
        self.out_path = "/dev/stdout"
        self.input_files = []
        self.infq_pair = []
        self.vcf_path = "../test/n2nl_199950-200130_snps.vcf.gz"
        self.bed_path = "../test/n2nl_region_199950-200130.bed"
        self.log_file = sys.stderr
        self.log_path = "/dev/stderr"
        self.bam_path = ""
        self.threads = 10

    def io(self):
        help_str = "Usage :\n"
        help_str += "Python3 MipAnalysis.py -i <path to bam> \n"
        help_str += "-h\t--help\t\tPrint detailed help information\n"
        if len(sys.argv) < 2:
            sys.stderr.write(help_str)
            sys.exit(-1)
        for i in range(100):
            help_str += "~"
        help_str += "\nFlags:\n"
        help_str += "-i\t--input      \t[path],[path]\tA comma separated list of paths to the sorted\n"
        help_str += "  \t             \t      \tindexed bams or fastqs with r1:r2,r1:r2...  file to analyse\n"
        help_str += "  \t             \t      \tnote to input fastqs you need to supply -R and -x\n"
        help_str += "-g\t--no-graphics\t[bool]\tIf this flag is included the program sill skip generating graphics\n"
        help_str += "-o\t--output     \t[path]\tPath to output information, default will print to stdout\n"
        help_str += "-b\t--out-bam    \t[path]\tPath to a bam file to put all novel reads in a region checked\n"
        help_str += "-v\t--vcf        \t[path]\tPath to a bgziped vcf with all variants you are interested in checking \n"
        help_str += "-d\t--bed        \t[path]\tPath to a bed file of the regions to check\n"
        help_str += "-R\t--mem-ref    \t[path]\tpath to align input fastqs to\n"
        help_str += "-x\t--bx-len     \t[int] \tthe lenght of the BX to strip from input fastq files\n"
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
            elif sys.argv[i] in ['-v', "--vcf"]:
                self.vcf_path = sys.argv[i + 1]
                i = i + 2
            elif sys.argv[i] in ['-l', "--log"]:
                self.log_path = sys.argv[i + 1]
                self.log_file = open(sys.argv[i + 1], "a")
                i = i + 2
            elif sys.argv[i] in ['-b', "--out-bam"]:
                self.bam_path = sys.argv[i + 1]
                i = i + 2
            elif sys.argv[i] in ['-R', "--mem-ref"]:
                self.bwa_ref = sys.argv[i + 1]
                self.trigger_split = True
                i = i + 2
            elif sys.argv[i] in ['-x', "--bx-len"]:
                self.bx_len = sys.argv[i + 1].split("-")
                self.trigger_split = True
                i = i + 2
            elif sys.argv[i] in ['-t', "--threads"]:
                self.threads = sys.argv[i + 1]
                i = i + 2
            elif sys.argv[i] in ['-d', "--bed"]:
                self.bed_path = sys.argv[i + 1]
                i = i + 2
            else:
                sys.stderr.write("Unrecognised input:\t" + sys.argv[i] + "\n")
                i = i + 1
        if self.trigger_split:
            if self.bwa_ref == "" or self.bx_len == -1 or ':' not in self.input_files[0]:
                print("Attempted to input raw files, without -x or -R\n \
                      OR you added those tags erroneously when trying to use pre-aligned bams")
                sys.exit(-3) ## crash becouse we are missing inputs we expect
            temp = self.input_files
            self.input_files = []
            for pair in temp:
                fq_pair = pair.split(":")
                self.infq_pair.append(fq_pair)
                # assume it is named "somthing.fastq" and make the bams be "somthing.bam"
                self.input_files.append(fq_pair[0].rstrip("fastq") + "bam")
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
        if "N" not in read.get_tag("BX"):
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

    return name, avMQ, avRQ, avRS, avCoverage, totBar


def main():
    arguments = ARGS()
    arguments.io()
    points_of_interest = parse_bed_vcf_pair(arguments.vcf_path, arguments.bed_path)
    make_header(points_of_interest, arguments.output_file)
    for file_index in range(len(arguments.input_files)):
        if arguments.trigger_split:
            make_bam(arguments, file_index)
        run_file(points_of_interest, arguments.input_files[file_index], arguments)



def combine_vars(raw_vars, current_mnp):
    next_mnp = []
    if len(raw_vars) == 0:
        return current_mnp
    else:
        for mnp in current_mnp:
            # add a pair for ref of the snp we are on
            next_mnp.append([mnp[0].lstrip('_')  + "_ref-" + raw_vars[0][0], mnp[1].copy()])
            next_mnp[-1][1][raw_vars[0][3]] = raw_vars[0][1]
            # add one for the varriant
            next_mnp.append([mnp[0].lstrip('_')  + "_var-" + raw_vars[0][0], mnp[1].copy()])
            next_mnp[-1][1][raw_vars[0][3]] = raw_vars[0][2]
        # now remove it from the list and try to repeate the process
        del raw_vars[0]
        return combine_vars(raw_vars, next_mnp)


def parse_bed_vcf_pair(vcf_path,bed_path):
    bed = open(bed_path, "r")
    vcf = pysam.VariantFile(vcf_path)
    regions_maped_to_vars = []

    for line in bed:
        # turn the bed into a list of [chr, start, stop]
        current_region = line.rstrip("\n").split("\t")
        variants = []
        #rip id,ref bp, alt bp, and pos from the vcf (fix pos to 0 based idx)
        for variant_record in vcf.fetch(contig = current_region[0], start = int(current_region[1]), end = int(current_region[2])):
            variants.append([variant_record.id, variant_record.ref, variant_record.alleles[1], variant_record.pos - 1 ])
        # combine it and format
        print(variants)
        regions_maped_to_vars.append([current_region,combine_vars(variants, [["",{}]])])
    return regions_maped_to_vars


def make_header(data_to_run, fl):
    # write out the standard portion:
    fl.write("File/Individual\tAverage Mapping Quality\tAverage Read Quality\tAverage Size of Reads\tAverage Coverage\tTotal Number of Barcodes")
    # make the dynamic header
    for region in data_to_run:
        fl.write("\t%s:%s-%s_average_coverage" % tuple(region[0]))
        fl.write("\t%s:%s-%s_total_barcodes" % tuple(region[0]))
        for mnp_record in region[1]:
            fl.write("\t%s_num_barcodes" % mnp_record[0])
    fl.write("\n")


def get_stats_on_region(region_to_run, bam):
    # make things more readable
    region_with_mnp = region_to_run[0]
    mnps = region_to_run[1]
    list_of_sets = [[set(), set()] for x in range(len(mnps))]
    if DEBUG > 0:
        print(region_with_mnp)
        print(mnps)
    for sequence in bam.fetch(contig = region_with_mnp[0], start = int(region_with_mnp[1]), end = int(region_with_mnp[2])):
        i = 0 # we need an iterator
        match_flag = False
        seq_str = sequence.query_sequence

        seq_map = sequence.get_reference_positions()
        rpos_seq = dict(zip(seq_map, seq_str))
        for mnp in mnps:
            for snp_pos in mnp[1].keys():
                if snp_pos in seq_map:
                    #print (mnp[0],mnp[1][snp_pos],rpos_seq[snp_pos], sequence.get_tag("BX"))
                    if mnp[1][snp_pos] == rpos_seq[snp_pos]:
                        # we have a match flip the flag and continue
                        match_flag = True
                    else:
                        # we dont match the motif
                        match_flag = False
                        break
                else:
                    break## if we don't have the bp for one just skip the sequence

            # when we match break before iterating so we can add properly
            if match_flag:
               break
            i = i + 1
        if sequence.get_tag("BX") == "CGTCGC-MIP-CG":
            print(mnp[0],mnps[i][0], rpos_seq[199991], rpos_seq[199996], len(seq_map), len(seq_str))
        if match_flag:
            match_flag = False # reset this
            try:
                list_of_sets[i][sequence.is_read1].add(sequence.get_tag("BX"))
            except:
                # if we can't get the BX tag throw error
                sys.stderr.write("Sequence is missing a BX tag, this shouldn't happen!!\n")
    return make_sure_r1_and_r2_agree(list_of_sets)


##TODO add disagreeing mip tags to an error file
def make_sure_r1_and_r2_agree(nr1_r2):
    list_of_number_of_barcodes = [0]

    for r1, r2 in nr1_r2:
        num_bx = len(r1.intersection(r2))
        list_of_number_of_barcodes.append(num_bx)
        list_of_number_of_barcodes[0] += num_bx
    return list_of_number_of_barcodes


def run_file(data_to_run, bam_file, user_args):
    if DEBUG > 0:
        print(bam_file)
    # open the bam file and get some summery statistics
    reader = pysam.AlignmentFile(bam_file, 'rb')
    user_args.output_file.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%8d" %
                                runAnalysis(reader, user_args, bam_file.rstrip(".bam")))
    # now run analysis on each of the regions of the bed
    for region in data_to_run:
        ##TODO get average coverage here
        user_args.output_file.write("\tNA")
        for value in get_stats_on_region(region, reader):
            user_args.output_file.write("\t%s" % value)
    user_args.output_file.write("\n")


###TODO add the bowtie2 portion to this !!!
def make_bam(arguments, fidx):
    os.system("FastqMipTag " + arguments.infq_pair[fidx][0] +
              " " + arguments.infq_pair[fidx][1] + " _tagged.fq "
              + str(arguments.bx_len[0]) + " " + arguments.bx_len[1]  )
    os.system("bwa mem -t " + str(arguments.threads) + " " +
              arguments.bwa_ref + " " +
              arguments.infq_pair[fidx][0] + "_tagged.fq " + ## read one tagged
              arguments.infq_pair[fidx][0] + "_tagged.fq " + ## read two tagged
              ">" + arguments.input_files[fidx].rstrip("bam")+"sam")
    os.system("samtools view -Sb@ "+arguments.threads + " " + arguments.input_files[fidx].rstrip("bam")+"sam" +
              " -o " + arguments.input_files[fidx].rstrip("bam") + "unsorted.bam")
    os.system("samtools sort -@ " + arguments.threads + " " +
              arguments.input_files[fidx].rstrip("bam") + "unsorted.bam" +
              " -o " + arguments.input_files[fidx])
    os.system("samtools index " + arguments.input_files[fidx])



if __name__ == "__main__":
    main() 
