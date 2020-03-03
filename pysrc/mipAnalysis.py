#!/usr/bin/env python3
# Name: Taylor Real & Nicholas Heyer (tdreal)
# Group Members: NOTCH2NL Group
# note the script is compatible with  3.5+

import sys
import pysam
import os
import matplotlib.pyplot as plt
from operator import add
DEBUG = 0
BITWISE_READ1 = 64
BITWISE_READ2 = 128

class ARGS:
    def __init__(self):
        self.make_graphs = True
        self.trigger_split = False
        self.bwa_ref = ""
        self.bow_ref = ""
        self.bx_len = [-1,-1]
        self.output_file = sys.stdout
        self.out_path = "/dev/stdout"
        self.input_files = []
        self.infq_pair = []
        self.vcf_path = ""
        self.bed_path = ""
        self.log_file = sys.stderr
        self.log_path = "/dev/stderr"
        self.bam_path = ""
        self.threads = "10"
        self.picard = "java -jar ~/bin/picard/build/libs/picard.jar"

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
        help_str += "-R\t--mem-ref    \t[path]\tpath to align input fastqs to it should be somthing.fa\n"
        help_str += "  \t--bow-ref    \t[path]\tif the bowtie2 index is in a different location put it here\n"
        help_str += "-x\t--bx-len     \t[int]-[int] \tthe lenght of the BX to strip from input fastq files\n"
        help_str += "-p\t--picard     \t[path]\tPath to the .jar file that contains picard default assumes picard was\n"
        help_str += "  \t             \t      \tinstalled in ~/bin/\n"
        i = 1
        while i < len(sys.argv):
            if sys.argv[i] in ['-i', "--input"]:
                self.input_files = sys.argv[i + 1].split(",")
                i = i + 2
            elif sys.argv[i] in ['-p', "--picard"]:
                self.picard = "java -jar " + sys.argv[i + 1]
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
            elif sys.argv[i] == "--bow-ref":
                self.bow_ref = sys.argv[i + 1]
                self.trigger_split = True
                i = i + 2
            elif sys.argv[i] in ['-x', "--bx-len"]:
                self.bx_len = sys.argv[i + 1].split("-")
                self.trigger_split = True
                i = i + 2
            elif sys.argv[i] in ['-t', "--threads"]:
                self.threads = str(sys.argv[i + 1])
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
            # by default assume bowtie2 reference is named the same and is in the same spot
            if not self.bow_ref:
                self.bow_ref = self.bwa_ref.rstrip("fasta").rstrip(".")
                print(self.bow_ref, self.bwa_ref)
            temp = self.input_files
            self.input_files = []
            for pair in temp:
                fq_pair = pair.split(":")
                self.infq_pair.append(fq_pair)
                # assume it is named "somthing.fastq" and make the bams be "somthing.bam"
                self.input_files.append(fq_pair[0].rstrip("fastq") + "bam")
        return True


def main():
    arguments = ARGS()
    arguments.io()
    points_of_interest = parse_bed_vcf_pair(arguments.vcf_path, arguments.bed_path)
    make_header(points_of_interest, arguments.output_file)
    if DEBUG > 2:
        print(points_of_interest)
        sys.exit(0)
    for file_index in range(len(arguments.input_files)):
        if arguments.trigger_split:
            make_bam(arguments, file_index)
        run_file(points_of_interest, arguments.input_files[file_index], arguments)


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
        if DEBUG > 1:
            print(variants)
        regions_maped_to_vars.append([current_region,combine_vars(variants, [["",{}]])])
    return regions_maped_to_vars


def combine_vars(raw_vars, current_mnp):
    next_mnp = []
    if len(raw_vars) == 0:
        return current_mnp
    else:
        tigger_insertion = 0
        for mnp in current_mnp:
            tigger_insertion = False #assume it is not an insertion
            #code indels
            if len(raw_vars[0][1]) > len(raw_vars[0][2]) :
                # refrence is greater then alt so it is a deletion
                raw_vars[0][2] = "D" * (len(raw_vars[0][1]) - len(raw_vars[0][2]))
                raw_vars[0][1] = raw_vars[0][1][1] ## code ref as the first bp in the deletion var
                raw_vars[0][3] += 1 ## increase our position so we are ] into the deletion instead of the pos just before
            elif len(raw_vars[0][1]) < len(raw_vars[0][2]) :
                # alt is greater then ref so it is an insertion
                tigger_insertion = True


            # add a pair for ref of the snp we are on
            next_mnp.append([mnp[0].lstrip('_')  + "_ref-" + raw_vars[0][0], mnp[1].copy()])
            if tigger_insertion:
                next_mnp[-1][1][str(raw_vars[0][3])] = raw_vars[0][1]
            else:
                next_mnp[-1][1][raw_vars[0][3]] = raw_vars[0][1]
            # add one for the varriant
            next_mnp.append([mnp[0].lstrip('_')  + "_var-" + raw_vars[0][0], mnp[1].copy()])
            if tigger_insertion:
                next_mnp[-1][1][str(raw_vars[0][3]) + "-I"+ str(len(raw_vars[0][2]) - len(raw_vars[0][1])) ] = raw_vars[0][2]
            else:
                next_mnp[-1][1][raw_vars[0][3]] = raw_vars[0][2]
        # now remove it from the list and try to repeate the process
        del raw_vars[0]
        return combine_vars(raw_vars, next_mnp)


def make_bam(arguments, fidx):
    # make the bx tags on each fastq
    os.system("FastqMipTag " + arguments.infq_pair[fidx][0] +
              " " + arguments.infq_pair[fidx][1] + " _tagged.fq "
              + str(arguments.bx_len[0]) + " " + arguments.bx_len[1] )
    # create an initial alignment with bwa adding the tags
    os.system("bwa mem -C -t " + arguments.threads + " " +
              arguments.bwa_ref + " " +
              arguments.infq_pair[fidx][0] + "_tagged.fq " + ## read one tagged
              arguments.infq_pair[fidx][1] + "_tagged.fq " + ## read two tagged
              ">" + "temp_file.sam")
    # compress to bam, sort, index, then remove reads that don't map to a mip position
    os.system("samtools view -Sb@ "+ arguments.threads + " temp_file.sam " +
              " -o " + arguments.input_files[fidx].rstrip("bam") + "unsorted_raw.bam")
    os.system("rm temp_file.sam &")
    os.system("samtools sort -@ " + arguments.threads + " " +
              arguments.input_files[fidx].rstrip("bam") + "unsorted_raw.bam" +
              " -o " + arguments.input_files[fidx].rstrip("bam") + "sorted_raw.bam")
    os.system("samtools index " + arguments.input_files[fidx].rstrip("bam") + "sorted_raw.bam")
    os.system("samtools view  -q 60 -b@ "+arguments.threads + " -L " + arguments.bed_path + " "
              + arguments.input_files[fidx].rstrip("bam") + "sorted_raw.bam -o "
              + arguments.input_files[fidx].rstrip("bam") + "hits.bam")
    # now unmap the bam into a ubam so we can use bowtie2 for more accurate mapping
    os.system(arguments.picard + " RevertSam I=" + arguments.input_files[fidx].rstrip("bam") + "hits.bam O="
              + arguments.input_files[fidx].rstrip(".bam") + "_unmapped.unsorted.bam")
    os.system("samtools sort -n@ " + arguments.threads + " "
              + arguments.input_files[fidx].rstrip(".bam") + "_unmapped.unsorted.bam -o"
              + arguments.input_files[fidx].rstrip(".bam") + "_unmapped.sorted.bam")
    os.system("samtools index -@ " + arguments.threads + " "
              + arguments.input_files[fidx].rstrip(".bam") + "_unmapped.sorted.bam")
    # now we will map with bowtie2, compress to a bam, sort and index the reads
    os.system("bowtie2 -x " + arguments.bow_ref + " --end-to-end --no-mixed --no-discordant --align-paired-reads --preserve-tags --threads "
              + arguments.threads + " -b " + arguments.input_files[fidx].rstrip(".bam") + "_unmapped.sorted.bam > "
              + "temp_file.sam")
    os.system("samtools view -q 30 -Sb@ " + arguments.threads + " temp_file.sam > "
              + arguments.input_files[fidx].rstrip("bam") + "unsorted.bam")
    os.system("rm temp_file.sam &")
    os.system("samtools sort -@ " + arguments.threads + " " + arguments.input_files[fidx].rstrip("bam") + "unsorted.bam -o "
              + arguments.input_files[fidx] )
    os.system("samtools index " + arguments.input_files[fidx])


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


def run_file(data_to_run, bam_file, user_args):
    if DEBUG > 0:
        print(bam_file)
    # open the bam file and get some summery statistics
    reader = pysam.AlignmentFile(bam_file, 'rb')
    user_args.output_file.write("%s\t%.4f\t%.4f\t%.4f\t%.4f\t%8d" %
                                runAnalysis(reader, user_args, bam_file.rstrip(".bam")))
    # now run analysis on each of the regions of the bed
    for region in data_to_run:
        if DEBUG > 0:
            print(region[0])
            user_args.log_file.write("Working on region==>\t%s:%s-%s\n" % (region[0][0],region[0][1],region[0][2]))
        # get a pileup of our region
        a,t,c,g = reader.count_coverage(contig = region[0][0],
                                        start = int(region[0][1]),
                                        stop = int(region[0][2]),
                                        quality_threshold = 0  )
        # add all them into a list of depths
        full_depth = list(map(add,map(add,a,t),map(add,c,g)))
        if DEBUG > 0:
            print("We have a total depth of %s" % sum(full_depth))
        # make an average and print it !
        user_args.output_file.write("\t%.2f" % (float(sum(full_depth))/(int(region[0][2])- int(region[0][1]))))
        for value in get_stats_on_region(region, reader, user_args):
            user_args.output_file.write("\t%s" % value)
    user_args.output_file.write("\n")


# Takes BAM file data using Pysam and outputs quality reports and graphs.
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


def get_stats_on_region(region_to_run, bam, user_args):
    # make things more readable
    region_with_mnp = region_to_run[0]
    mnps = region_to_run[1]
    list_of_sets = [[set(), set()] for x in range(len(mnps))]
    all_tags = set()
    use_tags = set()
    if DEBUG > 2:
        print(region_with_mnp)
        print(mnps)
    for sequence in bam.fetch(contig = region_with_mnp[0], start = int(region_with_mnp[1]), end = int(region_with_mnp[2])):
        # get the tag once
        try:
            s_tag = sequence.get_tag("BX")
            all_tags.add(s_tag)
        except:
            # if we can't get the BX tag throw error
            user_args.log_file.write("Sequence is missing a BX tag, this shouldn't happen!!\n")
            continue
        i = 0 # we need an iterator
        match_flag = False
        rpos_seq , seq_map = seq_to_usable(sequence, user_args)
        for mnp in mnps:
            for snp_pos in mnp[1].keys():
                if snp_pos in seq_map:
                    if DEBUG > 1:
                        print (mnp[0],mnp[1][snp_pos],rpos_seq[snp_pos])
                    if mnp[1][snp_pos] == rpos_seq[snp_pos]:
                        # we have a match flip the flag and continue
                        match_flag = True
                    else:
                        # we dont match the motif if any is wrong, so keep going
                        match_flag = False
                        break
                else:
                    break## if we don't have the bp for one just skip the sequence

            # when we match break before iterating so we can add properly
            if match_flag:
                use_tags.add(s_tag)
                break
            i = i + 1
        if match_flag:
            match_flag = False # reset this
            list_of_sets[i][sequence.is_read1].add(s_tag)

        if DEBUG > 0:
            print("Invalid Tags:\t%s\nAll Tags:\t%s" % (len(all_tags-use_tags),len(all_tags)) )
    return make_sure_r1_and_r2_agree(list_of_sets)



def seq_to_usable(seq, user_args):
    seq_str = seq.query_sequence
    cigar_info = seq.get_cigar_stats()
    seq_list = []
    previous_pos = '0'
    seq_index = 0
    seq_map = seq.get_reference_positions(full_length = True)
    new_map = []
    ins_len = 0
    for pos in seq_map:
        if pos is not None: # then we are not coding an insertion here sooo it is not needed to modify here
            ins_len = 0
            previous_pos = str(pos)
            new_map.append(str(pos))
            seq_list.append(seq_str[seq_index])
        elif ins_len == 0: # we need to start an insertion
            # make sure the ref will not map if we have an ins
            seq_list[-1] = "This is not a bp pair"
            new_map.append(previous_pos + "-I1")
            seq_list.append(seq_str[seq_index-1] + seq_str[seq_index])
            ins_len = 1
        elif ins_len > 0: # we need to extend our insertion
            ins_len += 1
            new_map[-1] = previous_pos + "-I" + str(ins_len)
            seq_list[-1] = seq_list[-1] + seq_str[seq_index]
        # increment the position on the mapped sequence
        seq_index += 1
    seq_dict =  dict(zip(seq_map, seq_list))
        #print(seq_str, "\n", seq_map)
    if cigar_info[1][2] > 0 : # if this means the CIGAR string informs at least 1 deletion, we need to add the refs and note them for later
        improper_pairs = [x[1]  for x in  seq.get_aligned_pairs() if x[0] is None]
        del_len_ds = [x[1] * 'D' for x in seq.cigartuples if x[0] == 2]
        if len(del_len_ds) != len(improper_pairs):
            user_args.log_file.write("CIGAR ref/seq mismatch on sequence, likely a malformed bam!!\n")
            sys.exit(-1)
        else:
            seq_dict.update(dict(zip(improper_pairs,del_len_ds))) ## add the dels
            new_map.extend(improper_pairs)
        if DEBUG > 2:
            print(seq_dict,improper_pairs,del_len_ds)

    return seq_dict , new_map


##TODO add disagreeing mip tags to an error file
def make_sure_r1_and_r2_agree(nr1_r2):
    list_of_number_of_barcodes = [0]
    true_sets = []
    # joins the r1 and r2 iff they both have this mnp
    for r1, r2 in nr1_r2:
        true_sets.append(r1.intersection(r2))
        # make sure no barcode overlaps
    for setp_idx in range(len(true_sets)):
        unq_set = true_sets[setp_idx].copy()
        for setc_idx in range(len(true_sets)):
            if setp_idx != setc_idx:
                unq_set.difference_update(true_sets[setc_idx])
        ubx = len(unq_set)
        list_of_number_of_barcodes[0] += ubx
        list_of_number_of_barcodes.append(ubx)

    return list_of_number_of_barcodes


if __name__ == "__main__":
    main() 
