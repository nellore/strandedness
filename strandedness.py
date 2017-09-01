#!/usr/bin/env python2

"""
Python 2.7 code for collecting p-values for strandedness of RNA-seq experiments

Required input:
    .csv File from SRA with SRA accession numbers to check.
    Reference genome for hisat2 alignment

Optional input:
    Path to fastq-dump
    Path to hisat2
    alpha, i.e. acceptable significance level for p-value
    Path to directory to store output files


Improvements to be made:
    - checking for paired-end experiment - is just the "PAIRED" tag OK?
    - How many useful reads are required?
"""

import argparse
import os
import random
import re
import subprocess as sp
from scipy import stats
import time

def extract_sra_data(csv_line):
    '''Returns required data for one SRA experiment.

    Input one line of a SRA "run info" csv file (string).

    Returns the SRA accession number, the total number of reads in the
    sequencing experiment, and whether the experiment layout was paired-end or
    single-end.
    '''
    expt_data = csv_line.split(',')
    accession_num = str(expt_data[0])
    total_reads = int(expt_data[3])
    # Is this "paired" label ever wrong?  Should I be re-checking this
    # with the alignment output?
    if expt_data[15] == 'PAIRED':
        paired_label = True
    else:
        paired_label = False
    return accession_num, total_reads, paired_label


def get_hisat_input(required, multiplier, total, fastq_path, acc, output,
                    pairedtag):
    '''Samples & downloads fastq reads, and prepares them for a hisat2 -12 run.

    Input:
        required: the target number of "useful"/junction reads desired (int)
        multiplier: a multiplier to correct for the fact that most reads are
            not junction reads, and some will not be downloaded due to low
            quality scores (int)
        total: the total number of spots in the SRA experiment available to
            sample (int)
        fastq_path: the path to fastq-dump (string)
        acc: the SRA accession number for the current experiment (string)
        output: the output path for writing out (string)
        paired_tag: whether the experiment is paired or not, only necessary for
            processing multiple reads per spot (true/false)

    A list of unique random numbers between 1 and the total number of spots is
    generated, then the reads at those spots are downloaded with fastq-dump.
    Each downloaded read is then formated appropriately for hisat2 with the -12
    option, which requires the following form: one read or read pair per line,
    tab separated as follows: read name, read 1 alignment, read 1 quality
    scores, then (if applicable) read 2 alignment and read 2 quality scores.

    Writes a file, SRAx_spots.txt, with the sampled spots used to get data.

    Returns all of the sampled reads formatted for batch alignment by hisat2,
    and the list of random spots (to be saved for later reference).
    '''
    dl_time_start = time.time()
    # # # 10,000 reads straight in a row.
    # spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
    # start_spot = random.randint(1, total-11000)
    # spot_file = open(spot_path, 'w')
    # spot_file.write('{}\n'.format(start_spot))
    # spot = str(start_spot)
    # stopspot = str(start_spot + 9999)
    # # If we have x unique spots, we will get <x reads out, since the
    # # -E quality filter rejects some reads before downloading them.
    # fastq = sp.check_output(['{}'.format(fastq_path), '-I', '-B',
    #                          '-W', '-E', '--split-spot',
    #                          '--skip-technical', '-N', spot, '-X',
    #                          stopspot, '-Z', acc])
    # lines = fastq.split('\n')
    # format_lines = []
    # read_format = []
    # for i, line in enumerate(lines, 1):
    #     if i % 2 == 0:
    #         format_lines += [line]
    #     if pairedtag:
    #         if i % 8 == 0:
    #             read_format.extend(['\t'.join(format_lines) + '\n'])
    #             format_lines = []
    #         if (i-1) % 8 == 0:
    #             format_lines += [line]
    #     else:
    #         if i % 4 == 0:
    #             read_format.extend(['\t'.join(format_lines) + '\n'])
    #             format_lines = []
    #         if (i-1) % 4 == 0:
    #             format_lines += [line]
    # read_input = ''.join(read_format)
    # hisat_formatted_input = read_input
    # spot_file.close()
    # dl_time_stop = time.time()
    # print 'download time was {}'.format(dl_time_stop - dl_time_start)
    # return hisat_formatted_input

    # # X bins, 10000/X from each
    required_spots = required * multiplier
    spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
    spot_file = open(spot_path, 'w')
    num_bins = 50
    bin_spots = required_spots/num_bins
    bin_size = total/num_bins
    bin = 1
    # bin_start = 1
    # bin_stop = bin_size
    read_format = []
    while bin <= num_bins:
        binstart = time.time()
        bin_start = (bin - 1) * bin_size + 1
        bin_stop = bin * bin_size
        start_spot = random.randint(bin_start, bin_stop - bin_spots)
        spot_file.write('{}\n'.format(start_spot))
        spot = str(start_spot)
        stop_spot = str(start_spot + bin_spots)
        # If we have x unique spots, we will get <x reads out, since the
        # -E quality filter rejects some reads before downloading them.
        fastq = sp.check_output(['{}'.format(fastq_path), '-I', '-B',
                                 '-W', '-E', '--split-spot',
                                 '--skip-technical', '-N', spot, '-X',
                                 stop_spot, '-Z', acc])
        lines = fastq.split('\n')
        format_lines = []
        for i, line in enumerate(lines, 1):
            if i % 2 == 0:
                format_lines += [line]
            if pairedtag:
                if i % 8 == 0:
                    read_format.extend(['\t'.join(format_lines) + '\n'])
                    format_lines = []
                if (i-1) % 8 == 0:
                    format_lines += [line]
            else:
                if i % 4 == 0:
                    read_format.extend(['\t'.join(format_lines) + '\n'])
                    format_lines = []
                if (i-1) % 4 == 0:
                    format_lines += [line]
        # bin_start += bin_size
        # bin_stop += bin_stop
        bin += 1
        binend = time.time()
        print 'bintime = {}'.format(binend-binstart)
    
    read_input = ''.join(read_format)
    hisat_formatted_input = read_input
    spot_file.close()
    dl_time_stop = time.time()
    print 'download time was {}'.format(dl_time_stop - dl_time_start)
    return hisat_formatted_input

#     # one spot at a time
#     required_spots = required * multiplier
#     spot_list = random.sample(range(total-1), required_spots)
#     spot_list.sort()
#     spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
#     spot_file = open(spot_path, 'w')
#     hisat_formatted_input = ''
#     for start_spot in spot_list:
#         spot_file.write('{}\n'.format(start_spot))
#         spot = str(start_spot)
#         # If we have x unique spots, we will get <x reads out, since the
#         # -E quality filter rejects some reads before downloading them.
#         fastq = sp.check_output(['{}'.format(fastq_path), '-I', '-B',
#                                  '-W', '--fasta', '--split-spot',
#                                  '--skip-technical', '-N', spot, '-X',
#                                  spot, '-Z', acc])
#         # did --fasta instead of -E quality filtering
#         lines = fastq.split('\n')
#         format_lines = lines[:1] + lines[1::2]
#         read_input = '\t'.join(format_lines) + '\n'
#         hisat_formatted_input += read_input
#     spot_file.close()
#     dl_time_stop = time.time()
#     print 'download time was {}'.format(dl_time_stop - dl_time_start)
#     return hisat_formatted_input


def align_sampled_reads(hisat2_path, reference_genome, reads):
    '''Returns SAM format reads aligned by hisat2.

    Input:
    hisat2_path: the path to hisat2 (string)
    reference_genome: the path to the reference genome to be used for
        alignment (string)
    reads: fastq reads correctly formatted for the hisat2 -12 option (string,
        one read per line tab separated, name seq qual if single-end or
        name seq qual seq qual if paired-end)

    Returns a list of aligned reads in SAM format.
    '''
    align_time_start = time.time()
    align_command = ('{h2} -k 1 --no-head --12 - -x {ref}'
                     ).format(h2=hisat2_path, ref=reference_genome)
    align_process = sp.Popen(align_command, stdin=sp.PIPE, stdout=sp.PIPE,
                             shell=True, executable='/bin/bash')
    align_output = align_process.communicate(input=reads)
    alignment = align_output[0].split('\n')[:-1]
    align_time_stop = time.time()
    print 'alignment time was {}'.format(align_time_stop - align_time_start)
    return alignment


def read_is_useful(previous_read, SAM_flag, has_XS_A_tag):
    '''Checks requirements for usability of the current read.

    Input
    previous read: True if the previous read was useful.
    SAM_flag (int)
    has_XS_A_tag: True if the alignment has this tag.

    If this is a paired read and its pair was a useful junction read, then this
    one can't be counted (to avoid "unfair" double counting).  If this read
    represents a non-primary alignment, or is not a junction read, then it is
    not useful.

    Returns True if the read is useful, otherwise returns False.
    '''
    secondary_read = SAM_flag & 256
    paired = SAM_flag & 1
    if has_XS_A_tag:
        if not secondary_read:
            if not (paired and previous_read):
                return True
    return False


def read_sense(SAM_flag, plus_or_minus):
    '''Checks a read's SAM flag and XS:A: tag, and determines its "direction".

    Input the SAM flag (int) and XS:A:? tag (string) from the aligned read.

    We have two states, arbitrarily called "sense" and "antisense," indicating
    whether all the first/second reads align with a gene or with its reverse
    complement, or whether this is random.

    Return the read's "sense"ness - if bit = 1, the read is "sense", otherwise
    it is "antisense."
    '''
    fwd_gene = plus_or_minus[0] == 'XS:A:+'
    paired = SAM_flag & 1
    rev_read = SAM_flag & 16
    first_read = SAM_flag & 64
    return (fwd_gene is not rev_read) is (first_read or not paired)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--srafile', '-s', required=True, help='File with '
                        'SRA accession numbers to be downloaded and '
                        'checked for strandedness.')
    parser.add_argument('--refgenome', '-r', required=True, help='Path to '
                        'reference genome to use for the aligner.')
    parser.add_argument('--fastqdumppath', '-f', default='fastq-dump',
                        help='specify the path for fastq-dump')
    parser.add_argument('--alignerpath', '-p', default='hisat2',
                        help='specify the path for hisat2')
    parser.add_argument('--alpha', '-a', type=float, default=0.05,
                        help='acceptable significance level of strandedness '
                        'call; must be between 0 and 1.')
    parser.add_argument('--outputpath', '-o', default='./', help='give path '
                        'for output files: sampled spots, aligned junction '
                        'reads, and SRA numbers with their p-values.')
    parser.add_argument('--requiredreads', '-n', type=int, default=50, help=''
                        'give the target number of useful/junction reads.')
    parser.add_argument('--multiplier', '-m', type=int, default=15, help='a '
                        'multiplier for generating the number of reads to '
                        'download from fastq-dump, to account for quality '
                        'filtering and for not all reads being useful.')

    args = parser.parse_args()
    alpha = args.alpha
    sra_file = args.srafile
    ref_genome = args.refgenome
    fastq_dump = args.fastqdumppath
    hisat2 = args.alignerpath
    out_path = args.outputpath
    required_reads = args.requiredreads
    read_multiplier = args.multiplier

    pval_file_path = os.path.join(out_path, 'SRA_pvals.txt')
    pval_file = open(pval_file_path, 'w')

    pval_list = []
    sra_array = open(sra_file, 'r')
    SRA_num = 0


    # sra_sample = random.sample(range(2499), 30)

    for i, experiment in enumerate(sra_array):
        if i == 0:
            continue
            
        sra_acc, num_spots, paired = extract_sra_data(experiment)

        # For a sample run: with a large .csv file, this will check X SRAs.
        if SRA_num >= 3:
            break
        SRA_num += 1

        reads_path = os.path.join(out_path, '{}_reads.sam'.format(sra_acc))
        reads_file = open(reads_path, 'w', 1)
        hisat_input = get_hisat_input(required_reads, read_multiplier,
                                      num_spots, fastq_dump, sra_acc, out_path,
                                      paired)
        aligned_reads = align_sampled_reads(hisat2, ref_genome, hisat_input)

        sense = 0
        checked_reads = 0
        useful = True
        last_read = not useful

        for entry in aligned_reads:
            read = entry.split()
            flag = int(read[1])
            XS_A_tag = re.findall('XS:A:[+-]', entry)
            if read_is_useful(last_read, flag, XS_A_tag):
                sense += read_sense(flag, XS_A_tag)
                checked_reads += 1
                last_read = useful
                reads_file.write('{}\n'.format(entry))
            else:
                last_read = not useful

        antisense = checked_reads - sense
        r = 0.5
        # 2-sided & symmetrical: same result whether we pick sense or antisense
        p_value = stats.binom_test(sense, checked_reads, r)
        print '\nthe SRA accession number is {}'.format(sra_acc)
        print 'number of sense reads is {}'.format(sense)
        print 'number of antisense reads is {}'.format(antisense)
        print 'total number of junction reads is {}'.format(checked_reads)
        print 'The unstranded p-value is {}'.format(p_value)
        print '\nmoving on to the next SRA accession number \n'

        pval_list.append([sra_acc, p_value])
        pval_file.write('{},{}\n'.format(sra_acc, p_value))
        reads_file.close()

    pval_file.close()
    sra_array.close()
