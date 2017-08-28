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


def extract_sra_data(csv_line):
    '''Returns required data for one SRA experiment.

    Input one line of a SRA "run info" csv file.

    Returns the SRA accession number, the total number of reads in the
    sequencing experiment, and whether the experiment layout was paired-end or
    single-end.
    '''
    expt_data = csv_line.split(',')
    accession_num = str(expt_data[0])
    total_reads = int(expt_data[3])
    # Is this "paired" label ever wrong?  Should I be re-checking this
    # with the alignment output?
    if csv_line[15] == 'PAIRED':
        paired_label = 1
    else:
        paired_label = 0
    return accession_num, total_reads, paired_label


def get_hisat_input(required, multiplier, total, fastq_path, acc, output):
    '''Samples & downloads fastq reads, and prepares them for a hisat2 -12 run.

    Input:
        - the target number of "useful"/junction reads desired
        - a multiplier to correct for the fact that most reads are not
            junction reads, and some will not be downloaded due to low
            quality scores
        - the total number of spots in the SRA experiment available to sample
        - the path to fastq-dump
        - the SRA accession number for the current experiment.
        - the output path for writing out.

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
    required_spots = required * multiplier
    spot_list = random.sample(range(total), required_spots)
    spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
    spot_file = open(spot_path, 'w')
    hisat_formatted_input = ''
    for start_spot in spot_list:
        spot_file.write('{}\n'.format(start_spot))
        spot = str(start_spot)
        # If we have x unique spots, we will get <x reads out, since the
        # -E quality filter rejects some reads before downloading them.
        fastq = sp.check_output(['{}'.format(fastq_path), '-I', '-B',
                                 '-W', '-E', '--split-spot',
                                 '--skip-technical', '-N', spot, '-X',
                                 spot, '-Z', acc])
        lines = fastq.split('\n')
        format_lines = lines[:1] + lines[1::2]
        read_input = '\t'.join(format_lines) + '\n'
        hisat_formatted_input += read_input
    spot_file.close()
    return hisat_formatted_input


def align_sampled_reads(hisat2_path, reference_genome, reads):
    '''Returns SAM format reads aligned by hisat2.

    Input: the path to hisat2, the path to the reference genome to be used for
    alignment, and fastq reads correctly formatted for the hisat2 -12 option.

    Returns a list of aligned reads in SAM format.
    '''
    align_command = ('{h2} -k 1 --no-head --12 - -x {ref}'
                     ).format(h2=hisat2_path, ref=reference_genome)
    align_process = sp.Popen(align_command, stdin=sp.PIPE, stdout=sp.PIPE,
                             shell=True, executable='/bin/bash')
    align_output = align_process.communicate(input=reads)
    alignment = align_output[0].split('\n')[:-1]
    return alignment


def read_is_useful(previous_read, SAM_flag, has_XS_A_tag):
    '''Checks requirements for usability of the current read.

    If this is a paired read and its pair was a useful junction read, then this
    one can't be counted (to avoid "unfair" double counting).  If this read
    represents a non-primary alignment, or is not a junction read, then it is
    not useful.

    Returns True if the read is useful, otherwise returns False.
    '''
    secondary_read = SAM_flag & 256
    paired = SAM_flag & 1
    if (paired and previous_read) or secondary_read or not has_XS_A_tag:
        return False
    else:
        return True


def read_sense(SAM_flag, plus_or_minus):
    '''Checks a read's SAM flag and XS:A: tag, and determines its "direction".

    Input the SAM flag and XS:A:? tag from the aligned read.

    We have two states, arbitrarily called "sense" and "antisense," indicating
    whether all the first/second reads align with a gene or with its reverse
    complement, or whether this is random.

    Return the read's "sense"ness - if bit = 1, the read is "sense", otherwise
    it is "antisense."
    '''
    fwd_gene = plus_or_minus[0] == 'XS:A:+'
    rev_read = SAM_flag & 16
    first_read = SAM_flag & 64
    second_read = SAM_flag & 128
    if (fwd_gene and rev_read) or (not fwd_gene and not rev_read):
        if first_read or not paired:
            bit = 1
        elif second_read:
            bit = 0
    elif (fwd_gene or rev_read):
        if first_read or not paired:
            bit = 0
        elif second_read:
            bit = 1
    return bit


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
    parser.add_argument('--requiredreads', '-n', default = 50, help='give the'
                        'target number of useful/junction reads.')

    args = parser.parse_args()
    alpha = args.alpha
    sra_file = args.srafile
    ref_genome = args.refgenome
    fastq_dump = args.fastqdumppath
    hisat2 = args.alignerpath
    out_path = args.outputpath
    required_reads = args.requiredreads

    pval_file_path = os.path.join(out_path, 'SRA_pvals.txt')
    pval_file = open(pval_file_path, 'w')
    pval_file.write('SRA accession number, p-value\n')

    read_multiplier = 15
    pval_list = []
    sra_array = open(sra_file, 'r')
    SRA_num = 0

    for experiment in sra_array:
        if experiment.startswith('S'):
            sra_acc, num_spots, paired = extract_sra_data(experiment)
        else:
            continue

        # For a sample run: with a large .csv file, this will check X SRAs.
        if SRA_num >= 3:
            break
        SRA_num += 1

        reads_path = os.path.join(out_path, '{}_reads.sam'.format(sra_acc))
        reads_file = open(reads_path, 'w', 0)
        hisat_input = get_hisat_input(required_reads, read_multiplier,
                                      num_spots, fastq_dump, sra_acc, out_path)
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
        pval_list.append([sra_acc, p_value])
        pval_file.write('{},{}\n'.format(sra_acc, p_value))
        reads_file.close()

    pval_file.close()
    sra_array.close()
