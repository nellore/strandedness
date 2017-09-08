#!/usr/bin/env python2

"""
Python 2.7 code for collecting p-values for strandedness of RNA-seq experiments

Required input:
    .csv File from SRA with SRA accession numbers to check.
    Reference genome for hisat2 alignment

Optional input:
    Path to fastq-dump
    Path to hisat2
    Path to directory to store output files
    Target number of useful reads expected per experiment.
    The multiplier to translate target number of useful end reads to number
        of reads required to download initially.
    Number of alignment attempts to try per experiment before moving on.
    Log level for logger.

Improvements to be made:
    - checking for paired-end experiment - is just the "PAIRED" tag OK?
    - How many useful reads are required?
"""

import argparse
import csv
from datetime import datetime
import gzip
from itertools import groupby
import logging
import os
import random
import re
import subprocess as sp
from scipy import stats


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
    dl_time_start = datetime.now()
    spot_path = os.path.join(output, '{}_spots.txt'.format(acc))
    with open(spot_path, 'w') as spot_file:
        required_spots = required * multiplier
        num_bins = 100
        bin_spots = required_spots/num_bins
        bin_size = total/num_bins
        bin = 1
        read_format = []
        while bin <= num_bins:
            bin_start = (bin - 1) * bin_size + 1
            bin_stop = bin * bin_size
            start_spot = random.randint(bin_start, bin_stop - bin_spots)
            spot_file.write('{}\n'.format(start_spot))
            spot = str(start_spot)
            stop_spot = str(start_spot + bin_spots)
            fastq = sp.check_output(['{}'.format(fastq_path), '-I', '-B',
                                     '-W', '-E', '--split-spot',
                                     '--skip-technical', '-N', spot, '-X',
                                     stop_spot, '-Z', acc])
            lines = fastq.split('\n')
            format_lines = []
            if pairedtag:
                for i, line in enumerate(lines, 1):
                    if i % 2 == 0:
                        format_lines += [line]
                    if i % 8 == 0:
                        read_format.extend(['\t'.join(format_lines) + '\n'])
                        format_lines = []
                    if (i-1) % 8 == 0:
                        format_lines += [line]
            else:
                for i, line in enumerate(lines, 1):
                    if i % 2 == 0:
                        format_lines += [line]
                    if i % 4 == 0:
                        read_format.extend(['\t'.join(format_lines) + '\n'])
                        format_lines = []
                    if (i-1) % 4 == 0:
                        format_lines += [line]
            bin += 1
        read_input = ''.join(read_format)
        hisat_formatted_input = read_input
        dl_time_stop = datetime.now()
        time_difference = dl_time_stop - dl_time_start
        elapsed_time = time_difference.total_seconds()
        logging.info('total download time was {} seconds'.format(elapsed_time))
        return hisat_formatted_input


def align_reads(hisat2_path, reference_genome, outpath, acc, reads):
    '''Returns SAM format reads aligned by hisat2.

    Input:
        hisat2_path: the path to hisat2 (string)
        reference_genome: the path to the reference genome to be used for
            alignment (string)
        reads: fastq reads correctly formatted for the hisat2 -12 option
            (string, one read per line tab separated, name seq qual if
            single-end or name seq qual seq qual if paired-end)
        temp_dir: dir in which to put temporary fastq

    Returns a list of aligned reads in SAM format.
    '''
    align_start = datetime.now()
    reads_path = os.path.join(outpath, '{}_reads.sam.gz'.format(acc))
    align_command = ('set -exo pipefail; {h2} --no-head --12 - -x {ref} | gzip'
                     ' > {file}').format(h2=hisat2_path, ref=reference_genome,
                                         file=reads_path)
    align_process = sp.Popen(align_command, stdin=sp.PIPE, stderr=sp.PIPE,
                             shell=True, executable='/bin/bash')
    align_process.communicate(input=reads)

    align_end = datetime.now()
    time_difference = align_start - align_end
    elapsed_time = time_difference.total_seconds()
    logging.info('total alignment time was {}'.format(elapsed_time))
    return reads_path, align_process.returncode


def read_is_useful(SAM_flag, has_XS_A_tag):
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
    if has_XS_A_tag:
        if not secondary_read:
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
    parser.add_argument('--sra-file', '-s', required=True, help='File with '
                        'SRA accession numbers to be downloaded and '
                        'checked for strandedness.')
    parser.add_argument('--ref-genome', '-r', required=True, help='Path to '
                        'reference genome to use for the aligner.')
    parser.add_argument('--fastq-dump-path', '-f', default='fastq-dump',
                        help='specify the path for fastq-dump')
    parser.add_argument('--aligner-path', '-p', default='hisat2',
                        help='specify the path for hisat2')
    parser.add_argument('--output-path', '-o', default='./', help='give path '
                        'for output files: sampled spots, aligned junction '
                        'reads, and SRA numbers with their p-values.')
    parser.add_argument('--required-reads', '-n', type=int, default=100, 
                        help='give the target number of useful reads.')
    parser.add_argument('--multiplier', '-m', type=int, default=10, help='a '
                        'multiplier for generating the number of reads to '
                        'download from fastq-dump, to account for quality '
                        'filtering and for not all reads being useful.')
    parser.add_argument('--max-attempts', '-a', type=int, default=3, help='the'
                        ' number of times to attempt hisat2 alignment one one '
                        'SRA accession number after alingment failure before '
                        'continuing.')
    parser.add_argument('--log-level', '-l', choices=['DEBUG', 'INFO', 'ERROR'
                                                      'WARNING', 'CRITICAL'],
                        default='INFO', help='choose what logging mode to run')

    args = parser.parse_args()
    sra_file = args.sra_file
    ref_genome = args.ref_genome
    fastq_dump = args.fastq_dump_path
    hisat2 = args.aligner_path
    out_path = args.output_path
    required_reads = args.required_reads
    read_multiplier = args.multiplier
    max_attempts = args.max_attempts
    log_mode = args.log_level

    time_now = str(datetime.now())
    log_file = os.path.join(out_path, '{}_log_file.txt'.format(time_now))
    logging.basicConfig(filename=log_file, level=log_mode)

    read_shuffle_seed = 1
    random.seed(read_shuffle_seed)
    SRA_num = 0
    name_tag = sra_file.split('/')[-1].split('.')[0]
    pv_path = os.path.join(out_path, 'SRA_pvals_{}.txt'.format(name_tag))
    with open(sra_file) as sra_array, open(pv_path, 'w', 0) as pval_file:
        sra_array.next()
        csv_reader = csv.reader(sra_array, delimiter=',', quotechar='"')
        for experiment in csv_reader:
            sra_acc, num_spots, paired = (experiment[0], int(experiment[3]),
                                          experiment[15] == 'PAIRED')
            hisat_input = get_hisat_input(required_reads, read_multiplier,
                                          num_spots, fastq_dump, sra_acc,
                                          out_path, paired)
            attempt = 0
            while attempt < max_attempts:
                reads_path, failure = align_reads(hisat2, ref_genome, out_path,
                                                  sra_acc, hisat_input)
                attempt += 1
                if not failure:
                    break
            else:
                logging.info('alignment failed for SRA {}'.format(sra_acc))
                continue
                
            sense = 0
            checked_reads = 0
            with gzip.open('{}'.format(reads_path), 'r') as aligned_reads:
                sort_by_names = lambda x: x.split('\t')[0]
                for key, alignments in groupby(aligned_reads, sort_by_names):
                    alignment_list = list(alignments)
                    random.shuffle(alignment_list)
                    for entry in alignment_list:
                        read = entry.split('\t')
                        flag = int(read[1])
                        XS_A_tag = re.findall('XS:A:[+-]', entry)
                        if read_is_useful(flag, XS_A_tag):
                            sense += read_sense(flag, XS_A_tag)
                            checked_reads += 1
                            break

            antisense = checked_reads - sense
            r = 0.5
            # 2-sided & symmetrical: same result whether we pick sense or
            # antisense
            p_value = stats.binom_test(sense, checked_reads, r)
            logging.info('The SRA accession number is {}'.format(sra_acc))
            logging.info('{} sense reads.'.format(sense))
            logging.info('{} antisense reads.'.format(antisense))
            logging.info('{} junction reads.'.format(checked_reads))
            logging.info('The unstranded p-value is {}\n'.format(p_value))

            print >>pval_file, '{},{}'.format(sra_acc, p_value)
            # pval_file.write('{},{}\n'.format(sra_acc, p_value))
