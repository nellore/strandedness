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
    - function definitions
    - How many reads are required?  Unknown.
"""

# import what needs to be imported.
import argparse
import numpy as np
import os
import random
import re
import subprocess as sp
from scipy import stats
import time


# Function definitions?


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
                        'for output files: junction reads and p-values.')

    args = parser.parse_args()
    alpha = args.alpha
    sra_file = args.srafile
    ref_genome = args.refgenome
    fastq_dump = args.fastqdumppath
    hisat2 = args.alignerpath
    out_path = args.outputpath

    pval_file_path = os.path.join(out_path, 'SRA_pvals.txt')
    pval_file = open(pval_file_path, 'w')

    # How many reads need to be sampled per experiment?
    # Currently set low to check for correct running.
    required_reads = 50
    pval_list = []
    sra_array = np.genfromtxt(sra_file, delimiter=',', dtype=None)
    SRA_num = 0
    for experiment in sra_array:
        start_time = time.time()
        if experiment[0][0] == 'S':
            sra_acc = str(experiment[0])
            num_spots = int(experiment[3])
        else:
            continue

        # For a sample run: with a large .csv file, this will check X SRAs.
        if SRA_num >= 5:
            break
        SRA_num += 1
        reads_path = os.path.join(out_path, '{}_reads.sam'.format(sra_acc))
        reads_file = open(reads_path, 'w')

        # Is the "paired" label ever wrong?  Should I be re-checking this
        # with the alignment output?
        paired = 0
        if experiment[15] == 'PAIRED':
            paired = 1
            print 'experiment is paired'

        # collect random spots - 15x required reads number
        required_spots = required_reads * 15
        spot_list = random.sample(range(num_spots), required_spots)

        dl_start = time.time()
        hisat_input = ''
        for start_spot in spot_list:
            spot = str(start_spot)

            # If we have x unique spots, we will get <x reads out, since the
            # -E quality filter rejects some reads before downloading them.
            fastq = sp.check_output(['{}'.format(fastq_dump), '-I', '-B',
                                     '-W', '-E', '--split-spot',
                                     '--skip-technical', '-N', spot, '-X',
                                     spot, '-Z', sra_acc])
            lines = fastq.split('\n')
            format_lines = lines[:1]+lines[1::2]
            read_input = '\t'.join(format_lines) + '\n'
            hisat_input += read_input

        dl_stop = time.time()
        print 'the time for read download was {}'.format(dl_stop-dl_start)

        align_command = ('{h2} -k 1 --no-head --12 - -x {ref}'
                         ).format(h2=hisat2, ref=ref_genome)
        align_process = sp.Popen(align_command, stdin=sp.PIPE, stdout=sp.PIPE,
                                 shell=True, executable='/bin/bash')
        align_output = align_process.communicate(input=hisat_input)
        aligned_reads = align_output[0].split('\n')

        sense = 0
        antisense = 0
        checked_reads = 0
        useful = 0
        for entry in aligned_reads[:-1]:
            # if the previous read was useful skip this one in case
            # it is its pair....
            if useful:
                useful = 0
                continue

            read = entry.split()
            flag = int(read[1])

            secondary = flag & 256
            if secondary:
                print 'not a primary read, go to next spot'
                continue

            check_sense = re.findall('XS:A:[+-]', entry)
            if not check_sense:
                print 'not a junction read, go to next spot'
                continue

            fwd_gene = check_sense[0] == 'XS:A:+'
            rev_read = flag & 16
            first_read = flag & 64
            second_read = flag & 128

            if (fwd_gene and rev_read) or (not fwd_gene and not rev_read):
                if first_read or not paired:
                    sense += 1
                elif second_read:
                    antisense += 1
            elif (fwd_gene or rev_read):
                if first_read or not paired:
                    antisense += 1
                elif second_read:
                    sense += 1

            checked_reads = sense + antisense
            print 'a useful read! {} so far'.format(checked_reads)
            useful = 1
            reads_file.write('{}\n'.format(entry))

        print '\nthe SRA accession number is {}'.format(sra_acc)
        print 'number of sense reads is {}'.format(sense)
        print 'number of antisense reads is {}'.format(antisense)
        print 'total number of junction reads is {}'.format(checked_reads)

        r = 0.5
        # 2-sided & symmetrical: same result whether we pick sense or antisense
        p_value = stats.binom_test(sense, checked_reads, r)

        print 'The unstranded p-value is {}'.format(p_value)
        print '\nmoving on to the next SRA accession number \n'

        pval_list.append((sra_acc, p_value))
        pval_file.write('{}, {}\n'.format(sra_acc, p_value))
        reads_file.close()
        end_time = time.time()
        print "elapsed time for this SRA was {}".format(end_time-start_time)

    pval_file.close()

    # Benjamini-Hochberg multiple-testing correction
    # Also available in a separate script, BH_correction.py
    m = len(pval_list)
    stranded_expts = []
    stranded_file = open('stranded_SRAs.txt', 'w')
    for k, expt in enumerate(pval_list, 1):
        if expt[1] <= alpha * k / m:
            stranded_expts.append(expt)
            print '{} expt: stranded with p value {}.'.format(expt[0], expt[1])
            stranded_file.write('{}\n'.format(expt))
        else:
            print 'BH correction complete!  Stranded list is finished.'
            break
