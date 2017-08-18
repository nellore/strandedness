#!/usr/bin/env python2

"""
Python 2.7 code for determining whether RNA-seq experiments are stranded or
unstranded, to within a given confidence..

Required input:
    .csv File from SRA with SRA accession numbers to check.
    Reference genome for hisat2 alignment

Optional input:
    Path to fastq-dump
    Path to hisat2
    Acceptable probability of correctness for strandedness call.
    Acceptable confidence of that probability.

Limitations:
    Many!!

Improvements to be made:
    - Improve hisat2 call - either take two reads as separate, unpaired reads
    - Improve checking for paired-end experiment - not just "PAIRED" tag
    - if the reads are paired, check both of them, not just one.
    - Put meat of code into function definitions
    - How many reads are required?  Unknown.
    - Statistical analysis of reads afterward.

"""

# import what needs to be imported.
import argparse
import math
import numpy as np
import random
import re
import subprocess as sp
from scipy import stats

# Function definitions?


if __name__ == '__main__':
    # arg parse section: for now, accept sam; later, sra numbers?
    parser = argparse.ArgumentParser(description="Determine sample "
                                                 "strandedness.")
    parser.add_argument('--srafile', '-s', required=True, help='File with '
                        'SRA accession numbers to be downloaded and '
                        'checked for strandedness.')
    parser.add_argument('--refgenome', '-r', required=True, help='Path to '
                        'reference genome to use for the aligner.')
    parser.add_argument('--fastqdumppath', '-f', default='fastq-dump',
                        help='specify the path for fastq-dump')
    parser.add_argument('--alignerpath', '-a', default='hisat2',
                        help='specify the path for hisat2')
    parser.add_argument('--probability', '-p', type=float, default=0.95,
                        help='acceptable probability of strandedness call; '
                             'must be between 0 and 1.')
    parser.add_argument('--confidence', '-c', type=float, default=0.95,
                        help='minimum confidence level of the call '
                             'probability; must be between 0 and 1.')

    args = parser.parse_args()
    probability = args.probability
    confidence = args.confidence
    sra_file = args.srafile
    ref_genome = args.refgenome
    fastq_dump = args.fastqdumppath
    hisat2 = args.alignerpath

    # arbitrary (??) alignment error rate for binomial distr.
    error_rate = 0.2

    # for now, arbitrary: how many reads need to be sampled per experiment
    # Currently set at a low number simply to check for correct running.
    required_reads = 10

    sra_array = np.genfromtxt(sra_file, delimiter=',', dtype=None)

    SRA_num = 0
    for line in sra_array:
        if line[0][0] == 'S':
            sra_acc = str(line[0])
            num_spots = int(line[3])
        else:
            continue

        # THE FOLLOWING IS JUST TO CREATE A SAMPLE RUN: with a large .csv file,
        # this will check X of the SRAs.
        if SRA_num > 1:
            break

        SRA_num += 1

        # Is the "paired" label ever wrong?  Should I be re-checking this
        # with the alignment output?
        paired = 0
        if line[15] == 'PAIRED':
            paired = 1
            print 'experiment is paired'

        # Counters
        sense = 0
        antisense = 0
        checked_reads = 0
        total_reads = 0

        while checked_reads < required_reads:
            print '\nnot enough reads for this SRA, starting at a random spot'
            useful = 0
            start_spot = random.randint(1, num_spots)
            while not useful:
                # get one read or read pair
                spot = str(start_spot)
                start_spot += 1
                print 'the current spot is {}'.format(spot)
                fastq = sp.check_output(['{}'.format(fastq_dump), '-I', '-B',
                                         '-W', '-E', '--split-spot',
                                         '--skip-technical', '-N', spot, '-X',
                                         spot, '-Z', sra_acc])
                lines = fastq.split('\n')
                format_lines = lines[:1]+lines[1::2]
                hisat_input = '\t'.join(format_lines)

                align_command = ('{h2} -k 1 --no-head --12 - -x {ref}'
                                 ).format(h2=hisat2, ref=ref_genome)

                align_process = sp.Popen('set -exo pipefail; ' + align_command,
                                         stdin=sp.PIPE, stdout=sp.PIPE,
                                         shell=True, executable='/bin/bash')

                align_output = align_process.communicate(input=hisat_input)

                entries = align_output[0].split('\n')
                for line in entries[:-1]:
                    print line
                    read = line.split()
                    total_reads += 1
                    print "{} reads have been checked".format(total_reads)
                    flag = int(read[1])

                    secondary = flag & 256
                    if secondary:
                        print 'not a primary read, go to next spot'
                        continue

                    check_sense = re.findall('XS:A:[+-]', line.upper())
                    if not check_sense:
                        print 'not a junction read, go to next spot'
                        continue

                    fwd_gene = check_sense[0] == 'XS:A:+'
                    rev_read = flag & 16
                    first_read = flag & 64
                    second_read = flag & 128

                    if ((fwd_gene and rev_read)
                        or (not fwd_gene and not rev_read)):
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
                    useful = 1
                    print 'a useful read! {} so far'.format(checked_reads)
                    print 'moving on to the next random spot.'
                    break

        print '\nfor SRA experiment number {}'.format(SRA_num)
        print 'the SRA accession number is {}'.format(sra_acc)
        print 'number of sense reads is {}'.format(sense)
        print 'number of antisense reads is {}'.format(antisense)
        print 'total number of checked reads is {}'.format(checked_reads)
        r = error_rate * 0.5
        # Now we have min_reads processed, X as sense and Y as antisense.
        # Check probability:
        errors = min(sense, antisense)
        stranded_prob = stats.binom.pmf(errors, checked_reads, r)

        r = 0.5
        unstranded_prob = stats.binom.pmf(errors, checked_reads, r)


        print 'number of "stranded errors" is {}'.format(errors)
        print ('probability of getting this number of sense and antisense'
               'reads if the expt is stranded is {}').format(stranded_prob)

        unstranded_errors = int(abs(sense-antisense)/2)
        print 'number of unstranded errors is {}'.format(unstranded_errors)
        print ('probability of getting this number of sense and antisense '
               'reads if the expt is unstranded is {}').format(unstranded_prob)

        print '\nmoving on to the next SRA accession number \n'
