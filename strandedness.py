#!/usr/bin/env python2
# coding=utf-8

"""
Python 2.7 code for determining whether an RNA-seq experiment is stranded or
unstranded, to within a given confidence..

Required input:
    SAM file to check.

Optional input:
    Acceptable probability of correctness for strandedness call.
    Acceptable confidence of that probability.

Limitations:
    Many!!

Improvements to be made:
    - Input = SRA accession numbers
    - Start fastq-dump of SRA at one random spot.
    - align each read as it's run many in a row.

"""

# import what needs to be imported.
import argparse
import math
import numpy as np
import random
import re
import subprocess
from scipy import stats


# function definitions?


if __name__ == '__main__':


    # arg parse section: for now, accept sam; later, sra numbers?
    parser = argparse.ArgumentParser(description="Determine sample "
                                                 "strandedness.")
    parser.add_argument('--samfile', '-s', required=True, help='SAM file of '
                        ' the experiment to be checked for strandedness.')
    parser.add_argument('--probability', '-p', type=float, default=0.95,
                        help='acceptable probability of strandedness call; '
                             'must be between 0 and 1.')
    parser.add_argument('--confidence', '-c', type=float, default=0.95,
                        help='minimum confidence level of the call '
                             'probability; must be between 0 and 1.')

    args = parser.parse_args()
    probability = args.probability
    confidence = args.confidence
    samfile = args.samfile

    if confidence <= 0 or confidence >= 1:
        print "Please retry with desired confidence level between 0 and 1."
        exit()

    if probability <= 0 or probability >= 1:
        print "Please retry with desired probability level between 0 and 1."
        exit()

    # Later: begin with sra accession numbers and run alignment first.

    # arbitrary (?) alignment error rate
    error_rate = 0.2

    # min_reads = 7000000
    min_reads = 1600
    # min_reads = 700

    # collect min_reads reads from sam file: samtools -s?  reservoir sample?
    print '\nSampling {} reads...'.format(min_reads)
    wc_output = subprocess.check_output(['wc', '-l', samfile])

    for s in wc_output.split():
        if s.isdigit():
            sam_length = float(s)

    ratio = min_reads/sam_length
    add_seed = ratio + random.randint(0, 1000000)
    sample_seed = '{}'.format(add_seed)
    sampled_reads = subprocess.check_output(['samtools', 'view', '-s',
                                             sample_seed, samfile])

    cutoff = sampled_reads[:-2]
    read_lines = cutoff.split('\n')
    #
    # read_lines = sampled_reads.split('\n')
    # number_reads = len(read_lines)
    # read_lines = read_lines[0:number_reads-2]

    # check number of sampled reads...
    print 'The length of the sam file is {} lines.'.format(int(sam_length))
    print 'The ratio should be {}.'.format(ratio)
    print 'The minimum number of reads should be {}.'.format(min_reads)
    print 'The number of sampled reads is {}.'.format(len(read_lines))

    # print "type of sampled_reads is {}.".format(type(sampled_reads))
    # print 'type of read_lines is {}.'.format(type(read_lines))
    # read = read_lines[0]
    # entries = read.split()
    # print 'qname is {}'.format(entries[0])
    # print 'FLAG is {}'.format(entries[1])
    # print 'Rname (ref seq name, or * for unmapped) is {}'.format(entries[2])
    # print 'position is {}'.format(entries[3])
    # print 'mapping quality is {} (255 = unavailable)'.format(entries[4])
    # print 'Cigar string is {}'.format(entries[5])
    # print 'next read is named {} (* indicates blank)'.format(entries[6])
    # print 'position of next read is {} unless 0'.format(entries[7])
    # print 'template length is {}'.format(entries[8])
    # print 'sequence is {}'.format(entries[9])
    # print ' quality is {}'.format(entries[10])

    # # Check sample for being paired end or single end
    print '\nChecking sample for paired- or single-endedness...'
    paired = 0
    read = read_lines[0]
    entries = read.split()
    flag = int(entries[1])
    if flag & 1:
        paired = 1
    if paired:
        print 'The sample is paired-end.'
    else:
        print 'The sample is single end.'


    print '\nChecking sampled reads: sense or antisense?'
    sense = 0
    antisense = 0
    checked_reads = 0
    total_reads = 0
    print_num = 2000

    for read in read_lines:
        entries = read.split()
        total_reads += 1
        flag = int(entries[1])

        secondary = flag & 256
        check_sense = re.findall('XS:A:[+-]', read.upper())
        if secondary or not check_sense:
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

    unstranded_errors = int(abs(sense - antisense) / 2)
    print 'number of unstranded errors is {}'.format(unstranded_errors)
    print ('probability of getting this number of sense and antisense '
           'reads if the expt is unstranded is {}').format(unstranded_prob)


    #
    # mean = stranded_errors/checked_reads
    # std = np.std(checked_reads)
    #
    # scaled_std = std/
    #
    # # Calculate confidence interval for current sample.
    # # conf_int = stats.norm.interval(confidence,loc=1,scale=std/math.sqrt(len(s)))
    # conf_int = stats.norm.interval(confidence,loc=mean,scale=std/math.sqrt(len(s)))
    #
