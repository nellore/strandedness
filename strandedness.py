#!/usr/bin/env python2

"""
Python 2.7 code for determining whether an RNA-seq experiment is stranded or
unstranded, to within a given confidence..

Required input:
    SAM file to check.

Optional input:
    Acceptable probability of correctness for strandness call.
    Acceptable confidence of that probability.

Limitations:
    Many!!

Improvements to be made:
    - Change to reservoir sampling instead of samtools?
        ->  Take fastq file input and run through an aligner.
            -> Take SRA accession numbers and download and run many in a row.

    - FIX THE SAMPLE SIZE ESTIMATE

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

    # Calculate for error_rate, minimum number of reads to give confidence
    #
    # COME BACK TO THIS: confusion on std_dev stuff.
    #       Standard deviation of read results, or of error?
    #
    print 'Checking for minimum number of required reads...'
    num_errors = int(math.ceil(1000*error_rate))
    sample_array = np.ones(1000)
    sample_array[:num_errors] = 0
    std_dev = np.std(sample_array)
    z = stats.norm.ppf(1 - ((1-confidence)/2))
    min_reads = int(math.ceil((z*std_dev / ((1-probability)/2)) ** 2))

    # # ALternate for now:
    # min_reads = 7000000
    # min_reads = 1600

    # collect min_reads reads from sam file: samtools -s?  reservoir sample?
    print '\nSampling reads...'
    wc_output = subprocess.check_output(['wc', '-l', samfile])

    for s in wc_output.split():
        if s.isdigit():
            sam_length = float(s)

    ratio = min_reads/sam_length
    add_seed = ratio + random.randint(0, 1000000)
    sample_seed = '{}'.format(add_seed)
    sampled_reads = subprocess.check_output(['samtools', 'view', '-s',
                                             sample_seed, samfile])
    read_lines = sampled_reads.split('\n')
    number_reads = len(read_lines)
    read_lines = read_lines[0:number_reads-2]

    # check number of sampled reads...
    print 'The length of the sam file is {} lines.'.format(int(sam_length))
    print 'The ratio should be {}.'.format(ratio)
    print 'The minimum number of reads should be {}.'.format(min_reads)
    print 'The number of sampled reads is {}.'.format(number_reads-1)

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
        print (flag & 1)
        paired += 1
    if paired:
        print 'The sample is paired-end.'
    else:
        print 'The sample is single end.'

    # Reads:
    # 	- Is it a primary alignment?  (position 2: not 256 in the bitmask

    print '\nChecking sampled reads: sense or antisense?'
    sense = 0
    antisense = 0
    checked_reads = 0
    if not paired:
        for read in read_lines:
            # Check for primary alignment?  Not done yet....
            check_sense = re.findall('XS:A:[+-]', read.upper())
            if not check_sense:
                continue
            elif check_sense[0] == 'XS:A:+':
                sense += 1
                checked_reads += 1
            elif check_sense[0] == 'XS:A:-':
                antisense += 1
                checked_reads += 1
            else:
                continue
    else:
        for read in read_lines:
            entries = read.split()
            flag = int(entries[1])
            check_sense = re.findall('XS:A:[+-]', read.upper())
            if not check_sense:
                continue
            elif check_sense[0] == 'XS:A:+':
                if flag & 64:
                    # indicates first read; arbitrarily call these forward.
                    sense += 1
                    checked_reads += 1
                elif flag & 128:
                    antisense += 1
                    checked_reads += 1
                else:
                    continue
            elif check_sense[0] == 'XS:A:-':
                if flag & 64:
                    # indicates first read; arbitrarily = main direction.
                    antisense += 1
                    checked_reads += 1
                elif flag & 128:
                    sense += 1
                    checked_reads += 1
                else:
                    continue

    print 'sense number is {}'.format(sense)
    print 'antisense number is {}'.format(antisense)
    print 'total number of checked reads is {}'.format(checked_reads)

    # # Now we have min_reads processed, X as sense and Y as antisense.
    # # Check probability:
    # stranded_errors = min(sense, antisense)
    # current_prob = stats.binom.pmf(stranded_errors, checked_reads, error_rate)
    #
    # mean = stranded_errors/checked_reads
    # std = np.std(sense_count)
    # scaled_std = np.std(sense_count)
    #
    # # Calculate confidence interval for current sample.
    # conf_int = stats.norm.interval(0.95,loc=1,scale=std/math.sqrt(len(s)))
