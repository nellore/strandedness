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
    Starting with pysam here...noted that this is not as good as using
        samtools as a subprocess, and/or other subprocess command line
        controls.

Improvements to be made:
    - Change to reservoir sampling instead of samtools?
        ->  Take fastq file input and run through an aligner.
            -> Take SRA accession numbers and download and run many in a row.

    - Use subprocess/samtools instead of Pysam

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
    sense = 0
    antisense = 0
    checked_reads = 0

    if confidence <= 0 or confidence >= 1:
        print "Please retry with desired confidence level between 0 and 1."
        exit()

    if probability <= 0 or confidence >= 1:
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

    num_errors = int(math.ceil(1000*error_rate))
    sample_array = np.ones(1000)
    sample_array[:num_errors] = 0
    std_dev = np.std(sample_array)
    z = stats.norm.ppf(1 - ((1-confidence)/2))
    min_reads = int(math.ceil((z*std_dev / ((1-probability)/2)) ** 2))
    # # OR alternately for now - set minimum number of reads.
    # min_reads = 50000

    # collect min_reads reads from sam file: samtools -s?  reservoir sample?
    wc_output = subprocess.check_output(['wc', '-l', samfile])

    for s in wc_output.split():
        if s.isdigit():
            sam_length = float(s)

    ratio = min_reads/sam_length
    add_seed = ratio + random.randint(0,1000)
    sample_seed = '{}'.format(add_seed)
    sampled_reads = subprocess.check_output(['samtools', 'view', '-s',
                                             sample_seed, samfile])

    # check length of sampled reads...
    print 'the total length of the sam file is {}'.format(sam_length)
    print 'the ratio should be {}'.format(ratio)
    print 'should be {} sampled reads'.format(min_reads)
    print 'number of sampled reads is {}'.format(len(sampled_reads.split('\n')))

    # # Check sample for being paired end or single end
    # paired_end = 


    # # Reads:
    # # Check paired end or single end? (position 2: 1 in the bitmask = paired or
    # # pysam is_paired)
    # # 	- Is it a primary alignment?  (position 2: not 256 in the bitmask or pysam
    # #       is_secondary)
    # # Does it have an extra XS:A string, and if so, + or -?
    # #   - If paired-end: is this the first or second pair?
    # #       (position 2: first = 64, second = 128 or pysam is_read1 or is_read2)
    # #   - If paired-end, check to make sure that the two pairs are on opposite
    # #       strands?  If read1, are 16 and 32 opposite?)
#
#
#     for read in sampled_reads:
#         if paired_end:
#             if pysam.is_secondary == 0:
#             if
#
#
#             checked_reads += 1
#         else:
#
#             check_sense = re.findall('XS:A:[+-]', line.upper())
#             if check_sense[0] == '+':
#                 sense += 1
#                 checked_reads += 1
#             elif check_sense[0] == '-':
#                 antisense += 1
#                 checked_reads += 1
#             else:
#                 continue
#
#
#     # Now we have min_reads processed, X as sense and Y as antisense.
#     # Check probability:
#     stranded_errors = min(sense, antisense)
#     current_prob = stats.binom.pmf(stranded_errors, checked_reads, error_rate)
#
#     mean = stranded_errors/checked_reads
#     std = np.std(sense_count)
#     scaled_std = np.std(sense_count)
#
#     # Calculate confidence interval for current sample.
#     conf_int = stats.norm.interval(0.95,loc=1,scale=std/math.sqrt(len(s)))


    # look up:
    # scipy binomial distribution
    # how to get confidence of the binomial distribution
    # gaussian distribution/central limit theorem/law of large numbers
    # how to get sigma from number of reads
    #