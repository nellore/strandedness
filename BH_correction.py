#!/usr/bin/env python2

"""
Python 2.7 code for performing a Benjamini-Hochberg procedure to correct for
multiple-testing false discovery rate.

Required input:
    File containing SRA (or other experiment ID numbers) and their p-values

Optional input:
    alpha value - desired significance level for discovery of strandedness/
        rejection of null hypothesis (unstrandedness)
    Path to output directory to store list of experiments declared "stranded"
"""

import argparse
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Determine sample '
                                                 'strandedness.')
    parser.add_argument('--pvaluefile', '-p', required=True, help='File with '
                        'SRA accession numbers and their p-values; generated'
                        'by strandedness.py.')
    parser.add_argument('--alpha', '-a', type=float, default=0.05,
                        help='acceptable significance level of strandedness '
                        'call; must be between 0 and 1.')
    parser.add_argument('--outputpath', '-o', default='./', help='give path '
                        'for output files: junction reads and p-values.')

    args = parser.parse_args()
    alpha = args.alpha
    pval_file_path = args.pvaluefile
    out_path = args.outputpath

    with open(pval_file_path, 'r') as pval_file:
        pval_list = []
        for line in pval_file:
            if not line == '\n':
                line_list = line.split(',')
                line_list[1] = line_list[1].strip('\n')
                pval_list.append(line_list)

    # print 'original list is {}'.format(pval_list)
    pval_list.sort(key=lambda line: line[1], reverse=True)
    # print 'sorted list is {}'.format(pval_list)

    discovery = 0
    m = len(pval_list)
    for k_rev, expt in enumerate(pval_list, 1):
        k = m - k_rev + 1
        p_value = float(expt[1])
        # print 'checking line {}'.format(k_rev)
        # print 'k is {}'.format(k)
        # print 'check is {}'.format(alpha * k / m)
        # print 'p value is {}\n'.format(p_value)
        if p_value <= alpha * k / m:
            discovery = 1
            largest_k = k
            print 'largest k is {}'.format(k)
            break
    if discovery:
        stranded_expts = list(reversed(pval_list))[:k]
        print 'BH correction complete!  Stranded experiments are these:'
        print stranded_expts

        stranded_file_path = os.path.join(out_path, 'stranded_SRAs.txt')

        with open(stranded_file_path, 'w') as stranded_file:
            for set in stranded_expts:
                stranded_file.write('{},{}\n'.format(set[0], set[1]))
    else:
        print 'BH correction complete! All experiments were unstranded.'
