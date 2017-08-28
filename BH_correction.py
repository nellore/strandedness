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
from operator import itemgetter
import os


def get_sorted_pvals(file_path):
    '''Returns a list of SRA accession numbers and their p-values from file.

    Input the path to a file with SRA numbers and p-values.

    Returns a list of data from the file, sorted from lowest to highest
    p-value.
    '''
    with open(file_path, 'r') as pval_file:
        list = []
        for line in pval_file:
            if not line == '\n':
                line_list = line.split(',')
                line_list[1] = line_list[1].strip('\n')
                list.append(line_list)
    list.sort(key=itemgetter(1))
    return list


def BH_correction_procedure(pvals, alpha):
    '''Performs Benjamini-Hochberg multiple testing correction procedure.

     Input list of SRA numbers and their p-values, sorted from smallest to
     largest p-value, and the acceptable significance level for declaring a
     discovery (alpha).

     Returns whether or not a discovery can be found, and if discovery exists,
     the list number that represents the last (highest p-value) discovery.
     '''
    m = len(pvals)
    for k_rev, expt in enumerate(list(reversed(pvals))):
        k = m - k_rev
        p_value = float(expt[1])
        if p_value <= alpha * k / m:
            return True, k
        elif k == 0:
            return False, 0


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

    pval_list = get_sorted_pvals(pval_file_path)
    discovery, largest_k = BH_correction_procedure(pval_list, alpha)

    if discovery:
        stranded_expts = pval_list[:largest_k]
        print 'BH correction complete!  Largest k is {}'.format(largest_k)
        print 'Stranded experiments are these:'
        print stranded_expts

        stranded_file_path = os.path.join(out_path, 'stranded_SRAs.txt')
        with open(stranded_file_path, 'w') as stranded_file:
            for set in stranded_expts:
                stranded_file.write('{},{}\n'.format(set[0], set[1]))
    else:
        print 'BH correction complete! All experiments were unstranded.'
