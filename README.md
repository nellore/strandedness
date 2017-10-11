# strandedness
For analyzing RNA-seq data for stranded or unstranded experimental protocol.

Current code allows user to run code with .csv file input; file is parsed for NIH Sequence Read Archive (SRA) accession numbers and paired- or single-end experiment protocol information.  Each SRA is then sampled: a fastq-dump batch download of sampled reads is performed, and alignment is performed on the entire set of reads.  Each read is tested for "usefulness" (i.e.: is it a junction read?).  At the end of each SRA, a binomial test/two-ended cumulative distribution function (CDF) is calculated to give a p value for the null hypothesis "experiment is unstranded" given the obtained number of sense/antisense reads.  Benjamini-Hochberg (BH) correction for multiple testing false discovery rate is performed either in the main script, or in the separate BH_correction.py script.

# Reproducing strandedness classifications of human SRA RNA-seq samples, each with >= 10M reads

1. Run
```
python X.py --arg1 /path/to/X
```
where `/path/to/X` is ...
