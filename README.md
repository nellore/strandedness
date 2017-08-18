# strandedness
For analyzing RNA-seq data for stranded or unstranded experimental protocol.

Current code allows user to run code with .csv file input; file is parsed for SRA accession numbers and paired- or single-end experiment protocol information.  Each SRA is then sampled: one read at a time is downloaded, aligned, and tested for usefulness.  At the end of each SRA, the binomial distribution function pmf is calculated for the obtained number of sense/antisense "errors."  i.e.: if the sample is stranded/unstranded what is the probability that we will have gotten X errors?
