# strandedness
For analyzing RNA-seq data for stranded or unstranded experimental protocol.

Current code allows user to run code with .sam file input, samples the file, parses the sampled reads, and calculated the binomial distribution function pmf for the resulting sense/antisense "errors."  i.e.: if the sample is stranded/unstranded what is the probability that we will have gotten X errors?
