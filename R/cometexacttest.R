comet_exact_test <- function(tbl, pvalthresh = 1.1) {
	N <- as.integer(sum(tbl))
	k <- as.integer(log2(length(tbl)))
	stopifnot(is.numeric(tbl), 2**k == length(tbl), is.integer(k), is.integer(N), is.numeric(pvalthresh))
	.Call('cometexacttest', k, N, as.integer(tbl), pvalthresh)
}