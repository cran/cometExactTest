\name{comet_exact_test}
\alias{comet_exact_test}
\title{CoMEt exact test for 2^k contingency tables}
\description{
  This function performs the CoMEt exact test for a given 2^k contingency table.
}
\usage{
comet_exact_test(tbl, pvalthresh=1.1)
}
\arguments{
  \item{tbl}{is a vector of 2^k non-negative integers (where k is the
number of dimensions in your contingency table).}

  \item{pvalthresh}{is a float specifying the threshold at which the
cometexacttest will halt execution if the P-value being computed is ever
more than the threshold.}

}

\examples{
  comet_exact_test(c(33, 10, 10, 1, 10, 0, 0, 1)) # 2^3 test => 0.02303503
}
