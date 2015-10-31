\name{comet_exact_test}
\alias{comet_exact_test}
\title{CoMEt exact test for 2^k contingency tables}
\description{
  This function performs the CoMEt exact test for a given 2^k contingency table.
}
\usage{
comet_exact_test(tbl, pvalthresh=1.1, mutmatplot=T)
}
\arguments{
  \item{tbl}{is a vector of 2^k non-negative integers (where k is the
number of dimensions in your contingency table). We represent each cell of the
contingency table as a binary string, and the order of values in the table
follows the Least Significant Bit for the binary representation of a given cell.
The binary representation for the 2x2 contingency table X for a pair of genes is
x00, x01, x10, x11, where x00 is the cell where neither gene is mutated, x01 is the
cell where the first gene is mutated and the second is not, etc. Take the
following 2x2x2 contingency table as an additional example, where the
ordering of the vector is x000, x001, x010, x011, x100, x101, x110, x111. 
  \if{html}{\figure{tablex.png}{options: width="35\%" alt="Figure: tablex.png"}}
  \if{latex}{\figure{tablex.pdf}{options: width=10cm}}
  }

  \item{pvalthresh}{is a float specifying the threshold at which the
cometexacttest will halt execution if the P-value being computed is ever
more than the threshold.}

  \item{mutmatplot}{is a boolean specifying to create mutation matrix of the
  input table or not. If \code{mutmatplot==T}, a file
  "cometExactTest_plot.pdf" will be created under your working directory.}

}

\examples{
  comet_exact_test(c(33, 10, 10, 1, 10, 0, 0, 1))  # 2^3 test => 0.02303503
}
