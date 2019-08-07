\name{cnaPanCO}
\alias{cnaPanCO}
\docType{data}
\title{
DNA CNV (germline) and CNA (tumor) data for the TCGA Colon Cancer Subjects
}
\description{
Per-Gene DNA Copy Number Variation (germline) and Copy Number Alteration
(tumor) data for the TCGA Colon Cancer Subjects
}
\usage{data("cnaPanCO")}
\format{
  A data frame with 6 subjects.  One non-responder and 5 matched responders
  \describe{
    \item{\code{CHROM}}{chromosome of gene}
    \item{\code{START}}{gene start position}
    \item{\code{STOP}}{gene stop position}
    \item{\code{Gene.Symbol}}{gene symbol, NCBI}
    \item{\code{TCGA-DM-A0XD}}{patient identifier}
 }
}
\details{

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(cnaPanCO)
str(cnaPanCO)
}
\keyword{datasets}
