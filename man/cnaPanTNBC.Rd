\name{cnaPanTNBC}
\alias{cnaPanTNBC}
\docType{data}
\title{
DNA CNV (germline) and CNA (tumor) data for TCGA samples
}
\description{
Per-Gene DNA Copy Number Variation (germline) and Copy Number Alteration
(tumor) data for TCGA samples
}
\usage{data("cnaPanTNBC")}
\format{
  A data frame with 77 subjects.  
  \describe{
    \item{\code{CHROM}}{chromosome of gene}
    \item{\code{START}}{gene start position}
    \item{\code{STOP}}{gene stop position}
    \item{\code{geneid}}{Gene ID}
    \item{\code{Gene.Symbol}}{gene symbol, NCBI}
    \item{\code{TCGA-A2-A04T}}{patient identifier}
 }
}
\details{
  %% more to come
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\examples{
data(cnaPanBC)
str(cnaPanBC)
}
\keyword{datasets}
