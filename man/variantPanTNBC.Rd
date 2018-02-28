\name{variantPanTNBC}
\alias{variantPanTNBC}
\docType{data}
\title{
Dataset with somatic and germline DNA variants for TCGA patients
}
\description{
Dataset with somatic and germline DNA variants for TCGA patients
patients
}
\usage{data("variantPanTNBC")}
\format{
  A data frame with the following columns:  
  \describe{
    \item{\code{Gene.Symbol}}{Gene Symbol}
    \item{\code{CHROM}}{Chromosome}
    \item{\code{POS}}{Start Position}
    \item{\code{Stop}}{Stop Position}
    \item{\code{Ref}}{Reference Allele}
    \item{\code{Alt}}{Alternate Allele}
    \item{\code{Effect}}{Effect}
    \item{\code{PatientID}}{Patient ID}
    \item{\code{V1B.refDP}}{Germline Reference allele depth}
    \item{\code{V1B.altDP}}{Germline Alternate allele depth}
    \item{\code{V1T.refDP}}{Tumor Reference allele depth}
    \item{\code{V1T.altDP}}{Tumor Alternate allele depth}
    \item{\code{AD}}{AD}
    \item{\code{GT}}{Genotype}
    \item{\code{SampleType}}{Sample Type \code{Blood} \code{Tissue}}
 }
}
\details{
Non-synonymous somatic variants, and high-impact germline variants.
}

\examples{
data(variantPanTNBC)
str(variantPanTNBC)
}
\keyword{datasets}
