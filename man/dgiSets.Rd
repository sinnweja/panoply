\name{dgiSets}
\alias{dgiSets}
\docType{data}
\title{
Druggable Genome Interaction Database for Cancer
}
\description{
  Druggable Genome Interaction Database annotation for Cancer Drugs and
  their targets
}
\usage{data("dgiSets")}
\format{
  The adjacency matrix \code{dgi.adj} has drugs (rows) and a 1 marked
   where the drug targets a gene (columns). The same information is
   stored as a list of gene names for each drug in \code{dgi.gs}. The
   data.frame \code{dgiDrug} has each drug in a row, with genes they
   target, and sources for which drug database contained it in these
   columns:
 %% A data frame with 6 subjects.  One non-responder and 5 matched responders
 \describe{ 
   
    \item{\code{Drug}}{Cancer drugs}
    \item{\code{GeneID}}{Gene Symbols of genes targeted by the drug}
    \item{\code{Synonym}}{Drug synonyms}
    \item{\code{DGIdb}}{Logical, curated by DGI-db?}
    \item{\code{DrugBankCurate}}{Logical, curated by Drug Bank list?}
  }
 }

\details{
  
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\keyword{datasets}
