\name{clinPanBC}
\alias{clinPanBC}
\docType{data}
\title{
 Clinical Data for Breast Cancer Subjects
}
\description{
 Clinical Data of Breast Cancer patients to be used in Panoply
}
\usage{data("clinPanBC")}
\format{
  A data frame with 6 subjects.  One non-responder and 5 matched responders
  \describe{
    \item{\code{PatientID}}{character, patient identifier}
    \item{\code{Clinical.Molecular.Subtype}}{character, Clinical Molecular Subtype
      \code{Triple Negative}}
    \item{\code{Breast.Nodal.CR.all}}{character, Complete Response \code{yes} \code{no}}
    \item{\code{T.stage}}{character, T stage \code{T2}}
    \item{\code{N.stage}}{character, N stage \code{N0}}
    \item{\code{Percentage.tumor.cells.baseline.internal}}{integer, Percentage
      Tumor Cells}
    \item{\code{Age.Group}}{ordered factor, Age in 10 year categories
      \code{60-69} < \code{70+}}
 }
}
\details{

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(clinPanBC)
str(clinPanBC)
}
\keyword{datasets}
