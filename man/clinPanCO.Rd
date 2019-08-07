\name{clinPanCO}
\alias{clinPanCO}
\docType{data}
\title{
 Clinical Data for TCGA Colon Cancer Patients
}
\description{
 Clinical Data of TCGA Colon Cancer patients to be used in Panoply
}
\usage{data("clinPanCO")}
\format{
  A data frame with 9 subjects.  5 cases (deceased) and 4 matched
  controls (alive)
  \describe{
    \item{\code{bcr_patient_barcode}}{character, patient identifier,
      TCGA barcode}
    \item{\code{bcr_patient_uuid}}{character, TCGA uuid}
    \item{\code{patient_id}}{character, simplified patient identifier}
    \item{\code{days_to_death}}{character, days to death}
    \item{\code{days_to_last_followup}}{character, days to last followup}
    \item{\code{race}}{character, race, \code{WHITE}}
    \item{\code{vital_status}}{character, vital status \code{Alive} \code{Dead}}
    \item{\code{gender}}{character, gender \code{FEMALE} \code{MALE}}
    \item{\code{kras_mutation_found}}{character, kras mutation present\code{YES}}
    \item{\code{surv}}{character, \code{alive} \code{recur}}
    \item{\code{CN.tumorpct}}{numeric, Tumor Percentage}
 }
}
\details{

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}

\examples{
data(clinPanCO)
str(clinPanCO)
}
\keyword{datasets}
