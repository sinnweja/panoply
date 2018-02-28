\name{panDrugGraph}
\alias{panDrugGraph}
\title{
  Drug Network Graph
}
\description{
Graph that shows the connection between the drugs, genes and pathways
for the drugs that target the most genes
}
\usage{
panDrugGraph(panDrug,usePathIDs = TRUE, ndrugs=8)
}
\arguments{
  \item{panDrug}{
    A data.frame of drug test results from panDrugSets
  }
  \item{usePathIDs}{
  If TRUE the shorter reactome pathway IDs will be used in the graph
  instead fo the complete name
  }
  \item{ndrugs}{
    Minimum number of Drugs to include in the graph
  }
 }
\details{
Graph where Drugs=Red Ellipse, Genes=Blue Ellipse, Pathways=Green Rectangles
}
\value{
  Nothing is returned
}
\examples{
data(clinPanTNBC)
data(genelistPan)
data(cnaPanTNBC)
data(gcPanTNBC)
data(gcinfoPan)
data(variantPanTNBC)
data(dgidbPan)
data(dgiSets)
data(reactome)
patient <- "TCGA-B6-A0IK"
match.index <- which(clinPanTNBC$Vital.Status=="Alive" & clinPanTNBC$PatientID!=patient & clinPanTNBC$age_at_diagnosis>50 & clinPanTNBC$N.stage == "N1" & clinPanTNBC$ndays > 1000)
ptmatch <- clinPanTNBC[match.index,"PatientID"]
drivGenes <- panGeneSets(caseid=patient, controlid=ptmatch, eventOnly=TRUE, variant=variantPanTNBC,
  cna=cnaPanTNBC, gcount=gcPanTNBC, tumorpct=0.3, tailEnd="upper", tailPct=0.1)
set.seed(1000)
drugResults <- panDrugSets(drivGenes, caseids=patient,
  controlids=ptmatch, gcount=gcPanTNBC,nsim=200, tailEnd="upper") 
drugResults<- drugResults[!is.na(drugResults$Pathways) & !is.na(drugResults$N.Cancer.Genes),]
panDrugGraph(drugResults)
}

