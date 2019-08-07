\name{panGeneGraph}
\alias{panGeneGraph}
\title{
  Creates reactome network graph for Genes connected to top Drugs
}
\description{
   Uses Reactome data to create graph that shows the connection between
   Driver Genes and Connected Genes for the top Drugs in panDrug
}
\usage{
panGeneGraph(panGene, panDrug, minTargets=2, minPathways=2, ndrugs=4, ndrivers=20, ...)
}
\arguments{
  \item{panGene}{
    A data.frame of drug test results from panGeneSets
  }
  \item{panDrug}{
    A data.frame of drug test results from panDrugSets
  }
  \item{minTargets}{
    Retain drug results for drugs that target at least minTargets genes
    in the cancer and network genes.
  }
  \item{minPathways}{
    Filter drug results to those that target genes in at least
    minPathways, after accounting for minPathPct and minPathSize.
  }
  \item{ndrugs}{
    Number of Drugs to use to find Cancer Driver and Network genes in graph
  }
  \item{ndrivers}{
    Number of unique Cancer Drive genes in graph
  }
  \item{\dots}{
    Dynamic parameter for the values of additional parameters for the
    graph plot. 
  }
}
\details{
  A graph that shows the Cancer Driver Network for the top Drugs in pan Drug.  Circles=Druggable,
  Ellipse=Expressed Druggable Drivers, Red=Cancer Driver Genes, Blue=Network Genes
}
\value{
  nothing is returned
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
drugResults <- panDrugSets(drivGenes, caseids=patient, controlids=ptmatch, gcount=gcPanTNBC,nsim=200, tailEnd="upper") 

panGeneGraph(drivGenes,drugResults)
}