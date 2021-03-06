\name{panDrugSets}
\alias{panDrugSets}
\alias{annotateDrugs}
\alias{panDrugStouffers}
\alias{panDrugConnect}
\alias{panallelPlot}
\alias{panallelPrepare}
\alias{findPathways}
\title{
  Tests on drugs for their effectiveness in targeting cancer genes in cases 
}
\description{
   Tests on drugs for their effectiveness in targeting activated cancer genes
   in cases against their matched controls.
}
\usage{
panDrugSets(panGene, caseids, controlids, gcount, minTargets=1, minPathPct=.05, minPathSize=8,
            minPathways=1, nsim=1000, tailEnd="both", gene.gs=NULL,
            gene.adj=NULL, drug.gs=NULL, drug.adj=NULL, 
            gageCompare = ifelse(length(caseids) > 1, "as.group", "unpaired"))
}
\arguments{
  \item{panGene}{
    A data.frame of drug test results from panGeneSets
  }
  \item{caseids}{
    identifiers of subjects case subjects that will match to variant
    column PatientID, and column name of cna and gcount.
  }
  \item{controlids}{
    identifiers of subjects in row names of variant and column names of gcount who are controls
  }
  \item{gcount}{
    Normalized gene expression counts for patients (columns) at gene
    symbols (rows). Gene symbols are the row names and should match
    dataset gcinfoPan gene symbols.
  }
  \item{minTargets}{
    Retain drug results for drugs that target at least minTargets genes
    in the cancer and network genes.
  }
  \item{minPathPct}{
   Of reactome pathways affected by gene networks targeted by the drug,
   only count those with total druggable by the drug in the pathway
   divided by the total pathway size >= minPathPct.
  }
  \item{minPathSize}{
    Of reactome pathways affected by gene networks targeted by the drug,
    only count the pathways that have total genes >= minPathSize.
  }
  \item{minPathways}{
    Filter drug results to those that target genes in at least
    minPathways, after accounting for minPathPct and minPathSize.
  }gene.gs=NULL, gene.adj=NULL, drug.gs=NULL, drug.adj=NULL,
  \item{drug.gs}{
    A list of genes, where the genes in each set are targeted by each
    drug in drug.adj.
  }
  \item{drug.adj}{
    Adjacency matrix (values of 1 or 0) for drug-gene interactions with the drug in the row
    directly impacting the gene in the column. The framework is equipped
    for coding for other interaction types of the drug and genes. 
  }
  \item{gene.gs}{
    List of genes per pathway/network to annotate the pathways targeted
    by the drug. Can be user-defined, but panoply contains reactome sets.
  }
  \item{gene.adj}{
    Adjacency matrix of genes connected via a graph of nodes and
    edges. Panoply contains reactome.adj, which is directed acyclic
    graph of Reactome nodes (genes) and edges (relationships). 
  }
  \item{tailEnd}{
    For Drug Network tests, perform differential expression tests that are one-sided ("upper" or
    "lower"), or two-sided ("both"). Also used with tailPct to include both,
    upper, or lower expression events to be used as cancer events.
  }
  \item{nsim}{
    Number of simulations to perform for evalulating significance of
    meta drug tests  
  }
  \item{gageCompare}{
    Character string for the R gage package to specify how to do
    within-network comparisons of the case(s) versus the controls across genes  
  }
}
\details{
  If tailEnd is upper, then testing for drugs that target over-expressed
  genes and gene networks in case(s) versus controls.
}
\value{
  A data.frame with the following columns
  \itemize{
  \item{Drug: Drug Name} 
  \item{N.Cancer.Genes: number of cancer genes that are targeted by the drug}  
  \item{Cancer.Genes: cancer genes that are targeted by the drug}
  \item{N.Network.Genes: number of network genes targeted by the drug} 
  \item{Network.Genes: network genes targeted by the drug}
  \item{N.Pathways: number of pathways with genes targeted by the drug
    (meeting minPathPct and minPathSize settings)} 
  \item{Pathways: pathway names of pathways with a gene targeted by the
    drug}
 %% remove ZTail.stat and Perm.pval
  \item{Network: Number of genes directly targeted by the drug}
  \item{DNT.pval: p-value comparing gene expression
    between caseids and controlids patients in the genes directly
    targeted by the drug}
  \item{DMT.Stat: Meta Z-statistic of combined gene networks for any gene
    targeted by the drug}
  \item{DMT.pval: p-value for DMT.Stat Z-statistic, based on multiplying
    s_i~N(0,1) by gene network test Z-statistics for i=1..nsim}  
  \item{Z.stat: Meta Z-statistic of combined gene networks for any gene
    targeted by the drug, weighted by number of genes targeted by the
    drug divided by gene network size}
  \item{PScore: A score used to sort the drug test results, which is a
    sum of the -log10(pval) of Network.pval and ZSim.pval}
 }
}
\references{
  %%  ~~ possibly secondary sources and usages ~~
  %% publication in progress
}
\examples{
data(clinPanTNBC)
data(genelistPan)
data(cnaPanTNBC)
data(gcPanTNBC)
data(gcinfoPan)
data(variantPanTNBC)
data(dgiSets)
data(reactome)
patient <- "TCGA-B6-A0IK"
match.index <- which(clinPanTNBC$Vital.Status=="Alive" & clinPanTNBC$PatientID!=patient & clinPanTNBC$age_at_diagnosis>50 & clinPanTNBC$N.stage == "N1" & clinPanTNBC$ndays > 1000)
ptmatch <- clinPanTNBC[match.index,"PatientID"]
drivGenes <- panGeneSets(caseid=patient, controlid=ptmatch, eventOnly=TRUE, variant=variantPanTNBC,
  cna=cnaPanTNBC, gcount=gcPanTNBC, tumorpct=0.3, tailEnd="upper", tailPct=0.1)
set.seed(1000)
drugResults <- panDrugSets(drivGenes, caseids=patient,
controlids=ptmatch, gcount=gcPanTNBC,nsim=200, tailEnd="upper",
drug.gs=dgi.gs, drug.adj=dgi.adj, gene.gs=reactome.gs, gene.adj=reactome.adj)

drugResults[1:10,!grepl("^Pathway", names(drugResults))]
}

