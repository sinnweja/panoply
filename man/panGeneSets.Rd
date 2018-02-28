\name{panGeneSets}
\alias{panGeneSets}
\alias{findDrivers}
\alias{filterNetworks}
\alias{drivDNA}
\alias{outRNA}
\alias{panConnect}
\alias{panGeneDruggable}
\title{
  Test gene networks for a set of cases versus a set of controls 
}
\description{
   Test differential gene expression in gene networks for a set of cases versus a set of controls, with
   gene networks defined by reactome directed graph for cancer genes,
   and genomic events determined from variant, cna, or gcount.
}
\usage{
panGeneSets(caseids, controlids, variant = NULL, cna = NULL,
           gcount, tumorpct = 0.5, tailPct = 0.1, tailEnd = "both",
           eventOnly = FALSE,  gene.adj = NULL, drug.adj = NULL,
           gageCompare = ifelse(length(caseids) > 1, "as.group", "unpaired"))
drivDNA(ids, variant = NULL, cna = NULL, tumorpct = 0.5, gene.adj=NULL)
outRNA(ids, gcount, tailPct = 0.1, tailEnd = "both")
}
\arguments{
  \item{caseids}{
    identifiers of subjects case subjects that will match to variant
    column PatientID, and column name of cna and gcount.
  }
  \item{controlids}{
    identifiers of subjects in row names of variant and column names of gcount who are controls
  }  
  \item{ids}{
    identifiers of subjects in row names of variant and column names of
    gcount and cna matrices
  }
  \item{variant}{
    data.frame for somatic and/or germline variants. Must contain
    columns CHROM, POS, Gene.Symbol, SampleType (Germline, Tumor) and
    PatientID that must match the caseid and controlid
  }
  \item{cna}{
    data.frame with log2 of copy number aberrations. Required columns
    CHROM, START, STOP, Gene.Symbol, and columns named to match caseid
    and controlid
  }
  \item{gcount}{
    Normalized gene expression counts for patients (columns) at gene
    symbols (rows). Gene symbols are the row names and should match
    dataset gcinfoPan gene symbols.
  }
  \item{tumorpct}{
    Vector of approximate tumor percentage for each patient samples given in
    the order they appear in cases and controls. If a single value is
    passed, it is used for all patients.  Copy number mutation events are
    called based on their approximate tumor percentage, transformed to
    the log2-ratio copy number call.
  }
  \item{tailPct}{
    percentile used to define the most extreme RNA expression events as
    cancer events per patient; only used when eventOnly=TRUE. Values 0 <
    tailPct < .5 allowed.
  }
  \item{tailEnd}{
    perform differential expression tests that are one-sided ("upper" or
    "lower"), or two-sided ("both"). Also used with tailPct to include both,
    upper, or lower expression events to be used as cancer events.
  }
  \item{gene.adj}{
    Adjacency matrix of genes connected via a graph of nodes and
    edges. Panoply contains reactome.adj, which is directed acyclic
    graph of Reactome nodes (genes) and edges (relationships). 
  }
  \item{drug.adj}{
    Adjacency matrix (values of 1 or 0) for drug-gene interactions with the drug in the row
    directly impacting the gene in the column. The framework is equipped
    for coding for other interaction types of the drug and genes. 
  } 
  \item{gageCompare}{
    Character string telling the gage package, which performs the gene
    network tests, to perform tests of 1-vs-M ("unpaired") or N-vs-M ("as.group")
  }
  \item{eventOnly}{
    Logical, test only gene networks that have a cancer event (variant,
    cna, and rna-expression outside tailPct)).
  }
}
\details{
  If eventOnly=FALSE, then test all gene networks for differential
  expression between cases and controls. Otherwise, cancer events are
  determined per gene from variant (any variant in the gene for the
  case(s)), cna (copy gain or loss for case(s) given tumorpct), and
  gcount (outlying expression from tailPct and tailEnd settings) and
  only those gene networks with a cancer event are tested. 
}
\value{
  A data.frame with the following columns
  \itemize{
  \item{Cancer.Gene: genes whose networks were tested}
  \item{Network: total connections of the cacer gene to other genes
    via reactome adjacency matrix}
  \item{Network.pval: p-value comparing gene expression
    between caseids and controlids patients}
  \item{Network.mean: Z-statistic of differential expression between
    caseids and controlids}
  \item{Druggable: Is the cancer gene druggable? 1=TRUE}
  \item{Total.Druggable: Number of the connected network genes that are druggable}
  \item{NetGenes: Network genes that are over- or under-expressed}
  \item{NetDrugs: Drug names that would target the genes in Total.Druggable}
  \item{Drugs: Drug names that would target the cancer gene}
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
data(dgidbPan)
data(dgiSets)
data(reactome)

patient <- "TCGA-B6-A0IK"
match.index <- which(clinPanTNBC$Vital.Status=="Alive" & clinPanTNBC$PatientID!=patient & clinPanTNBC$age_at_diagnosis>50 & clinPanTNBC$N.stage == "N1" & clinPanTNBC$ndays > 1000)
ptmatch <- clinPanTNBC[match.index,"PatientID"]
drivGenes <- panGeneSets(caseid=patient, controlid=ptmatch, eventOnly=TRUE, variant=variantPanTNBC,
   cna=cnaPanTNBC, gcount=gcPanTNBC, tumorpct=0.3, tailEnd="upper", tailPct=0.1)
drivGenes[1:10,1:6]
}

