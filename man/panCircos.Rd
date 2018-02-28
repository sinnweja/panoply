\name{panCircos}
\alias{panCircos}
\title{
  plot method for panPath/panDruggable objects
}
\description{
  Create Circos Plot of copy number, RNA expression events, germline and
  somatic variants. By default, plots DNA mutation (variant/cna) and RNA expression
  events and shows connections between events in top pathways with a driver gene
  event. Allows either germline/somatic variants, cna, or all three.
}
\usage{
panCircos(panGene, panDrug, caseids, variant = NULL, 
    cna = NULL, gcount, gcinfo, tumorpct = 0.5, 
    tailPct = 0.05, tailEnd = "upper", minTargets = 1, minPathPct = 0.05, 
    minPathSize = 8, minPathways = 1, ...) 
}
\arguments{
\item{panGene}{
    A data.frame of drug test results from panGeneSets
}
\item{panDrug}{
    A data.frame of results from panDrugSets
}
\item{caseids}{
  identifiers of subjects that will match to variant
    column PatientID, and column name of cna and gcount.
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
}
}
\value{
  Nothing is returned
}
\details{
 circos plot of the driver gene events in outer bands, and gene
 expression events in the inner bands, colored blue if they are over- or
 under-expressed, and if one of the key drugs target them. Across the
 center circle, we show genes that are connected as the top cancer
 drivers connected to their druggable in-network genes. 
}
% docclass is function
% Converted by Sd2Rd version 43267.
