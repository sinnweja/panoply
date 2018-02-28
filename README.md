# PANOPLY: *P*recision C*a*ncer Ge*no*mic Re*p*ort: Single Samp*l*e Inventor*y*

Panoply is an R package that is part of a workflow to create reports to match cancer drugs to a patient's tumor using genomic data. The method can be applied to multiple cancer types.

# Datasets

We have downloaded publically available data from The Cancer Genome Atlas (TCGA) for Colon Cancer (COAD) and Triple Negative Breast Cancer (TNBC). We curated a set of poor responders to be matched to  set of long-term responders to the same drug, for patients matched on clinical features.

We have also provided drug and gene annotation datasets that represent key cancer driver genes and the genes that are linked to them via public annotation sources (i.e. Reactome and Tthe Druggable Genome Atlas Database, or DGIdb)

# Functions 

The primary functions unique to panoply are *panGeneSets* and *panDrugSets*.  The *panGeneSets* function tests for differentially-expressed cancer networks in a set of cases (non-responders) versus controls (responders).  The *panDrugSets* function provides two tests: DNT is the Drug Network Test, which is a test of a differentially-expressed sets of genes that are directly linked to the drug via drug-gene annotation; DMT is the Drug Meta Test, which is a combined, or meta, test of all of the gene sets with any gene targeted by a the drug, weighted by the inverse of the size of the gene network. 
 
# Vignettes

We provide example vignettes for each of the COAD and TNBC cancer datasets, which are in Sweave format.  We also provide a vignette for how we use multiple drug-gene annotation sources to build the drug-gene datasets within the vignette titled *buildDrugSets*.