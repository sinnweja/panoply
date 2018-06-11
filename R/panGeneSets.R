
## Purpose: connetivity tests for driver genes for cases versus controls and carriers vs non-carriers
## Author: Jason Sinnwell, Rani Kalari, Kevin Thompson
## Created: 8/2015
## Updated: 6/2016 add gene-expr outliers to drivers. Add panDrugConnect to test drug networks rather than just reactant adjacency

## caseid: subject ids of cases that are rownames of mutMat and outMat. same length
## controlid: subject ids of controls
## variant: CHROM, POS, Gene.Symbol, PatientID, SampleType
## cna: CHROM, START, STOP, Gene.Symbol, Cytoband, cn-for patients
## gcount: rows=Gene.Symbol, cols=genecounts for patients

## mutMat: mutation matrix; rows=subjects, columns=genes, colnames are genes
## outMat: rnaseq outlier matrix  rows=subjects, columns=genes, colnames are genes
## adjMat: adjacency matrix: rows and columns are genes, with gene names as dimnames

panGeneSets <- function(caseids, controlids, variant=NULL, cna=NULL, gcount, tumorpct=0.5, tailPct=0.1, tailEnd="both", eventOnly=FALSE, gene.adj=NULL, drug.adj=NULL, gageCompare=ifelse(length(caseids)>1,"as.group","unpaired")) {

  if(is.null(gene.adj)) {
    data(reactome)
    gene.adj <- reactome.adj
  }
  if(is.null(drug.adj)) {
    data(dgiSets)
    drug.adj <- dgi.adj
  }
  ## call panConnect, then panGeneDruggable
  pconn <- panConnect(caseids, controlids, variant=variant, cna=cna, gcount=gcount, tumorpct=tumorpct, tailPct=tailPct, tailEnd=tailEnd, eventOnly=eventOnly, gene.adj=gene.adj, gageCompare=gageCompare)
  pgene <- panGeneDruggable(pconn, gcount=gcount, caseids=caseids, tailPct=tailPct, tailEnd=tailEnd, gene.adj=gene.adj, drug.adj=drug.adj)
  return(pgene)
}
  
drivDNA <- function(ids, variant=NULL, cna=NULL, tumorpct=0.5, gene.adj=NULL) {
  
  if(length(tumorpct) != length(ids)) {
  ##  warning("using first element of tumorpct for all subjects")
    tumorpct <- rep(tumorpct[1], length(ids))
  }
  if(is.null(gene.adj)) {
    # loads reactome.adj
    data(reactome)
    gene.adj <- reactome.adj
  }
  drivMat <- matrix(0, nrow=length(ids), ncol=nrow(gene.adj)) #
  dimnames(drivMat) <- list(ids,rownames(gene.adj))

  if(!is.null(cna)) {
    cnamat <- t(cna[,ids,drop=FALSE])
    colnames(cnamat) <- cna[,'Gene.Symbol']  
    addcnagenes <- colnames(drivMat)[!(colnames(drivMat) %in% colnames(cnamat))]
    if(length(addcnagenes)>0) {
       oldcnancol <- ncol(cnamat)
       cnamat <- cbind(cnamat, matrix(0, nrow=nrow(cnamat), ncol=length(addcnagenes)))
       colnames(cnamat) <- c(colnames(cnamat)[1:oldcnancol], addcnagenes)
       cnamat <- cnamat[,order(colnames(cnamat))]
    }
    ## we use a flexible cutoff for cna events based on the tumor percentage
    ## estimated from ordered cna plots, or alt allele frequency distributions
    ##pct <- 0.4
    ##log2((2*(1-pct)+4*pct)/2)  ## gain 2
    ##log2((2*(1-pct)+3*pct)/2) ## gain 1
    ##log2((2*(1-pct)+2*pct)/2)  # neutral
    ##log2((2*(1-pct)+1*pct)/2) ## loss 1
    ##log2((2*(1-pct)+0*pct)/2)  ## loss 2 
    
    cnamat <- cnamat[ids, colnames(drivMat)] #colnames(cnamat) %in% colnames(drivMat)]
    cnacut1 <- matrix(log2((2*(1-tumorpct)+1*tumorpct)/2),
                      nrow=nrow(cnamat), ncol=ncol(cnamat), byrow=FALSE)
    cnacut3 <- matrix(log2((2*(1-tumorpct)+3*tumorpct)/2),
                      nrow=nrow(cnamat), ncol=ncol(cnamat), byrow=FALSE)
    ## or quick and dirty would be:
    ##  cnacut1 <- 0.4
    ##  cnacut3 <- 1.5
  
    drivMat <- drivMat + 1*(cnamat < cnacut1) + 1*(cnamat > cnacut3)
  }
  
  if(!is.null(variant)) {
    variantscores <- data.frame(CHROM=tapply(variant$CHROM, variant$Gene.Symbol, function(x) x[1]),
                 GENE=tapply(variant$Gene.Symbol, variant$Gene.Symbol, function(x) x[1]))
  
    for(subj in ids) { #   #unique(variant$Subject)) {
      variantscores[,subj] <- 0
      vscore <- tapply(variant[variant$PatientID==subj,"GT"],
               variant[variant$PatientID==subj, "Gene.Symbol"],
               FUN=function(x) { alleles <- unlist(strsplit(c(x), split="/")); sum(alleles!="0")})
      ## count of non-zero (non-REF) alleles per gene    
      variantscores[match(names(vscore), variantscores$GENE),subj] <- unlist(vscore)
    }
    
    variantmat <- t(variantscores[,3:ncol(variantscores)])
    colnames(variantmat) <- variantscores$GENE
    variantmat <- variantmat[order(rownames(variantmat)),order(colnames(variantmat))]
    addvariantcols <- colnames(drivMat)[!(colnames(drivMat) %in% colnames(variantmat))]
    oldvariantncol <- ncol(variantmat)
    variantmat <- cbind(variantmat, matrix(0, nrow=nrow(variantmat), ncol=length(addvariantcols)))
    colnames(variantmat) <- c(colnames(variantmat)[1:oldvariantncol], addvariantcols)
    variantmat <- variantmat[,colnames(variantmat) %in% colnames(drivMat)]
    variantmat <- variantmat[ids,order(colnames(variantmat))]
    ## check: table(colnames(variantmat) == colnames(drivMat))
    drivMat <- drivMat + variantmat[,colnames(drivMat)]    
  }
 
  ## drivMat <- variantmat + 1*(cnamat < cnacut1) + 1*(cnamat > cnacut3)
  ## feature to add: require dna-repair genes to need 2 "hits", oncogenes 1
  ## data(genelistPan)
  if(0) {   ## disable 2-hit coding for TS, 1-hit for OG
    codeVec <- with(genelistPan[genelistPan$AllLists==TRUE,],ifelse(grepl("OG",HC.class),1,2))
    names(codeVec) <- genelistPan[genelistPan$AllLists==TRUE, "GeneSymbol"]
    codeMat <- t(matrix(codeVec[match( colnames(drivMat), names(codeVec))], nrow=ncol(drivMat), ncol=nrow(drivMat)))
    colnames(codeMat) <- names(codeVec)[match( colnames(drivMat), names(codeVec))]
  }
  drivMat <- ifelse(drivMat >= 1, 1, 0)    ## if enable 2-hit coding, drivMat >= codeMat, 1, 0
  ## simpler:  drivMat <- ifelse(drivMat >= 1, 1, 0)
  ## complex:  drivMat <- ifelse(drivMat >= codeMat, 1, 0)
  return(drivMat)
}

outRNA <- function(ids, gcount, tailPct=0.1, tailEnd="both") {
  tailEnd <- casefold(tailEnd, upper=FALSE)
  if(!(tailEnd %in% c("upper","lower","both"))) {
    warning("tailEnd argument should be upper, lower, or both, set to upper\n")
  }
  if(tailPct > .5 | tailPct <= 0) {
    warning("tailPct argument should be < 0.5 and >0\n")
  }
  gcount <- gcount[order(rownames(gcount)),]
  if(tailEnd=="upper") {
    cutoff.subj <- apply(gcount, 2, function(x) { quantile(x, prob=1-tailPct)})
   
    outMat <- t(1*(gcount >= matrix(cutoff.subj, nrow=nrow(gcount), ncol=ncol(gcount), byrow=TRUE)))
  }
  if(tailEnd=="lower") {
    cutoff.subj <- apply(gcount, 2, function(x) { quantile(x, prob=tailPct)})
    outMat <- t(1*(gcount <= matrix(cutoff.subj, nrow=nrow(gcount), ncol=ncol(gcount), byrow=TRUE)))
  }
  if(tailEnd=="both") {
    cutoff.lower <- apply(gcount, 2, function(x) { quantile(x, prob=tailPct/2)})
    cutoff.upper <- apply(gcount, 2, function(x) { quantile(x, prob=1-tailPct/2)})
    outMat <- t(1*(gcount <= matrix(cutoff.lower, nrow=nrow(gcount), ncol=ncol(gcount), byrow=TRUE) |
                   gcount >= matrix(cutoff.upper, nrow=nrow(gcount), ncol=ncol(gcount), byrow=TRUE)))
  }
 
  #table(colSums(outMat))
  outMat <- outMat[ids,order(colnames(outMat)),drop=FALSE]
  return(outMat)
}

  
## result: data.frame with gene, carrierlist, ncarriers, node.degree,
##         prop.pval.carrier, prop.pval.case,    :: prop.test p-value of events/node.degree in carrier/case
##         caseMut.cov.events, caseWT.cov.events,
##         controlMut.cov.events, controlWT.cov.events,
##         ttest.pval.carrier,  :: two-sided t-test of connectivity of MUTation carriers vs non-carriers(WT)

panConnect <- function(caseids, controlids, variant=NULL, cna=NULL, gcount, tumorpct=0.5, tailPct=0.1, tailEnd="both", eventOnly=FALSE, gene.adj=NULL, gageCompare="unpaired") {

  if(is.null(gene.adj)) {
    data(reactome)
    gene.adj <- reactome.adj
  }
  
  tailEnd <- casefold(tailEnd)
  if(!(tailEnd %in% c("upper","lower","both"))) {
    warning("invalid tailEnd, setting to 'both'.")
  }
  if(length(tumorpct) != length(c(caseids,controlids))) {
    ##  warning("using first element of tumorpct for all subjects")
    tumorpct <- rep(tumorpct[1], length(caseids) + length(controlids))
  }
  outMat <- outRNA(ids=c(caseids, controlids), gcount, tailPct=tailPct, tailEnd=tailEnd)
  drivMat <- drivDNA(ids=c(caseids, controlids), variant, cna, gene.adj=gene.adj, tumorpct=tumorpct)
  if(eventOnly) {
    drivRNA <- colnames(drivMat)[colnames(drivMat) %in% colnames(outMat)]
    drivMat[1,drivRNA] <- drivMat[1,drivRNA] + outMat[1,drivRNA]
    drivMat <- drivMat[,drivMat[1,]>0]   
  }
  
#  gene.adj <- gene.adj[!duplicated(gene.adj),]
#  data(reactome)
  drivers <- data.frame(Cancer.Gene=colnames(drivMat), Network=0,
     carrierlist=apply(drivMat,2, function(x) paste(names(x)[x>0 & names(x) %in% c(caseids, controlids)], collapse=",")), 
     Network.pval=1, Network.mean=0, pval.WT=1, pval.MUT=1, 
     caseMUT.events="",caseWT.events="",
     controlMUT.events="", controlWT.events="")
  warnset <- options()$warn
  options(warn = -1) 
  for(k in 1:nrow(drivers)) {
    gconnect <- colnames(gene.adj)[which(gene.adj[drivers[k,1],]>0)]
    gconnect <- gconnect[gconnect %in% rownames(gcount)]
    
    kcase <- which(colnames(gcount) %in% caseids)
    kcontrol <- which(colnames(gcount) %in% controlids)
#browser()                   
    gage.direct <- gage(gcount[,c(kcontrol,kcase)], gsets = list(gconnect),
                        ref = 1:length(kcontrol), same.dir=TRUE,set.size=c(1,500),
                        samp = (1+length(kcontrol)):(length(kcontrol)+length(kcase)),
                        saaTest=gs.zTest, compare =  gageCompare) #"unpaired")
    gage.both <- gage(gcount[,c(kcontrol,kcase)], gsets = list(gconnect),
                        ref = 1:length(kcontrol), same.dir=FALSE,set.size=c(1,500),
                        samp = (1+length(kcontrol)):(length(kcontrol)+length(kcase)),
                      saaTest=gs.zTest, compare = gageCompare) #as.group or "unpaired")
    if(tailEnd=="upper") {    
      drivers$Network.pval[k] <- gage.direct$greater[1,'p.val']
      drivers$Network.mean[k] <- gage.direct$greater[1,'stat.mean']
    }
    if(tailEnd=="lower") {
      drivers$Network.pval[k] <- gage.direct$less[1,'p.val']
      drivers$Network.mean[k] <- gage.direct$less[1,'stat.mean']
    }
    if(tailEnd=="both") {
      drivers$Network.pval[k] <- gage.both$greater[1,'p.val']
      drivers$Network.mean[k] <- gage.both$greater[1,'stat.mean']
    }

    #drivers$Network.pval[k] <- gage.direct$greater[1,'p.val']
    #drivers$Network.pval.both[k] <- gage.both$greater[1,'p.val']
    ## drivers$pval.less[k] <- gage.direct$greater[1,'p.val']
    ## gagedf <- data.frame(Cancer.Gene=rownames(gage.direct$greater), p.gage=gage.direct$greater[,'p.val'])
   
    carriers <- strsplit(drivers$carrierlist[k], split=",")[[1]]
    caseMut.ce <- rowSums(outMat[rownames(outMat) %in% carriers & rownames(outMat) %in% caseids,colnames(outMat) %in% gconnect,drop=FALSE])
    caseWT.ce <- rowSums(outMat[!(rownames(outMat) %in% carriers) & rownames(outMat) %in% caseids,colnames(outMat) %in% gconnect,drop=FALSE])
    controlMut.ce <- rowSums(outMat[rownames(outMat) %in% carriers & rownames(outMat) %in% controlids,colnames(outMat) %in% gconnect,drop=FALSE])
    controlWT.ce <- rowSums(outMat[!(rownames(outMat) %in% carriers) & rownames(outMat) %in% controlids,colnames(outMat) %in% gconnect,drop=FALSE])
    ##  drivers$mut.events[k] <- sum(c(caseMut.ce, controlMut.ce))
    drivers$Network[k] <- length(gconnect)
    ##  if( length(c(caseMut.ce,controlMut.ce)) > 1 & length(c(caseWT.ce,controlWT.ce)) > 1) {
    
    if(drivers$Network[k] > 0 & sum(caseMut.ce) > 0 & sum(controlMut.ce) > 0) {
      drivers$pval.MUT[k] <- prop.test(c(sum(caseMut.ce),sum(controlMut.ce)),c(length(caseMut.ce), length(controlMut.ce))*as.numeric(drivers$Network[k]),
                                       alternative="greater")$p.value
    }
    if(drivers$Network[k] > 0 & sum(caseMut.ce > 0) &  sum(controlWT.ce)>0) {
      drivers$pval.WT[k] <- prop.test(c(sum(caseMut.ce),sum(controlWT.ce)),c(length(caseMut.ce), length(controlWT.ce))*as.numeric(drivers$Network[k]),
                                      alternative="greater")$p.value
    }
  
    drivers$caseMUT.events[k] <- paste(caseMut.ce, collapse=",")
    drivers$caseWT.events[k] <- paste(caseWT.ce, collapse=",")
    drivers$controlMUT.events[k] <- paste(controlMut.ce, collapse=",")
    drivers$controlWT.events[k] <- paste(controlWT.ce, collapse=",")   
  } ## end for()
  options(warn=warnset) 
  drivers <- drivers[order(drivers$Network.pval,decreasing=FALSE),-c(grep("WT",colnames(drivers)),grep("MUT",colnames(drivers)),grep("both$",colnames(drivers)),grep("carrierlist",colnames(drivers)))]
  return(drivers)
}

## formerly known as panDruggable
panGeneDruggable <- function(pconn, gcount, caseids, tailPct=0.1, tailEnd="both", gene.adj=NULL, drug.adj=NULL) {
  ## pconn: result from panPath
  ## gcount: data.frame of genecounts
  ## caseids: case ids to subset network genes to those that are outliers
#  data(dgidbPan)
  if(is.null(gene.adj)) {
    data(reactome)
    gene.adj <- reactome.adj
  }
  if(is.null(drug.adj)) {
    data(dgiSets)
    drug.adj <- dgi.adj
  }
  pdrug <- pconn
  pdrug$Druggable <- ifelse(pdrug$Cancer.Gene %in% colnames(dgi.adj), 1, 0) ##dgidbPan$Gene, 1, 0)
#  pdrug$ConnDruggable <- sapply(strsplit(pdrug$Case.PathGenes,split=","), function(x) sum(x %in% dgidbPan$Gene))
  pdrug$Network.Druggable <- pdrug$Total.Druggable <- 0
#  pdrug$TotDruggable <- pdrug$DrivDruggable + pdrug$ConnDruggable
#  pdrug <- merge(pdrug, dgidbPan[,c("Gene","Drugs")], by.x="Cancer.Gene",by.y="Gene",all.y=FALSE, all.x=TRUE)
#  pdrug <-pdrug[nchar(pdrug$caseMUT.event)>0,]
  pdrug$Drugs <- pdrug$NetDrugs <- pdrug$NetGenes <- ""
  outMat <- outRNA(ids=caseids, gcount, tailPct=tailPct, tailEnd=tailEnd)
  for(k in 1:nrow(pdrug)) {
    pdrug[k,"Drugs"] <- paste(rownames(dgi.adj)[rowSums(dgi.adj[,which(colnames(dgi.adj) %in% pdrug[k,"Cancer.Gene"]),drop=FALSE])>0],collapse="/")
  
    ## get drugs for connected genes
    conngenes <- c(pdrug$Cancer.Gene[k],colnames(gene.adj)[gene.adj[pdrug$Cancer.Gene[k],]>0])
    conngenes <- conngenes[conngenes %in% colnames(outMat)[colSums(outMat)>0]]
    pdrug$NetGenes[k] <- paste(conngenes,collapse=",")
    conndrugs <- paste(rownames(dgi.adj)[rowSums(dgi.adj[,which(colnames(dgi.adj) %in% conngenes),drop=FALSE])>0],collapse="/")
          ##paste(dgidbPan$Drugs[dgidbPan$Gene %in% conngenes],sep="/")    
    pdrug$Network.Druggable[k] <- sum(colnames(dgi.adj)[colSums(dgi.adj)>0] %in% conngenes)
      ##  sum(dgidbPan$Gene %in% conngenes)
      ##  pdrug$ConnGenes[k] <- ifelse(!is.null(conngenes),conngenes,"")
    pdrug$NetDrugs[k] <- ifelse(!is.null(conndrugs), conndrugs, "")
    pdrug$Total.Druggable[k] <- pdrug$Druggable[k] + pdrug$Network.Druggable[k]
  }
  return(pdrug[with(pdrug, order(pmin(Network.pval))),-grep("Network.Druggable",names(pdrug))])
                          #rowMeans(cbind(log10(pval.WT),log10(pval.MUT),log10(pval.gage))),decreasing=FALSE)),])
}

