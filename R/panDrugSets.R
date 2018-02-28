
## Purpose: connectivity tests for driver genes for cases versus controls and carriers vs non-carriers
## Author: Jason Sinnwell, Rani Kalari, Kevin Thompson
## Created: 8/2015
## Last Updated: 2/6/2018
panDrugSets <- function(panGene, caseids, controlids, gcount, minTargets=1, minPathPct=.05, minPathSize=8,
             minPathways=1, nsim=1000, tailEnd="both", gene.gs=NULL, gene.adj=NULL, drug.gs=NULL, drug.adj=NULL, 
             gageCompare = ifelse(length(caseids) > 1, "as.group", "unpaired")) {

  ## if drug and gene network data/sets not given, use data from package
   if(is.null(drug.gs) | is.null(drug.adj)) {
     warning("drug.gs or drug.adj not found, loading dgidb data...")
     data(dgiSets)
     drug.gs <- dgi.gs
     drug.adj <- dgi.adj
   }
  
   if(is.null(gene.gs) | is.null(gene.adj)) {
     warning("gene.gs or gene.adj not found, loading reactome data...")
     data(reactome)
     gene.adj <- reactome.adj
     gene.gs <- reactome.gs
   }

   ## make sure drug sources agree
   drug.adj <- drug.adj[rownames(drug.adj) %in% names(drug.gs),]
   drug.gs <- drug.gs[rownames(drug.adj)]
   uDrugs <- names(drug.gs)
   
   ## older function, still calculates ZTail.pval, and annotates pathways
   drugScores <- panDrugStouffers(panGene, nsim=nsim, tailEnd=tailEnd, minTargets=minTargets,
         minPathPct=minPathPct, minPathSize=minPathSize, drug.gs=drug.gs, drug.adj=drug.adj, gene.gs=gene.gs, gene.adj=gene.adj)
 
  drugConnect <- panDrugConnect(caseids=caseids, controlids=controlids, gcount=gcount, gageCompare=gageCompare, tailEnd=tailEnd, drug.adj=drug.adj)
  drugTable <- merge(drugScores, drugConnect,by.x="Drug", by.y="Drug", all.x=TRUE, all.y=TRUE)
  ## Function to get ZSim where Z.stat compared to Z.stat.g*,
  ##     using Z.stat.g* = si x Zg i in 1..nsim
  colnms <- c("Drug", "N.Cancer.Genes", "Cancer.Genes", "N.Network.Genes", "Network.Genes", "Network","Network.pval","Z.Meta",
              "ZSim.pval", "PScore", "N.Pathways", "Pathways")
  ## order by sum of log10(p) of two p-values
  drugTable$PScore <- -1*log10(ifelse(drugTable$ZSim.pval==0,1/(2*nsim), drugTable$ZSim.pval)) - 1*log10(drugTable$Network.pval)
  return(drugTable[with(drugTable, order(PScore,decreasing=TRUE)),colnms])
}

panDrugStouffers <- function(panGene, nsim=1000, tailEnd="both", minTargets=1,
                             minPathPct=.05, minPathSize=5, drug.gs, drug.adj, gene.gs, gene.adj) {
  uDrugs <- names(drug.gs)
  #uDrugs <- unique(c(unlist(strsplit(panGene$Drugs, split="/")),unlist(strsplit(panGene$NetDrugs, split="/"))))
  #uDrugs <- uDrugs[!is.na(uDrugs) & nchar(uDrugs)>1]
  drugsdf <- data.frame(Drug=uDrugs, N.Cancer.Genes=0,Cancer.Genes="", N.Network.Genes=0, Network.Genes="",
                        Z.Meta=0, ZSim.pval=1, N.Pathways=0, Pathways="")

  ## for permutation-based p-values, permute the gage network p-values
  panGene$Network.mean.std <- (panGene$Network.mean - mean(panGene$Network.mean,na.rm=TRUE))/sd(panGene$Network.mean,na.rm=TRUE)
  panGene$Network.mean.std[is.na(panGene$Network.mean.std)] <- 0
  zMeans.perm <- matrix(sample(panGene$Network.mean.std, nsim*nrow(panGene), replace=TRUE),
                     nrow=nrow(panGene))
  
  zMeans.sim <- matrix(panGene$Network.mean.std*rnorm(nrow(panGene)*nsim), nrow=nrow(panGene), ncol= nsim)
  zTest.sim <- matrix(0, nrow=nrow(drugsdf), ncol=nsim)
  ## taking out ZPerm p-value
  #zTest.perm <- matrix(0, nrow=nrow(drugsdf), ncol=nsim)
  for(dr in 1:nrow(drugsdf)) {   
    drividx <- grep(drugsdf[dr,"Drug"], panGene$Drugs)
    drugsdf[dr,"N.Cancer.Genes"] <- length(unique(panGene$Cancer.Gene[drividx]))
    drugsdf[dr,"Cancer.Genes"] <- paste(unique(panGene$Cancer.Gene[drividx]),sep=",",collapse=",")

    conngenes <- unique(unlist(strsplit(panGene[grep(drugsdf[dr,"Drug"], panGene$NetDrugs),"NetGenes"],split=",")))
 ##browser()
    if(length(conngenes)>0) {
      conngenes <- conngenes[conngenes %in% names(which(abs(drug.adj[dr, ])>0))]
      drugsdf[dr,"N.Network.Genes"] <- length(conngenes)
      drugsdf[dr,"Network.Genes"] <- paste(conngenes, collapse=",")
    }    
    gset <- unique(c(panGene$Cancer.Gene[drividx], conngenes))
    gmatch <- sapply(gene.gs, function(x) sum(gset %in% x))
    gtotal <- sapply(gene.gs, length)
    gpct <- signif(gmatch/gtotal,2)

    drpathidx <- which(gpct >= minPathPct & gtotal >= minPathSize)
    if(length(drpathidx)>0) {    ##   any(gpct >= minPathPct & gtotal >= minPathSize)) {
      drugsdf$N.Pathways[dr] <- length(drpathidx)
      drugsdf$Pathways[dr] <- paste(names(gene.gs)[drpathidx], collapse="; ")                                         
    }    

    
    drugName <- drugsdf[dr,"Drug"]
    drugTargets <- drug.gs[[drugName]][which(drug.gs[[drugName]] %in% colnames(gene.adj))]
    drugTargetNetworks <- rownames(gene.adj)[rowSums(gene.adj[,drugTargets,drop=FALSE])>0,drop=FALSE]
    dTN.ngenes <- rowSums(gene.adj[drugTargetNetworks,drugTargets,drop=FALSE])
    dTN.size <- rowSums(gene.adj[drugTargetNetworks,,drop=FALSE])
    dTN.size <- dTN.size[names(dTN.size) %in% panGene$Cancer.Gene]
    dTN.ngenes <- dTN.ngenes[names(dTN.ngenes) %in% panGene$Cancer.Gene]
    dTN.idx <- match(names(dTN.size), panGene$Cancer.Gene)
    wts <- sqrt(dTN.ngenes/dTN.size)
    if(length(dTN.idx)>0) {
      drugsdf[dr,"Z.Meta"] <- sum(panGene[dTN.idx,"Network.mean.std"]*wts)/sum(wts^2)
      zTest.sim[dr,] <- apply(zMeans.sim[dTN.idx,,drop=FALSE], 2, function(x) sum(x*wts)/sum(wts^2))
#      zTest.perm[dr,] <- apply(zMeans.perm[dTN.idx,,drop=FALSE], 2, function(x) sum(x*wts)/sum(wts^2))
    }
  }
  
  if(tailEnd=="upper") {
   # drugsdf[,"ZPerm.pval"] <- rowSums(zTest.perm >= drugsdf[,"Z.Meta"])/nsim
    drugsdf[,"ZSim.pval"] <- rowSums(zTest.sim >= drugsdf[,"Z.Meta"])/nsim
   }
  if(tailEnd=="lower") {
   # drugsdf[,"ZPerm.pval"] <- rowSums(zTest.perm <= drugsdf[,"Z.Meta"])/nsim
    drugsdf[,"ZSim.pval"] <- rowSums(zTest.sim <= drugsdf[,"Z.Meta"])/nsim
  }
  if(tailEnd=="both") {
   #  drugsdf[,"ZPerm.pval"] <- rowSums(zTest.perm^2 >= drugsdf[,"Z.Meta"]^2)/nsim
     drugsdf[,"ZSim.pval"] <- rowSums(zTest.sim^2 >= drugsdf[,"Z.Meta"]^2)/nsim  
   }
  return(drugsdf)
}


  
## Purpose: similar steps to run gage as done in panConnect, but on gene networks that are targets per drug
panDrugConnect <- function(caseids, controlids, gcount, tailEnd="both", gageCompare=ifelse(length(caseids)>1,"as.group","unpaired"), drug.adj) {

  tailEnd <- casefold(tailEnd)
  if(!(tailEnd %in% c("upper","lower","both"))) {
    warning("invalid tailEnd, setting to 'both'.")
  }
   
  drivers <- data.frame(Drug=rownames(drug.adj), Network=0, Network.pval=1) ## Network.pval.both=1)
     
  warnset <- options()$warn
  options(warn = -1)

  for(k in 1:nrow(drivers)) {

    gconnect <- colnames(drug.adj)[which(drug.adj[drivers[k,1],]>0)]
    gconnect <- gconnect[gconnect %in% rownames(gcount)]
    kcase <- which(colnames(gcount) %in% caseids)
    kcontrol <- which(colnames(gcount) %in% controlids)
                  
    gage.direct <- gage(gcount[,c(kcontrol,kcase)], gsets = list(gconnect),
                        ref = 1:length(kcontrol), same.dir=TRUE,set.size=c(1,500),
                        samp = (1+length(kcontrol)):(length(kcontrol)+length(kcase)),
                        saaTest=gs.zTest,compare =gageCompare)
    gage.both <- gage(gcount[,c(kcontrol,kcase)], gsets = list(gconnect),
                        ref = 1:length(kcontrol), same.dir=FALSE,set.size=c(1,500),
                        samp = (1+length(kcontrol)):(length(kcontrol)+length(kcase)),
                      saaTest=gs.zTest, compare = gageCompare) #"unpaired")
    if(tailEnd=="upper") {    
      drivers$Network.pval[k] <- gage.direct$greater[1,'p.val']
    }
    if(tailEnd=="lower") {
      drivers$Network.pval[k] <- gage.direct$less[1,'p.val']
    }
    if(tailEnd=="both") {
      drivers$Network.pval[k] <- gage.both$greater[1,'p.val']
    }   
   
    ##  drivers$mut.events[k] <- sum(c(caseMut.ce, controlMut.ce))
    drivers$Network[k] <- length(gconnect)
   
  } # end for()
  options(warn=warnset)
  drivers <- drivers[order(drivers$Network.pval,decreasing=FALSE),]
  return(drivers)
}


