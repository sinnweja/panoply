
## Purpose: post-processing functions to summarize panPath and panDruggable results in tables.
## Author: Jason Sinnwell, Rani Kalari, Kevin Thompson
## Created: 2/2016

## Both these functions work on either panPath or panDruggable objects

## filter panDruggable results by quantiles of p-values, and the max times a gene or pathway shows up to show more results.
filterNetworks <- function(panDrug, netpmax=.2, connQuant=1, minDrug=1, cancerDrug=TRUE) {
 
  fpath <- panDrug[!is.na(panDrug$Network) & panDrug$Network>0,]

  if(any(grepl("Network.pval", names(fpath)))) {
     fpath <- fpath[with(fpath, !is.na(Network.pval) & Network.pval <= netpmax),]
                   # |  with(fpath, !is.na(Network.pval.both) & Network.pval.both <= netpmax),]
                         #< quantile(p.match.greater,gageQuant,na.rm=TRUE)),] 
  }

  ## filter by connQuant  
  ## fpath <- fpath[with(fpath,pval.MUT <= quantile(pval.MUT,connQuant,na.rm=TRUE) |
  ##                pval.WT <= quantile(pval.WT,connQuant,na.rm=TRUE)),] 
  
  if(cancerDrug & any(grepl("Drug", colnames(fpath)))) {
    drugpath <- data.frame()
    drugpath <- fpath[fpath$Druggable > 0,]
    fpath <- drugpath
  }
  if(any(grepl("Drug", colnames(fpath)))) {
    drugpath <- data.frame()
    drugpath <- fpath[fpath$Total.Druggable >= minDrug,]
    fpath <- drugpath
  }  

  ## order by pvals  (sum of 1og10(p)
  fpath <- fpath[order(with(fpath, pmin(Network.pval)), decreasing=FALSE),]

  return(fpath)
}

findPathways <- function(drTable, minTargets=1, minPathPct=.05, source="reactome") {

  ## only reactome for now... expand to kegg
  if(source=="reactome") {
    data(reactome)
    pathdf <- data.frame(Drug="",Pathway="", N.Genes=0,Druggable.Pct=0)
    for(k in 1:nrow(drTable)) {
      gset <- unique(unlist(strsplit(paste(drTable[k,"Cancer.Genes"],drTable[k,"Network.Genes"],collapse=",",sep=""),split=",")))
      gmatch <- sapply(reactome.gs, function(x) sum(gset %in% x))
      gtotal <- sapply(reactome.gs, length)
      gpct <- signif(gmatch/gtotal,2)
      if(any(gpct > minPathPct)) { 
        pathdf <- rbind.data.frame(pathdf, data.frame(Drug=drTable[k,"Drug"],
                                                    Pathway=names(gpct)[gpct > minPathPct],
                                                    N.Genes=gtotal[gpct > minPathPct],
                                                    Druggable.Pct=gpct[gpct > minPathPct]))
      }
    }
  }
  ## rm first row, and ensure minPct
  pathdf <- pathdf[pathdf$Druggable.Pct>minPathPct,]
  drTable <- merge(drTable, pathdf, by.x="Drug",by.y="Drug",all.x=TRUE)
  
  return(drTable[order(drTable$Druggable.Pct,decreasing=TRUE),])
}


## Return key info on the variant and/or cna that are driving expression events in pathways
findDrivers <- function(ppath, caseid, controlid, variant=NULL, cna=NULL, esnv=NULL, 
                        gcount=NULL, rppa=NULL, tumorpct=0.5, max.control=length(controlid),tailEnd="upper",tailPct=.2) {
  data(gcinfoPan)
  vardf <- cnadf <- data.frame(CHROM=character(), POS=numeric(), Gene.Symbol = character(),
                               Case.Type = character(), Control.Count=numeric())

  returncols <- names(vardf)[-2]
  if(!is.null(variant)) {
    vargenes <- unique(variant$Gene.Symbol[variant$Gene.Symbol %in% ppath$Cancer.Gene &
                                           variant$PatientID %in% caseid])    
    for(dr in vargenes) {
      vardf <- rbind.data.frame(vardf,
          data.frame(CHROM=variant[variant$PatientID %in% caseid & variant$Gene.Symbol==dr,"CHROM"],
                     POS=variant[variant$PatientID %in% caseid & variant$Gene.Symbol==dr,"POS"],
                     Gene.Symbol=dr,  Case.Type = paste0("DNA:",
             paste(variant[variant$PatientID %in% caseid & variant$Gene.Symbol==dr,"SampleType"],collapse=",")),
          Control.Count = length(unique(variant[variant$PatientID %in% controlid &
            variant$Gene.Symbol==dr,"PatientID"]))))
    }
    vardf <- vardf[vardf$Control.Count <= max.control,]
    if(!is.null(esnv)) {
      vardf <- merge(vardf, esnv[,c("POS","GENE","AD")], by.x=c("POS","Gene.Symbol"),
                     by.y=c("POS","GENE"), all.y=FALSE, all.x=TRUE)      
    }   
  }
  if(!is.null(cna)) {
    cnaDriv <- drivDNA(ids=c(caseid, controlid),cna=cna, variant=NULL,tumorpct=tumorpct)
    if(length(caseid)>1){      
      cnagenes <- colnames(cnaDriv)[colSums(cnaDriv[caseid,])>0 & colnames(cnaDriv) %in% ppath$Cancer.Gene]
    } else cnagenes <- colnames(cnaDriv)[cnaDriv[caseid,]>0 & colnames(cnaDriv) %in% ppath$Cancer.Gene]
    for(dr in cnagenes) {
      ##  print(dr)
      cnadf <- rbind.data.frame(cnadf, data.frame(CHROM=cna[cna$Gene.Symbol==dr,"CHROM"],
                                         POS=cna[cna$Gene.Symbol==dr,"START"],Gene.Symbol=dr,
         Case.Type = paste0("CN:", ifelse(rowMeans(cna[cna$Gene.Symbol==dr,caseid,drop=FALSE])>0,"Amp", "Del")),
         Control.Count = sum(cnaDriv[controlid,dr]>0)))
    }
    cnadf <- cnadf[cnadf$Control.Count <= max.control,]
    if(!is.null(esnv)) {      
      cnadf <- merge(cnadf, esnv[,c("GENE","AD")], by.x="Gene.Symbol", by.y="GENE", all.y=FALSE, all.x=TRUE)
    }
  }

  
  drivdf <- merge(vardf, cnadf, all.x=TRUE, all.y=TRUE)
  names(drivdf) <- gsub("AD","RNA.AD", names(drivdf))
  drivdf$CHROM <- gsub("^chr","", drivdf$CHROM)
  if(any(grepl("AD$", colnames(drivdf)))) {
    returncols <- c(returncols, "RNA.AD")
    drivdf[is.na(drivdf$RNA.AD),"RNA.AD"] <- " "
  }

  if(!is.null(gcount)) {
    gcount <- gcount[ppath$Cancer.Gene[ppath$Cancer.Gene %in% rownames(gcount)],]
    if(tailEnd=="both"){   
      gcdf <- data.frame(Gene.Symbol=rownames(gcount),Expressed=ifelse(rowMeans(gcount[,caseid, drop=FALSE]) > apply(gcount[,controlid,drop=FALSE], 1, quantile, prob=1-tailPct/2),"+",ifelse(rowMeans(gcount[,caseid, drop=FALSE]) < apply(gcount[,controlid,drop=FALSE], 1, quantile, prob=tailPct/2),"-",NA),GE.Rank=rank(gcount[,caseid])))

    }
    if(tailEnd=="upper"){   
      gcdf <- data.frame(Gene.Symbol=rownames(gcount),Expressed=ifelse(rowMeans(gcount[,caseid, drop=FALSE]) > apply(gcount[,controlid,drop=FALSE], 1, quantile, prob=1-tailPct),"+",NA),GE.Rank=rank(gcount[,caseid]))
    }
    gcdf <- merge(gcdf,gcinfoPan,by="Gene.Symbol")
    gcdf.keep <- which(gcdf$Gene.Symbol %in% drivdf$Gene.Symbol|!is.na(gcdf$Expressed))
    
    gcdf.genes <- gcdf[gcdf.keep,"Gene.Symbol"]
    drivdf <- merge(drivdf,gcdf[gcdf.keep,c("Gene.Symbol","Expressed","GE.Rank")],by="Gene.Symbol",keep=T,all=T)
    drivdf[is.na(drivdf$POS),c("POS","CHROM")] <- c(gcinfoPan[match(drivdf$Gene.Symbol[is.na(drivdf$POS)],gcinfoPan$Gene.Symbol),c("START","CHROM")])
                                        #    drivdf$Expressed <- ifelse(rowMeans(gcount[,caseid, drop=FALSE]) > rowMeans(gcount[,controlid, drop=FALSE]), "+", "-")
                                        #    drivdf$Expressed <- ifelse(apply(gcount[,caseid, drop=FALSE], 1, min) > apply(gcount[,controlid, drop=FALSE], 1, max), "++", drivdf$Expressed)
                                        #    drivdf$Expressed <- ifelse(apply(gcount[,caseid, drop=FALSE], 1, max) < apply(gcount[,controlid, drop=FALSE], 1, min), "--", drivdf$Expressed)
                                        #    returncols <- c(returncols, "Expressed")
    drivdf <- drivdf[order(drivdf$GE.Rank, decreasing=TRUE),]
    returncols <- c(returncols, "Expressed")
  }

  if(!is.null(rppa)) {   
    data(RPPAinfoPan)
    prot.idx <- which(!grepl("_p",RPPAinfoPan$Antiboty.Name))
    phos.idx <- grepl("_p",RPPAinfoPan$Antiboty.Name)
    Prot.diff <- apply(rppa[c(caseid, controlid),prot.idx],2,function(x) 
                   x[1]-mean(x[2:length(x)], na.rm=TRUE))
    Phos.diff <- apply(rppa[c(caseid, controlid),phos.idx],2,function(x) 
                   x[1]-mean(x[2:length(x)], na.rm=TRUE))

    RPPAdf <- data.frame(GeneSymbol=RPPAinfoPan$GeneSymbol[prot.idx],Prot.Diff=Prot.diff)
    RPPAdf <- merge(RPPAdf, data.frame(GeneSymbol=RPPAinfoPan$GeneSymbol[phos.idx],Phos.Diff=Phos.diff), by.x="GeneSymbol", by.y="GeneSymbol", all.x=TRUE, all.y=TRUE)
    RPPAdf <- RPPAdf[order(RPPAdf$GeneSymbol,decreasing=T),]    
    drivdf <- merge(drivdf, RPPAdf, by.x="Gene.Symbol",by.y="GeneSymbol",all.y=FALSE,all.x=TRUE)
    drivdf <- drivdf[order(drivdf$Prot.Diff, decreasing=TRUE,na.last=TRUE),]
    returncols <- c(returncols,"Prot.Diff","Phos.Diff")
  }
  
  drivdf <- drivdf[order(drivdf$CHROM, drivdf$Gene.Symbol, decreasing=FALSE),returncols]
  return(drivdf)
}
