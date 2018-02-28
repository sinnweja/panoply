
## Purpose:plot method for panPath/panDruggable objects
## Authors: Erin Carlson and Jason Sinnwell, with much help from Xiaojia Tang
## ppath: panPath result
## ids: the caseids, or individual id use in panPath
## variant: germline and somatic variants in a data.frame
## cna: copy number aberration data.frame
## gcount: normalized gene expression data.frame
## gcinfo: gene info for gene expression genes
## controls: tell to get Cont.GenesPath from ppath object rather than Case

panCircos <- function(panGene, panDrug, caseids, variant=NULL, cna=NULL, gcount, gcinfo, tumorpct=0.5, tailPct=0.05, tailEnd="upper", minTargets=1, minPathPct=.05, minPathSize=8, minPathways=1, ...) {

  requireNamespace("RColorBrewer")
  nolegend=FALSE

  ## top drugs
  drTable <- panDrug

  drTable <- drTable[!duplicated(drTable[,c("Cancer.Genes","Network.Genes")]),]
  drTable <- drTable[1:min(4, nrow(drTable)),]
  drTable$AllGenes <- gsub(",$","",gsub("^,","",paste(drTable$Cancer.Genes, drTable$Network.Genes, sep=",")))
  ## from data given, vectors to say what is represented in regions
  legendvec <- character()
  pchvec <- numeric()
  ltyvec <- numeric()
  lwdvec <- numeric()
  colvec <- character()
  
  ## CNVs
  if(!is.null(cna) && nrow(cna) > 1) {
    legendvec <- c(legendvec,"CN-Amp","CN-Del")
    pchvec <- c(pchvec,NA,NA)
    ltyvec <- c(ltyvec,1,1)
    lwdvec <- c(lwdvec,3,3)
    colvec <- c(colvec,"red","green")
#    tumorpct <- .5
    cnacut1 <- log2((2*(1-tumorpct)+1*tumorpct)/2)
    cnacut3 <- log2((2*(1-tumorpct)+3*tumorpct)/2)
    cnidx <- which(apply(cna[,caseids,drop=FALSE] <= cnacut1, 1, any) | apply(cna[,caseids,drop=FALSE] >= cnacut3, 1, any))
    if(length(cnidx)==0) {
      ## put at least the most extreme CN region in
      cnidx <- which(abs(cna[,caseids]) == max(abs(cna[,caseids])))
    }
    cncols <- c("CHROM","START","STOP","Gene.Symbol",caseids)
    cngene <- cna[cnidx,cncols]
    cngene <- cngene[,c("CHROM","START","STOP","Gene.Symbol",caseids)]
    cngene <- cngene[!is.na(cngene$CHROM),]
    ## bed format with scores for Circos
    ## cnbed <-cngene[,c("CHROM","START","STOP","Gene.Symbol",caseids)]
 
    cngene$gene_chr <- paste("chr",cngene$CHROM,sep="")
    cngene$Driv.score <- ifelse(rowMeans(cngene[,caseids,drop=FALSE])>0, 1, -1) *
      ifelse(cngene$Gene.Symbol %in% panGene$Cancer.Gene, 2, 1)    
  }
  
  if(!is.null(variant) && nrow(variant)>0) { 
  ###germline and somatic variants
    
    ## could also subset to those in panGene$Cancer.Gene
    if(sum(variant$SampleType=="Blood")>0) {
      legendvec <- c(legendvec,"Germline")
      pchvec <- c(pchvec,18)
      ltyvec <- c(ltyvec,0)
      lwdvec <- c(lwdvec,0)
      colvec <- c(colvec,"blue")
      germdf <- variant[variant$SampleType=="Blood" & variant$PatientID %in% caseids,c("CHROM","POS","Effect","Gene.Symbol")]
      germdf <- germdf[!is.na(germdf$CHROM),]
      germdf$Driv.score <- ifelse(germdf$Gene.Symbol %in% panGene$Cancer.Gene, 1, 2)
    }
    if(sum(variant$SampleType=="Tissue")>0) {
      legendvec <- c(legendvec,"Somatic")
      pchvec <- c(pchvec,18)
      ltyvec <- c(ltyvec,0)
      lwdvec <- c(lwdvec,0)
      colvec <- c(colvec,"orange")
      somdf <- variant[variant$SampleType=="Tissue" & variant$PatientID %in% caseids,c("CHROM","POS","Effect","Gene.Symbol")]
      somdf <- somdf[!is.na(somdf$CHROM),]
      somdf$Driv.score <- ifelse(somdf$Gene.Symbol %in% panGene$Cancer.Gene, 1, 2) 
    }
  }
  if(ncol(gcount)>1) {
 ###   legendvec <- c(legendvec,"Diff-Expr")
    pchvec <- c(pchvec,NA)   
    outMat <- outRNA(ids=caseids, gcount, tailPct=tailPct, tailEnd=tailEnd)
    gcount <- gcount[rownames(gcount) %in% colnames(outMat)[colSums(outMat[caseids,,drop=FALSE]) >= 1],]
    gcgene <- data.frame(gene=rownames(gcount), expr=rowMeans(gcount[,caseids,drop=FALSE]))
    drGeneList <- strsplit(drTable$AllGenes,split=",")
    names(drGeneList) <- drTable$Drug
    ###    gcgene$Conn.score <- ifelse(gcgene$gene %in% dgidbPan$Gene, 2, 1)
    gcgene$gcColor <- "white"
  
    drugBlues <- brewer.pal(6,"Blues")[6:3]
    for(k in length(drGeneList):1) {
      gcgene$gcColor[gcgene$gene %in% drGeneList[[k]]] <- drugBlues[k]
    } 
    ## table(gcgene$Conn.score)
   
    if(sum(gcgene$gcColor != "white")>0) {
      gcgene <- gcgene[gcgene$gcColor != "white",]
    }
    gcgene <- merge(gcgene, gcinfo[,c("Gene.Symbol","CHROM","START")],
                    by.x="gene",by.y="Gene.Symbol", all.x=TRUE, all.y=FALSE)
    
    gcgene <- gcgene[!is.na(gcgene$CHROM),]
    gcgene$gene_chr <- paste("chr",gcgene$CHROM,sep="")  
  }

  par(xpd=TRUE, mar=c(0,1,0,1), ...)  #, mar=c(0,0,0,2)  # add at the beginning 
  circos.initializeWithIdeogram()

# copy number
  circos.trackPlotRegion(ylim = c(-1, 1), bg.border = NA, bg.col ="#EEEEEE",
                       track.height = 0.15, panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         ylim = get.cell.meta.data("ylim")
                         for(i in -2:2) {
                           circos.lines(xlim, c(i, i)/2, col = "#999999", lwd = 0.2)
                         }
                         xrange = get.cell.meta.data("xrange")
                       })
  maxcn <- ceiling(max(abs(cngene[,caseids]),na.rm=TRUE))
  ## add scale to plot
  for(i in -1:1)  {
    circos.text(0,i,labels = maxcn*i,sector.index = "chr1",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxcn*i,sector.index = "chr5",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxcn*i,sector.index = "chr9",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxcn*i,sector.index = "chr15",col="#777777",cex=0.75)
  }
  ## plot copy number genes
#  colcn <- sapply(rowMeans(cngene[,caseids,drop=FALSE]),function(x) if(x>0) "red" else "green" )
  for(i in 1:nrow(cngene)) {
    ## gene_start, gene_stio is the x coordinates,
    ## rep(EX**,2) is the y coordinates, gene_chr is the sector index, track by default. 
    ## lwd is match.score, col=red for positive, green for negative. lty is by defaul
    for(k in 1:length(caseids)) {
      circos.lines(unlist(cngene[i,2:3]), rep(cngene[i,caseids[k]]/maxcn,2),
                   sector.index = cngene[i,"gene_chr"],lty=1,lwd=abs(cngene[i,caseids[k]])*10,col=ifelse(cngene[i,caseids[k]]>0, "red","green")) 
    }
  }
  
  ## germbed and sombed
  if(exists("germdf") | exists("somdf")) {
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, bg.col ="#EEEEEE",
                         track.height = 0.1, panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           for(i in 0:1) {
                             circos.lines(xlim, c(i, i)/2, col = "#999999", lwd = 0.2)
                           }
                           xrange = get.cell.meta.data("xrange")
                         })

    if(exists("germdf") && nrow(germdf) > 0) {
      for(i in 1:nrow(germdf)[1]){
        ## pos is the x coordinates, Match.score is the y coordinates,
        ## chr is the sector index, track by default. 
        ## pch is defaut, col=green for germline. cex is effect.score
        circos.points(germdf$POS[i], .8, #germdf$Driv.score[i]/max(germdf$Driv.score),
           sector.index = germdf$CHROM[i],col="blue",pch=18,cex=0.75*germdf$Driv.score[i])
      }
    }
    if(exists("somdf") && nrow(somdf) > 0){
      for(i in 1:nrow(somdf)) {
        
        ##pos is the x coordinates, Match.score is the y coordinates,
        ## chr is the sector index, track by default. 
        ##pch is defaut, col=red for germline. cex is effect.score
        circos.points(somdf$POS[i], .4, #somdf$Driv.score[i]/max(somdf$Driv.score),
           sector.index = somdf$CHROM[i],col="orange", pch=18,cex=0.75*somdf$Driv.score[i]) 
      }
    }
  } 
  ## gene expression z-score
  circos.trackPlotRegion(ylim = c(0,1), bg.border = NA, bg.col ="#EEEEEE",
                         track.height = 0.15, panel.fun = function(x, y) {
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           for(i in 0:2) {
                             circos.lines(xlim, c(i, i)/2, col = "#999999", lwd = 0.2)
                           }
                           xrange = get.cell.meta.data("xrange")
                         })

  maxgc <- 1    ### max(gcgene$Conn.score)
  for(i in 0:1)  {
    circos.text(0,i,labels = maxgc*i,sector.index = "chr1",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxgc*i,sector.index = "chr5",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxgc*i,sector.index = "chr9",col="#777777",cex=0.75)
    circos.text(0,i,labels = maxgc*i,sector.index = "chr15",col="#777777",cex=0.75)
  }
 
###  colgc <- rep("pink3",nrow(gcgene)) ## if not all over-expressed, do this: sapply(gcgene$zscore,function(x) if(x>0) "red" else "green" )
  for(i in 1:nrow(gcgene)){
       ### color by drug
       circos.lines(x=as.numeric(rep(gcgene$START[i],2)), y=c(0,1), ### gcgene$Conn.score[i]/maxgc),
                   sector.index = gcgene$gene_chr[i],lty=1,lwd=10,col=gcgene$gcColor[i])
   ###  if(gcgene$Conn.score[i] > 0) {   
      ## gene_start is the x coordinates, c(0,zscore) is the y coordinates,
      ## gene_chr is the sector index, track by default. 
      ## lwd is 1, col=red for positive, green for negative. lty is by default      
   ###   circos.lines(x=as.numeric(rep(gcgene$START[i],2)), y=c(0,gcgene$Conn.score[i]/maxgc),
   ###                sector.index = gcgene$gene_chr[i],lty=1,lwd=3,col=colgc[i])
      
    ###}
  }
  ## lines connecting driver genes with connected genes
  panGene <- panGene[order(panGene$Total.Druggable,decreasing=TRUE),]
  drivgenesplot <- panGene[1:min(nrow(panGene),4),"Cancer.Gene"]
#  drivgenesplot <- which(!duplicated(panGene$Cancer.Gene))
#  drivgenesplot <- drivgenesplot[1:min(c(length(drivgenesplot),4))]
  netPurples <- brewer.pal(6,"Purples")[-c(1,2)]
#  clinecols <- c("limegreen","purple3","skyblue2","darkorange2")
    ## c('hotpink','pink','pink3','deeppink')
  netlty <- c(1:4)
  netlwd <- seq(1,2,by=.25)
  panGene <- panGene[panGene$Cancer.Gene %in% drivgenesplot,]
   
  for(i in 1:nrow(panGene)) {
    ## cat("\n i=", i, ";j=")
      conngenes <- colnames(reactome.adj)[reactome.adj[panGene$Cancer.Gene[i],]>0]
      conngenes <- conngenes[conngenes %in% unlist(strsplit( panDrug[,"Network.Genes"], split=","))] #dgidbPan$Gene]
      ##strsplit(panGene$Cont.PathGenes[i],",")[[1]]
      ## } else  {
      ##   conngenes <- strsplit(panGene$Case.PathGenes[i],",")[[1]]
      ## }
      if(length(conngenes)>0 & any(conngenes %in% gcinfo$Gene.Symbol)) {
        conngenes.info <- gcinfo[gcinfo$Gene.Symbol %in% conngenes,,drop=FALSE]
        drivgene.info <- gcinfo[gcinfo$Gene.Symbol==panGene$Cancer.Gene[i],,drop=FALSE]
        for(j in 1:nrow(conngenes.info)){
          circos.link(sector.index1 = paste("chr",drivgene.info$CHROM,sep=""),point1 = as.numeric(drivgene.info$START),
            sector.index2 = paste("chr",conngenes.info$CHROM[j],sep=""),point2 = as.numeric(conngenes.info$START[j]),
            col=netPurples[which(drivgene.info$Gene.Symbol==panGene$Cancer.Gene)],
            lwd=netlwd[which(drivgene.info$Gene.Symbol==panGene$Cancer.Gene)],
            lty=netlty[which(drivgene.info$Gene.Symbol==panGene$Cancer.Gene)])
        }
      } # if length(conngenes)
    } # for nrow panGene
  if(!nolegend) {
    legend("bottomright",legend=legendvec, pch=pchvec, lty=ltyvec, lwd=lwdvec,col=colvec,
           ##legend=c("CN-Amp","CN-Deletion","Somatic","Germline","GE:Over"),
           ##pch=c(NA,NA,18,18,NA,NA),lty=c(1,1,0,0,1),lwd=c(3,3,0,0,2),
           ##col=c("red","green","orange","blue","red"),
           bty="n",inset=c(0,0))
    
    legend("bottomleft", legend=c("Genes", panGene$Cancer.Gene),
           col=c("white", netPurples),
           lty=c(1,netlty), lwd=c(1,netlwd), bty="n", inset=c(0,0))
    
    legend("topleft", legend=c("Drug Targets", names(drGeneList)), col=c("white",drugBlues), lty=1, lwd=10, bty="n", inset=c(0,0),cex=.8)
  }
  invisible()
}

