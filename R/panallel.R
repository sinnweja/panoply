
## Make parallel coordinates plot to show gene expressions that go into why a drug is significant.
## drugNames takes the top 4 drugs, or a vector of drug names from the user.
panallelPrepare <- function(panGene, panDrug, caseids, controlids, gcount)
  {
                     
  ### data frame with one row for each drug/driver combo  
  
  glist <- dgi.gs[unique(panDrug$Drug)]
  mat1 <- NULL

  for(di in unique(panDrug$Drug)) {
    mat1 <- rbind(mat1,cbind(rep(di,length(glist[[di]])),glist[[di]]))
    mat1 <- mat1[mat1[,2] %in% colnames(reactome.adj) & mat1[,2] %in% rownames(gcount),,drop=FALSE]
  }
  colnames(mat1) <- c("Drug","Cancer.Gene")
  mat1 <- mat1[!duplicated(mat1),,drop=FALSE]
 #browser()
  mat2 <- data.frame()
  for(i in 1:nrow(mat1)){
    ##  connlist <- colnames(reactome.adj)[reactome.adj[dg,]>0 ]
    drivlist <- rownames(reactome.adj)[reactome.adj[,mat1[i,2]]>0]
    drivlist <- drivlist[drivlist %in% panGene$Cancer.Gene]
    if(length(drivlist)>0) {
      for(k in 1:length(drivlist)) {
      ## for checking  cat("i=",i, ", driver", drivlist[k], "\n")
        connlist <- colnames(reactome.adj)[reactome.adj[drivlist[k],]>0]
        ##  connlist <- colnames(reactome.adj)[reactome.adj[panGene$Cancer.Gene,mat1[i,2]]>0 ]
        connlist <- connlist[connlist %in% rownames(gcount)]
        if(length(connlist)>0){
          mat2 <- rbind(mat2,cbind(rep(mat1[i,1],length(connlist)),
                                   rep(drivlist[k],length(connlist)),connlist))
        }
      }
    }
  }  
 
  colnames(mat2) <- c("Drug","Cancer.Gene","Network.Gene")

  mat2 <- as.data.frame(mat2)
#  mat <- data.frame(mat2,row.names=NULL)
##  mat$Drug.PScore <- panDrug[match(mat$Drug,panDrug$Drug),"PScore"]
  minzperm <- min(panDrug$ZPerm.pval[panDrug$ZPerm.pval>0])
  mat2$Drug.log10ZPerm <- -log10(pmax(panDrug[match(mat2[,"Drug"],panDrug$Drug),"ZPerm.pval"],minzperm/2))
  mat2$Cancer.Gene.NetMean <- panGene[match(mat2$Cancer.Gene,panGene$Cancer.Gene),"Network.mean"]
  mat2$Network.Gene.count <- apply(gcount[mat2$Network.Gene,caseids,drop=FALSE], 1, median)
  mat2$Cancer.Gene.NetMeanP <- -log10(panGene[match(mat2$Cancer.Gene,panGene$Cancer.Gene),"Network.pval"])
  panDrug$Network.pval <- ifelse(is.na(panDrug$Network.pval), 1.01, panDrug$Network.pval)
  mat2$Drug.log10NetP <- -log10(panDrug[match(mat2$Drug,panDrug$Drug),"Network.pval"])
 # mat$Drug.log10NetP <- ifelse(mat$Cancer.Gene==mat$Network.Gene,mat$Drug.log10NetP,NA)
  return(mat2)
  
}

panallelPlot <- function(mat2,drugNames){
  tmpgenes <- unlist(dgi.gs[drugNames])
  tmpmat2 <- rbind(mat2[mat2$Drug==drugNames & !(mat2$Network.Gene %in% tmpgenes),],
                   mat2[mat2$Drug==drugNames & mat2$Network.Gene %in% tmpgenes,])
  mat2 <- mat2[!duplicated(mat2[,c(2,3,4,7)]),]
  mat2 <- rbind(mat2,tmpmat2)

  plot.col=ifelse(mat2$Drug==drugNames,"black","grey")
  plot.lty=ifelse(mat2$Drug==drugNames,1,2)
  plot.col=ifelse(plot.col=="black" & mat2$Network.Gene %in% tmpgenes,"red",plot.col)
  #table(plot.col)
  colnames(mat2)[c(8,6,5,4)] <- c("log(Drug.NetP)","Gene.Expr","Driver.Network.DE","log(Drug.PermP)")

  parcoord(mat2[,c(8,6,5,4)],col=plot.col,lty=plot.lty,var.label=T,xaxt='n',cex=0.8)
#  axis(1,at=c(1,2,3,4),labels=c("Log10(Drug.Net.Pval)","Gene.Expr","Driver.Network.DE","Log10(Drug.ZPerm.Pval"))
  title(drugNames)

}

## test code:
#unique(mat$Drug[mat[,4]>2.4])

#pdf("testpanollel.pdf")
#panollel(drivGenes, drugTable, patient, ptmatch, gcPanBC,drugName="OLAPARIB")
#panollel(drivGenes, drugTable, patient, ptmatch, gcPanBC,drugName="4SC-202")
#panollel(drivGenes, drugTable, patient, ptmatch, gcPanBC,drugName="DANUSERTIB")
#panollel(drivGenes, drugTable, patient, ptmatch, gcPanBC,drugName="BARASERTIB")
#dev.off()
#"4SC-202"
#"DANUSERTIB"
#"BARASERTIB"
