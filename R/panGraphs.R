

## Graph networks of pan-connect plot
panGeneGraph <- function(panGene, panDrug, minTargets=2, minPathways=2, ndrugs=4, ndrivers=20, ...) {
  ## g count
  require(Rgraphviz)
  data(react.graph)
 
  ## load("/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/applications/Rlib/panoply/data/graphNEL_ad.RData")

  drTable <- panDrug[with(panDrug,!is.na(N.Pathways) & N.Pathways >= minPathways 
                     & !is.na(N.Network.Genes) & N.Network.Genes >= minTargets),]
  drTable <- drTable[!duplicated(drTable[,c("Cancer.Genes","Network.Genes")]),]
  drTable <- drTable[1:max(2,min(nrow(drTable),ndrugs)),]
  ### drTable <- panDrugSets(panGene, caseids, controlids, gcount=gcount, minTargets=minTargets, minPathPct=minPathPct, minPathSize=minPathSize, minPathways=minPathways, nperm=500)
  ##  drTable <- panDrugTable(pconn, minTargets=minTargets, minPathPct=minPathPct, minPathSize=minPathSize)
  drTable$AllGenes <- gsub(",$","",gsub("^,","",paste(drTable$Cancer.Genes, drTable$Network.Genes, sep=",")))
  drGeneList <- strsplit(drTable$AllGenes,split=",")
  names(drGeneList) <- drTable$Drug
  drugGenes <- unique(unlist(strsplit(drTable$AllGenes,split=",")))

  zdr <- panGene$Cancer.Gene[1:ndrivers]
  zexpr <- drugGenes[drugGenes %in% drTable$Network.Genes]
  zconn <- unlist(strsplit(panGene$NetGenes,split=","))
  ##  zconn <- colnames(reactome.adj)[colSums(reactome.adj[zdr,])>0]
  ##  zconn <- zconn[zconn %in% dgidbPan$Gene]
  ##  zconn.out <- outRNA(ids=id, gcount, tailPct=tailPct, tailEnd=tailEnd)
  ##  zconn <- zconn[zconn %in% colnames(zconn.out)[zconn.out[1,]>0]]
  zconn <- zconn[zconn %in% drugGenes]

  ## driver List  for sub graph selection
  ##z2=c('MED12','CDK8','LCK','PPARG','LYN','GALNT12','SYK','PTEN','GNAS','RAN','FGR','PRSS1','XPO1','PTPN11','HSD3B1','CEBPA','CD274','PDCD1LG2','PTPRJ','AKT2')
  ##  zvec <- unique(c(zdr,zconn))
  
  sg1 <- subGraph(unique(c(zdr,zconn)), react.graph)
  ##  nodes(sg1) <- paste0(nodes(sg1), "*")
  sgAttrs <- list()

  zdr.target <- panGene$Cancer.Gene[panGene$Cancer.Gene %in%  zexpr]
 
  ## rn = rownames(reactome.adj) more precisely rownames(reactome.adj)[apply(reactome.adj,1,sum)!=0]
  sgAttrs$fillcolor=rep(rgb(128/256,177/256,211/256), length(nodes(sg1)))  # blue
  names(sgAttrs$fillcolor)=nodes(sg1)
  ##,ifelse(nodes(sg1) %in% zexpr, "*", ""))
  sgAttrs$fillcolor[intersect(zdr, nodes(sg1))]=rgb(251/256,128/256,114/256) #red
  sgAttrs$shape=rep('box', length(nodes(sg1)))   # box for all
  names(sgAttrs$shape)=nodes(sg1)
  sgAttrs$shape[unique(c(zconn,zdr.target))]='circle'  # for netgenes (druggable only)
  sgAttrs$shape[intersect(zdr.target,nodes(sg1))]='ellipse' # for druggable expressed drivers
  
  sgAttrs$fontsize=rep(25, length(nodes(sg1)))
  names(sgAttrs$fontsize)=nodes(sg1)
  sgAttrs$height=rep(1.5, length(nodes(sg1)))
  names(sgAttrs$height)=nodes(sg1)
  sgAttrs$width=rep(2.2, length(nodes(sg1)))
  names(sgAttrs$width)=nodes(sg1)
  
  ## options are: circle (the default), ellipse, plaintext, and box
  ##png('DriverNetwork.png', height=6, width=6, units='in', res=300)
  plot(sg1, nodeAttrs=sgAttrs,
       main="Cancer Driver Network:\n red=drivers, blue=network \n circles=druggable, ellipse=expressed druggable drivers")
}


panDrugGraph <- function(panDrug,  usePathIDs=TRUE, ndrugs=8)  {
 
  panDrug <- panDrug[1:min(ndrugs, nrow(panDrug)),]
  drugs <- panDrug$Drug
  genes <- unique(c(unlist(strsplit(panDrug$Cancer.Genes,split=",")),unlist(strsplit(panDrug$Network.Genes,split=","))))
  pathfull <- pathways <- unique(unlist(strsplit(panDrug$Pathways, split="; ")))
   
  if(usePathIDs) {
##    panDrug$path <- sapply(pathways,function(x) unlist(strsplit(x,":"))[[1]])
    pathways <- unlist(sapply(strsplit(pathways, split=": "), function(x) x[1]))
    pathways <- gsub("R-","",pathways)
  } else {
    pathways <- unlist(sapply(strsplit(pathways, split=": "), function(x) x[2]))
    pathways <- gsub("Negative", "Negt", gsub("Signaling by", "Sig", gsub("regulation of", "reg",pathways)))
    pathways <- substring(pathways,1,13)
  }
##  pathways <- unique(panDrug$path)
#  reg <- new("graphNEL",nodes=c(drugs,genes,pathways),edgemode="directed")
  reg <- new("graphNEL", nodes = unique(c(drugs, genes, pathways)),edgemode="directed")

  for(i in c(1:nrow(panDrug))){
     tmpdrug <- panDrug[i,"Drug"]
     tmpgenes <- unique(c(unlist(strsplit(panDrug[i,"Cancer.Genes"],",")),unlist(strsplit(panDrug[i,"Network.Genes"],split=","))))
     tmppath <- strsplit(panDrug[i,"Pathways"], split="; ")[[1]]
         
     for(j in c(1:length(tmpgenes))){
       reg <- addEdge(tmpdrug,tmpgenes[j],reg,1)
       
       for(kpath in pathways[pathfull %in% tmppath]) {
         reg <- addEdge(tmpgenes[j],kpath,reg,1)
       }
     }
   }
  set.seed(100) 
  sgAttrs <- list()
  sgAttrs$fontsize=c(rep(35, length(drugs)),rep(35,length(genes)),rep(25,length(pathways)))
  names(sgAttrs$fontsize)=nodes(reg)
  sgAttrs$height=c(rep(3, length(drugs)),rep(3,length(genes)),rep(3,length(pathways)))
  names(sgAttrs$height)=nodes(reg)
  sgAttrs$width=c(rep(5, length(drugs)),rep(5,length(genes)),rep(5,length(pathways)))
  names(sgAttrs$width)=nodes(reg)
  sgAttrs$shape=c(rep("ellipse",length(drugs)),rep("ellipse",length(genes)),rep("rectangle",length(pathways)))
 
  
  names(sgAttrs$shape)=nodes(reg)
  sgAttrs$fillcolor=c(rep("red",length(drugs)),rep("dodger blue",length(genes)),rep("green",length(pathways)))
  names(sgAttrs$fillcolor)=nodes(reg)


#  width1 <- c(rep(1, length(drugs)),rep(1,length(genes)),rep(1.2,length(pathways)))
  width2 <- c(rep(1.5, length(drugs)),rep(1,length(genes)),rep(1.9,length(pathways)))

  sgAttrs$width=width2
  names(sgAttrs$width)=nodes(reg)
  plot(reg, nodeAttrs=sgAttrs)
}
