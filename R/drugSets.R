
if(0) {
file1 <- "/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/dgiAntiNeo1.json"
library(RJSONIO)
dgiAntiNeo1 <- fromJSON(paste(readLines(file1), collapse=""))

dgiAntiNeo1 <- do.call(rbind,sapply(dgiAntiNeo1$matchedTerms,function(x) cbind(x$geneName,matrix(unlist(x$interactions),byrow=T,ncol=4))[,c(1,2,4,5)]))[,1:3]
colnames(dgiAntiNeo1)<-c("Gene","Drug","interactionType")

file2 <- "/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/dgiAntiNeo2.json"
dgiAntiNeo2 <- fromJSON(paste(readLines(file2), collapse=""))

dgiAntiNeo2 <- do.call(rbind,sapply(dgiAntiNeo2$matchedTerms,function(x) cbind(x$geneName,matrix(unlist(x$interactions),byrow=T,ncol=4))[,c(1,2,4,5)]))[,1:3]
colnames(dgiAntiNeo2)<-c("Gene","Drug","interactionType")

file3 <- "/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/dgiAntiNeo3.json"
dgiAntiNeo3 <- fromJSON(paste(readLines(file3), collapse=""))
dgiAntiNeo3 <- do.call(rbind,sapply(dgiAntiNeo3$matchedTerms,function(x) cbind(x$geneName,matrix(unlist(x$interactions),byrow=T,ncol=4))[,c(1,2,4,5)]))[,1:3]
colnames(dgiAntiNeo3)<-c("Gene","Drug","interactionType")

dgiAntiNeo <- rbind.data.frame(dgiAntiNeo1,dgiAntiNeo2,dgiAntiNeo3)

dim(dgiAntiNeo)
length(unique(dgiAntiNeo[,"Gene"]))
length(unique(dgiAntiNeo[,"Drug"]))
dgiAntiNeo <- dgiAntiNeo[!duplicated(dgiAntiNeo[,1:2]),]
dim(dgiAntiNeo)
dgiAntiNeo[,"Drug"] <- gsub("ENDOSTATIN (84-114)-NH2 (JKC367)","ENDOSTATIN",dgiAntiNeo[,"Drug"],fixed=TRUE)
dgidb <- dgiAntiNeo
dgidb$Source <- "DGIdb"
names(dgidb) <- gsub("interactionType", "type", names(dgidb))
## Drug Bank curated by Bioinformatcs paper authors

dbankLines <- readLines("/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/DrugGeneCurated.txt")
dbankStr <- strsplit(dbankLines, split="\t")
dbank <- do.call("rbind.data.frame", dbankStr[2:length(dbankStr)])
names(dbank) <- dbankStr[[1]]
dim(dbank)
colnames(dbank)

uniprot <- readLines("/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/mapped.results.txt")
length(uniprot)
nrow(dbank)
uniprot[1:10]
dbank$GeneID <- uniprot[2:length(uniprot)]
dbank[1:5,]
ugenes.db <- unique(unlist(strsplit(dbank[,"GeneID"], split=";")))
#write.table(dbank, sep="\t", quote=FALSE, file="/data2/bsi/tertiary/Weinshilboum_Richard_weinsh/s112047.beauty/analyses/manuscripts/panoply/paperFinal/panoplyDev/buildData/bionfPaperDrugs.tsv")

## Merge DGI Anti-Neoplastic and Curated DrugBank together

dbank$DRUG <- casefold(dbank$Drug_name, upper=TRUE)
table(unique(dgidb$Drug) %in% unique(dbank$DRUG))
table( unique(dbank$DRUG) %in% unique(dgidb$Drug) )
udrugs.dgi <- unique(c(dgidb$Drug,dbank$DRUG))
length(udrugs.dgi)
udrugs.dgi <- udrugs.dgi[!(grepl("\\[",udrugs.dgi) | grepl("\\{",udrugs.dgi) | grepl("\\(",udrugs.dgi))]
length(udrugs.dgi)
dsynonyms <- unique(unlist(strsplit(casefold(dbank$Drug_Synonyms,upper=TRUE), split="; ")))
table((unique(dgidb$Drug) %in% dsynonyms))
table((unique(dbank$DRUG) %in% dsynonyms))
length(dsynonyms)
#udrugs.dgi <- udrugs.dgi[!(udrugs.dgi %in% dsynonyms)]
length(udrugs.dgi)
glist <- strsplit(dbank$GeneID, split=";")
dbankdf <- data.frame(Drug=NULL, Gene=NULL, type=NULL, Source=NULL)
for(k in 1:nrow(dbank) ) {
  if(length(glist[[k]])>0) {
    dbankdf <- rbind.data.frame(dbankdf, data.frame(Drug=dbank$DRUG[k], Gene=glist[[k]], type="n/a", Source=dbank[k,"Annotation From"]))
  }
}

## check
#dgiAntiNeo[dgiAntiNeo$Gene %in% c("ERBB2","BRCA1","ATM","TP53","ATR","PTEN"),]
dgidb[dgidb$Gene %in% c("ERBB2","BRCA1","TP53","ATR","PTEN"),]
dgidb[dgidb[,3] %in% c("agonist"),]

dbankdf[dbankdf$Gene %in% c("ERBB2","BRCA1","TP53","ATR","PTEN"),]


drugdbPan <- rbind.data.frame(dgidb, dbankdf)
table(duplicated(drugdbPan[,1:2]))

annoDrugs <- annotateDrugs(drugdbPan)
annoDrugs <- annotateDrugs(dgiDrug, byDrug=TRUE)
dsets <- annoDrugs[[1]]
dadj <- annoDrugs[[2]]
table(rowSums(dadj))
table(colSums(dadj))
table(sapply(dsets, length))
table(rowSums(dgi.adj))
table(colSums(dgi.adj))
}

annotateDrugs <- function(drugPan, byDrug=FALSE) {

  ## Make drug-gene list and adjacency matrix for use in panDrugSets
  
  if(byDrug) {
    ## for drug/gene in this format (by drug)
    ##         Drug                              GeneID Synonyms DGIdb DBankCurate
    ## ELLAGIC ACID            PRKCA;SQLE;PRKCB;SYK;CA9           TRUE       FALSE
    ## BRYOSTATIN-1 PRKCA;PRKCB;PRKCQ;PRKCD;PRKCG;PRKCE           TRUE       FALSE
    allgenes <- unique(unlist(strsplit(drugPan[,2], split=";")))
    drug.adj <- matrix(0, nrow=length(unique(drugPan[,1])), ncol=length(allgenes))
    dimnames(drug.adj) <- list(unique(drugPan[,1]), allgenes)
    drug.gs <- list()
    for(drugk in unique(drugPan[,1])) {
      k <- which(drugPan[,1] %in% drugk)
      drug.gs[[drugk]] <- unique(unlist(strsplit(drugPan[k,2],split=";")))
      gidx <- which(colnames(drug.adj) %in% drug.gs[[drugPan[k,"Drug"]]])
      drug.adj[drugk,gidx] <- 1
    }
  } else {
    ## this format: 
    ##       Gene         Drug                  type Source
    ## PRKCA ELLAGIC ACID inhibitor,competitive  DGIdb
    ## PRKCA BRYOSTATIN-1                   n/a  DGIdb
    ## PRKCA   SOPHORETIN             inhibitor  DGIdb
    allgenes <- unique(drugPan[,1])
    alldrugs <- unique(drugPan[,2])
    drug.adj <- matrix(0, nrow=length(alldrugs), ncol=length(allgenes))
    dimnames(drug.adj) <- list(alldrugs, allgenes)
    drug.gs <- list()
    for(drugk in alldrugs) {
      k <- which(drugPan[,2] %in% drugk)
      drug.gs[[drugk]] <- unique(drugPan[k,1])
      gidx <- which(colnames(drug.adj) %in% drug.gs[[drugk]])
      drug.adj[drugk,gidx] <- 1
    }    
  }
  return(list(drug.gs, drug.adj))
}





