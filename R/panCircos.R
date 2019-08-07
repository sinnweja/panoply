panCircos <- function (panGene, panDrug, caseids, variant = NULL, cna = NULL, 
    gcount, gcinfo, tumorpct = 0.5, tailPct = 0.05, tailEnd = "upper", 
    minTargets = 1, minPathPct = 0.05, minPathSize = 8, minPathways = 1, 
    ...) 
{
    requireNamespace("RColorBrewer")
    nolegend = FALSE
    drTable <- panDrug
    drTable <- drTable[!duplicated(drTable[, c("Cancer.Genes", 
        "Network.Genes")]), ]
    drTable <- drTable[1:min(4, nrow(drTable)), ]
    drTable$AllGenes <- gsub(",$", "", gsub("^,", "", paste(drTable$Cancer.Genes, 
        drTable$Network.Genes, sep = ",")))
    legendvec <- character()
    pchvec <- numeric()
    ltyvec <- numeric()
    lwdvec <- numeric()
    colvec <- character()
    if (!is.null(cna) && nrow(cna) > 1) {
        legendvec <- c(legendvec, "CN-Amp", "CN-Del")
        pchvec <- c(pchvec, NA, NA)
        ltyvec <- c(ltyvec, 1, 1)
        lwdvec <- c(lwdvec, 3, 3)
        colvec <- c(colvec, "red", "green")
        cnacut1 <- log2((2 * (1 - tumorpct) + 1 * tumorpct)/2)
        cnacut3 <- log2((2 * (1 - tumorpct) + 3 * tumorpct)/2)
        cnidx <- which(apply(cna[, caseids, drop = FALSE] <= 
            cnacut1, 1, any) | apply(cna[, caseids, drop = FALSE] >= 
            cnacut3, 1, any))
        if (length(cnidx) == 0) {
            cnidx <- which(abs(cna[, caseids]) == max(abs(cna[, 
                caseids])))
        }
        cncols <- c("CHROM", "START", "STOP", "Gene.Symbol", 
            caseids)
        cngene <- cna[cnidx, cncols]
        cngene <- cngene[, c("CHROM", "START", "STOP", "Gene.Symbol", 
            caseids)]
        cngene <- cngene[!is.na(cngene$CHROM), ]
        cngene$gene_chr <- paste("chr", cngene$CHROM, sep = "")
        cngene$Driv.score <- ifelse(rowMeans(cngene[, caseids, 
            drop = FALSE]) > 0, 1, -1) * ifelse(cngene$Gene.Symbol %in% 
            panGene$Cancer.Gene, 2, 1)
    }
    if (!is.null(variant) && nrow(variant) > 0) {
        if (sum(variant$SampleType == "Germline") > 0) {
            legendvec <- c(legendvec, "Germline")
            pchvec <- c(pchvec, 18)
            ltyvec <- c(ltyvec, 0)
            lwdvec <- c(lwdvec, 0)
            colvec <- c(colvec, "sienna")
            germdf <- variant[variant$SampleType == "Germline" & 
                variant$PatientID %in% caseids, c("CHROM", "POS", 
                "Gene.Symbol")]
            germdf <- germdf[!is.na(germdf$CHROM), ]
            germdf$Driv.score <- ifelse(germdf$Gene.Symbol %in% 
                panGene$Cancer.Gene, 1, 2)
        }
        if (sum(variant$SampleType == "Tumor") > 0) {
            legendvec <- c(legendvec, "Somatic")
            pchvec <- c(pchvec, 18)
            ltyvec <- c(ltyvec, 0)
            lwdvec <- c(lwdvec, 0)
            colvec <- c(colvec, "orange")
            somdf <- variant[variant$SampleType == "Tumor" & 
                variant$PatientID %in% caseids, c("CHROM", "POS", 
                "Gene.Symbol")]
            somdf <- somdf[!is.na(somdf$CHROM), ]
            somdf$Driv.score <- ifelse(somdf$Gene.Symbol %in% 
                panGene$Cancer.Gene, 1, 2)
        }
    }
    if (ncol(gcount) > 1) {
        pchvec <- c(pchvec, NA)
        outMat <- outRNA(ids = caseids, gcount, tailPct = tailPct, 
            tailEnd = tailEnd)
        gcount <- gcount[rownames(gcount) %in% colnames(outMat)[colSums(outMat[caseids, 
            , drop = FALSE]) >= 1], ]
        gcgene <- data.frame(gene = rownames(gcount), expr = rowMeans(gcount[, 
            caseids, drop = FALSE]))
        drGeneList <- strsplit(drTable$AllGenes, split = ",")
        names(drGeneList) <- drTable$Drug
        gcgene$gcColor <- "white"
        drugBlues <- brewer.pal(6, "Blues")[6:3]
        for (k in length(drGeneList):1) {
            gcgene$gcColor[gcgene$gene %in% drGeneList[[k]]] <- drugBlues[k]
        }
        if (sum(gcgene$gcColor != "white") > 0) {
            gcgene <- gcgene[gcgene$gcColor != "white", ]
        }
        gcgene <- merge(gcgene, gcinfo[, c("Gene.Symbol", "CHROM", 
            "START")], by.x = "gene", by.y = "Gene.Symbol", all.x = TRUE, 
            all.y = FALSE)
        gcgene <- gcgene[!is.na(gcgene$CHROM), ]
        gcgene$gene_chr <- paste("chr", gcgene$CHROM, sep = "")
    }
    par(xpd = TRUE, mar = c(0, 1, 0, 1), ...)
    circos.initializeWithIdeogram()
    circos.trackPlotRegion(ylim = c(-1, 1), bg.border = NA, bg.col = "#EEEEEE", 
        track.height = 0.15, panel.fun = function(x, y) {
            xlim = get.cell.meta.data("xlim")
            ylim = get.cell.meta.data("ylim")
            for (i in -2:2) {
                circos.lines(xlim, c(i, i)/2, col = "#999999", 
                  lwd = 0.2)
            }
            xrange = get.cell.meta.data("xrange")
        })
    maxcn <- ceiling(max(abs(cngene[, caseids]), na.rm = TRUE))
    for (i in -1:1) {
        circos.text(0, i, labels = maxcn * i, sector.index = "chr1", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxcn * i, sector.index = "chr5", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxcn * i, sector.index = "chr9", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxcn * i, sector.index = "chr15", 
            col = "#777777", cex = 0.75)
    }
    for (i in 1:nrow(cngene)) {
        for (k in 1:length(caseids)) {
            circos.lines(unlist(cngene[i, 2:3]), rep(cngene[i, 
                caseids[k]]/maxcn, 2), sector.index = cngene[i, 
                "gene_chr"], lty = 1, lwd = abs(cngene[i, caseids[k]]) * 
                10, col = ifelse(cngene[i, caseids[k]] > 0, "red", 
                "green"))
        }
    }
    if (exists("germdf") | exists("somdf")) {
        circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, 
            bg.col = "#EEEEEE", track.height = 0.1, panel.fun = function(x, 
                y) {
                xlim = get.cell.meta.data("xlim")
                ylim = get.cell.meta.data("ylim")
                for (i in 0:1) {
                  circos.lines(xlim, c(i, i)/2, col = "#999999", 
                    lwd = 0.2)
                }
                xrange = get.cell.meta.data("xrange")
            })
        if (exists("germdf") && nrow(germdf) > 0) {
            for (i in 1:nrow(germdf)[1]) {
                circos.points(germdf$POS[i], 0.8, sector.index = germdf$CHROM[i], 
                  col = "sienna", pch = 18, cex = 0.75 * germdf$Driv.score[i])
            }
        }
        if (exists("somdf") && nrow(somdf) > 0) {
            for (i in 1:nrow(somdf)) {
                circos.points(somdf$POS[i], 0.4, sector.index = somdf$CHROM[i], 
                  col = "orange", pch = 18, cex = 0.75 * somdf$Driv.score[i])
            }
        }
    }
    circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, bg.col = "#EEEEEE", 
        track.height = 0.15, panel.fun = function(x, y) {
            xlim = get.cell.meta.data("xlim")
            ylim = get.cell.meta.data("ylim")
            for (i in 0:2) {
                circos.lines(xlim, c(i, i)/2, col = "#999999", 
                  lwd = 0.2)
            }
            xrange = get.cell.meta.data("xrange")
        })
    maxgc <- 1
    for (i in 0:1) {
        circos.text(0, i, labels = maxgc * i, sector.index = "chr1", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxgc * i, sector.index = "chr5", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxgc * i, sector.index = "chr9", 
            col = "#777777", cex = 0.75)
        circos.text(0, i, labels = maxgc * i, sector.index = "chr15", 
            col = "#777777", cex = 0.75)
    }
    for (i in 1:nrow(gcgene)) {
        circos.lines(x = as.numeric(rep(gcgene$START[i], 2)), 
            y = c(0, 1), sector.index = gcgene$gene_chr[i], lty = 1, 
            lwd = 10, col = gcgene$gcColor[i])
    }
    panGene <- panGene[order(panGene$Total.Druggable, decreasing = TRUE), 
        ]
    drivgenesplot <- panGene[1:min(nrow(panGene), 4), "Cancer.Gene"]
    netPurples <- brewer.pal(6, "Purples")[-c(1, 2)]
    netlty <- c(1:4)
    netlwd <- seq(1, 2, by = 0.25)
    panGene <- panGene[panGene$Cancer.Gene %in% drivgenesplot, 
        ]
    for (i in 1:nrow(panGene)) {
        conngenes <- colnames(reactome.adj)[reactome.adj[panGene$Cancer.Gene[i], 
            ] > 0]
        conngenes <- conngenes[conngenes %in% unlist(strsplit(panDrug[, 
            "Network.Genes"], split = ","))]
        if (length(conngenes) > 0 & any(conngenes %in% gcinfo$Gene.Symbol)) {
            conngenes.info <- gcinfo[gcinfo$Gene.Symbol %in% 
                conngenes, , drop = FALSE]
            drivgene.info <- gcinfo[gcinfo$Gene.Symbol == panGene$Cancer.Gene[i], 
                , drop = FALSE]
            for (j in 1:nrow(conngenes.info)) {
                circos.link(sector.index1 = paste("chr", drivgene.info$CHROM, 
                  sep = ""), point1 = as.numeric(drivgene.info$START), 
                  sector.index2 = paste("chr", conngenes.info$CHROM[j], 
                    sep = ""), point2 = as.numeric(conngenes.info$START[j]), 
                  col = netPurples[which(drivgene.info$Gene.Symbol == 
                    panGene$Cancer.Gene)], lwd = netlwd[which(drivgene.info$Gene.Symbol == 
                    panGene$Cancer.Gene)], lty = netlty[which(drivgene.info$Gene.Symbol == 
                    panGene$Cancer.Gene)])
            }
        }
    }
    if (!nolegend) {
        legend("bottomright", legend = legendvec, pch = pchvec, 
            lty = ltyvec, lwd = lwdvec, col = colvec, bty = "n", 
            inset = c(0, 0))
        legend("bottomleft", legend = c("Genes", panGene$Cancer.Gene), 
            col = c("white", netPurples), lty = c(1, netlty), 
            lwd = c(1, netlwd), bty = "n", inset = c(0, 0))
        legend("topleft", legend = c("Drug Targets", names(drGeneList)), 
            col = c("white", drugBlues), lty = 1, lwd = 10, bty = "n", 
            inset = c(0, 0), cex = 0.8)
    }
    invisible()
}
