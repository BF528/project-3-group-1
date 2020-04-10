par(mar=c(4,4,2,4), mfrow = c(3,1))


CLOTRIMAZOLE_limma_results <- read.csv("~/BF528/Project 3/CLOTRIMAZOLE_limma_results.csv")
significant_de_clot <- CLOTRIMAZOLE_limma_results[CLOTRIMAZOLE_limma_results$adj.P.Val < 0.05,]
hist(significant_de_clot$logFC, main = "A. Histogram of significant DEGs CAR/PXR", xlab = "number of DEGs", breaks=100, col='blue')

CHLOROFORM_limma_results <- read.csv("~/BF528/Project 3/CHLOROFORM_limma_results.csv", stringsAsFactors=FALSE)
significant_de_chl <- CHLOROFORM_limma_results[CHLOROFORM_limma_results$adj.P.Val < 0.05,]
hist(significant_de_chl$logFC, main = "B. Histogram of significant DEGs Cytotoxic", xlab = "number of DEGs", breaks=100, col='blue')

METHYLCHOLANTHRENE_limma_results <- read.csv("~/BF528/Project 3/3-METHYLCHOLANTHRENE_limma_results.csv", stringsAsFactors=FALSE)
significant_de_3ME <- METHYLCHOLANTHRENE_limma_results[METHYLCHOLANTHRENE_limma_results$adj.P.Val < 0.05,]
hist(significant_de_3ME$logFC, main = "C. Histogram of significant DEGs Ahr", xlab = "number of DEGs", breaks=100, col='blue')

par(mar=c(4,4,2,4), mfrow = c(3,1))

library(EnhancedVolcano)
row.names(CLOTRIMAZOLE_limma_results) <- CLOTRIMAZOLE_limma_results$X
CLOTRIMAZOLE_limma_results <- subset(CLOTRIMAZOLE_limma_results,select= c("logFC","P.Value","adj.P.Val"))


EnhancedVolcano(CLOTRIMAZOLE_limma_results,
                lab = rownames(CLOTRIMAZOLE_limma_results),
                x = 'logFC',
                y = 'adj.P.Val',
                title = "A. CAR/PXR",
                pointSize=2,
                pCutoff = 0.05,
                xlim = c(-5,8),
                widthConnectors = 0.5,colAlpha = 0.5,labSize = 3.0)

row.names(CHLOROFORM_limma_results) <- CHLOROFORM_limma_results$X
CHLOROFORM_limma_results <- subset(CHLOROFORM_limma_results,select= c("logFC","P.Value","adj.P.Val"))


EnhancedVolcano(CHLOROFORM_limma_results,
                lab = rownames(CHLOROFORM_limma_results),
                x = 'logFC',
                y = 'adj.P.Val',
                title = "B. Cytotoxic",
                pointSize=2,
                pCutoff = 0.05,
                xlim = c(-5,8),
                widthConnectors = 0.5,colAlpha = 0.5,labSize = 3.0)

row.names(METHYLCHOLANTHRENE_limma_results) <- METHYLCHOLANTHRENE_limma_results$X
METHYLCHOLANTHRENE_limma_results <- subset(METHYLCHOLANTHRENE_limma_results,select= c("logFC","P.Value","adj.P.Val"))


EnhancedVolcano(METHYLCHOLANTHRENE_limma_results,
                lab = rownames(METHYLCHOLANTHRENE_limma_results),
                x = 'logFC',
                y = 'adj.P.Val',
                title = "C. Ahr",
                pointSize=2,
                pCutoff = 0.05,
                xlim = c(-5,8),
                widthConnectors = 0.5,colAlpha = 0.5,labSize = 3.0)

