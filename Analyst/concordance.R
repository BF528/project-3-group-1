CLOTRIMAZOLE_limma_results <- read.csv("~/BF528/Project 3/CLOTRIMAZOLE_limma_results.csv", stringsAsFactors = F)
deseq_cytotoxic_results <- read.csv("~/BF528/Project 3/deseq_results_cytotoxic.csv", stringsAsFactors = F)
mapping <- read.csv("~/BF528/Project 3/refseq_affy_map.csv", stringsAsFactors = F)

significant_de_mic <- CLOTRIMAZOLE_limma_results[CLOTRIMAZOLE_limma_results$P.Value < 0.05 ,]
significant_de_RNA <- deseq_cytotoxic_results[deseq_cytotoxic_results$pvalue < 0.05 ,]
mapping <- mapping[!duplicated(mapping),]

names(significant_de_mic)[names(significant_de_mic) == "X"] <- "PROBEID"
names(significant_de_RNA)[names(significant_de_RNA) == "X"] <- "REFSEQ"

sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.clotrimazole = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))
RNA.DEG.clotrimazole = nrow(significant_de_RNA)
MIC.DEG.clotrimazole = nrow(significant_de_mic)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ABOVE- MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr > median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean > median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.abovemedian.clotrimazole = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))

#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
#BELOW MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr < median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean < median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.belowmedian.clotrimazole = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))

#----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

METHYLCHOLANTHRENE_limma_results <- read.csv("~/BF528/Project 3/3-METHYLCHOLANTHRENE_limma_results.csv", stringsAsFactors = F)
deseq_ahr_results <- read.csv("~/BF528/Project 3/deseq_results_ahr.csv", stringsAsFactors = F)

significant_de_mic <- METHYLCHOLANTHRENE_limma_results[METHYLCHOLANTHRENE_limma_results$P.Value < 0.05,]
significant_de_RNA <- deseq_ahr_results[deseq_ahr_results$pvalue < 0.05,]

names(significant_de_mic)[names(significant_de_mic) == "X"] <- "PROBEID"
names(significant_de_RNA)[names(significant_de_RNA) == "X"] <- "REFSEQ"

sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.METHYLCHOLANTHRENE = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))
RNA.DEG.METHYLCHOLANTHRENE = nrow(significant_de_RNA)
MIC.DEG.METHYLCHOLANTHRENE = nrow(significant_de_mic)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ABOVE- MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr > median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean > median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.abovemedian.methylcholanthrene = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))

#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
#BELOW MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr < median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean < median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.belowmedian.methylcholanthrene = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))




#------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------

CHLOROFORM_limma_results <- read.csv("~/BF528/Project 3/CHLOROFORM_limma_results.csv", stringsAsFactors = F)
deseq_car_results <- read.csv("~/BF528/Project 3/deseq_CAR_PXR_results.csv", stringsAsFactors = F)

significant_de_mic <- CHLOROFORM_limma_results[CHLOROFORM_limma_results$P.Value < 0.05,]
significant_de_RNA <- deseq_car_results[deseq_car_results$pvalue < 0.05,]

names(significant_de_mic)[names(significant_de_mic) == "X"] <- "PROBEID"
names(significant_de_RNA)[names(significant_de_RNA) == "X"] <- "REFSEQ"

sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.CHLOROFORM = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))
RNA.DEG.CHLOROFORM = nrow(significant_de_RNA)
MIC.DEG.CHLOROFORM = nrow(significant_de_mic)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ABOVE- MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr > median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean > median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.abovemedian.chloroform = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))

#--------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------
#BELOW MEDIAN

significant_de_RNA <- na.omit(significant_de_RNA)

significant_de_mic <- significant_de_mic[significant_de_mic$AveExpr < median(significant_de_mic$AveExpr) ,]
significant_de_RNA <- significant_de_RNA[significant_de_RNA$baseMean < median(significant_de_RNA$baseMean) ,]
sig_de_rna_sym <- merge(significant_de_RNA, mapping, on = "REFSEQ")

merged <- merge(sig_de_rna_sym , significant_de_mic , on = c("PROBEID","REFSEQ"))

distinct.merged <- merged[!duplicated(merged$SYMBOL),]
distinct.merged <- subset(distinct.merged, select = c("REFSEQ", "PROBEID", "SYMBOL", "log2FoldChange","logFC"))

mic.mapping<-mapping[,2:3]
rna.mapping <- subset(mapping, select = c(1,3))

significant_de_mic <- merge(significant_de_mic, mic.mapping, on="PROBEID")
significant_de_RNA <- merge(significant_de_RNA, rna.mapping, on="REFSEQ")

significant_de_mic <- significant_de_mic[!duplicated(significant_de_mic$SYMBOL),] 
significant_de_RNA <- significant_de_RNA[!duplicated(significant_de_RNA$SYMBOL),] 
significant_de_RNA <- na.omit(significant_de_RNA)


concordance.belowmedian.chloroform = (2 * nrow(distinct.merged)) / (nrow(significant_de_mic) + nrow(significant_de_RNA))




#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
bar.rna <- as.data.frame(rbind(c(concordance.clotrimazole,RNA.DEG.clotrimazole), c(concordance.CHLOROFORM,RNA.DEG.CHLOROFORM), c(concordance.METHYLCHOLANTHRENE,RNA.DEG.METHYLCHOLANTHRENE)))
bar.mic <- as.data.frame(rbind(c(concordance.clotrimazole,MIC.DEG.clotrimazole), c(concordance.CHLOROFORM,MIC.DEG.CHLOROFORM), c(concordance.METHYLCHOLANTHRENE,MIC.DEG.METHYLCHOLANTHRENE)))

par(mfrow=c(1,2))
fit <- lm(V1~V2, data=bar.rna)
plot(bar.rna$V2,bar.rna$V1,pch=19, main ="A. Concordance vs Treatment effect", xlab = "number of DEGs from RNA-seq", ylab = "Concordance", xlim = c(1000,4000), ylim = c(0.15,0.45))
abline(fit, lty = "dotted")
text(V1~V2, labels = c("CLO","CHL","3ME"),data=bar.rna, cex=0.5, font =4 , pos= 4)

fit <- lm(V1~V2, data=bar.mic)
plot(bar.mic$V2,bar.mic$V1,pch=19, main ="B. Concordance vs DEGs of microarray", xlab = "number of DEGs from microarray", ylab = "Concordance", xlim = c(1000,7500), ylim = c(0.15,0.45))
abline(fit, lty = "dotted")
text(V1~V2, labels = c("CLO","CHL","3ME"),data=bar.mic, cex=0.5, font =4 , pos= 4)


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

clotrimazole <- c(concordance.abovemedian.clotrimazole, concordance.clotrimazole, concordance.belowmedian.clotrimazole)
me <- c(concordance.abovemedian.methylcholanthrene, concordance.METHYLCHOLANTHRENE, concordance.belowmedian.methylcholanthrene)
chloroform <- c(concordance.abovemedian.chloroform, concordance.CHLOROFORM, concordance.belowmedian.chloroform)

bar <- as.data.frame(rbind(clotrimazole,me,chloroform))
row.names(bar) <- c("CAR/PXR","Ahr","Cytotoxic")
colnames(bar) <-c("Above Median", "All", "Below Median")

bar <- t(bar)


barplot(bar, main="Concordance for 3 MOAs",
        xlab="MOA", col=c("darkblue","red","white"),
        legend = TRUE, args.legend = list( x = "topright", cex= 0.6), ylim = c(0,0.60),xlim= c(0,15) ,beside=TRUE)


