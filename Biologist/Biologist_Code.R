setwd("~/Desktop/BF528/")


norm_count_ahr <- read.csv("deseq_norm_counts_ahr.csv", header = TRUE)
norm_count_car_pxr <- read.csv("deseq_CAR_PXR_norm_counts.csv", header = TRUE)
norm_count_cytotoxic <- read.csv("deseq_norm_counts_cytotoxic.csv", header = TRUE)

combined <- as.data.frame(c(norm_count_ahr,norm_count_car_pxr,norm_count_cytotoxic))


#Reading in the 100 top expressed genes for each MOA

top100_ahr <- read.csv("top100_ahr.csv", header = TRUE)
top100_car_pxr <- read.csv("top100_car_pxr.csv", header = TRUE)
top100_cytotoxin <- read.csv("top100_cytotoxic.csv", header= TRUE)

#Renaming the first column to be ID to be able to merge 
colnames(norm_count_ahr)[1] <- "ID"
colnames(norm_count_car_pxr )[1] <- "ID"
colnames(norm_count_cytotoxic )[1] <- "ID"

#subsetting them
top100_ahr_ID <- as.list(levels(top100_ahr[,1]))
top100_car_pxr_ID <- as.list(levels(top100_car_pxr[,1]))
top100_cytotoxin_ID <- as.list(levels(top100_cytotoxin[,1]))

#merging all 3 ID lists into 1
x <- append(top100_ahr_ID,top100_car_pxr_ID)
x <- append(x,top100_cytotoxin_ID)

merged_MOA<- merge(norm_count_ahr[,1:4],norm_count_car_pxr[,1:4], by = "ID") 
merged_MOA <- merge(merged_MOA,norm_count_cytotoxic[,1:4], by = "ID") 

rownames(merged_MOA) <- merged_MOA[, 1]
merged_MOA <- as.matrix(merged_MOA[,-1])

#subsetting the counts of only the above IDs
y <- merged_MOA[row.names(merged_MOA) %in% x,]





#heatmap
heatmap(as.matrix(y))



  heatmap_data <- data.frame("3_METHYLCHOLANTHRENE_1" = norm_count_ahr$SRR1177997,
                   "3_METHYLCHOLANTHRENE_2" = norm_count_ahr$SRR1177999,
                   "3_METHYLCHOLANTHRENE_3" = norm_count_ahr$SRR1178002,
                   "CAR_PXR_1" = norm_count_car_pxr$SRR1177997,
                   "CAR_PXR_2" = norm_count_car_pxr$SRR1177999,
                      "CAR_PXR_3" = norm_count_car_pxr$SRR1178002,
                    "Cytotoxic_1" = norm_count_cytotoxic$SRR1177987,
                   "Cytotoxic_2" = norm_count_cytotoxic$SRR1177988,
                   "Cytotoxic_3" = norm_count_cytotoxic$SRR1177989
)

heatmap(t(heatmap_data), Rowv = NA, scale = 'col')









