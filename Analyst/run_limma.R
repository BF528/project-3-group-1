library(limma)

# sample info dataframe with array_id and chemical columns
samples <- read.csv('~/BF528/Project 3/group_1_mic_info.csv',as.is=TRUE)

samples.1 <- samples[samples$chemical == "3-METHYLCHOLANTHRENE",]
samples.2 <- samples[samples$chemical == "CHLOROFORM",]
samples.3 <- samples[samples$chemical == "CLOTRIMAZOLE",]

samples.1 <- rbind(samples.1, samples[(samples$chemical == 'Control' & samples$vehicle == 'CMC_.5_%'),])
samples.2 <- rbind(samples.2, samples[(samples$chemical == 'Control' & samples$vehicle == 'CORN_OIL_100_%'),])
samples.3 <- rbind(samples.3, samples[(samples$chemical == 'Control' & samples$vehicle == 'CORN_OIL_100_%'),])

# the full RMA normalized matrix of all experiments
rma <- read.table('~/BF528/Project 3/liver-normalization-rma.txt',
  sep='\t',
  as.is=TRUE,
  header=TRUE,
  row.names=1,
)

# subset the full expression matrix to just those in this comparison
rma.subset.1 <- rma[paste0('X',samples.1$array_id)]
rma.subset.2 <- rma[paste0('X',samples.2$array_id)]
rma.subset.3 <- rma[paste0('X',samples.3$array_id)]



# construct a design matrix modeling treatment vs control for use by limma
design.1 <- model.matrix(
  ~factor(
    samples.1$chemical,
    levels=c('Control','3-METHYLCHOLANTHRENE')
  )
)
colnames(design.1) <- c('Intercept','3-METHYLCHOLANTHRENE')

# run limma
fit <- lmFit(rma.subset.1, design.1)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.1), adjust='BH')

# write out the results to file
write.csv(t,'~/BF528/Project 3/3-METHYLCHOLANTHRENE_limma_results.csv')
#-------------------------------------------------------------------------------------------
design.2 <- model.matrix(
  ~factor(
    samples.2$chemical,
    levels=c('Control','CHLOROFORM')
  )
)
colnames(design.2) <- c('Intercept','CHLOROFORM')

# run limma
fit <- lmFit(rma.subset.2, design.2)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.2), adjust='BH')

# write out the results to file
write.csv(t,'~/BF528/Project 3/CHLOROFORM_limma_results.csv')


#-------------------------------------------------------------------------------------------
design.3 <- model.matrix(
  ~factor(
    samples.3$chemical,
    levels=c('Control','CLOTRIMAZOLE')
  )
)
colnames(design.3) <- c('Intercept','CLOTRIMAZOLE')

# run limma
fit <- lmFit(rma.subset.3, design.3)
fit <- eBayes(fit)
t <- topTable(fit, coef=2, n=nrow(rma.subset.3), adjust='BH')

# write out the results to file
write.csv(t,'~/BF528/Project 3/CLOTRIMAZOLE_limma_results.csv')
