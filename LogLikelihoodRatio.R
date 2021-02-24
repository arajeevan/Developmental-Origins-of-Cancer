library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(preprocessCore)
library(data.table)


outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
setwd(outputDir)

logLikDf <- data.frame()
tissues <- c("Liver", "Kidney", "Brain", "Ovary", "Testis")

# Getting warning for Brain and Ovary: glm.fit: fitted probabilities numerically 0 or 1 occurred 

for (tissue in tissues){

  sampleTable <- read.table(file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
                            header=TRUE, sep="\t", row.names=1)
  adultSamples <- sampleTable[ (sampleTable$stage == "Adult") ,] %>% .$samples
  logFCDf <- read.table(file=paste(analyzedDataDir, paste(tissue, "GeneGroups.tsv", sep=""), sep="/"), header=TRUE, sep="\t")
  dds <- readRDS(file=paste(analyzedDataDir, paste(tissue, "dds.RData", sep=""), sep="/") )
  dds <- estimateSizeFactors(dds)
  ddsCounts <- counts(dds, normalized=TRUE)
  colnames(ddsCounts) <- dds$samples
  normAdultCounts <- ddsCounts[, adultSamples] 
  meanAdultCol <- rowMeans(normAdultCounts)
  
  foldChanges <- merge(logFCDf, meanAdultCol, by.x=0, by.y=0)
  rownames(foldChanges) <- foldChanges$Row.names
  names(foldChanges)[names(foldChanges) == "y"] <- "meanAdultCount"
  
  foldChanges$FA <- ifelse(foldChanges$groupFA == "UP", 1, 0)
  foldChanges$CA <- ifelse(foldChanges$groupCA == "UP", 1, 0)
  
  # Binomial
  nested <- glm(CA ~ meanAdultCount,
                data=foldChanges, family='binomial')
  baseline <- glm(CA ~ FA,
                  data=foldChanges, family='binomial')
  complex <- glm(CA ~ meanAdultCount + FA,
                 data=foldChanges, family='binomial')
  
  (A <- logLik(nested))
  (B <- logLik(complex))
  (C <- logLik(baseline))
  (teststat <- -2 * (as.numeric(A)-as.numeric(B)))
  (p.val <- pchisq(teststat, df = 1, lower.tail = FALSE))
  
  logLikRow <- c(tissue, teststat, p.val)
  logLikDf <- rbind(logLikRow, logLikDf)
  
}


colnames(logLikDf) <- c("tissue", "testStat", "pVal")
write.table(logLikDf,
            file = paste(analyzedDataDir, paste("AllTissuesLogLikRatio.tsv", sep=""), sep="/"),
            sep = "\t", row.names = T, quote=F)


