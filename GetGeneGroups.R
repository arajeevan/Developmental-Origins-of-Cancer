library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(preprocessCore)
library(ggplot2)


outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
dir.create(plotsDir)
setwd(outputDir)

tissues <- c("Liver", "Kidney", "Brain", "Ovary", "Testis")

for (tissue in tissues){
  logFCTable <- read.table(file=paste(analyzedDataDir, paste(tissue,"LogFC", sep=""), sep="/"), header=TRUE, sep="\t", row.names=1)
  sampleTable <- read.table(file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
                            header=TRUE, sep="\t", row.names=1)
  
  logFCDf <- logFCTable %>% na.omit
  logFCDf_significant <- logFCDf[logFCDf[, "padj_FA"] < 0.05,] %>% .[.[, "padj_CA"] < 0.05,]
  
  logFCDf_significant$groupFA <- ifelse(logFCDf_significant$logFC_FA > 2, "DOWN", "UP")
  logFCDf_significant$groupCA <- ifelse(logFCDf_significant$logFC_CA < -2, "DOWN", "UP")
  
  logFCDf_A <- logFCDf_significant[ (logFCDf_significant$logFC_FA > 2) | (logFCDf_significant$logFC_FA < -2),]
  logFCDf_B <- logFCDf_A[ (logFCDf_A$logFC_CA > 2) | (logFCDf_A$logFC_CA < -2),]
  
  logFCDf_B$group <- paste(logFCDf_B$groupFA, logFCDf_B$groupCA, sep="")
  print(table(logFCDf_B$group))
  
  write.table(logFCDf_B,
              file = paste(analyzedDataDir, paste(tissue, "GeneGroups.tsv", sep=""), sep="/"),
              sep = "\t", row.names = T, quote=F)
  
  countsDf <- read.table(file = paste(analyzedDataDir, (paste(tissue, "NormalizedCounts", sep="")), sep="/"),
  header=TRUE, sep="\t")
  
  fetalSamples <- sampleTable[ (sampleTable$stage == "Fetal") ,] %>% .$samples
  adultSamples <- sampleTable[ (sampleTable$stage == "Adult") ,] %>% .$samples
  cancerSamples <- sampleTable[ (sampleTable$stage == "Cancer") ,] %>% .$samples
  
  fetalExp <- countsDf[ , fetalSamples]
  fetalExp$geneID <- rownames(fetalExp)
  adultExp <- countsDf[ , adultSamples]
  adultExp$geneID <- rownames(adultExp)
  cancerExp <- countsDf[ , cancerSamples]
  cancerExp$geneID <- rownames(cancerExp)
  
  groups <- logFCDf_B$group %>% unique
  
  for (groupName in groups){
    groupGenes <- logFCDf_B[ (logFCDf_B$group == groupName), ] %>% rownames
    fetalGroupExp <- fetalExp[ groupGenes, ] %>%
      dplyr::select(-c(geneID)) %>% unlist %>% as.data.frame
    fetalGroupExp$stage <- "Fetal"
    rownames(fetalGroupExp) <- NULL
    adultGroupExp <- adultExp[ groupGenes, ] %>%
      dplyr::select(-c(geneID)) %>% unlist %>% as.data.frame
    adultGroupExp$stage <- "Adult"
    rownames(adultGroupExp) <- NULL
    cancerGroupExp <- cancerExp[ groupGenes, ] %>%
      dplyr::select(-c(geneID)) %>% unlist %>% as.data.frame
    cancerGroupExp$stage <- "Cancer"
    rownames(cancerGroupExp) <- NULL
    
    groupExpDf <- rbind(fetalGroupExp, adultGroupExp, cancerGroupExp)
    colnames(groupExpDf) <- c("expression", "stage")
    
    groupExpDf$stage <- factor(groupExpDf$stage , levels=c("Fetal", "Adult", "Cancer"))
    
    ggplot(data=groupExpDf, aes(x=stage, y=expression)) +
      geom_boxplot() +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = c(0,2500)) +
      ggtitle(paste(tissue, groupName,"Gene Expression",sep=" "))
    ggsave(paste(plotsDir, paste(tissue, groupName,"GeneExpressionCheck.png", sep=""), sep="/"))
  
  }
  
}

