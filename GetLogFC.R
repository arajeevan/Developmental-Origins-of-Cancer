library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(DESeq2)
library(preprocessCore)


# Load counts and sample tables from CombineKallistoSamples.R pipeline



outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
dir.create(plotsDir)
setwd(outputDir)

tissues <- c("Liver", "Kidney", "Brain", "Ovary", "Testis")

for (tissue in tissues){
  tissue <- "Testis"
  countsDf <- read.table(file=paste(analyzedDataDir, paste(tissue, "Counts.tsv", sep=""), sep="/"),
                         header=TRUE, sep="\t", row.names=1) 
  sampleTable <- read.table(file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
                            header=TRUE, sep="\t", row.names=1)
  
  if(tissue=="Testis"){
    countsDf <- read.table(file=paste(analyzedDataDir, paste(tissue, "Counts.tsv", sep=""), sep="/"),
                           header=TRUE, sep="\t", row.names=1) %>% dplyr::select(-c("Testis_84days_Rep2"))
    sampleTable <- sampleTable[!(sampleTable$samples == "Testis_84days_Rep2"),]
    row.names(sampleTable) <- NULL
  }
  
  sampleTable <- sampleTable %>% 
    mutate(
      source = case_when(
        (source == "Fetal") ~ "paper",
        (source == "Adult") ~ "GTex",
        (source == "Cancer") ~ "TCGA"
      ) 
    )
  
  
  
  # Run DESeq on all samples
  
  ddsFile <- paste(analyzedDataDir, paste(tissue, "dds.RData", sep=""), sep="/")
  
  if(!file.exists(ddsFile)){
    dds <- DESeqDataSetFromMatrix(countData = countsDf, colData=sampleTable, design = ~stage)
    dds <- DESeq(dds)
    beep()
    saveRDS(dds, file=ddsFile)
  }
  
  dds <- readRDS(file=ddsFile)
  dds <- estimateSizeFactors(dds)
  ddsCounts <- counts(dds, normalized=TRUE)
  colnames(ddsCounts) <- dds$samples
  write.table(ddsCounts,
              file = paste(analyzedDataDir, (paste(tissue, "NormalizedCounts", sep="")), sep="/"),
              sep = "\t", row.names = T, quote=F)
  
  
  res_FA <- results(dds, name = "stage_Fetal_vs_Adult")
  res_CA <- results(dds, name = "stage_Cancer_vs_Adult")
  
  res_All <- cbind(res_FA$log2FoldChange, res_FA$padj, res_CA$log2FoldChange, res_CA$padj)
  rownames(res_All) <- rownames(res_FA)
  colnames(res_All) <- c("logFC_FA", "padj_FA", "logFC_CA", "padj_CA")
  
  write.table(res_All, file=paste(analyzedDataDir, paste(tissue,"LogFC", sep=""), sep="/"),
              quote=F, sep="\t")
  write.table(sampleTable, file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
              quote=F, sep="\t")
}


