library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(preprocessCore)
library(tidyverse)
library(beepr)

# Input: Kallisto transcript read counts files for each tissue (fetal, adult, and cancer)

# 1) Create Sample Table From Data Columns

outputDir <- "~/Developmental-Origins-of-Cancer"
rawDataDir <- paste(outputDir, "data", sep="/")
setwd(outputDir)

tissue = "Kidney"
cancerType <- "KIRC"

analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
dir.create(analyzedDataDir)

tx2gene_file <- paste(outputDir,"mart_export.txt", sep="/") # use new reference annotation file
tx2gene <- read_tsv(tx2gene_file) %>% dplyr::select(c("Transcript_ID_version", "Gene_ID"))



countsFile <- paste(analyzedDataDir, paste(tissue, "Counts.tsv", sep=""), sep="/")
samplesFile <- paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/")

if(!file.exists(countsFile)){

  cancerStageFile <- paste(paste(rawDataDir, "Cancer", "TumorStages", sep="/"), cancerType, sep="")    #Download from TCGA
  cancerStages <- read_tsv(file=cancerStageFile)
  cancerStages$samples <- str_replace_all(cancerStages$case_submitter_id, "-", ".")    #match ID name with data columns
  cancerStages <- cancerStages[duplicated(cancerStages),]
  stages <- c("Fetal", "Adult", "Cancer")
  add_stages <- "No"    #Yes only for cancer (for cancer stages)
  sampleTable <- data.frame(matrix(vector(), 0,4))    #create empty vector for sample table
  colnames(sampleTable) <- c('samples', 'tumor_stage', 'tissue', 'stage')
  for (stage in stages){
    subtissue <- tissue
    if (stage == "Cancer"){
      subtissue <- cancerType
      add_stages <- "Yes"
    }
    stageDir <- paste(rawDataDir, stage, sep="/")
    rawFiles <- list.files(stageDir)
    
    tissueFilePath <- paste(stageDir, "KallistoRawCountsForTranscriptExpressionIn", sep="/") %>%
      paste(.,subtissue, sep="")
    countsDf <- read.table(file=tissueFilePath, header=TRUE, row.names=1)
    
    # Take the antilog of adult and cancer counts values
    if(stage=="Adult" | stage=="Cancer"){ 
      countsDf <- (2^countsDf)-1
    }
  
    head(countsDf)
    samples <- colnames(countsDf)

    if(stage == "Cancer"){
      samples <- samples %>% grep('[.]+01', ., value=TRUE)
      countsDf <- countsDf %>% dplyr::select(all_of(samples))
      samples <- samples %>% str_remove("[.]+01")
    }
    
    stageSampleTable = as.data.frame(samples)
    rownames(stageSampleTable) <- samples
    
    if (stage == "Fetal"){
      stageSampleTable$days <- sapply(strsplit(stageSampleTable$samples, "_"), `[`, 2) %>% strsplit("days")
      stageSampleTable$tumor_stage <- "NA"
    }
  
    if (stage == "Adult"){
      stageSampleTable$days <- "NA"
      stageSampleTable$tumor_stage <- "NA"
    }
    
    if (stage == "Cancer"){
      stageSampleTable$days <- "NA"
      stageSampleTable <- left_join(stageSampleTable, cancerStages, by="samples") %>% dplyr::select(-c(case_submitter_id))
      rownames(stageSampleTable) <- stageSampleTable$samples
    }
    
    
    stageSampleTable$tissue <- tissue
    stageSampleTable$stage <- stage
    sampleTable <- rbind(sampleTable, stageSampleTable)
  
    
    if (stage=="Fetal"){
      fetalCounts <- countsDf
    }
    if (stage=="Adult"){
      adultCounts <- countsDf
    }
    if (stage=="Cancer"){
      cancerCounts <- countsDf
    }
    
  }


  # 2) Merge Fetal, Adult, and Cancer Data
  tempCountsDf <- merge(fetalCounts, adultCounts, by=0, all=FALSE)
  rownames(tempCountsDf) <- tempCountsDf$Row.names
  
  allCountsDf <- merge(tempCountsDf, cancerCounts, by=0, all=FALSE) 
  rownames(allCountsDf) <- allCountsDf$Row.names
  allCountsDf <- allCountsDf %>% dplyr::select(-c(Row.names))
  # rownames(allCountsDf) <- sapply(strsplit(x=rownames(allCountsDf), split=".", fixed = TRUE), `[`, 1)
  allCountsDf$TXNAME <- rownames(allCountsDf)
  joinTxGenes <- inner_join(tx2gene, allCountsDf, by=c("Transcript_ID_version" = "TXNAME")) %>% dplyr::select(-c(Transcript_ID_version))
  allSamples <- colnames(joinTxGenes %>% dplyr::select(-c(Gene_ID)))
  groupGenes <- aggregate(. ~ Gene_ID, data=joinTxGenes, FUN=sum)
  rownames(groupGenes) <- groupGenes$Gene_ID
  countDataDf <- groupGenes %>% dplyr::select(-c(Gene_ID)) %>% round
  sampleTable <- apply(sampleTable, 2, as.character) %>% as.data.frame
  sampleTable$source <- sampleTable$stage
  sampleTablePaper <- sampleTable[sampleTable$source == "Fetal",]
  sampleTablePaper$days <- sampleTablePaper$days %>% as.numeric
  
  sampleTablePaper <- sampleTablePaper %>% 
    mutate(
      stage = case_when(
        (days<270) ~ "Fetal",
        (days >= 270 & days < 6840) ~ "Child",
        (days>=6840) ~ "Adult"
      ) 
    )
  
  sampleTableNotPaper <- sampleTable[!sampleTable$source == "Fetal",]
  finalSampleTable <- rbind(sampleTablePaper, sampleTableNotPaper)

  
 
  # 3) Save count dataframe and samples dataframe
  write.table(countDataDf, file=countsFile, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)
  write.table(finalSampleTable, file=samplesFile, quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

}
beep()




