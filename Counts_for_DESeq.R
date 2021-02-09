library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(tximport)
library(DESeq2)
library(preprocessCore)
library(tidyverse)
library(DESeq)


rawDataDir <- "~/DevelopmentalGenesAnalysis/SalmonOutput"

tx2gene_file <- paste(rawDataDir,"tx2gene.tsv", sep="/")

if(!file.exists(tx2gene_file)) {
  
  # Build tx2gene conversion table from salmon index file and run tximport
  
  txdb <- makeTxDbFromGFF(file="~/salmon_partial_sa_index/hg38.gtf") # remove pseudogenes!
  k <- keys(txdb, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
  head(tx2gene)
  
  # Remove pseudogenes
  pseudogene_table <- read.table(file.path("~/Liver_Extracted_Reads/hg38_transcript_biotype.tsv.gz"), header = TRUE)
  pseudogenes <- pseudogene_table[str_detect(pseudogene_table$transcript_biotype, "unprocessed_pseudogene"), ]
  names(pseudogenes)[names(pseudogenes) == "ensembl_transcript_id"] <- "TXNAME"
  tx2gene_rm <- anti_join(tx2gene, pseudogenes, by= c("TXNAME"))
  
  write.table(tx2gene, file=paste(rawDataDir,"tx2gene.tsv", sep="/"), quote=F, row.names=F, sep="\t")
  
}

tx2gene <- read_tsv(tx2gene_file)


tissue <- "Kidney"
cancerType <- "KIRC"

tissueDir <- paste(rawDataDir, tissue, sep="/")
tissueDataDir <- paste(rawDataDir, tissue, "RawData", sep="/")

ddsFile <- paste(tissueDir, "dds.RData", sep="/")

if(!file.exists(ddsFile)){
  rawFiles <- list.files(tissueDataDir)
  
  # Extract count files
  
  for (countsFile in rawFiles) {
    newFile <- paste(tissueDataDir,
                     str_remove(countsFile, '.tar') %>% str_remove(., '.gz'),
                     "RNAseqAnalysis", tissue, "Measurements/quant.sf",
                     sep="/")
    
    
    if(!file.exists(newFile)) {
      print(paste("Extracting file:", newFile, sep = " "))
      untar(paste(tissueDataDir, countsFile, sep="/"),
            exdir = paste(tissueDataDir, str_remove(countsFile, '.tar') %>% str_remove(., '.gz'), sep="/"))
    }
    else(print("File extracted"))
    
  }
  
  
  # Create table with stages
  
  samples <- data.frame(filename= character(0), sample= character(0), tissue= character(0), days_pc= numeric(0))
  
  for (filename in rawFiles) {
    
    filename_short <- str_remove(filename, '.tar') %>% 
      str_remove(., '.gz')
    stage <- str_remove(filename, '.tar') %>% 
      str_remove(., '.gz') %>% 
      str_remove(., tissue) %>% 
      str_remove(., 'SAMEA') %>% 
      substring(., 8)
    
    stage_num <- as.numeric(gsub(pattern = "[[:alpha:]]", replacement = "", stage))
    stage_unit <- gsub(pattern = "[[:digit:]]", replacement = "", stage)
    
    if(stage_unit=="ypb") {
      days <- 270 + (365 * stage_num)
    }
    
    if(stage_unit=="dpb") {
      days <- 270 + stage_num
    }
    
    if(stage_unit=="mpb") {
      days <- 270 + (30 * stage_num)
    }
    
    if(stage_unit=="w") {
      days <- 7 * stage_num
    }
    
    if(stage_unit=="wpc") {
      days <- 7 * stage_num
    }
    
    if(stage_unit=="dpc") {
      days <- stage_num
    }
    
    counts_file <- paste(tissueDataDir,
                         filename_short,
                         "RNAseqAnalysis", tissue, "Measurements/quant.sf",
                         sep="/")
    counts_file
    
    sample_row <- c(counts_file, filename_short, tissue, days)
    samples[nrow(samples)+1, ] <- sample_row
    
  }
  
  UniqueSamples <- unique(samples)
  
  all_counts_files <- as.vector(UniqueSamples$filename)
  
  # Keep all samples with quant.sf files
  all_quantFiles <- c()
  for (quantFile in all_counts_files) {
    if(file.exists(quantFile)) {
      all_quantFiles <- c(all_quantFiles, quantFile)
    }
  }
  
  UniqueSamplesQuant <- UniqueSamples[UniqueSamples$filename %in% all_quantFiles ,]
  

  
  
  # Add cancer datasets to samples list
  
  cancer_dir <- paste(tissueDir, paste(cancerType, "_HTSeq", sep=""), sep="/")
  
  sampleTable <- read.table(paste(tissueDir, paste(cancerType, "_SampleTable.tsv", sep=""), sep="/"), header=T, sep="\t")# %>% dplyr::select(FileID, FileName, SampleType)
  sampleTable$filename <- paste(sampleTable$FileID, sampleTable$FileName, sep="/")
  sampleTableCancer <- sampleTable %>% 
    mutate(
      stage = case_when(
        SampleType == "Primary Tumor" ~ "cancer",
        SampleType == "Solid Tissue Normal" ~ "adult",
        SampleType == "Recurrent Tumor" ~ "cancer"
      ) 
    )
  
  
  sampleTableCancer$tissue <- tissue
  sampleTableCancer$sample <- sampleTableCancer$FileID
  sampleTableCancer$days_pc <- -1
  
  AllSamples <- rbind(UniqueSamplesQuant, dplyr::select(sampleTableCancer, c(filename, tissue, sample, days_pc)))
  
  cancerFilePaths <- list.files(cancer_dir, recursive=T) %>% grep(".htseq.counts.gz", ., value=T)
  geneCol <- read_tsv(paste(cancer_dir, cancerFilePaths[1], sep="/"),
                      col_names=c("gene", "NumReads"))
  geneCol$geneID <- gsub(pattern="\\..*",
                         replacement="",
                         x=geneCol$gene)
  
  cancerCounts <- geneCol %>% dplyr::select(geneID)
  
  for (cancerFile in cancerFilePaths){
    filePath <- paste(cancer_dir, cancerFile, sep="/")
    cancerSample <- gsub( "/.*$", "", cancerFile)
    sampleCounts <- read_tsv(filePath,
                             col_names=c("gene", cancerSample))
    sampleCounts$geneID <- gsub(pattern="\\..*",
                                replacement="",
                                x=sampleCounts$gene)
    sampleDf <- sampleCounts %>% dplyr::select(-gene)
    cancerCounts <- inner_join(cancerCounts, sampleDf, by="geneID")
    
    
  }
  
  cancerCountsDf <- as.data.frame(cancerCounts)
  row.names(cancerCountsDf) <- cancerCountsDf$geneID 
  cancer_counts_df <-cancerCountsDf %>% dplyr::select(-geneID)
  
  
  
  AllSamples$days_pc <- as.numeric(AllSamples$days_pc)
  
  sample_table <- AllSamples %>% 
    mutate(
      stage = case_when(
        days_pc == -1 ~ "cancer",
        days_pc < 270 & days_pc >= 0 ~ "fetal",
        days_pc >= 270 & days_pc < 6840 ~ "child",
        days_pc >= 6840 ~ "adult"
      ))
  sample_table
  rownames(sample_table) <- sample_table$sample
  
  write.table(sample_table, file=paste(tissueDataDir,"Sample_Files_Table", sep="/"), quote=F, row.names=F, sep="\t")
  
  
  # Tximport counts data
  
  
  
  # Run DESeq to get log2FC
  
  txi <- tximport(all_quantFiles, type="salmon", tx2gene=tx2gene, countsFromAbundance = "scaledTPM", ignoreTxVersion = TRUE)
  counts_df <- txi$counts
  
  sample_list <- sapply(strsplit(all_quantFiles, "/", fixed = TRUE), "[", 6)
  colnames(counts_df) <- sample_list
  
  SalmonCountsMatrix <- data.matrix(counts_df) %>% round
  HTSeqCountsMatrix <- data.matrix(cancer_counts_df) %>% round
  
  commonGenes <- intersect(rownames(SalmonCountsMatrix), rownames(HTSeqCountsMatrix))
  SalmonCountsMatrix_sub <- SalmonCountsMatrix[commonGenes, ]
  HTSeqCountsMatrix_sub <- HTSeqCountsMatrix[commonGenes, ]
  
  AllCountsMatrix <- cbind(SalmonCountsMatrix_sub, HTSeqCountsMatrix_sub)
  
  SamplesOrder <- rownames(sample_table)
  AllCountsMatrixOrdered <- AllCountsMatrix[, SamplesOrder]
  
  allRPKM <- read.table(file="~/Enhancer_Data/Human.RPKM.txt", header=TRUE, sep=" ", row.names=1)
  tissueRPKM <- allRPKM %>% dplyr::select(starts_with(tissue))
  
  dds <- DESeqDataSetFromMatrix(countData = AllCountsMatrixOrdered,
                                colData = sample_table,
                                design = ~stage)
  dds <- DESeq(dds)
  
  saveRDS(dds, file=paste(tissueDir, "dds.RData", sep="/"))
  
}
  

dds <- readRDS(file=ddsFile)
sample_table <- read_tsv(file=paste(tissueDataDir,"Sample_Files_Table", sep="/"))
View(sample_table)

res_FA <- results(dds, name="stage_fetal_vs_adult")
res_CA <- results(dds, name="stage_cancer_vs_adult")
# fold change log2(fetal/adult)
# fold change log2(cancer/adult)
write.table(res_FA, file = paste(tissueDir, "log2FoldChange_FetalAdult.csv", sep="/"), sep = "\t", row.names = T, quote=F)
write.table(res_CA, file = paste(tissueDir, "log2FoldChange_CancerAdult.csv", sep="/"), sep = "\t", row.names = T, quote=F)

res_FA_df <- na.omit(as.data.frame(res_FA))
res_FA_df_significant <- res_FA_df[res_FA_df[, "padj"] < 0.05,]
nrow(res_FA)
nrow(res_FA_df_significant)

res_CA_df <- na.omit(as.data.frame(res_CA))
res_CA_df_significant <- res_CA_df[res_CA_df[, "padj"] < 0.05,]
nrow(res_CA)
nrow(res_CA_df_significant)

res_FA_df_significant$groupFA <- ifelse(res_FA_df_significant$log2FoldChange > 2, "DOWN", "UP")
res_CA_df_significant$groupCA <- ifelse(res_CA_df_significant$log2FoldChange < -2, "DOWN", "UP")

res_FA_filtered <- res_FA_df_significant[ (res_FA_df_significant$log2FoldChange >2) | (res_FA_df_significant$log2FoldChange < -2),]
res_CA_filtered <- res_CA_df_significant[ (res_CA_df_significant$log2FoldChange >2) | (res_CA_df_significant$log2FoldChange < -2),]




FC_groups <- merge(dplyr::select(res_FA_filtered, groupFA), dplyr::select(res_CA_filtered, groupCA), by.x = 0, by.y = 0)
FC_groups$group <- paste(FC_groups$groupFA, FC_groups$groupCA, sep="_")
print(table(FC_groups$group))

write.table(FC_groups, file = paste(tissueDir, "GeneGroups.tsv", sep="/"), sep = "\t", row.names = F, quote=F)


# vst <- vst(dds, blind=TRUE)
# nVariable <- 1000
# pcaPlot <- DESeq2::plotPCA(vst, intgroup='stage', ntop=nVariable)
# pcaPlot + ggtitle(paste(tissue, "PCA, Top", nVariable, "Genes", sep=" "))

sample_table
samplesDf <- as.data.frame(sample_table)
rownames(samplesDf) <- samplesDf$sample

fetalSamples <- rownames(samplesDf[samplesDf$stage == 'fetal',])
adultSamples <- rownames(samplesDf[samplesDf$stage == 'adult',])

# Get DESeq counts df

dds <- estimateSizeFactors(dds)
ddsCounts <- counts(dds, normalized=TRUE)
ddsCountsFetal <- ddsCounts[, fetalSamples]
ddsCountsAdult <- ddsCounts[, adultSamples]

write.table(ddsCountsAdult, file = paste(tissueDir, "DENormalizedAdultCounts.tsv", sep="/"), sep = "\t", row.names = T, quote=F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





UP_CA <- res_CA_filtered[res_CA_filtered$groupCA == 'UP',] %>% rownames
DOWN_CA <- res_CA_filtered[res_CA_filtered$groupCA == 'DOWN',] %>% rownames

devUP <- ddsCountsFetal[UP_CA,]
devDOWN <- ddsCountsFetal[DOWN_CA,]

# get stage max for all UP genes

maxStage <- colnames(devUP)[apply(devUP, 1, which.max)] %>% as.data.frame
colnames(maxStage) <- 'maxStage'
maxStage$age <- sub(paste(".*", tissue, sep=""), "", maxStage$maxStage)
rownames(maxStage) <- rownames(devUP)


# get stage min for all DOWN genes

minStage <- colnames(devDOWN)[apply(devDOWN, 1, which.min)] %>% as.data.frame
colnames(minStage) <- 'minStage'
minStage$age <- sub(paste(".*", tissue, sep=""), "", minStage$minStage)
rownames(minStage) <- rownames(devDOWN)


# get stage max and min for all genes

maxStageAll <- colnames(ddsCountsFetal)[apply(ddsCountsFetal, 1, which.max)] %>% as.data.frame
colnames(maxStageAll) <- 'maxStageAll'
maxStageAll$age <- sub(paste(".*", tissue, sep=""), "", maxStageAll$maxStageAll)
rownames(maxStageAll) <- rownames(ddsCountsFetal)

minStageAll <- colnames(ddsCountsFetal)[apply(ddsCountsFetal, 1, which.min)] %>% as.data.frame
colnames(minStageAll) <- 'minStageAll'
minStageAll$age <- sub(paste(".*", tissue, sep=""), "", minStageAll$minStageAll)
rownames(minStageAll) <- rownames(ddsCountsFetal)



maxUPDf <- as.data.frame(sort(table(maxStage$age)), decreasing=FALSE)
colnames(maxUPDf) <- c("stage", "maxUP")
minDOWNDf <- as.data.frame(sort(table(minStage$age)), decreasing=FALSE)
colnames(minDOWNDf) <- c("stage", "minDOWN")
maxAllDf <- as.data.frame(sort(table(maxStageAll$age)), decreasing=FALSE)
colnames(maxAllDf) <- c("stage", "maxAll")
minAllDf <- as.data.frame(sort(table(minStageAll$age)), decreasing=FALSE)
colnames(minAllDf) <- c("stage", "minAll")

maxStagesDf <- inner_join(maxUPDf, maxAllDf, by='stage')
maxStagesDf$normMaxStage <- maxStagesDf$maxUP / maxStagesDf$maxAll

minStagesDf <- inner_join(minDOWNDf, minAllDf, by='stage')
minStagesDf$normMinStage <- minStagesDf$minDOWN / minStagesDf$minAll

stageOrder <- c('4w', 'CS13', '5w', 'CS17', 'CS18', 'CS19', '8w', 'CS20', 'CS21', 'CS22', 'CS23', '9w', '10w', '11w', '12w', '13w', '16w', '18w', '19w', '20w') 
maxStagesDf$stage <- factor(maxStagesDf$stage,levels = stageOrder)
minStagesDf$stage <- factor(minStagesDf$stage,levels = stageOrder)


UP <- ggplot(maxStagesDf, aes(x=stage, y=normMaxStage))
UP + geom_col() + ggtitle(paste("UP Genes Adult to Cancer", tissue, sep=" ")) + xlab("Developmental Stage Max Expression") +
  ylab("Number of Genes Normalized") + theme_classic()

DOWN <- ggplot(minStagesDf, aes(x=stage, y=normMinStage))
DOWN + geom_col() + ggtitle(paste("DOWN Genes Adult to Cancer", tissue, sep=" ")) + xlab("Developmental Stage Min Expression") +
  ylab("Number of Genes Normalized") + theme_classic()




