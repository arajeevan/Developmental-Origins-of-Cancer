library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(DESeq2)
library(preprocessCore)
library(ggplot2)

tissue <- "Testis"

outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
dir.create(plotsDir)
setwd(outputDir)


logFCTable <- read.table(file=paste(analyzedDataDir, paste(tissue,"LogFC", sep=""), sep="/"), header=TRUE, sep="\t", row.names=1)
dds <- readRDS(file=paste(analyzedDataDir, paste(tissue, "dds.RData", sep=""), sep="/") )
sampleTable <- read.table(file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
                          header=TRUE, sep="\t", row.names=1)

# Get top 2000 positive log fc in fetal vs adult, do prcomp with just those genes

nsub_val <- 2000
DEFetalGenes <- logFCTable[order(-logFCTable$logFC_FA), ] %>% .[1:(2*nsub_val),] %>% rownames
ddsDEF <- dds[rownames(dds) %in% DEFetalGenes,]
dim(ddsDEF)

# For vst use top 1000 most variable genes

vst <- vst(ddsDEF, blind=T, nsub=nsub_val)
vstMat <- assay(vst)
pca <- prcomp(t(vstMat))
pca_df <- cbind(sampleTable, pca$x)

# Group cancer samples into early and late stage tumors
pca_df$tumor_stage_grouped <- pca_df$tumor_stage %>% str_remove("stage ")
pca_df$tumor_stage_grouped <- pca_df$tumor_stage_grouped %>% gsub("[abc]", "", .)

pca_df <- pca_df %>% 
  mutate(
    tumor_general_stage = case_when(
      (tumor_stage_grouped == "i" | tumor_stage_grouped == "ii") ~ "early",
      (tumor_stage_grouped == "iii" | tumor_stage_grouped == "iv") ~ "late"
    ) 
  )

# Plot Fetal, Child, Adult, and Cancer Stages PCA Plot (top 1000 DE)
ggplot(pca_df) + geom_point(aes(x=PC1, y=PC2, color = stage)) +
  ggtitle(paste("PCA Gene Expression for ", tissue, sep=""))
plotfilename <- paste(plotsDir, paste(tissue,"DevStagePCA.png", sep=""), sep="/")
ggsave(plotfilename)

if (length(pca_df$tumor_stage %>% unique) > 2){
  
  # Plot Tumor Stage PCA (all stages and substages)
  
  ggplot(pca_df) + geom_point(aes(x=PC1, y=PC2, color = tumor_stage)) +
    ggtitle(paste("PCA Gene Expression for ", tissue, sep=""))
  ggsave(paste(plotsDir, paste(tissue,"TumorStagePCA.png", sep=""), sep="/"))

  # Plot Tumor Stage PCA Grouped (i, ii, iii, iv)
  
  ggplot(pca_df) + geom_point(aes(x=PC1, y=PC2, color = tumor_stage_grouped)) +
    ggtitle(paste("PCA Gene Expression for ", tissue, sep=""))
  ggsave(paste(plotsDir, paste(tissue,"TumorStageGroupedPCA.png", sep=""), sep="/"))
  
  
  # Plot Tumor Stage PCA Grouped (Early/Late)
  
  ggplot(pca_df) + geom_point(aes(x=PC1, y=PC2, color = tumor_general_stage)) +
    ggtitle(paste("PCA Gene Expression for ", tissue, sep=""))
  ggsave(paste(plotsDir, paste(tissue,"TumorGeneralStagePCA.png", sep=""), sep="/"))
  
}






# Fetal sample distances boxplot

top50PCs <- pca_df %>% dplyr::select( head(grep("PC", colnames(pca_df)), 50) )
dist_df <- as.matrix(dist(top50PCs, method="euclidean"))
colnames(dist_df) <- pca_df$stage
rownames(dist_df) <- pca_df$stage

fetalAdult_distances <- dist_df[grepl("Adult*", rownames(dist_df)), ] %>%
  .[ , (grepl("Fetal*", colnames(.)))]

fetalCancer_distances <- dist_df[grepl("Cancer*", rownames(dist_df)), ] %>%
  .[ , (grepl("Fetal*", colnames(.)))]

FAdistDf <- as.vector(fetalAdult_distances) %>% as.data.frame
FAdistDf$samples <- "Fetal Adult"
colnames(FAdistDf) <- c("distance", "samples")
FCdistDf <- as.vector(fetalCancer_distances) %>% as.data.frame
FCdistDf$samples <- "Fetal Cancer"
colnames(FCdistDf) <- c("distance", "samples")

distancesDf <- rbind(FAdistDf, FCdistDf)

give.n <- function(x){
  return(c(y = quantile(x, 0.64)[[1]], label = length(x)))
}

faDists <- distancesDf[distancesDf$samples == "Fetal Adult",] %>% .$distance
fcDists <- distancesDf[distancesDf$samples == "Fetal Cancer",] %>% .$distance

pval_label <- paste("p=", wilcox.test(x=faDists,
                                      y=fcDists,
                                      alternative = "greater") %>% .[[3]], sep="")
print(pval_label)

labelHeight <- distancesDf$distance %>% max
ggplot(distancesDf, aes(x=samples, y=distance, color=samples)) +
  geom_boxplot() +
  stat_summary(fun.data = give.n, geom = "text", size=3) +
  ggtitle(paste("Sample Distances ", tissue, sep="")) +
  ylab("50-PCA Euclidean Distance") +
  xlab("") +
  theme_classic() +
  annotate("text", label=pval_label, x=1.5, y=labelHeight)
ggsave(paste(plotsDir, paste(tissue,"FetalSampleDistancesBoxplot.png", sep=""), sep="/"))


# Tumor stage distance boxplot

if (length(pca_df$tumor_stage %>% unique) > 2){
  
  top50PCs <- pca_df %>% dplyr::select( head(grep("PC", colnames(pca_df)), 50) )
  dist_df <- as.matrix(dist(top50PCs, method="euclidean"))
  colnames(dist_df) <- pca_df$samples
  rownames(dist_df) <- pca_df$samples
  
  sampleNames1 <- subset(pca_df, pca_df[,'tumor_general_stage'] == "early") %>% .[,'samples']
  sampleNames2 <- subset(pca_df, pca_df[,'tumor_general_stage'] == "late") %>% .[,'samples']
  sampleNames3 <- subset(pca_df, pca_df[,'stage'] == "Fetal") %>% .[,'samples']
  
  fetalTo1_distances <- dist_df[rownames(dist_df) %in% sampleNames1,] %>%
    .[, colnames(dist_df) %in% sampleNames3]
  
  fetalTo2_distances <- dist_df[rownames(dist_df) %in% sampleNames2,] %>%
    .[, colnames(dist_df) %in% sampleNames3]
  
  F1distDf <- as.vector(fetalTo1_distances) %>% as.data.frame
  F1distDf$samples <- "Fetal / Early stage"
  colnames(F1distDf) <- c("distance", "samples")
  F2distDf <- as.vector(fetalTo2_distances) %>% as.data.frame
  F2distDf$samples <- "Fetal / Late stage"
  colnames(F2distDf) <- c("distance", "samples")
  
  distancesDf <- rbind(F1distDf, F2distDf)
  
  give.n <- function(x){
    return(c(y = quantile(x, 0.64)[[1]], label = length(x)))
  }
  
  f1Dists <- distancesDf[distancesDf$samples == "Fetal / Early stage",] %>% .$distance
  f2Dists <- distancesDf[distancesDf$samples == "Fetal / Late stage",] %>% .$distance
  
  pval_label <- paste("p=", wilcox.test(x=f1Dists,
                                        y=f2Dists,
                                        alternative = "greater") %>% .[[3]], sep="")
  
  labelHeight <- distancesDf$distance %>% max
  
  ggplot(distancesDf, aes(x=samples, y=distance, color=samples)) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", size=3) +
    ggtitle(paste("Sample Distances for ", tissue, sep="")) +
    ylab("50-PCA Euclidean Distance") +
    xlab("") +
    theme_classic() +
    scale_color_manual(values=c("green", "purple")) +
    annotate("text", label=pval_label, x=1.5, y=labelHeight)
  ggsave(paste(plotsDir, paste(tissue,"SampleDistancesBoxplotEarlyLateStage.png", sep=""), sep="/"))
  
}



