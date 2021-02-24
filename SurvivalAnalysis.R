library(dplyr); library(magrittr); 
library(preprocessCore); library(reshape2);
library(clusterProfiler); library(glue)
library(TCGAbiolinks); library(ggplot2)
library(survival); library(OneR)
library(stringr)

# Set tissue and cancer type
outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
setwd(outputDir)

cancerType <- "LIHC"
tissue <- "Liver"

sampleTable <- read.table(file=paste(analyzedDataDir, paste(tissue, "Samples.tsv", sep=""), sep="/"),
                          header=TRUE, sep="\t", row.names=1)
adultSamples <- sampleTable[ (sampleTable$stage == "Adult") ,] %>% .$samples
cancerSamples <- sampleTable[ (sampleTable$stage == "Cancer") ,] %>% .$samples

# Load clinical data from TCGA for survival stats
clinicalData <- read.table("TCGA-CDR-Table", sep = "\t", header = T)
clinicalData <- clinicalData[clinicalData$type == cancerType, ]

# Load gene names from Gencode
gencodeDf <- read.table("/Users/rajeevana2/gencode_genes_hg38", sep="\t", header= T)
colnames(gencodeDf) <- c("chr", "start", "end", "GeneName", "GeneID")
geneGroups <- read.table(file=paste(analyzedDataDir, paste(tissue, "GeneGroups.tsv", sep=""), sep="/"), header=TRUE, sep="\t") %>%
  dplyr::select(c(group))
geneGroups$GeneID <- rownames(geneGroups)

nonGroupGenes <- gencodeDf[!gencodeDf$GeneID %in% geneGroups$GeneID,]
randomGroup <- sample_n(nonGroupGenes, 1000) %>% dplyr::select(c(GeneID))
randomGroup$group <- "CTRL"

GeneDf <- inner_join(gencodeDf, geneGroups %>% rbind(., randomGroup), by=c("GeneID"))
normalizedCounts <- read.table(file=paste(analyzedDataDir, (paste(tissue, "NormalizedCounts", sep="")), sep="/"), header=TRUE, sep="\t")


# Load cancer RNAseq normalized counts df
cancerExpression <- normalizedCounts[ , grepl( "TCGA" , names( normalizedCounts ) ) ]
colnames(cancerExpression) <- substr(colnames(cancerExpression), 0, 12) %>% gsub("\\.", "-", .)

cancerExpressionLog <- cancerExpression %>% log1p()


# Function to calculate hazard ratio

survivalFunction <- function(clinicalData, geneExpression, test, binStep, binned) {
  
  geneExpression = unlist(geneExpression)
  dataToTest = merge(clinicalData, geneExpression, by.x = "bcr_patient_barcode", by.y = 0)
  names(dataToTest)[names(dataToTest) == "y"] <- "geneExpression"
  dataToTest$age = scale(dataToTest$age,center = T,scale = T)
  head(dataToTest)
  
  if (test == "logRank") {
    fit <- NA; sd <- NA; groupMedians <- c(NA, NA)
    medianLevel <- median(dataToTest$geneExpression, na.rm = T, na.action = na.pass) # median gene expression of one gene for all patients
    dataToTest$BinnedExpression <- ifelse(dataToTest$geneExpression >= medianLevel, "HighExpression", "LowExpression")
    groupMedians <- dataToTest %>% group_by(BinnedExpression) %>% summarise(median = median(geneExpression, na.rm = T, na.action = na.pass))
    groupMedians <- na.omit(groupMedians)
    groupMedians <- groupMedians[ ,'median'] %>% unlist() %>% set_names(unlist(groupMedians[ ,'BinnedExpression']))
    res <- list(fit = fit, sd = sd, groupMedians = groupMedians)

    fit <- survfit(Surv(OS.time,O) ~ BinnedExpression, data = dataToTest)
    sd <- survdiff(Surv(OS.time,O) ~ BinnedExpression, data = dataToTest)
    res <- list(fit = fit, sd = sd, groupMedians = groupMedians)

    return(res)
  }
  
  if (test == "cox") {
    HR <- NA; Pval <- NA;
    cox.base = coxph(Surv(OS.time,OS) ~ geneExpression + age, data = dataToTest) # Add age because it might be a confounding variable
    HR = summary(cox.base)$coefficients['geneExpression',2]
    Pval = summary(cox.base)$coefficients['geneExpression',5]
    res <- list(HR = HR, Pval = Pval)
    return(res)
  }
}


# Prepare inputs for survival function and create df of hazard ratios

# # Expression (not logged)
# survResults <- list()
# for (i in 1:nrow(GeneDf)) {
#   Gene <- GeneDf$GeneID[i]
#   GeneIndex <- which(rownames(cancerExpression) == Gene)
#   
#   res <- NA
#   if (length(GeneIndex) == 1) {
#     expressionLevel <- cancerExpression[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
#     res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionLevel, test = "cox", binStep = 10, binned = F)
#   }
#   survResults[[i]] <- res
# 
# }
# names(survResults) <- GeneDf$GeneID %>% as.character() # can save this as .rda file for later use
# 
# 
# EffectSize <- lapply(survResults, function(x) x[1]) %>% unlist()
# Pvalue <- lapply(survResults, function(x) x[2]) %>% unlist()
# survivalVals <- data.frame(GeneID = names(survResults), EffectSize, Pvalue) %>% set_rownames(NULL)




# # Expression (logged)
# survResultsLog <- list()
# for (i in 1:nrow(GeneDf)) {
#   Gene <- GeneDf$GeneID[i]
#   GeneIndex <- which(rownames(cancerExpressionLog) == Gene)
# 
#   res <- NA
#   if (length(GeneIndex) == 1) {
#     expressionLevel <- cancerExpressionLog[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
#     res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionLevel, test = "cox", binStep = 10, binned = F)
#   }
#   survResultsLog[[i]] <- res
# 
# }
# names(survResultsLog) <- GeneDf$GeneID %>% as.character() # can save this as .rda file for later use
# 
# EffectSizeLog <- lapply(survResultsLog, function(x) x[1]) %>% unlist()
# PvalueLog <- lapply(survResultsLog, function(x) x[2]) %>% unlist()
# survivalValsLog <- data.frame(GeneID = names(survResultsLog), EffectSizeLog, PvalueLog) %>% set_rownames(NULL)




############

# Grouped expression means (not logged)
geneGroupDf <- GeneDf [ (GeneDf$group != "CTRL"),]
groupNames <- geneGroupDf$group %>% unique
survResultsMeans <- list()
for (groupName in groupNames){
  print(groupName)
  Genes <- geneGroupDf[(geneGroupDf$group == groupName) ,] %>% .$GeneID
  GeneIndex <- cancerExpression[(rownames(cancerExpression) %in% Genes) ,]
  
  res <- NA
  if (length(GeneIndex) > 0) {
    expressionMeans <- colMeans(GeneIndex)
    # expressionLevel <- cancerExpression[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
    res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionMeans, test = "cox", binStep = 10, binned = F)
  }
  survResultsMeans[[groupName]] <- res
  
}

names(survResultsMeans) <- groupNames %>% as.character() # can save this as .rda file for later use


EffectSize <- lapply(survResultsMeans, function(x) x[1]) %>% unlist()
Pvalue <- lapply(survResultsMeans, function(x) x[2]) %>% unlist()
survivalValsMeans <- data.frame(geneGroup = names(survResultsMeans), EffectSize, Pvalue) %>% set_rownames(NULL)

############

############

# Grouped expression medians (not logged)
geneGroupDf <- GeneDf [ (GeneDf$group != "CTRL"),]
groupNames <- geneGroupDf$group %>% unique
survResultsMedians <- list()
for (groupName in groupNames){
  print(groupName)
  Genes <- geneGroupDf[(geneGroupDf$group == groupName) ,] %>% .$GeneID
  GeneIndex <- cancerExpression[(rownames(cancerExpression) %in% Genes) ,]
  
  res <- NA
  if (length(GeneIndex) > 0) {
    expressionMedians <- apply(GeneIndex,2,median)
    # expressionLevel <- cancerExpression[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
    res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionMedians, test = "cox", binStep = 10, binned = F)
  }
  survResultsMedians[[groupName]] <- res
  
}

names(survResultsMedians) <- groupNames %>% as.character() # can save this as .rda file for later use


EffectSize <- lapply(survResultsMedians, function(x) x[1]) %>% unlist()
Pvalue <- lapply(survResultsMedians, function(x) x[2]) %>% unlist()
survivalValsMedians <- data.frame(geneGroup = names(survResultsMedians), EffectSize, Pvalue) %>% set_rownames(NULL)

############




# Expression (not logged)
survResults <- list()
for (i in 1:nrow(GeneDf)) {
  Gene <- GeneDf$GeneID[i]
  GeneIndex <- which(rownames(cancerExpression) == Gene)
  
  res <- NA
  if (length(GeneIndex) == 1) {
    expressionLevel <- cancerExpression[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
    res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionLevel, test = "cox", binStep = 10, binned = F)
  }
  survResults[[i]] <- res
  
}
names(survResults) <- GeneDf$GeneID %>% as.character() # can save this as .rda file for later use


EffectSize <- lapply(survResults, function(x) x[1]) %>% unlist()
Pvalue <- lapply(survResults, function(x) x[2]) %>% unlist()
survivalVals <- data.frame(GeneID = names(survResults), EffectSize, Pvalue) %>% set_rownames(NULL)

# Expression (logged)
survResultsLog <- list()
for (i in 1:nrow(GeneDf)) {
  Gene <- GeneDf$GeneID[i]
  GeneIndex <- which(rownames(cancerExpressionLog) == Gene)

  res <- NA
  if (length(GeneIndex) == 1) {
    expressionLevel <- cancerExpressionLog[GeneIndex, ] #One gene's expression in all cancer samples (1D vector)
    res <- survivalFunction(clinicalData = clinicalData, geneExpression = expressionLevel, test = "cox", binStep = 10, binned = F)
  }
  survResultsLog[[i]] <- res

}
names(survResultsLog) <- GeneDf$GeneID %>% as.character() # can save this as .rda file for later use

EffectSizeLog <- lapply(survResultsLog, function(x) x[1]) %>% unlist()
PvalueLog <- lapply(survResultsLog, function(x) x[2]) %>% unlist()
survivalValsLog <- data.frame(GeneID = names(survResultsLog), EffectSizeLog, PvalueLog) %>% set_rownames(NULL)



# rerun cox regression with mean and median expression of all genes in group
# if significant, do log rank (survival curves of those two groups with top and bottom 25% mean expression in group)



# survivalResults <- merge(survivalResults, GeneDf, by = "GeneID")
survivalValsCombined <- inner_join(survivalVals, survivalValsLog, by = "GeneID")
survivalResults <- merge(survivalValsCombined, GeneDf, by="GeneID")

# cor.test(survivalResults$EffectSize, survivalResults$EffectSizeLog, method="spearman")
# 
# ggplot(survivalResults, aes(x=EffectSize)) +
#   geom_density() +
#   coord_cartesian(xlim=c(0.5,1.5))






forPlot <- na.omit(survivalResults)
forPlot$FDR <- p.adjust(forPlot$Pvalue, method = "fdr")
# forPlot <- forPlot %>% na.omit() %>% .[.$FDR < 0.20, ]

groupNames <- forPlot$group %>% unique()
fracTable <- data.frame()
for (groupName in groupNames){
  groupVals <- forPlot[ (forPlot$group == groupName) ,] %>% .[(.$FDR <0.20),]
  groupCount <- forPlot[ (forPlot$group == groupName) ,] %>% nrow
  # print(groupName)
  # print(paste("Frac Greater than 1:", ( nrow(groupVals[(groupVals$EffectSize > 1),]) / groupCount) , sep=" "))
  # print(paste("Frac Less than 1:", ( nrow(groupVals[(groupVals$EffectSize < 1),]) / groupCount), sep=" "))
  
  fracGreater <- ((nrow(groupVals[(groupVals$EffectSize > 1),])) / groupCount ) %>% round(., digits=3)
  fracLess <- ((nrow(groupVals[(groupVals$EffectSize < 1),])) / groupCount ) %>% round(., digits=3)
  fracRow <- c(groupName, fracGreater, fracLess)
  fracTable <- rbind(fracTable, fracRow)
}

colnames(fracTable) <- c("Group", "Frac. >1", "Frac <1")
fracTable <- fracTable[(order(fracTable$`Frac. >1`, decreasing=TRUE)) ,]




survivalValsMeans %>%
  kbl(caption=paste(tissue, "Hazard Ratio Mean Group Expression", sep=" ")) %>%
  kable_classic(full_width = F, html_font = "Cambria")

survivalValsMedians %>%
  kbl(caption=paste(tissue, "Hazard Ratio Median Group Expression", sep=" ")) %>%
  kable_classic(full_width = F, html_font = "Cambria")

fracTable  %>%
  kbl(caption=paste(tissue, "HR Fractions", sep=" "), row.names=FALSE) %>%
  kable_classic(full_width = F, html_font = "Cambria")




# for each group, count fraction that are significant and >1, fraction that are sig. and <1


# give.n <- function(x){
#   return(c(y = 0, label = length(x))) 
# }
# 
# ggplot(data = forPlot, aes(x = group,  y = EffectSizeLog)) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_summary(fun.data = give.n, geom = "text", fun = median,
#                position = position_dodge(width = 0.75)) +
#   coord_cartesian(ylim = c(0,3)) + ###set the y axis limits if needed
#   ggtitle(paste(tissue, "Hazard Ratio Effect Size (Log Expression)", sep=" ")) +
#   ylab("Effect Size") +
#   xlab("Gene Group")
# 
# ggsave(filename = paste(plotsDir, paste(tissue,"LogExpHazardRatioBoxplot.png", sep=""), sep="/"))
# 
# ggplot(data = forPlot, aes(x = group,  y = EffectSize)) +
#   geom_boxplot(outlier.shape = NA) +
#   stat_summary(fun.data = give.n, geom = "text", fun = median,
#                position = position_dodge(width = 0.75)) +
#   coord_cartesian(ylim = c(0,3)) + ###set the y axis limits if needed
#   ggtitle(paste(tissue, "Hazard Ratio Effect Size", sep=" ")) +
#   ylab("Effect Size") +
#   xlab("Gene Group")
# 
# ggsave(filename = paste(plotsDir, paste(tissue,"HazardRatioBoxplot.png", sep=""), sep="/"))
