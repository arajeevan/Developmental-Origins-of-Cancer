library(readr)
library(dplyr)
library(stringr)
library(GenomicFeatures)
library(preprocessCore)
library(ggplot2)
library(fgsea)


outputDir <- "~/Developmental-Origins-of-Cancer"
analyzedDataDir <- paste(outputDir, "AnalyzedData", sep="/")
plotsDir <- paste(outputDir, "Plots", sep="/")
dir.create(plotsDir)
setwd(outputDir)

tissues <- c("Liver", "Kidney", "Brain", "Ovary", "Testis")

for (tissue in tissues){
  
  genesTable <- read.table(file = paste(analyzedDataDir, paste(tissue, "GeneGroups.tsv", sep=""), sep="/"),
                           header=TRUE)
  genesTable$gene_id <- rownames(genesTable)
  countsDf <- read.table(file = paste(analyzedDataDir, (paste(tissue, "NormalizedCounts", sep="")), sep="/"),
                         header=TRUE, sep="\t")
  
  
  # genesTable <- read.table(paste("~/developmentalRawCounts/AnalyzedData/", tissue, "GeneGroups", sep=""), header=TRUE)
  # enhancersDf <- read.table(paste("~/developmentalRawCounts/AnalyzedData/", tissue, "EnhancerGroups", sep=""), header=TRUE)
  # countsDf <- read.table(paste("~/developmentalRawCounts/AnalyzedData/", tissue, "Counts.tsv", sep=""), header=TRUE)
  
  
  geneConv <- read.table("~/gencode_genes_hg38", header=TRUE)
  genesDf <- left_join(x=genesTable, y=geneConv, by='gene_id')
  groups <- genesDf$group %>% unique
  
  for(groupName in groups){
    groupCount <- genesDf[ (genesDf$group == groupName) ,] %>% nrow
    print(groupName)
    print(groupCount)
  }
  
  
  
  countsDf$sumCol <- rowSums(countsDf)
  ensIDs <- subset(countsDf, countsDf['sumCol'] > 0) %>% rownames %>% as.data.frame 
  colnames(ensIDs) <- "gene_id"
  allGenes <- left_join(x=ensIDs, y=geneConv, by='gene_id') %>% .[,c('gene_name')]
  
  c2Pathways <- gmtPathways("c2.all.v7.2.symbols.gmt")
  c4Pathways <- gmtPathways("c4.all.v7.2.symbols.gmt")
  c6Pathways <- gmtPathways("c6.all.v7.2.symbols.gmt")
  
  concatPathways <- c(c2Pathways, c4Pathways, c6Pathways)
  pathwayList <- names(concatPathways)
  
  getGroupName <- function(group){
    if(group == "UPUP"){
      newGroup <- "Cancer-high"
    }
    if(group == "UPDOWN"){
      newGroup <- "Re-repressed"
    }
    if(group == "DOWNUP"){
      newGroup <- "Re-activated"
    }
    if(group == "DOWNDOWN"){
      newGroup <- "Development-high"
    }
    return(newGroup)
  }
  
  # Re-activated (fg) and Cancer-high (bg)
  
  reactivated <- subset(genesDf, genesDf['group'] == "DOWNUP") %>% .[,c('gene_name')]
  cancerhigh <- subset(genesDf, genesDf['group'] == "UPUP") %>% .[,c('gene_name')]
  
  count = 1
  
  gseaDf <- data.frame(pathway=as.Date(character()),
                       oddsRatio=double(),
                       pVal=double(),
                       stringsAsFactors=FALSE)
  
  for(pathway in concatPathways){
    pathwayName <- pathwayList[[count]]
    quad1 <- reactivated[reactivated %in% pathway] %>% length #in reactivated, in pathway
    quad2 <- cancerhigh[cancerhigh %in% pathway] %>% length #in cancer high, in pathway
    quad3 <- reactivated[!reactivated %in% pathway] %>% length #in reactivated, not in pathway
    quad4 <- cancerhigh[!cancerhigh %in% pathway] %>% length #in cancer high, not in pathway
    
    fisherTable <- matrix(c(quad1, quad3, quad2, quad4), nrow = 2,
                          dimnames = list(c("InPathway", "NotInPathway"),
                                          c("Reactivated", "CancerHigh")))
    fisherTest <- fisher.test(fisherTable, alternative = "greater")
    oddsRatio<- fisherTest[[3]][[1]]
    pVal <- fisherTest[[1]][[1]]
    pathwayRow <- c(pathwayName, oddsRatio, pVal)
    gseaDf <- rbind(gseaDf, pathwayRow)
    count <- count + 1
    
    
    # table = [ [a, b], [ c, d] ]
    # OR = (a*d) / (b*c)
    # if b or c is 0, then OR will be INF, if a or d is 0, then OR will be 0 (add 1 to everything)
    
    print(pathwayName)
    print(fisherTable)
    print(oddsRatio)
    
  }
  
  
  colnames(gseaDf) <- c("pathway", "oddsRatio", "pVal")
  gseaDf$padj <- p.adjust(gseaDf$pVal %>% as.numeric, method="fdr")
  rankedGseaDf <- subset(gseaDf, gseaDf['padj'] < 0.05) %>%
    .[order(.$oddsRatio, decreasing=TRUE),]
  row.names(rankedGseaDf) <- NULL
  topRankedGsea <- rankedGseaDf[1:10,] %>% na.omit
  topRankedGsea$oddsRatioRounded <- round(x=topRankedGsea[,'oddsRatio'] %>%
                                            as.numeric, digits=3)
  
  groupGseaPlot <- ggplot(topRankedGsea, aes(x=reorder(pathway, oddsRatioRounded), y=oddsRatioRounded)) +
    geom_col() +
    coord_flip() +
    ggtitle(paste(tissue, "Gene Set Enrichment:", "Reactivated/Cancer-high", sep=" ")) +
    ylab("odds ratio") +
    xlab("") +
    theme(plot.title = element_text(size = 10))
  
  plotfilename <- paste(plotsDir, paste(tissue,"GeneSetEnrichmentPlot", "Reactivated", ".png", sep=""), sep="/")
  ggsave(plotfilename, width = 10, height = 6, units = "in")
  
  
  
  
  
  
  # Re-repressed (fg) and Development-high (bg)
  
  rerepressed <- subset(genesDf, genesDf['group'] == "UPDOWN") %>% .[,c('gene_name')]
  devhigh <- subset(genesDf, genesDf['group'] == "DOWNDOWN") %>% .[,c('gene_name')]
  
  count = 1
  
  gseaDf <- data.frame(pathway=as.Date(character()),
                       oddsRatio=double(),
                       pVal=double(),
                       stringsAsFactors=FALSE)
  
  for(pathway in concatPathways){
    pathwayName <- pathwayList[[count]]
    
    quad1 <- rerepressed[rerepressed %in% pathway] %>% length
    quad2 <- devhigh[devhigh %in% pathway] %>% length
    quad3 <- rerepressed[!rerepressed %in% pathway] %>% length
    quad4 <- devhigh[!devhigh %in% pathway] %>% length
    
    fisherTable <- matrix(c(quad1, quad3, quad2, quad4), nrow = 2,
                          dimnames = list(inPathway = c("Reactivated", "CancerHigh"),
                                          notInPathway = c("Reactivated", "CancerHigh")))
    fisherTest <- fisher.test(fisherTable, alternative = "greater")
    oddsRatio<- fisherTest[[3]][[1]]
    pVal <- fisherTest[[1]][[1]]
    pathwayRow <- c(pathwayName, oddsRatio, pVal)
    gseaDf <- rbind(gseaDf, pathwayRow)
    count <- count + 1
    
    
    
  }
  
  
  colnames(gseaDf) <- c("pathway", "oddsRatio", "pVal")
  gseaDf$padj <- p.adjust(gseaDf$pVal %>% as.numeric, method="fdr")
  rankedGseaDf <- subset(gseaDf, gseaDf['padj'] < 0.05) %>%
    .[order(.$oddsRatio, decreasing=TRUE),]
  row.names(rankedGseaDf) <- NULL
  topRankedGsea <- rankedGseaDf[1:10,] %>% na.omit
  topRankedGsea$oddsRatioRounded <- round(x=topRankedGsea[,'oddsRatio'] %>%
                                            as.numeric, digits=3)
  
  groupGseaPlot <- ggplot(topRankedGsea, aes(x=reorder(pathway, oddsRatioRounded), y=oddsRatioRounded)) +
    geom_col() +
    coord_flip() +
    ggtitle(paste(tissue, "Gene Set Enrichment:", "Rerepressed/Development-high", sep=" ")) +
    ylab("odds ratio") +
    xlab("") +
    theme(plot.title = element_text(size = 10))
  
  plotfilename <- paste(plotsDir, paste(tissue,"GeneSetEnrichmentPlot", "Rerepressed", ".png", sep=""), sep="/")
  ggsave(plotfilename, width = 10, height = 6, units = "in")
  
  
  
  
  # Motif Tables
  
  # for(group in groups){
  #   groupName <- getGroupName(group)
  #   groupDf <- read.table(paste("~/developmentalRawCounts/AnalyzedData/Motifs/", tissue, "_Enr_TFs_Table_", group, sep=""), header=TRUE)
  #   groupDf$oddsRatioRounded <- round(x=groupDf[,'oddsRatio'] %>%
  #                                       as.numeric, digits=3)
  #   groupDfTop <- groupDf[1:10,] %>% drop_na
  #   
  #   groupMotifsPlot <- ggplot(groupDfTop, aes(x=reorder(motif, oddsRatio), y=oddsRatioRounded)) +
  #     geom_col() +
  #     coord_flip() +
  #     ggtitle(paste(tissue, "Motif Enrichment:", groupName, sep=" ")) +
  #     ylab("odds ratio") +
  #     xlab("") +
  #     theme(plot.title = element_text(size = 10))
  #   
  #   plotfilename <- paste(plotsDir, paste(tissue,"TFEnrichmentPlot", group, ".png", sep=""), sep="/")
  #   ggsave(plotfilename, width = 10, height = 6, units = "in")
  
}


