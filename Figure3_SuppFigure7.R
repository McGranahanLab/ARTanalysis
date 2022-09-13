#####################################################################
## Figure 3. Gene based analyses (& Supplementary Figure 7)
#####################################################################
# written by Haoran Zhai (haoran.zhai.17@ucl.ac.uk) and run in R version 4.0.2

# Description:
# Scripts to create Figure 5 and Supplementary Figure 12 in the manuscript named "Replication timing alterations impact mutation acquisition during tumour evolution".
# Data accessibility statement can be found in the manuscript.

#libraries
options(stringsAsFactors = F)
library(CNTools)
library(TCGA2STAT)
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ComplexUpset)
library(limma)
library(edgeR) 
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gtable)
library(ggbeeswarm)
library(ggrepel)
library(ggalluvial)
library(readxl)
library(writexl)
library(tidyverse)
library(dplyr)

#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

cancer.types <- c("BRCA", "LUAD", "LUSC")
cancer_cellLines <- list('BRCA' = c('SK-BR3', 'MCF-7', 'T47D', 'MDA453'),
                         'LUAD' = c('H1650', 'H1792', 'H2009', 'A549'),
                         'LUSC' = c('H520', 'H2170', 'SW900'))
normal_cellLines <- c('BRCA' = 'HMEC', 'LUAD' = 'T2P', 'LUSC' = 'HBEC3')
chr_to_use <- paste0('chr', c(1:22))

### Figure 5C: gene density ===========
# Load overlaping.ART regions:
overlap.ART_list <- lapply(cancer.types, function(x){
  print(x)
  normal <- normal_cellLines[[x]]
  data <- readRDS(paste0(data_dir, x, '_overlapping_ARTregions.rds')) # Overlapping ART regions identified in Figure 2E.
  return(data)
})
names(overlap.ART_list) <- cancer.types

# find % of genes per cancer: 
overlap.geneProp.ART_list <- lapply(cancer.types, function(cancer){
  print(cancer)
  overlappingART <- overlap.ART_list[[cancer]]
  #compare fraction of genes within 50kb windows across different replication timings
  #get position of all genes
  genesPos_gr       <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  overlappingART_gr <- GRanges(seqnames = overlappingART$chr, IRanges(start = overlappingART$start + 1, end = overlappingART$stop))
  overlap           <- findOverlaps(genesPos_gr, overlappingART_gr)
  
  #calculate fraction of gene overlap
  intersects        <- pintersect(overlappingART_gr[subjectHits(overlap)], genesPos_gr[queryHits(overlap)])
  intersects$index  <- subjectHits(overlap)
  intersects_perBin <- lapply(unique(intersects$index), function(x){
    tmp<-intersects[intersects$index == x] %>% as.data.frame()
    data <- GenomicRanges::reduce(intersects[intersects$index == x], ignore.strand = T)
    dt <- data.frame(index = x, overlap = sum(width(data)))
  })
  intersects_perBin <- Reduce(rbind, intersects_perBin)
  
  overlappingART$geneOverlap <- 0
  overlappingART$geneOverlap[intersects_perBin$index] <- intersects_perBin$overlap
  overlappingART$geneOverlapFrac <- overlappingART$geneOverlap / (overlappingART$stop - overlappingART$start)
  overlappingART <- overlappingART %>% mutate(cancerType = cancer)
  return(overlappingART)
})
names(overlap.geneProp.ART_list) <- cancer.types

# read data & make plots to show % of genes: 
plot_data_raw <- data.frame()
for (i in names(overlap.geneProp.ART_list)) {
  dt <- overlap.geneProp.ART_list[[i]] %>% 
    dplyr::select(timing, geneOverlapFrac, class, cancerType)
  plot_data_raw <- rbind(plot_data_raw, dt)
}
plot_data_raw[1:2,]

# annot.col:
annot.col <- list(RT = c(brewer.pal(11, 'PiYG')[c(2,9,3,10)]))
names(annot.col$RT) <- c("early","later","earlier", "late")

# plot:
plot_data <- plot_data_raw %>% filter(class %in% c('not_altered', 'shared', 'recurrent'))
comparisons <- list(c('early', 'later'),c('earlier', 'late')) 

p.violin_gene.density <- plot_data %>% 
  ggplot(aes(x=factor(timing,levels=c(names(annot.col$RT))), y = geneOverlapFrac, 
             fill=factor(timing,levels=names(annot.col$RT)))) + 
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_fill_manual(name = 'Replication Timing', values = annot.col$RT) +
  scale_colour_manual(name = 'Replication Timing', values = annot.col$RT) +
  stat_compare_means(method = 'wilcox', label = 'p.signif',
                     comparisons = comparisons) +
  facet_grid(.~factor(cancerType, levels = cancer.types)) +
  xlab('') + ylab('Gene proprotion within 50kb bins') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),
        axis.text.y = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text = element_text(size = 10, color = 'black'),
        legend.title = element_text(size = 12, color = 'black', face = 'bold'),
        plot.title = element_text(size = 15, color = 'black', face = 'bold'),
        legend.position = 'none',
        strip.text.x = element_text(size = 12),
        panel.grid = element_blank())
p.violin_gene.density

pdf(paste0(output_dir, 'geneProp.RT_violin.per.cancer.pdf'), width = 7, height = 5, useDingbats = F)
print(p.violin_gene.density)
dev.off()

#####################################################################
### Supplementary Figure 12A: Proportions of genes with consistent RT or shared/recurrent ART (in BRCA, LUAD and LUSC) ======
all.cellLines <- data.frame(cellLine = c(unlist(cancer_cellLines), normal_cellLines)) %>%
  rownames_to_column(var = 'cancer.type') %>% 
  mutate(cancer.type = substr(cancer.type, 1,4)) %>%
  mutate(tissue.type = if_else(cellLine %in% normal_cellLines, 'normal', 'cancer'))
all.cellLines[1:4,]

# load ARTextreme.genes:
ARTgenes.extreme_list <- readRDS(paste0(data_dir, 'ARTgenes.per.cell_list.RDS'))

# shared ARTgenes:
ARTgenes.consist_type_list <- lapply(cancer.types, function(x){
  print(x)
  # genes with shared/recurrent ART or consistent RT in BRCA, LUAD and LUSC (processed in Figure 2E)
  data <- fread( paste0(data_dir, x, "_sharedARTgenes.txt") ) %>%
    unique
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

# shared art.genes in 3 cancer:
{
  ARTgenes.consist_df <- ARTgenes.consist_type_list[[1]] %>% 
    dplyr::select('gene_name','chr','start','stop','timing') #%>% rename(BRCA = timing) #
  colnames(ARTgenes.consist_df)[5] <- cancer.types[1]
  ARTgenes.consist_df[1:2,]
  
  for (i in cancer.types[2:3]) {
    df <- ARTgenes.consist_type_list[[i]] %>% dplyr::select('gene_name','timing')
    colnames(df)[2] <- i
    ARTgenes.consist_df <- ARTgenes.consist_df %>% 
      full_join(df, by=c('gene_name'))
  }
  ARTgenes.consist_df[1:2,]
  
  count_ART <- apply(ARTgenes.consist_df[,c(5:7)], 1, function(x) length(x[x%in%c('earlier','later')]))
  count_ART %>% table()
  
  timings <- apply(ARTgenes.consist_df[,c(5:7)], 1, function(x){
    if( length(unique(x))==1 ){
      y <- unique(x)
    }else{
      y <- 'discordant'
    }
    return(y)
  })
  ARTgenes.consist_df$timing <- timings
  ARTgenes.consist_df$count_ART <- count_ART
  
  keep <- apply(ARTgenes.consist_df[,c(5:7)], 1, function(x) !anyNA(x) )
  genes.consistCRT <- ARTgenes.consist_df[keep, ]%>%
    filter(count_ART==0 & timing!='discordant')
  
}
genes.sharedART <- ARTgenes.consist_df%>%filter(count_ART==3)

## make plots 
plot.gene.frac_list <- lapply(cancer.types, function(x){
  print(x)
  genes.all <- ARTgenes.extreme_list[[x]]%>%pull(gene_name) %>% unique()
  shared.art <- ARTgenes.consist_type_list[[x]]
  genes.df <- data.frame(gene_name = genes.all) %>%
    left_join(shared.art[,c('gene_name','timing')]) %>%
    mutate(timing=if_else(is.na(timing), 'unique', timing))
  genes.df$timing %>% unique()
  df <- table(genes.df$timing) %>% as.data.frame() %>% rename(timing=Var1) %>%
    mutate(Frac = Freq/length(unique(genes.df$gene_name))) %>%
    mutate(cancer.type = x)
  df
})
plot.gene.frac <- Reduce(rbind, plot.gene.frac_list) %>%
  mutate(timing=as.character(timing)) %>%
  mutate(text=paste0(round(Frac*100, digits = 1), '%', ' (N=',Freq,')'))

# annot.color:
annot.color <- list(ART=c('early'=brewer.pal(n = 11, 'PiYG')[2],
                          'earlier' = brewer.pal(n = 11, 'PiYG')[3], 
                          'later' = brewer.pal(n = 11, 'PiYG')[9],
                          'late'=brewer.pal(n = 11, 'PiYG')[10],
                          'unique'='grey'),
                    CAG=c("OG"=brewer.pal(9,'Set1')[1],
                          "TSG"=brewer.pal(9,'Set1')[2],
                          "driver.gene"=brewer.pal(9,'Set1')[3]))

# plots:
plot.gene.frac[1:4,]
p.frac <- plot.gene.frac %>%
  ggplot(aes(x=factor(cancer.type,levels = cancer.types), y=Frac, 
             fill=factor(timing,levels = names(annot.color$ART)))) + 
  geom_bar(stat = 'identity', position = 'stack', colour = 'black') + 
  geom_text(aes(label=text), position = position_stack(vjust = 0.5),size = 3,color='white') +
  scale_fill_manual(name = 'Shared.Timing', values = annot.color$ART,
                    labels=c('early'="Early_consistent",'earlier'="Earlier_shared/recurrent",
                             'later'="Later_shared/recurrent", 'late'="Late_consistent",
                             'unique'='Unique_per.cell.line')) +
  ylab(paste0('% Genes with replication timing')) +
  theme_bw() + theme_hz +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 0), #axis.ticks.x = element_blank(), 
        panel.grid = element_blank() )
p.frac

pdf(paste0(output_dir,'Frac.sharedART.genes_per.cancer','.pdf'), 
    width = 7, height = 6)
print(p.frac)
dev.off()

#####################################################################
### Figure 3B: Cancer genes (and Supplementary Figure 7B: Oncogenes and Tumour suppressor genes) ===========
# Load overlaping.ART genes:
overlap.ARTgenes_list <- lapply(cancer.types, function(x){
  print(x)
  normal <- normal_cellLines[[x]]
  data <- readRDS(paste0(data_dir, x, '_overlapping_ARTgenes.rds')) %>% # genes with recurrent/shared/unique ART or consistent RT (identified in Figure 2E)
    mutate(cancerType = x)
  return(data)
})
names(overlap.ARTgenes_list) <- cancer.types
overlap.ARTgenes_list[[1]][1:2,]

overlap.ART_cancer.genes_list <- lapply(cancer.types, function(cancer){
  overlap.ART <- overlap.ARTgenes_list[[cancer]] %>% filter(class %in% c('shared', 'recurrent'))
  overlap.ART[1:2,] 
  # "onco","tsg","cancerType_cancerGene":
  column.used <- c("onco","tsg","cancerType_cancerGene",'cancerGene','essential')
  dt.genes_list <- lapply(column.used, function(column){
    print(column)
    art.dt <- overlap.ART
    art.dt[1:2,]
    colnames(art.dt)[match(column, colnames(art.dt))] <- 'column.used'
    gene.dt <- art.dt %>% group_by(timing) %>%
      summarise(value = sum(column.used), total = n(),
                genes = paste(gene_name[column.used][order(gene_name[column.used],decreasing = F)], 
                              collapse = '\n')) %>%
      mutate(fraction = value / total, geneType=column, cancerType = cancer)
    gene.dt
  })
  data.genes <- Reduce(rbind, dt.genes_list)
  data.genes
  out <- data.genes
  return(out)
})
names(overlap.ART_cancer.genes_list) <- cancer.types

# make plot:
plot_data_raw <- Reduce(rbind, overlap.ART_cancer.genes_list)
plot_data <- plot_data_raw %>%
  mutate(label=if_else(geneType%in%c("onco","tsg","cancerType_cancerGene"), genes, NULL)) %>%
  filter(geneType %in% c("cancerType_cancerGene", "onco", "tsg")) %>%
  mutate(fraction = fraction,
         geneType=if_else(geneType=="cancerType_cancerGene", 'Cancer genes',
                          if_else(geneType=="onco", 'Oncogenes',
                                  if_else(geneType=="tsg", 'Tumour suppressor genes','Other'))))

gene.types <- c("Cancer genes", "Oncogenes", "Tumour suppressor genes")

# annot.col:
annot.col <- list(RT = c(brewer.pal(11, 'PiYG')[c(3,9)]))
names(annot.col$RT) <- c("earlier", "later")

# geneType:
{
  plot_data.cancer <- plot_data %>% filter(geneType=="Cancer genes")
  plot_data.cancer[1:2,]
  p.cancer_genes <- plot_data.cancer %>%
    ggplot(aes(x = factor(cancerType,levels = c(cancer.types)), 
               y = fraction, fill = factor(timing, levels = names(annot.col$RT)))) + 
    geom_bar(stat = 'identity', position = position_dodge()) + 
    geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = 1.1, size = 3) +
    scale_fill_manual(name = 'Altered.RT', values = annot.col$RT) +
    xlab('') + ylab('Proportion of genes') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    theme_bw() + theme_hz +
    theme(axis.text.x = element_text(angle = 0), strip.text.x = element_text(size=12),
          panel.grid = element_blank())
  p.cancer_genes
  
  plot_data.CAG <- plot_data %>% filter(geneType%in%c("Oncogenes","Tumour suppressor genes"))
  plot_data.CAG[1:2,]
  p.CAG <- plot_data.CAG %>%
    ggplot(aes(x = factor(cancerType,levels = c(cancer.types)), 
               y = fraction, fill = factor(timing, levels = names(annot.col$RT)))) + 
    geom_bar(stat = 'identity', position = position_dodge()) + 
    geom_text(aes(label = label), position = position_dodge(width = 0.9), vjust = 0.9, size=4) +
    scale_fill_manual(name = 'Altered.RT', values = annot.col$RT) +
    facet_grid(.~factor(geneType, levels=c(gene.types))) +
    xlab('') + ylab('Proportion of genes') +
    scale_y_continuous(expand = c(0,0,0.05,0.001)) +
    theme_bw() + theme_hz +
    theme(axis.text.x = element_text(angle = 0), strip.text.x = element_text(size=12),
          panel.grid = element_blank())
  p.CAG
  
}
p.cancer_genes
p.CAG


pdf(paste0(output_dir, 'bar_frac_cancerGenes.pdf'), width = 7, height = 4, useDingbats = F)
print(p.cancer_genes)
dev.off()

pdf(paste0(output_dir, 'bar_frac_OG.TSG.pdf'), width = 12, height = 9, useDingbats = F)
print(p.CAG)
dev.off()

#####################################################################
### Figure 3C: Differential gene analyses among genes with ART, followed by the scripts for bootstrapping analyses ================
## Step 1. Download and process TCGA data: ========
# function to download count data (getTCGA function)
func.tcga_download_hz <- function(cancer_type="LUSC", tumour.idx = paste0(0, c(1:9)),
                                  data.type = c("count", "RPKM")[1],
                                  normal.idx = paste0(1, c(0:9)),
                                  control.idx = paste0(2, c(0:9))){
  library(TCGA2STAT)
  data <- getTCGA(disease=cancer_type, data.type="RNASeq2", 
                  type= data.type , clinical = TRUE)
  # extract count data:
  count_data <- as.data.frame(data$dat)
  colnames(count_data) <- substr(colnames(count_data), 1, 15)
  colnames(count_data) <- paste(cancer_type, colnames(count_data), sep="-")
  
  # sample.data:
  sample.df <- data.frame(sampleID = colnames(count_data), cancerType = cancer_type) %>%
    mutate(tissueID = matrix(unlist(strsplit(.$sampleID, split='-')), byrow = T, ncol = 5)[,5]) %>%
    mutate(tissueType = if_else(.$tissueID %in% tumour.idx, 'tumour',
                                if_else(.$tissueID %in% normal.idx, 'normal',
                                        if_else(.$tissueID %in% control.idx, 'control', 'other')))) %>%
    mutate(PatientID = substr(.$sampleID, 1, 17))
  
  # survival data
  surv <- data$merged.dat %>% dplyr::select(bcr, status, OS)
  colnames(surv)[1] <- "PatientID"
  #convert survival data to years
  surv$OS <- surv$OS/365
  #censor at 5 years
  censor_patients <- dplyr::filter(surv, OS > 5 & status==1) %>% 
    dplyr::select(PatientID) 
  censor_patients
  surv$status_5yrs <- surv$status
  surv$status_5yrs[surv$PatientID %in% censor_patients$PatientID] <- 0
  surv$OS_5yrs <- surv$OS
  surv$OS_5yrs[which(surv$OS > 5)] <- 5
  #add cancer-type to patient id's
  surv$PatientID <- paste(cancer_type, surv$PatientID, sep="-")
  
  # matched normal and tumour: a subset
  tum.norm.rnaseq2 <- TumorNormalMatch(data)
  colnames(tum.norm.rnaseq2$primary.tumor) <- paste(cancer_type, colnames(tum.norm.rnaseq2$primary.tumor), sep="-")
  # tum.norm.rnaseq2$primary.tumor[1:2, 1:4]
  colnames(tum.norm.rnaseq2$normal) <- paste(cancer_type, colnames(tum.norm.rnaseq2$normal), sep="-")
  
  #output
  ls <- list(count_data, sample.df, surv, tum.norm.rnaseq2)
  names(ls) <- c('RNASeq2_exprSet','sample.df', 'survival', 'TumorNormalMatch')
  return(ls)
}

## then run&save one by one:
TCGA.download <- lapply(c('LUAD','LUSC','BRCA'), function(cancer_type){
  print(cancer_type)
  data.type.x <- c("count", "RPKM")[1]
  tcga <- func.tcga_download_hz(cancer_type = cancer_type, 
                                data.type = data.type.x)
  saveRDS(object = tcga, 
          file = paste0(data_dir, cancer_type,'_RNAseq2.',data.type.x,'_clin.survival_',Sys.Date(), ".rds"))
  data.type.x <- c("count", "RPKM")[2]
  tcga <- func.tcga_download_hz(cancer_type = cancer_type, 
                                data.type = data.type.x)
  saveRDS(object = tcga, 
          file = paste0(data_dir, cancer_type,'_RNAseq2.',data.type.x,'_clin.survival_',Sys.Date(), ".rds"))
  return(tcga)
})
names(TCGA.download) <- c('LUAD','LUSC','BRCA')

## Step 2. Run differential gene analyses: ========
exprCounts_list <- lapply(cancer.types, function(x){
  dt <- TCGA.download[[x]]
})
names(exprCounts_list) <- cancer.types

# filter BRCA --> only keep lobular and ductal BRCA:
subtype.brca <- read_xlsx(paste0(data_dir, 'SubtypeBRCA_TCGA_Thennavan_CellGenomics2021/1-s2.0-S2666979X21000835-mmc2_Histologic.BRCA.xlsx')) %>%
  as.data.frame()
subtype.brca.filter <- subtype.brca %>% 
  filter(`2016 Histology Annotations` %in% c('Invasive ductal carcinoma', 'Invasive lobular carcinoma')) %>%
  mutate(PatientID = paste0('BRCA-', substr(CLID, start = 1, stop = 12)), .after='CLID')
subtype.brca.pts <- subtype.brca.filter$PatientID %>% unique()
subtype.brca.CLID <- subtype.brca.filter$CLID %>% unique()

# functions for DE analyses:
DEG_DESeq2_function <- function(gene_counts_df, sample_id_df, geneList,
                                group.levels=c('tumour','normal'),
                                Cutoff.per.count = 1, Cutoff.fract = 0.2,
                                FC_cutoff= 2, padj_cutoff= 0.05,
                                qvalue_cutoff = 0.01){
  # 1-group.list:
  group_list <- data.frame(sample.id = colnames(gene_counts_df)) %>%
    left_join(sample_id_df, by=c('sample.id'='sampleID'))
  group_list <- factor(sample_id_df$tissueType, levels = group.levels)
  group_list <- relevel(group_list, ref = group.levels[2])
  
  # 2-countData as.matrix:
  countData <- as.matrix(round(gene_counts_df))
  colData <- data.frame(row.names = colnames(gene_counts_df),
                        group_list = group_list)
  library(DESeq2)
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(gene_counts_df)),
                                colData = colData,
                                design = ~ group_list)
  
  # 3-Pre-filtering: 
  keep <- rowSums(counts(dds) > Cutoff.per.count) >= (ncol(gene_counts_df) * Cutoff.fract)
  dds <- dds[keep, ]
  dds <- DESeq(dds) 
  
  # 5-normalise counts:
  # 1-cor() between samples: using normalized_counts 
  normalized_counts <- DESeq2::counts(dds, normalized = FALSE) 
  normalized_counts_mad <- apply(normalized_counts, 1, stats::mad) # calculate mad: numeric
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  
  # normalise counts: by log transformation
  vst_d <- vst(dds, blind=FALSE) 
  
  ## extract exprData: assay()
  exprMatrix_vst <- assay(vst_d) # then order by normalized_counts_mad
  exprMatrix_vst <- exprMatrix_vst[order(normalized_counts_mad, decreasing=T), ]
  
  ### extract DE results: 
  res <- results(dds, contrast=c("group_list", group.levels[1], group.levels[2]) )
  
  # to see mean() of each group_list:
  base_group.level.1 <- counts(dds, normalized=TRUE)[, colData(dds)$group_list == group.levels[1]]
  base_group.level.2 <- counts(dds, normalized=TRUE)[, colData(dds)$group_list == group.levels[2]]
  
  ## group.levels[1] = "Tumour"
  if (is.vector(base_group.level.1)){
    baseMean_group.level.1 <- as.data.frame(base_group.level.1)
  } else {
    baseMean_group.level.1 <- as.data.frame(rowMeans(base_group.level.1)) # mean in all samples
  }
  colnames(baseMean_group.level.1) <- 'baseMean_group.level.1'
  
  ## group.levels[2] = "Normal"
  if (is.vector(base_group.level.2)){
    baseMean_group.level.2 <- as.data.frame(base_group.level.2)
  } else {
    baseMean_group.level.2 <- as.data.frame(rowMeans(base_group.level.2))
  }
  colnames(baseMean_group.level.2) <- 'baseMean_group.level.2'
  
  ##combine: 
  res.comb <- cbind(baseMean_group.level.1, 
                    baseMean_group.level.2,
                    as.data.frame(res)) %>% 
    rownames_to_column(var = 'gene_id')
  colnames(res.comb)[grep('group.level.', colnames(res.comb))] <- c(paste0('baseMean_group.', group.levels))
  
  # rank by padj: reorder by the smallest p value (the largest DE)
  DEG_DEseq2 <- res.comb[order(res.comb$padj, decreasing = F), ] %>% na.omit() %>%
    mutate(change = if_else(.$padj < padj_cutoff & abs(.$log2FoldChange) > log2(FC_cutoff), ##
                            if_else(.$log2FoldChange > log2(FC_cutoff),'Up','Down'), 'Stable'))
  sample_id_df[1:2,]
  DEG_DEseq2$change %>% table()
  tab <- xtabs( ~ tissueType, data = sample_id_df)
  export <- list(DEG_DEseq2, as.data.frame(exprMatrix_vst), tab, sample_id_df)
  names(export) <- c('DESeq2_DEG.df', 'DESeq2_exprMatrix_vst', 'DESeq2_Cohort_tab',
                     'clin_df')
  return(export)
}

DEG_limma.voom_function <- function(gene_counts_df, sample_id_df,
                                    geneList,
                                    group.levels=c('tumour','normal'),
                                    Cutoff.per.count = 1, Cutoff.fract = 0.2,
                                    FC_cutoff= 2, padj_cutoff=0.05){
  library(limma)
  library(edgeR)
  # 1-group.list:
  sample_id_df[1:2,]
  gene_counts_df[1:2, 1:4]
  group_list <- data.frame(sample.id = colnames(gene_counts_df)) %>%
    left_join(sample_id_df, by=c('sample.id'='sampleID'))
  group_list[1:4,]
  group_list <- factor(sample_id_df$tissueType, levels = group.levels)
  group_list <- relevel(group_list, ref = group.levels[2])
  
  # 2-enable comparison between groups: model.matrix()
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  head(design)
  
  # 3.1-exprSet_counts:
  exprSet <- as.matrix(round(gene_counts_df)) 
  rownames(design)=colnames(exprSet)
  exprSet[1:2,1:4]
  # 3.2-Pre-filtering: Cutoff.per.count=1, Cutoff.fract= 20% (0.2) [HR added; consistent with DESeq2]
  keep <- rowSums(exprSet>Cutoff.per.count) >= (ncol(gene_counts_df) * Cutoff.fract)
  exprSet <- exprSet[keep, ]
  
  # 4-contrast.matrix:
  contrast.matrix <- limma::makeContrasts(contrasts = paste0(group.levels[1],'-',group.levels[2]),
                                          levels = design)
  # 5-DGEList: normalize
  dge <- DGEList(counts=exprSet)
  
  # Filter out low expression genes for the voom fit
  library(edgeR)
  keep.exprs <- filterByExpr(dge, group = group_list, min.count = 30)
  dge <- dge[keep.exprs,, keep.lib.sizes = FALSE]
  
  # Normalise by RNA content using TMMs 
  dge <- calcNormFactors(dge, method = "TMM") # Carlos
  
  # 6-voom:
  v <- voom(dge, design, plot=FALSE) 
  fit <- lmFit(v, design)
  
  fit2_contrasts <- contrasts.fit(fit, contrast.matrix)
  fit3_robust <- eBayes(fit2_contrasts, robust = TRUE)
  DEG_limma_voom = topTable(fit3_robust,sort.by = "P", coef=colnames(contrast.matrix), 
                            n=Inf) %>% na.omit(.) # colname: keep to use P.Value
  
  # 7-define change of Expr:
  DEG_limma_voom[1:2,]
  DEG_limma_voom <- DEG_limma_voom %>% 
    mutate(change = if_else(.$adj.P.Val < padj_cutoff & abs(.$logFC) > log2(FC_cutoff),
                            if_else(.$logFC > log2(FC_cutoff),'Up','Down'), 'Stable'))
  DEG_limma_voom$change %>% table()
  DEG_limma_voom$adj.P.Val %>% range()
  DEG_limma_voom$P.Value %>% range()
  sample_id_df$tissueType %>% table()
  tab <- xtabs( ~ tissueType, data = sample_id_df)
  
  export <- list(DEG_limma_voom, as.data.frame(exprSet), tab, sample_id_df)
  names(export) <- c('limma.voom_DEG.df', 'limma.voom_exprMatrix', 'limma.voom_Cohort_tab',
                     'clin_df')
  return(export)
}

draw_volcano <- function(method = c('limma', 'DESeq2'),
                         DES_output, change.factor=c("Up","Down","Stable"),
                         sample.cohort,
                         FC_cutoff = 2, 
                         padj_cutoff=0.05,
                         cancer_type='LUAD'){
  theme_hz <- theme(axis.text = element_text(size = 12, color = 'black'),
                    axis.title = element_text(size = 14, color = 'black'),
                    legend.text = element_text(size = 10, color = 'black'),
                    legend.title = element_text(size = 12, color = 'black', face = 'bold'),
                    plot.title = element_text(size = 15, color = 'black', face = 'bold'))
  
  if(method=='DESeq2'){
    DES_output[1:2,]
    g <- DES_output %>%
      ggplot(aes(x=log2FoldChange, y= -log10(padj), 
                 colour=factor(change, levels = change.factor))) +
      geom_point(alpha=0.5, size=1.75) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
      geom_vline(xintercept = c(-log2(FC_cutoff), log2(FC_cutoff)), linetype = "dashed") +
      scale_colour_manual(name='GeneExpression', values = c('Up'='#AE3121','Down'="#3182BD", 'Stable'='lightgrey')) +
      xlab(paste0("Log2 fold change (FC_cutoff=", (FC_cutoff),
                  '; p.adj_cutoff=', padj_cutoff, ')')) + 
      ylab("-log10 (p.adj value)") +
      ggtitle(paste0('DEGs in ',cancer_type,' (Normal=', sample.cohort[1],'; Tumour=', sample.cohort[2],')'), 
              subtitle = paste0('(Up genes=', nrow(DES_output[DES_output$change =='Up',]),
                                '; Down genes=', nrow(DES_output[DES_output$change =='Down',]),
                                ')(method=', method,')'))+
      theme_bw()+ theme_hz
  }else if(method=='limma'){
    DES_output[1:2,]
    
    g <- DES_output %>%
      ggplot(aes(x=logFC, y= -log10(adj.P.Val), 
                 colour=factor(change, levels = change.factor))) +
      geom_point(alpha=0.5, size=1.75) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
      geom_vline(xintercept = c(-log2(FC_cutoff), log2(FC_cutoff)), linetype = "dashed") +
      scale_colour_manual(name='GeneExpression', values = c('Up'='#AE3121','Down'="#3182BD", 'Stable'='lightgrey')) +
      xlab(paste0("Log fold change (FC_cutoff=", FC_cutoff,
                  '; p.adj_cutoff=', padj_cutoff, ')')) + 
      ylab("-log10 (p.adj value)") +
      ggtitle(paste0('DEGs in ',cancer_type,' (Normal=', sample.cohort[1],'; Tumour=', sample.cohort[2],')'), 
              subtitle = paste0('(Up genes=', nrow(DES_output[DES_output$change =='Up',]),
                                '; Down genes=', nrow(DES_output[DES_output$change =='Down',]),
                                ')(method=', method,')')) +
      theme_bw()+ theme_hz
    
  }
  
  print(g)
}

# run DE analyses: 
DE.outputs_list <- lapply(cancer.types, function(ca){
  print(ca)
  exprCounts <- exprCounts_list[[ca]]
  # define cohort for DE: 
  patients <- colnames(exprCounts[["TumorNormalMatch"]][["primary.tumor"]])
  patients[1:10]
  
  if(filter.BRCA==TRUE & ca=='BRCA'){
    patients <- patients[patients %in% subtype.brca.pts]
  }
  
  sample_id_df <- exprCounts$sample.df %>% dplyr::filter(.$PatientID %in% patients)
  samples.ID <- unique(sample_id_df$sampleID)
  exprCounts_df <- exprCounts$RNASeq2_count[, colnames(exprCounts$RNASeq2_count)%in%
                                              samples.ID]
  geneList_chr <- rownames(exprCounts_df)
  
  # run DE functions:
  {
    DEG.DESeq2_list <- DEG_DESeq2_function(gene_counts_df = exprCounts_df,
                                           sample_id_df = sample_id_df,
                                           geneList = geneList_chr,
                                           group.levels=c('tumour','normal'),
                                           Cutoff.per.count = 1, 
                                           Cutoff.fract = 0.2, 
                                           FC_cutoff = 2, padj_cutoff = 0.05,
                                           qvalue_cutoff = 0.01)
    # draw_volcano:
    p.vol.deseq <- draw_volcano(method = c('DESeq2','limma')[1],
                                DES_output = DEG.DESeq2_list[[1]],
                                sample.cohort = DEG.DESeq2_list[[3]],
                                change.factor=c("Up","Down","Stable"),
                                FC_cutoff= 2, padj_cutoff=0.05,
                                cancer_type= ca)
    p.vol.deseq
    
    # limma:
    DEG.limma_list <- DEG_limma.voom_function(gene_counts_df = exprCounts_df ,
                                              sample_id_df = sample_id_df,
                                              geneList = geneList_chr,
                                              group.levels=c('tumour','normal'),
                                              Cutoff.per.count = 1, 
                                              Cutoff.fract = 0.2, 
                                              FC_cutoff= 2, padj_cutoff=0.05)
    # draw_volcano:
    p.vol.limma <- draw_volcano(method = c('DESeq2','limma')[2],
                                DES_output = DEG.limma_list[[1]], 
                                sample.cohort = DEG.limma_list[[3]],
                                change.factor=c("Up","Down","Stable"),
                                FC_cutoff= 2, padj_cutoff=0.05,
                                cancer_type= ca)
    p.vol.limma
    
    ## save data: 
    DESeq2_out <- DEG.DESeq2_list[[1]] 
    colnames(DESeq2_out)[2:length(colnames(DESeq2_out))] <- paste('DESeq2', colnames(DESeq2_out)[2:length(colnames(DESeq2_out))], sep = '_')
    
    DEG.limma_list[[1]][1:2,]
    limma_out <- DEG.limma_list[[1]] %>%
      rownames_to_column(var = 'gene_id') ## NOTE: limma--> use logFC, not log2FC!!!
    colnames(limma_out)[2:length(colnames(limma_out))] <- paste('limma', colnames(limma_out)[2:length(colnames(limma_out))], sep = '_')
    
    DEG_output.comb <- DESeq2_out %>%
      left_join(limma_out, by=c('gene_id')) %>%
      mutate(change = if_else(.$DESeq2_change==.$limma_change, .$DESeq2_change, 'discordant'))
  }
  
  # export data & vol plots:
  list <- list(DEG.DESeq2_list, DEG.limma_list, DEG_output.comb, sample_id_df)
  names(list) <- c('DEG.DESeq2_list', 'DEG.limma_list', 'DEG_output.comb', 'sample_id_df')
  
  subset.title <- paste0(paste0('matched.', substr(use.matched, 1,3), '_by.', matched.type, '_', ca),
                         '.filter.', substr('FALSE', 1,3))
  saveRDS(list, file = paste0(output_dir,'step12_DEGs_DESeq2.Limma_',subset.title, '.RDS') )
  
  pdf(paste0(output_dir, 'volcano_DEGs.TCGA_',subset.title, '.pdf'), width = 8, height = 5)
  plot(p.vol.deseq)
  plot(p.vol.limma)
  dev.off()
  
  # return:
  return(list)
})
names(DE.outputs_list) <- cancer.types

## Step 3. Run bootstrapping analyses using log2.fold.change of each gene in BRCA, LUAD and LUSC: ======
DEG.outputs.TCGA_list <- DE.outputs_list
names(DEG.outputs.TCGA_list) <- cancer.types

DE.package <- c('DESeq2', 'limma')[1]
FC_cutoff <- 2 
P.value_cutoff <- 0.05
DEG.TCGA.list.used <- names(DEG.outputs.TCGA_list[[1]])[3]

DEG.TCGA_filter <- lapply(cancer.types, function(x){
  sub.deg <- DEG.outputs.TCGA_list[[x]][[ DEG.TCGA.list.used ]] %>%
    dplyr::select(c('gene_id', grep(DE.package, colnames(.), value = T) )) %>%
    dplyr::select(c('gene_id', grep('log|value|adj', colnames(.), value = T) ))
  sub.deg[1:2,]
  colnames(sub.deg) <- c('gene_id','log2.FC','p.value','p.adj')
  sub.deg <- sub.deg %>%
    mutate(DEG = if_else(.$p.adj < P.value_cutoff & abs(.$log2.FC) >= log2(FC_cutoff),
                         if_else(.$log2.FC >= log2(FC_cutoff),'Up','Down'),'Stable')) %>%
    dplyr::select(c('gene_id', 'DEG'))
  colnames(sub.deg)[2] <- paste0(x, '_DEG')
  return(sub.deg)
})
DEG.TCGA_filter <- Reduce(full_join, DEG.TCGA_filter) %>% dplyr::rename(gene = gene_id)

# load ART.genes:
ARTgenes.extreme_list <- readRDS(paste0(data_dir, 'ARTgenes.per.cell_list.RDS'))

# shared ART:
ARTgenes.consist_type_list <- lapply(cancer.types, function(ca){
  print(ca)
  # genes with shared ART or consistent ART in BRCA, LUAD, LUSC (saved in data_dir)
  data <- fread(paste0(data_dir, ca,"_sharedARTgenes.txt"))
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

## bootstrapping function:
boot.log2FC_per.cancer <- function(genes_test.ART=ARTgenes,
                                   ref.genes.df = ref.genes.df,
                                   y.axis = c('log2.FC')[1],
                                   column.ART = c('timing')[1],
                                   column.ref = c('refTiming')[1],
                                   niter = 10000,
                                   timing.ART = c('earlier','later'),
                                   timing.ref = c('early','late'),
                                   replace_or.not = c(TRUE, FALSE)[1],
                                   remove.ART.in.ref = c(TRUE, FALSE)[2]
){
  # process data:
  colnames(genes_test.ART)[match(column.ART, colnames(genes_test.ART))] <- 'column.ART'
  colnames(ref.genes.df)[match(column.ref, colnames(ref.genes.df))] <- 'column.ref'
  
  colnames(genes_test.ART)[match(y.axis, colnames(genes_test.ART))] <- 'y.axis'
  colnames(ref.genes.df)[match(y.axis, colnames(ref.genes.df))] <- 'y.axis'
  
  # bootstrapping: ===> Use cancer log2FC!
  boot.per.cancer_list <- lapply(timing.ART, function(t){
    # t=timing.ART[1]
    print(t)
    ref.timing <- if_else(t=='earlier', 'late', 'early')
    set.seed(123)
    test.genes <- genes_test.ART %>% filter(column.ART==t)
    if(remove.ART.in.ref==TRUE){
      ref.genes <- ref.genes.df %>% filter(column.ref==ref.timing) %>%
        filter(!gene_name%in%unique(test.genes$gene_name))
    }else if(remove.ART.in.ref==FALSE){
      ref.genes <- ref.genes.df %>% filter(column.ref==ref.timing)
    }
    n.size <- nrow(test.genes)
    cancer.type <- test.genes$cancer.type[1]
    
    test.genes[1:2,]
    ref.genes[1:2,]
    freq.original <- mean(test.genes$y.axis)
    freq.random_ls <- lapply(1:niter, function(i){
      index <- sample(1:nrow(ref.genes), size = n.size, replace = replace_or.not) # !!
      bg.sub <- ref.genes[index, ]
      random.mean.log2FC <- mean(bg.sub$y.axis)
      frac <- data.frame( matrix(random.mean.log2FC, nrow = 1, dimnames = list(i, cancer.type)) )
      return(frac)
    })
    freq.random <- Reduce(rbind, freq.random_ls)
    freq.random[1:6,]
    mean_frac <- colMeans(freq.random)
    mean_frac
    quantiles_boot <- apply(freq.random, 2, function(x) quantile(x, c(0.025, 0.975)))
    quantiles_boot
    # combine df:
    out <- data.frame(cancer.type=cancer.type, timing=t, refTiming=ref.timing,
                      comparisons = paste(t, ref.timing, sep = ':'),
                      n.test.genes = n.size,
                      mean_original = mean(test.genes$y.axis),
                      sd.original = sd(mean(test.genes$y.axis)),
                      median.original = median(test.genes$y.axis),
                      mean_iter = as.numeric(mean_frac), 
                      lowerCI = quantiles_boot[1,], 
                      upperCI = quantiles_boot[2,])
    out
  })
  combine <- Reduce(rbind, boot.per.cancer_list)
  return(combine)
}

# define expressed genes in TCGA:
TCGA.tumour_expressed.genes <- lapply(c('LUAD','LUSC','BRCA'), function(cancer_type){
  print(cancer_type)
  exprCounts <- TCGA.download[[cancer_type]]
  names(exprCounts)
  sample.df <- exprCounts[['sample.df']]
  sample.df[1:2,]
  if(cancer_type=='BRCA'){
    patients <- subtype.brca.filter$PatientID %>% unique()
  }else{
    patients <- sample.df$PatientID %>% unique()
  }
  sample.df <- sample.df %>% filter(PatientID%in%patients)
  
  samp.normal <- sample.df%>%filter(tissueType!='tumour')%>%pull(sampleID)%>%unique()
  samp.tumour <- sample.df%>%filter(tissueType=='tumour')%>%
    pull(sampleID)%>%unique()

  expr.all <- exprCounts[["RNASeq2_count"]]
  expr.tumour <- expr.all %>% dplyr::select(samp.tumour)

  # filtering: Cutoff.per.count=1, Cutoff.fract= 20% (0.2)
  Cutoff.per.count <- 1
  Cutoff.fract <- 0.2
  # only in tumour:
  keep <- rowSums(expr.tumour > Cutoff.per.count) >= (ncol(expr.tumour)*Cutoff.fract)
  keep.expr.tumour <- expr.tumour[keep, ]
  keep <- rowSums(expr.all > Cutoff.per.count) >= (ncol(expr.all)*Cutoff.fract)
  keep.expr.all <- expr.all[keep, ]

  genes <- rownames(keep.expr.tumour)
  return(genes)
})
names(TCGA.tumour_expressed.genes) <- c('LUAD','LUSC','BRCA')

# run bootstrapping of log2FC:
ART.column.used <- c('timing')
data.sets <- c('TCGA')
DE.package <- c('DESeq2')
iteration <- 100000
y.axis_used <- c('log2.FC')

ARTgenes.DEG_comb.list <- lapply(data.sets, function(data.set){
  print(data.set)
  
  to.comb_DEG.df <- DEG.TCGA_filter
  DEG.out_filter.list <- list(BRCA = DEG.outputs.TCGA_list[[3]][['DEG_output.comb']],
                              LUAD = DEG.outputs.TCGA_list[[1]][['DEG_output.comb']],
                              LUSC = DEG.outputs.TCGA_list[[2]][['DEG_output.comb']])
  cancer.types <- names(DEG.out_filter.list)
  
  #### Process DEG results with ART: (define up/down/stable) 
  DEG.out_filter_sub <- lapply(cancer.types, function(x){
    FC_cutoff <- 2 
    P.value_cutoff <- 0.05
    print(paste0('DEG:', x))
    sub.deg <- DEG.out_filter.list[[x]]%>%
      dplyr::select(c('gene_id', grep(DE.package, colnames(.), value = T, ignore.case = T) )) %>%
      dplyr::select(c('gene_id', grep('log|value|adj', colnames(.), value = T, ignore.case = T) ))
    colnames(sub.deg) <- c('gene_id','log2.FC','p.value','p.adj')
    sub.deg <- sub.deg %>% 
      mutate(change_to.use = if_else(.$p.adj < P.value_cutoff & abs(.$log2.FC) >= log2(FC_cutoff),
                                     if_else(.$log2.FC >= log2(FC_cutoff),'Up','Down'),'Stable'),
             cancer.type = x) %>%
      rename(gene = gene_id)
    return(sub.deg)
  })
  names(DEG.out_filter_sub) <- cancer.types
  
  comb.ARTgenes.DEG_list <- lapply(cancer.types, function(cancer.type){
    print(paste0(cancer.type))
    cell.info <- cell.line_list[[cancer.type]]
    tissue.type <- if_else(cancer.type=='BRCA','breast','lung')
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[cancer.type]]
    
    # comb DEG & ART.genes data:
    ARTgenes <- ARTgenes.consist_type_list[[cancer.type]][,-c('cancerType')] %>%
      filter(gene_name %in% sub.expressed.TCGA) %>% 
      mutate(refTiming = if_else(timing%in%c('later','early'), 'early',
                                 if_else(timing%in%c('late','earlier'), 'late','other'))) %>%
      left_join(DEG.out_filter_sub[[cancer.type]], by=c('gene_name'='gene')) %>%
      dplyr::filter(!is.na(change_to.use)) %>%
      mutate(timing_DEG = paste(timing, change_to.use, sep='_'))
    ARTgenes[1:2,]
    return(ARTgenes)
  })
  names(comb.ARTgenes.DEG_list) <- cancer.types
  
  #final output:
  return(comb.ARTgenes.DEG_list)
})
names(ARTgenes.DEG_comb.list) <- data.sets

## run boot:
ref.used <- c('vs.allRT_')[1] 

boot.logFC_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  cancers <- names(ARTgenes.DEG_comb.list[[data.set]])
  
  to.comb_DEG.df <- DEG.TCGA_filter
  DEG.out_filter.list <- list(BRCA = DEG.outputs.TCGA_list[[3]][['DEG_output.comb']],
                              LUAD = DEG.outputs.TCGA_list[[1]][['DEG_output.comb']],
                              LUSC = DEG.outputs.TCGA_list[[2]][['DEG_output.comb']])
  
  DEG.out_filter_sub <- lapply(cancers, function(x){
    FC_cutoff <- 2 
    P.value_cutoff <- 0.05
    print(paste0('DEG:', x))
    if(data.set=="TRACERx"){
      sub.deg <- DEG.out_filter.list[[x]]%>%filter(coef=="Tumour - Normal")
    }else{
      sub.deg <- DEG.out_filter.list[[x]]%>%
        dplyr::select(c('gene_id', grep(DE.package, colnames(.), value = T, ignore.case = T) ))
    }
    sub.deg <- sub.deg%>%
      dplyr::select(c('gene_id', grep('log|value|adj', colnames(.), value = T, ignore.case = T) ))
    sub.deg[1:2,]
    colnames(sub.deg) <- c('gene_id','log2.FC','p.value','p.adj')
    sub.deg <- sub.deg %>% 
      mutate(change_to.use = if_else(.$p.adj < P.value_cutoff & abs(.$log2.FC) >= log2(FC_cutoff),
                                     if_else(.$log2.FC >= log2(FC_cutoff),'Up','Down'),'Stable'),
             cancer.type = x) %>%
      rename(gene = gene_id)
    return(sub.deg)
  })
  names(DEG.out_filter_sub) <- cancers
  DEG.out_filter_sub[[1]][1:2,]
  
  # run boot of logFC:
  boot.list <- lapply(cancers, function(ca){
    # ca=cancers[1]
    print(ca)
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[ca]]
    ARTgenes <- ARTgenes.DEG_comb.list[[data.set]][[ca]]
    ARTgenes[1:2,]
    
    ref.genes.df <- ARTgenes.extreme_list[[ca]] %>% 
      dplyr::select(gene_name,chr,start,stop,normal_l2r) %>% distinct() %>%
      mutate(refTiming = if_else(normal_l2r>0, 'early', 'late')) %>%
      left_join(DEG.out_filter_sub[[ca]], by=c('gene_name'='gene')) %>%
      filter(gene_name %in% sub.expressed.TCGA) %>%
      filter(!is.na(log2.FC))
    
    # run boot:
    boot.out <- boot.log2FC_per.cancer(genes_test.ART=ARTgenes,
                                       ref.genes.df = ref.genes.df,
                                       y.axis = c('log2.FC')[1],
                                       column.ART = c('timing')[1],
                                       column.ref = c('refTiming')[1],
                                       niter = iteration,
                                       timing.ART = c('earlier','later'),
                                       timing.ref = c('early','late'),
                                       replace_or.not = c(TRUE, FALSE)[1],
                                       remove.ART.in.ref = c(TRUE, FALSE)[2])
    return(boot.out)
  })
  # out:
  out <- Reduce(rbind, boot.list)
  return(out)
})
names(boot.logFC_ART_list) <- data.sets

# count % DEG among earlier/later:
frac.DEG_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  cancers <- names(ARTgenes.DEG_comb.list[[data.set]])
  # frac of DEG per cancer:
  frac.list <- lapply(cancers, function(ca){
    print(ca)
    ARTgenes <- ARTgenes.DEG_comb.list[[data.set]][[ca]]
    ARTgenes[1:2,]
    sum <- ARTgenes %>% group_by(timing) %>% summarise(total=n())
    sum
    frac <- ARTgenes %>% group_by(timing, change_to.use) %>%
      summarise(count=n()) %>% left_join(sum) %>%
      mutate(frac=count/total) %>% mutate(cancer.type=ca)
    frac
    return(frac)
  })
  # out:
  out <- Reduce(rbind, frac.list)
  return(out)
})
names(frac.DEG_ART_list) <- data.sets

# load boot.logFC data: for plotting 
# annot.color:
{
  color.key <- c(brewer.pal(9,'PiYG')[c(1,2, 8,9)])
  names(color.key) <- c("early", "earlier","later","late")
  color.key
  
  col.DEG <- c(rev(brewer.pal(8,'Set2')[c(1:2)]), 'darkgrey')
  names(col.DEG) <- c('Up','Down','Stable')
  
  theme_hz <- theme(axis.text.x = element_text(angle=90, hjust = 0.5,size = 12, color = 'black'),
                    axis.text.y = element_text(size = 12, color = 'black'),
                    axis.title = element_text(size = 14, color = 'black'),
                    legend.text = element_text(size = 10, color = 'black'),
                    legend.title = element_text(size = 12, color = 'black', face = 'bold'),
                    plot.title = element_text(size = 15, color = 'black', face = 'bold'))
}

data.set <- c("TCGA")
cancer.types <- c("BRCA", "LUAD", "LUSC")

plot.boot.logFC <- boot.logFC_ART_list[[data.set]]

p_boot_by.cancer2 <- plot.boot.logFC %>% 
  ggplot(aes(x=factor(timing, levels = (names(color.key))),# factor(cancer.type, levels = cancer.types), 
             color=timing)) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, colour=factor(refTiming, levels = names(color.key))), width = 0.1) +
  geom_point(aes(y = mean_original, colour = factor(timing, levels = names(color.key)),
                 shape = paste0('Observed(mean)_ART')), size=5) +
  geom_point(aes(y = mean_iter, colour=factor(refTiming,levels=names(color.key)), 
                 shape = paste0('Random(mean)_Reference') ), size=3) +
  geom_text(aes(y=mean_original, label=paste0('n=', n.test.genes)), 
            nudge_y=-0.03, stat = 'identity',check_overlap=T) +
  scale_colour_manual(name = paste0('rep.Timing'), values = c(color.key)) +
  scale_shape_manual(name =paste0('Data'), values = c(18,16)) +
  facet_grid(~factor(cancer.type, levels = cancer.types),scales="free") +
  xlab('') + ylab('Log2.Fold.Change in tumour') +
  theme_bw() + theme_hz +
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
        panel.grid = element_blank(),
        strip.text = element_text(size=12))
p_boot_by.cancer2

pdf(paste0(data.dir, 'p.boot.wide_logFC_DEG.matched_',use.matched, '_', ref.used, data.set,'.pdf'), 
    width = 12, height = 4)
print(p_boot_by.cancer2)
dev.off()

# pie: earlier/later per.cancer
plot.frac.DEG <- frac.DEG_ART_list[[data.set]] %>%
  mutate(change_to.use=factor(change_to.use,levels = c('Up','Down','Stable'))) %>%
  mutate(timing=factor(timing,levels = names(color.key)))%>%
  mutate(cancer.type=factor(cancer.type,levels = c(cancer.types)))

plot.frac.DEG <- plot.frac.DEG[order(plot.frac.DEG$cancer.type,
                                     plot.frac.DEG$timing,
                                     plot.frac.DEG$change_to.use),] %>%
  group_by(cancer.type, timing) %>%
  mutate(cumulative = cumsum(frac), midpoint = cumulative-frac/2)

p.pie_list <- lapply(cancer.types, function(x){
  sub.p.pie <- lapply(c('earlier','later'), function(i){
    sub <- plot.frac.DEG %>% 
      filter(timing %in% i) %>% filter(cancer.type==x) %>%
      ggplot(aes(x = 1, weights = frac, fill=factor(change_to.use,levels = names(col.DEG)))) +
      geom_bar(width = 1, position = position_stack(reverse=T), colour = 'white') +
      coord_polar(theta = "y") +
      geom_text(aes(x=c(1, 1.2, 0.9), y = midpoint,label=paste0(round(frac*100, digits=1),'%\n(n=', count,')')),
                stat = 'identity', size=4) +
      scale_fill_manual(name = paste0('Diff.Expr'), values = c(col.DEG)) +
      xlab('') + ylab('') + ggtitle(paste0(i, ' in ', x)) +
      theme_void() +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) 
    return(sub)
  })
  p <- sub.p.pie[[1]]+sub.p.pie[[2]]
  p
})
names(p.pie_list) <- cancer.types

pdf(paste0(output_dir, 'p.pie_DEG.in.ART_', data.set,'.pdf'), 
    width = 8, height = 5)
for (i in names(p.pie_list)) {
  print(p.pie_list[[i]])
}
dev.off()

#####################################################################
### Figure 3D: bootstrapping of CN/Ploidy ========
#load shared ART genes
ARTgenes.consist_type_list <- lapply(cancer.types, function(ca){
  print(ca)
  # genes with shared ART or consistent ART in BRCA, LUAD, LUSC (saved in data_dir)
  data <- fread(paste0(data_dir, ca,"_sharedARTgenes.txt"))
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

#reptTiming of genes
repTiming_genes <-  readRDS(paste0(data_dir, 'cohort_meanL2R_genes.rds'))

#calculated copy number difference to ploidy per gene for BRCA, LUAD and LUSC from ascat data provided by TCGA
#--> only lobular and ductal samples used for BRCA
gene.CN_TCGA.list <- list('BRCA' = readRDS(paste0(data_dir, 'BRCAsubset_mean_cnDiff_perGene.rds')),
                          'LUAD' = readRDS(paste0(data_dir, 'LUAD_mean_cnDiff_perGene.rds')),
                          'LUSC' = readRDS(paste0(data_dir, 'LUSC_mean_cnDiff_perGene.rds')))

# run bootstrapping tests:
iteration <- 100000
y.axis_used <- c('mean_diff')
data.sets <- c('TCGA')

ARTgenes.CNA_comb.list <- lapply(data.sets, function(data.set){
  print(data.set)
  gene.CN_used.list <- gene.CN_TCGA.list
  cancer.types <- names(gene.CN_used.list)
  
  #### Process CNA results with ART: 
  comb.ARTgenes.CNA_list <- lapply(cancer.types, function(cancer.type){
    print(paste0(cancer.type))
    cell.info <- cell.line_list[[cancer.type]]
    tissue.type <- if_else(cancer.type=='BRCA','breast','lung')
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[cancer.type]]
    sub.CNA <- gene.CN_used.list[[cancer.type]]
    sub.CNA[1:2,]
    ARTgenes.consist_type_list[[cancer.type]][1:2,]
    # comb CNA & ART.genes data:
    ARTgenes <- ARTgenes.consist_type_list[[cancer.type]][,-c('cancerType')] %>%
      mutate(refTiming = if_else(timing%in%c('later','early'), 'early',
                                 if_else(timing%in%c('late','earlier'), 'late','other'))) %>%
      left_join(sub.CNA, by=c('gene_name', 'chr',  'start',   'stop')) %>%
      mutate(cancer.type = cancer.type)
    ARTgenes[1:2,]
    return(ARTgenes)
  })
  names(comb.ARTgenes.CNA_list) <- cancer.types
  
  #final output:
  return(comb.ARTgenes.CNA_list)
})
names(ARTgenes.CNA_comb.list) <- data.sets

## run bootstrapping:
ref.used <- c('vs.allRT_')
genes.used <- c('expressed')[1]

boot.CN.ploidy_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  gene.CN_used.list <- gene.CN_TCGA.list
  cancers <- names(gene.CN_used.list)
  cancers
  
  #### Process CNA results with ART: 
  # run boot of CN.ploidy:
  boot.list <- lapply(cancers, function(ca){
    print(ca)
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[ca]]
    ARTgenes <- ARTgenes.CNA_comb.list[[data.set]][[ca]] %>% filter(!is.na(mean_diff))
    sub.CNA <- gene.CN_used.list[[ca]]
    
    ARTgenes.extreme_list[[ca]][1:2,]
    ref.genes.df <- ARTgenes.extreme_list[[ca]] %>% 
      dplyr::select(gene_name,chr,start,stop,normal_l2r) %>% distinct() %>%
      mutate(refTiming = if_else(normal_l2r>0, 'early', 'late')) %>%
      left_join(sub.CNA, by=c('gene_name', 'chr', 'start', 'stop'))
    ref.genes.df[1:2,]
    
    if(genes.used=="expressed"){
      ref.genes.df <- ref.genes.df %>% filter(gene_name %in% sub.expressed.TCGA)
      ARTgenes <- ARTgenes %>% filter(gene_name %in% sub.expressed.TCGA)
    }
    
    # run boot:
    boot.out <- boot.log2FC_per.cancer(genes_test.ART=ARTgenes,
                                       ref.genes.df = ref.genes.df,
                                       y.axis = c('mean_diff')[1],
                                       column.ART = c('timing')[1],
                                       column.ref = c('refTiming')[1],
                                       niter = iteration,
                                       timing.ART = c('earlier','later'),
                                       timing.ref = c('early','late'),
                                       replace_or.not = c(TRUE, FALSE)[1],
                                       remove.ART.in.ref = c(TRUE, FALSE)[2])
    return(boot.out)
  })
  # out:
  out <- Reduce(rbind, boot.list)
  return(out)
})
names(boot.CN.ploidy_ART_list) <- data.sets

# Make plots:
data.set <- c("TCGA")
cancer.types <- c("BRCA", "LUAD", "LUSC")
plot.boot.CN.ploidy <- boot.CN.ploidy_ART_list[[data.set]]

p_boot_by.cancer2 <- plot.boot.CN.ploidy %>% 
  ggplot(aes(x=factor(timing, levels = (names(color.key))),# factor(cancer.type, levels = cancer.types), 
             color=timing)) +
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, colour=factor(refTiming, levels = names(color.key))), width = 0.1) +
  geom_point(aes(y = mean_original, colour = factor(timing, levels = names(color.key)),
                 shape = paste0('Observed(mean)_ART')), size=5) +
  geom_point(aes(y = mean_iter, colour=factor(refTiming,levels=names(color.key)), 
                 shape = paste0('Random(mean)_Reference') ), size=3) +
  geom_text(aes(y=mean_original, label=paste0('n=', n.test.genes)), 
            nudge_y=-0.004, stat = 'identity',check_overlap=T) +
  scale_colour_manual(name = paste0('rep.Timing'), values = c(color.key)) +
  scale_shape_manual(name =paste0('Data'), values = c(18,16)) +
  facet_grid(~factor(cancer.type, levels = cancer.types),scales="free") +
  xlab('') + ylab('Mean CN/Ploidy in tumour') +
  theme_bw() + theme_hz+
  theme(axis.text.x = element_text(angle=0,hjust=0.5,vjust=0.5),
        panel.grid = element_blank(),
        strip.text = element_text(size=12))

p_boot_by.cancer2

pdf(paste0(output_dir,'p.boot.wide_CN.ploidy-1_', genes.used,'_',ref.used, data.set,'_BRCAsubtype.pdf'), 
    width = 12, height = 4)
print(p_boot_by.cancer2)
dev.off()

#####################################################################
### Figure 3E: conserved RT (in all 10 normal) ========
tissue_info <- read.table(paste0(data_dir, 'tissueInfo_cellLines_20210309.tsv'), header = T, sep = '\t')
tissue_info$cellLine <- sub('A549$', 'A549encode', tissue_info$cellLine)
tissue_info$cellLine <- sub('A549rep', 'A549', tissue_info$cellLine)

#load log2-ratios
log2ratio_df        <- readRDS(paste0(data_dir,'cohort_50kb_l2r.rds'))
log2ratio_df        <- log2ratio_df[,grep('LTX', colnames(log2ratio_df), invert = T)]
normal_cellLines    <- colnames(log2ratio_df)[colnames(log2ratio_df) %in% tissue_info$cellLine[tissue_info$tissueType == 'Normal']]
normal_log2ratio_df <- log2ratio_df[,c('chr', 'start', 'stop', normal_cellLines)]

#load log2-ratios per gene
log2ratio_genes        <- readRDS(paste0(data_dir, 'cohort_meanL2R_genes.rds') )
log2ratio_genes        <- log2ratio_genes[,grep('LTX', colnames(log2ratio_genes), invert = T)]
normal_cellLines       <- colnames(log2ratio_genes)[colnames(log2ratio_genes) %in% tissue_info$cellLine[tissue_info$tissueType == 'Normal']]
normal_log2ratio_genes <- log2ratio_genes[,c('chr', 'start', 'stop', 'gene_name', normal_cellLines)]

#identify consereved replication timing regions (with NA values allowed)
normal_log2ratio_df$repTiming           <- 'unknown'
ww_nonConserved                         <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(sign(x[!is.na(x)]) > 0) > 0 & sum(sign(x[!is.na(x)]) < 0) > 0)
normal_log2ratio_df$repTiming[ww_nonConserved] <- 'non conserved'
ww_early                                <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) > 0))
normal_log2ratio_df$repTiming[ww_early] <- 'early'
ww_late                                 <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) < 0))
normal_log2ratio_df$repTiming[ww_late]  <- 'late'

#compare consereved RT in IN-HOUSE and ENCODE data
normal_ENCODE <- normal_cellLines[normal_cellLines %in% tissue_info$cellLine[tissue_info$RepliSeq_dataset == 'ENCODE']]
ENCODE_normal <- normal_log2ratio_df[, c('chr', 'start', 'stop', normal_ENCODE)]
ENCODE_normal$repTiming           <- 'unknown'
ww_nonConserved                   <- apply(ENCODE_normal[,normal_ENCODE], 1, function(x) sum(sign(x[!is.na(x)]) > 0) > 0 & sum(sign(x[!is.na(x)]) < 0) > 0)
ENCODE_normal$repTiming[ww_nonConserved] <- 'non conserved'
ww_early                          <- apply(ENCODE_normal[,normal_ENCODE], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) > 0))
ENCODE_normal$repTiming[ww_early] <- 'early'
ww_late                           <- apply(ENCODE_normal[,normal_ENCODE], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) < 0))
ENCODE_normal$repTiming[ww_late]  <- 'late'

normal_INHOUSE <- normal_cellLines[normal_cellLines %in% tissue_info$cellLine[tissue_info$RepliSeq_dataset == 'Haoran']]
INHOUSE_normal <- normal_log2ratio_df[, c('chr', 'start', 'stop', normal_INHOUSE)]
INHOUSE_normal$repTiming           <- 'unknown'
ww_nonConserved                    <- apply(INHOUSE_normal[,normal_INHOUSE], 1, function(x) sum(sign(x[!is.na(x)]) > 0) > 0 & sum(sign(x[!is.na(x)]) < 0) > 0)
INHOUSE_normal$repTiming[ww_nonConserved] <- 'non conserved'
ww_early                           <- apply(INHOUSE_normal[,normal_INHOUSE], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) > 0))
INHOUSE_normal$repTiming[ww_early] <- 'early'
ww_late                            <- apply(INHOUSE_normal[,normal_INHOUSE], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) < 0))
INHOUSE_normal$repTiming[ww_late]  <- 'late'

exclude_unknowns <- ENCODE_normal$repTiming == 'unknown' | INHOUSE_normal$repTiming == 'unknown'
ENCODE_normal    <- ENCODE_normal[!exclude_unknowns,]
INHOUSE_normal   <- INHOUSE_normal[!exclude_unknowns,]

plot_data <- ENCODE_normal[,c('chr', 'start', 'repTiming')] %>%
  left_join(INHOUSE_normal[,c('chr', 'start', 'repTiming')], by = c('chr', 'start'))
colnames(plot_data) <- c('chr', 'start', 'ENCODE', 'INHOUSE')
plot_data <- plot_data %>%
  group_by(ENCODE, INHOUSE) %>%
  summarise(count = n(),
            freq = count / nrow(ENCODE_normal),
            perc = freq * 100)

pdf(paste0(output_dir, 'alluvial_conservedRT_ENCODE_INHOUSE.pdf'), width = 5, height = 3)
ggplot(plot_data, aes(y = perc, axis1 = INHOUSE, axis2 = ENCODE)) +
  geom_alluvium(aes(fill = ENCODE), width = 1/12) +
  geom_stratum(width = 1/12, aes(fill = ENCODE), color = "black") +
  scale_fill_manual(name = '', values = c('early' = '#c51b7d', 'late' = '#4d9221', 'non conserved' = '#bababa')) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),colour = c('black', 'white', 'white', 'black', 'white', 'white')) +
  scale_x_discrete(limits = c("ENCODE", "In-HOUSE"), expand = c(.05, .05)) +
  theme_classic() + theme(legend.position = 'none') +
  scale_y_continuous(expand = c(0,0)) +
  coord_flip() +
  ylab('% genome') +
  ggtitle("Conserved RT in normal cell-lines")
dev.off()

#####################################################################
### Figure 3F: ART per cancer type being conserved RT or not (SuppFig 7C: conserved RT in cancer ========
## Conserved genes in all normal 
normal_log2ratio_df[1:2, ]
normal_log2ratio_df$repTiming <- 'unknown'
ww_nonConserved <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(sign(x[!is.na(x)]) > 0) > 0 & sum(sign(x[!is.na(x)]) < 0) > 0)
normal_log2ratio_df$repTiming[ww_nonConserved] <- 'non.conserved'

ww_early <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) > 0))
normal_log2ratio_df$repTiming[ww_early] <- 'early'
ww_late <- apply(normal_log2ratio_df[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) < 0))
normal_log2ratio_df$repTiming[ww_late] <- 'late'

# annot.col:
annot.col <- list(conserved.RT = c(brewer.pal(9, 'PiYG')[c(1,9)], 'darkgrey'), 
                  RT = c(brewer.pal(11, 'PiYG')[c(2,3,9,10)], brewer.pal(11, 'RdGy')[9]) )
names(annot.col$conserved.RT) <- c("early", "late", "non.conserved") 
names(annot.col$RT) <- c("early","earlier", "later", "late", "unique ART")

# make plots:
plot_data.conserved <- normal_log2ratio_df %>% filter(repTiming != 'unknown') %>%
  group_by(repTiming) %>%
  summarise(count = n()) %>%
  mutate(fraction = count / sum(count) * 100) %>%
  mutate(cumulative = cumsum(fraction),
         midpoint = cumulative - fraction / 2,
         label = paste0(round(fraction,1), '%'))

# shared ART:
cancer.types <- c("BRCA", "LUAD", "LUSC")
ARTgenes.shared_list <- lapply(cancer.types, function(ca){
  print(ca)
  data <- fread(paste0(data_dir, ca,"_sharedARTgenes.txt"))
  return(data)
})
names(ARTgenes.shared_list) <- cancer.types

ARTregions.shared_list <- lapply(cancer.types, function(ca){
  print(ca)
  data <- fread(paste0(data_dir, ca,"_sharedARTregions.txt"))
  return(data)
})
names(ARTregions.shared_list) <- cancer.types

# process data:
compare_repTiming_raw <- normal_log2ratio_df[,c('chr', 'start', 'stop', 'repTiming')]
compare_repTiming_raw[1:2,]

compare_repTiming <- compare_repTiming_raw %>% 
  filter(repTiming != 'unknown') %>%
  left_join(ARTregions.shared_list[[1]][,1:4], by = colnames(ARTregions.shared_list[[1]])[1:3]) %>% 
  left_join(ARTregions.shared_list[[2]][,1:4], by = colnames(ARTregions.shared_list[[2]])[1:3]) %>% 
  left_join(ARTregions.shared_list[[3]][,1:4], by = colnames(ARTregions.shared_list[[3]])[1:3])
compare_repTiming[1:2,]
colnames(compare_repTiming) <- c('chr', 'start', 'stop', 'consistent', cancer.types)

# make plot: among conserved.RT
plot_data <- compare_repTiming %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop', 'consistent')) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(value=if_else(is.na(value), 'unique ART', value)) %>% 
  rename(timing=value, cancer.type=variable)

plot_data.used <- plot_data %>% 
  group_by(consistent, cancer.type, timing) %>%
  summarise(count = n()) %>% as.data.frame() %>%
  group_by(consistent, cancer.type) %>% 
  mutate(total = sum(count), fraction = count/total) %>%
  mutate(consistent=if_else(consistent=="non.conserved", consistent, 
                            paste0('Conserved ', consistent, ' in all normal')))

# plotting:
annot.labels <- c("early"='EarlyN+T', "earlier"='LateN-to-EarlyT', 
                  "later"='EarlyN-to-LateT',"late"='LateN+T', "unique ART"='Unique ART')

p_in.conserved.RT <- plot_data.used %>% filter(consistent!="non.conserved") %>%
  ggplot(aes(x=factor(cancer.type,levels = cancer.types), y = fraction, 
             fill = factor(timing,levels = names(annot.col$RT)))) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label=paste0(round(fraction*100, 1),'%'), colour=timing),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = 'repTiming in cancer', values = annot.col$RT,
                    labels=c(annot.labels)) +
  scale_color_manual(name='repTiming in cancer', 
                     values = c("early"='grey', "earlier"='black', "later"='black',"late"='grey', "unique ART"='grey'),
                     labels=c(annot.labels)) +
  facet_grid(. ~ consistent) +
  scale_y_continuous(expand = c(0,0.02)) +
  xlab('') + ylab('% timings in conserved RT regions') +
  theme_bw() + theme_hz +
  theme(axis.text.x = element_text(angle = 0), strip.text.x = element_text(size=12),
        legend.position = 'bottom',
        panel.grid = element_blank())
p_in.conserved.RT

# make plot: among ART per.cancer
compare_repTiming[1:2,]
plot_data2 <- compare_repTiming %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop', 'consistent')) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(value=if_else(is.na(value), 'unique ART', value)) %>% 
  rename(timing=value, cancer.type=variable)
plot_data2[1:4,]

plot_data.used2 <- plot_data2 %>% 
  group_by(consistent, cancer.type, timing) %>%
  summarise(count = n()) %>% as.data.frame() %>%
  group_by(cancer.type, timing) %>% 
  mutate(total = sum(count), fraction = count/total)
plot_data.used2[1:4,]

p_in.ART <- plot_data.used2 %>% filter(timing %in% c("earlier", "later") ) %>%
  ggplot(aes(x=factor(cancer.type,levels = cancer.types), y = fraction, 
             fill = factor(consistent,levels = names(annot.col$conserved.RT)))) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label=paste0(round(fraction*100, 1),'%'), colour=consistent),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(name = 'repTiming in cancer', values = annot.col$conserved.RT,
                    labels = c("early"='Conserved early', "late"='Conserved late', "non.conserved"='Non-conserved')) +
  scale_color_manual(name='repTiming in cancer', 
                     values = c("early"='grey', "late"='grey', "non.conserved"='black'),
                     labels = c("early"='Conserved early', "late"='Conserved late', "non.conserved"='Non-conserved')) +
  facet_grid(. ~ factor(timing, levels = c('earlier', 'later'))) +
  scale_y_continuous(expand = c(0,0.02)) +
  xlab('') + ylab('% timings in conserved RT regions') +
  theme_bw() + theme_hz +
  theme(axis.text.x = element_text(angle = 0), strip.text.x = element_text(size=12),
        legend.position = 'bottom',
        panel.grid = element_blank())
p_in.ART

# Figure 3F:
pdf(paste0(output_dir, 'bar_frac.conserRT_ART.pdf'), width=6.5, height = 5)
print(p_in.ART)
dev.off()

# SuppFig 7C: 
pdf(paste0(output_dir, 'bar_frac.ART_conserved.RT.pdf'), width=7, height = 5)
print(p_in.conserved.RT)
dev.off()

#####################################################################
### SuppFig 7D: essential genes being conserved early RT ========
#classify genes as essential or not
lung_essential   <- as.data.frame(readxl::read_xlsx(paste0(data_dir, '20220203_Essential_lung.xlsx')))
lung_essential   <- unique(lung_essential$`Url Label`)
breast_essential <- as.data.frame(readxl::read_xlsx(paste0(data_dir, '20220203_Essential_breast.xlsx')))
breast_essential <- unique(breast_essential$`Url Label`)

normal_log2ratio_genes$lung_essential   <- ifelse(normal_log2ratio_genes$gene_name %in% lung_essential, 'essential', 'non essential')
normal_log2ratio_genes$breast_essential <- ifelse(normal_log2ratio_genes$gene_name %in% breast_essential, 'essential', 'non essential')

plot_data <- reshape2::melt(normal_log2ratio_genes[,c('chr', 'start', 'stop', 'gene_name', 'repTiming', 'lung_essential', 'breast_essential')], id.vars = c('chr', 'start', 'stop', 'gene_name', 'repTiming'))
plot_data <- plot_data[plot_data$value %in% 'essential',]
plot_data <- plot_data %>%
  filter(repTiming != 'unknown') %>%
  group_by(repTiming, variable, value) %>%
  summarise(count = n()) %>%
  group_by(variable, value) %>%
  mutate(total = sum(count),
         fraction = count / total)
plot_data$variable <- sub('_essential', '', plot_data$variable)

pdf(paste0(output_dir, 'bar_fractionEssentialGenes_conservedNormalRT.pdf'), width = 4, height = 3)
ggplot(plot_data, aes(x = variable, y = fraction, fill = repTiming)) + 
  geom_bar(stat = 'identity', colour = 'white') + 
  scale_fill_manual(name = 'Conserved Timing', values = c('early' = '#c51b7d', 'late' = '#4d9221', 'non conserved' = '#bababa')) +
  scale_y_continuous(expand = c(0,0)) +
  xlab('') + ylab('Proportion Essential Genes') +
  theme_bw()
dev.off()



















