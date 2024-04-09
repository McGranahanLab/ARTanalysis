
#####################################################################
## Figure 6 & Figure 2F. Gene based analyses (& Supplementary Figure 11)
#####################################################################
# written by Haoran Zhai (haoran.zhai.17@ucl.ac.uk) and run in R version 4.0.2

# Description:
# Scripts to create Figure 6 and Supplementary Figure 11 in the manuscript named "Replication timing alterations impact mutation acquisition during tumour evolution".
# Data accessibility statement can be found in the manuscript.

#libraries
options(stringsAsFactors = F)
library(CNTools)
library(TCGA2STAT)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
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
library(tibble)
library(data.table)
library(fst)
library(readxl)
library(writexl)
library(tidyverse)
library(dplyr)

# Parameters =======
data_dir   <- './data/' #set full path to the directory where the data for this analysis has been saved
output_dir <- './output/' #set full path to the directory where the results for this analysis should be saved

cancer.types <- c("BRCA", "LUAD")
cancer_cellLines <- list('BRCA' = c('SK-BR3', 'MCF-7', 'T47D', 'MDA453'),
                         'LUAD' = c('H1650', 'H1792', 'H2009', 'A549'))
normal_cellLines <- c('BRCA' = 'HMEC', 'LUAD' = 'T2P')
chr_to_use <- paste0('chr', c(1:22))

###############################################################################################################################
# Figure 2F:
###############################################################################################################################

overlap.ART_list <- readRDS(paste0(data_dir,'overlappingART.rds'))

# find gene.prop per cancer: 
overlap.ART_list[[1]][1:2,]
overlap.geneProp.ART_list <- lapply(cancer.types, function(cancer){
  # cancer=cancer.types[1]
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

saveRDS(overlap.geneProp.ART_list, file = paste0(output_dir, 'overlapping_ARTregions_per.cancer.rds'))

# read data & make gene.prop plot: --------
overlap.geneProp.ART_list <- readRDS(paste0(output_dir, 'overlapping_ARTregions_per.cancer.rds'))

plot_data_raw <- data.frame()
for (i in names(overlap.geneProp.ART_list)) {
  dt <- overlap.geneProp.ART_list[[i]] %>% 
    dplyr::select(timing, geneOverlapFrac, class, cancerType)
  plot_data_raw <- rbind(plot_data_raw, dt)
}
plot_data_raw[1:2,]

# annot.col:
{
  display.brewer.all()
  display.brewer.pal(11, 'PiYG')
  annot.col <- list(RT = c(brewer.pal(11, 'PiYG')[c(2,9,3,10)]))
  names(annot.col$RT) <- c("early","later","earlier", "late")
  
  theme_hz <- theme(axis.text.x = element_text(angle=90, hjust = 0.5,size = 12, color = 'black'),
                    axis.text.y = element_text(size = 12, color = 'black'),
                    axis.title = element_text(size = 14, color = 'black'),
                    legend.text = element_text(size = 10, color = 'black'),
                    legend.title = element_text(size = 12, color = 'black', face = 'bold'),
                    plot.title = element_text(size = 15, color = 'black', face = 'bold'))
}

# plot:
plot_data <- plot_data_raw %>% 
  filter(class %in% c('not_altered', 'shared', 'recurrent'))%>% 
  filter(cancerType%in%c('BRCA','LUAD'))

n.gene <- xtabs(~timing+cancerType, plot_data) %>% as.data.frame()
n.gene[1:2,]

comparisons <- list(c('early', 'later'),c('earlier', 'later'),c('earlier', 'late'),c('early', 'late')) #
stat.test <- c("t.test",'wilcox')[1]

p.violin_gene.density <- plot_data %>% 
  ggplot(aes(x=factor(timing,levels=c(names(annot.col$RT))), y = geneOverlapFrac, 
             fill=factor(timing,levels=names(annot.col$RT)))) + 
  # geom_quasirandom(aes(colour=factor(timing,levels=names(annot.col$RT))), alpha = 0.5) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, alpha = 0.2) +
  scale_fill_manual(name = 'Replication Timing', values = annot.col$RT) +
  scale_colour_manual(name = 'Replication Timing', values = annot.col$RT) +
  stat_compare_means(method = stat.test, label = 'p.format',
                     comparisons = comparisons) +
  scale_x_discrete(labels = c(paste0('EarlyN+T\n(N=)'), 
                              paste0('EarlyN-to-LateT\n(N=)'), 
                              paste0('LateN-to-EarlyT\n(N=)'), 
                              paste0('LateN+T\n(N=)') )) +
  facet_grid(.~factor(cancerType, levels = cancer.types)) +
  xlab('') + ylab('Gene proprotion within 50kb bins') +
  theme_bw() + theme_hz +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12),legend.position = 'none',
        strip.text.x = element_text(size = 12),
        panel.grid = element_blank()) #  
p.violin_gene.density

# Perform a stat test between groups
stat.test1 <- compare_means(geneOverlapFrac ~ timing, 
                            data = plot_data[plot_data$cancerType==cancer.types[1],], 
                            method = stat.test) %>% mutate(cancertype =cancer.types[1])
stat.test1
stat.test2 <- compare_means(geneOverlapFrac ~ timing, 
                            data = plot_data[plot_data$cancerType==cancer.types[2],], 
                            method = stat.test) %>% mutate(cancertype =cancer.types[2])
stat.test2
stat.test <- rbind(stat.test1, stat.test2)

pdf(paste0(output_dir, 'geneProp.RT_violin.per.cancer.pdf'), width = 7, height = 5, useDingbats = F)
print(p.violin_gene.density)
dev.off()

write.csv(stat.test, paste0(output_dir, 'stat.test_gene.proport.csv'))


###############################################################################################################################
# Figure 6A:
###############################################################################################################################

## Download & Process TCGA
#####################################
library('CNTools')
library(TCGA2STAT)
library(tidyverse)

cancer.types <- c('BRCA','LUAD')

## Function ========================
# function to download count data (getTCGA function)
func.tcga_download_hz <- function(cancer_type="LUAD", tumour.idx = paste0(0, c(1:9)),
                                  data.type = c("count", "RPKM")[1],
                                  normal.idx = paste0(1, c(0:9)),
                                  control.idx = paste0(2, c(0:9))){
  # use TCGA2STAT function to download count data, with 15-digit TCGA barcodes
  # help('getTCGA')
  # For RNA-Sequencing data from the second analysis pipeline (RNASeqV2), 
  # our package provides importation of data from Illumina HiSeq only as this is the most common profiling method. 
  # The values returned are the RSEM values for the genes, 
  # which usually contains 20501 genes, with RSEM values ranging from 0 to 10^6.
  library(TCGA2STAT)
  data <- getTCGA(disease=cancer_type, data.type="RNASeq2", 
                  type= data.type , clinical = TRUE)
  # extract count data:
  count_data <- as.data.frame(data$dat)
  colnames(count_data) <- substr(colnames(count_data), 1, 15)
  colnames(count_data) <- paste(cancer_type, colnames(count_data), sep="-")
  # # sample.data:
  # # # specify tumour/normal barcodes: [4]
  # # Tumor types range from 01 - 09, normal types from 10 - 19 and 
  # # control samples from 20 - 29. # See Code Tables Report for a complete list of sample codes
  count_data[1:2, 1:10]
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
TCGA.download <- lapply(cancer.types, function(cancer_type){
  print(cancer_type)
  data.type.x <- c("count", "RPKM")[1]
  tcga <- func.tcga_download_hz(cancer_type = cancer_type, 
                                data.type = data.type.x)
  saveRDS(object = tcga, 
          file = paste0(output_dir, cancer_type,'_RNAseq2.',data.type.x,'_clin.survival',".rds"))
  #
  data.type.y <- c("count", "RPKM")[2]
  tcga <- func.tcga_download_hz(cancer_type = cancer_type, 
                                data.type = data.type.y)
  saveRDS(object = tcga, 
          file = paste0(output_dir, cancer_type,'_RNAseq2.',data.type.y,'_clin.survival',".rds"))
  
})

# expressed genes in each cancer type in TCGA ======
# only keep ductal and lobular BRCA:====
# using Data.S1 in the published paper here with histological annotation of BRCA in TCGA: 
# https://www.cell.com/cell-genomics/fulltext/S2666-979X(21)00083-5#supplementaryMaterial
subtype.brca <- read_xlsx(paste0(data_dir,'mmc2_Histologic.BRCA.inTCGA.xlsx')) %>%
  as.data.frame()
subtype.brca[1:2, 1:5]
{
  subtype.brca$`2016 Histology Annotations` %>% table()
  subtype.brca.filter <- subtype.brca %>% 
    filter(`2016 Histology Annotations` %in% c('Invasive ductal carcinoma', 'Invasive lobular carcinoma')) %>%
    mutate(PatientID = paste0('BRCA-', substr(CLID, start = 1, stop = 12)), .after='CLID')
  subtype.brca.pts <- subtype.brca.filter$PatientID %>% unique()
  BRCA.patients <- subtype.brca.pts
  subtype.brca.CLID <- subtype.brca.filter$CLID %>% unique()
}

TCGA.tumour_expressed.genes <- lapply(cancer.types, function(cancer_type){
  print(cancer_type)
  exprCounts <- readRDS(paste0(output_dir, cancer_type, '_RNAseq2.count_clin.survival.rds'))
  names(exprCounts)
  sample.df <- exprCounts[['sample.df']]
  sample.df[1:2,]
  if(cancer_type=='BRCA'){
    patients <- subtype.brca.filter$PatientID %>% unique()
  }else{
    patients <- sample.df$PatientID %>% unique()
  }
  sample.df <- sample.df %>% filter(PatientID%in%patients)
  sample.df$tissueType%>% table()
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
  # tumour and normal: 
  keep <- rowSums(expr.all > Cutoff.per.count) >= (ncol(expr.all)*Cutoff.fract)
  keep.expr.all <- expr.all[keep, ]
  
  genes <- rownames(keep.expr.tumour)
  return(genes)
})
names(TCGA.tumour_expressed.genes) <- cancer.types

saveRDS(TCGA.tumour_expressed.genes, file = paste0(output_dir,'TCGA.tumour_expressed.genes.RDS'))

##================================================================
## DEG analyses using TCGA data: 
##================================================================

suppressMessages(library(tximport))
suppressMessages(library(optparse))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(dplyr))
suppressMessages(library(EnhancedVolcano))
library(ggplot2)
library(ComplexHeatmap)
library(gplots) 

# Functions =========
DEG_DESeq2_function <- function(gene_counts_df, sample_id_df, geneList,
                                group.levels=c('tumour','normal'),
                                Cutoff.per.count = 1, Cutoff.fract = 0.2, # counts > 1 in 20% patients
                                FC_cutoff= 2, padj_cutoff= 0.05,
                                qvalue_cutoff = 0.01){
  ## start the function:
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
  
  # 3-Pre-filtering: Cutoff.per.count=1, Cutoff.fract= 20% (0.2)
  keep <- rowSums(counts(dds) > Cutoff.per.count) >= (ncol(gene_counts_df) * Cutoff.fract)
  dds <- dds[keep, ]
  # 4-DESeq() processing: !
  dds <- DESeq(dds) 
  
  # 5-normalise counts:
  # 1-cor() between samples: using normalized_counts 
  normalized_counts <- DESeq2::counts(dds, normalized = FALSE) 
  ## rank by mad values:
  normalized_counts_mad <- apply(normalized_counts, 1, stats::mad) 
  normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
  
  # normalise counts: by log transformation
  vst_d <- vst(dds, blind=FALSE) 
  
  ## extract exprData: assay()
  exprMatrix_vst <- assay(vst_d) # then order by normalized_counts_mad
  exprMatrix_vst <- exprMatrix_vst[order(normalized_counts_mad, decreasing=T), ]
  
  ### extract DE results: results();  group.levels <- c("Tumour" "Normal")
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
  res.comb[1:2,]
  
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
                                    Cutoff.per.count = 1, Cutoff.fract = 0.2, # counts > 1 in 20% patients
                                    FC_cutoff= 2, padj_cutoff=0.05){
  library(limma)
  library(edgeR) #DGEList
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
  # 3.2-Pre-filtering: 
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
  
  # 6-voom --> then multiple fits (voom--lmFit--contrasts.fit--eBayes)
  v <- voom(dge, design, plot=FALSE) 
  fit <- lmFit(v, design)
  
  fit2_contrasts <- contrasts.fit(fit, contrast.matrix)
  fit3_robust <- eBayes(fit2_contrasts, robust = TRUE)
  DEG_limma_voom = topTable(fit3_robust,sort.by = "P", coef=colnames(contrast.matrix), 
                            n=Inf) %>% na.omit(.) 
  
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

## run DEG analyses:
# Load data:
cancer.types
exprCounts_list <- lapply(cancer.types, function(x){
  dt <- readRDS(paste0(output_dir, x, '_RNAseq2.count_clin.survival.rds'))
})
names(exprCounts_list) <- cancer.types
tumour.expr <- exprCounts_list[[1]][["TumorNormalMatch"]][["primary.tumor"]]

# sub-cohort in TCGA with matched exome sequencing data and passing QC =======
CN.cohort_list <- readRDS(paste0(data_dir, 'CN.cohort_sample.list.RDS'))

# per cancer: ===========
# run DE analyses: --> keep mean.counts in tumour or in normal (DESeq2 & limma)
CN.overlapped <- c(TRUE)[1]
output.dir_base <- paste0(data_dir, '_CN.overlapped/')
use.matched <- c(TRUE)[1]
matched.type <- c('sample')
filter.BRCA <- c(TRUE)[1] 

DE.outputs_list <- lapply(cancer.types, function(ca){
  print(ca)
  exprCounts <- exprCounts_list[[ca]]
  cn.samples <- paste0(ca, '-', unique(CN.cohort_list[[ca]]$sample))
  # define cohort for DE: ======
  if(use.matched==FALSE){
    tmp.output_dir <- paste0(output.dir_base,'_nonMatched/')
    list.files(tmp.output_dir)
    if(!dir.exists(tmp.output_dir)){dir.create(tmp.output_dir)}
    
    patients <- substr(colnames(exprCounts[["RNASeq2_count"]]), 1, 17) %>%
      unique()
    patients[1:10]
    subtype.brca.pts[1:10]
    
    if(filter.BRCA==TRUE & ca=='BRCA'){
      patients <- patients[patients %in% subtype.brca.pts] 
    }
    
    if(CN.overlapped==TRUE){
      patients <- patients[patients %in% cn.samples] 
    }
    
    length(patients)
    sample_id_df <- exprCounts$sample.df %>% filter(PatientID %in% patients)
    exprCounts_df <- exprCounts$RNASeq2_count[,unique(sample_id_df$sampleID)]
    exprCounts_df[1:4, 1:4]
    geneList_chr <- rownames(exprCounts_df)
  }else if(matched.type=='sample'){
    tmp.output_dir <- paste0(output.dir_base,'_matched.samples/')
    if(!dir.exists(tmp.output_dir)){dir.create(tmp.output_dir)}
    patients <- colnames(exprCounts[["TumorNormalMatch"]][["primary.tumor"]])
    
    if(filter.BRCA==TRUE & ca=='BRCA'){
      patients <- patients[patients %in% subtype.brca.pts] # 112--88
    }
    
    if(CN.overlapped==TRUE){
      patients <- patients[patients %in% cn.samples] 
    }
    
    sample_id_df <- exprCounts$sample.df %>% dplyr::filter(.$PatientID %in% patients)
    samples.ID <- unique(sample_id_df$sampleID)
    exprCounts_df <- exprCounts$RNASeq2_count[, colnames(exprCounts$RNASeq2_count)%in%
                                                samples.ID]
    exprCounts_df[1:4, 1:4]
    geneList_chr <- rownames(exprCounts_df)
  }
  
  # run DE: ======
  {
    DEG.DESeq2_list <- DEG_DESeq2_function(gene_counts_df = exprCounts_df,
                                           sample_id_df = sample_id_df,
                                           geneList = geneList_chr,
                                           group.levels=c('tumour','normal'),
                                           Cutoff.per.count = 1, 
                                           Cutoff.fract = 0.2, # if=0.2: counts > 1 in 20% patients
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
                                              Cutoff.fract = 0.2, # counts > 1 in 20% patients
                                              FC_cutoff= 2, padj_cutoff=0.05)
    
    DEG.limma_list[[1]][1:2,]
    DEG.limma_list[[1]]$logFC %>% range()
    DEG.limma_list[[1]]$change %>% table()
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
  
  if(filter.BRCA==TRUE&ca=='BRCA'){
    subset.title <- paste0(paste0('matched.', substr(use.matched, 1,3), '_by.', matched.type, '_', ca),
                           '.filter.', substr(filter.BRCA, 1,3))
  }else{
    subset.title <- paste0(paste0('matched.', substr(use.matched, 1,3), '_by.', matched.type, '_', ca),
                           '.filter.', substr('FALSE', 1,3))
  }
  subset.title
  saveRDS(list, file = paste0(output_dir,'DEGs_DESeq2.Limma_',subset.title, '.RDS') )
  
  patients[1:10]
  write.csv(sample_id_df, paste0(output_dir,'sample.df_',ca,'_matched.', substr(use.matched, 1,3),'.csv'))
  
  # pdf(paste0(output_dir, 'volcano_DEGs.TCGA_',subset.title, '.pdf'), width = 8, height = 5)
  # plot(p.vol.deseq)
  # plot(p.vol.limma)
  # dev.off()
  
  # return:
  return(list)
})
names(DE.outputs_list) <- cancer.types

################################
## run boostrapping analyses:
################################
data.dir <- paste0(output.dir_base, "_matched.samples/")
list.files(data.dir)

DEG.outputs.TCGA_list <- lapply(cancer.types, function(ca){
  print(ca)
  if(filter.BRCA==TRUE&ca=='BRCA'){
    subset.title <- paste0(paste0('DEGs_DESeq2.Limma_matched.', substr(use.matched, 1,3), '_by.', matched.type, '_', ca),
                           '.filter.', substr(filter.BRCA, 1,3))
  }else{
    subset.title <- paste0(paste0('DEGs_DESeq2.Limma_matched.', substr(use.matched, 1,3), '_by.', matched.type, '_', ca),
                           '.filter.', substr('FALSE', 1,3))
  }
  subset.title
  grep(subset.title, list.files(data.dir), value = T)
  dt <- read_rds( paste0(data.dir, grep(subset.title, list.files(data.dir), value = T)) )
  return(dt)
})
names(DEG.outputs.TCGA_list) <- cancer.types

DE.package <- c('DESeq2')
FC_cutoff <- 2 #(not log2.FC_cutoff) same with Carlos: who used 1 for log2(FC) in his codes when running Tx421
P.value_cutoff <- 0.05
DEG.TCGA.list.used <- names(DEG.outputs.TCGA_list[[1]])[3]
DEG.TCGA.list.used

DEG.TCGA_filter <- lapply(cancer.types, function(x){
  sub.deg <- DEG.outputs.TCGA_list[[x]][[ DEG.TCGA.list.used ]] %>%
    dplyr::select(c('gene_id', grep(DE.package, colnames(.), value = T) )) %>%
    dplyr::select(c('gene_id', grep('log|value|adj', colnames(.), value = T) ))
  colnames(sub.deg) <- c('gene_id','log2.FC','p.value','p.adj')
  sub.deg <- sub.deg %>%
    mutate(DEG = if_else(.$p.adj < P.value_cutoff & abs(.$log2.FC) >= log2(FC_cutoff),
                         if_else(.$log2.FC >= log2(FC_cutoff),'Up','Down'),'Stable')) %>%
    dplyr::select(c('gene_id', 'DEG'))
  colnames(sub.deg)[2] <- paste0(x, '_DEG')
  return(sub.deg)
})
DEG.TCGA_filter <- Reduce(full_join, DEG.TCGA_filter) %>% dplyr::rename(gene = gene_id)

## ----------------------------------------------------
## Cell lines
{
  cell.line_list <- lapply(cancer.types, function(cancer_type){
    print(cancer_type)
    if(cancer_type == 'LUAD'){
      reference <- c('T2P')
      cancer.cells <- c('H1650', 'H1792','H2009',"A549")
      all.cells <- c('TT1','IMR90', 'H1650', 'H1792','H2009',"A549","LTX1000", "LTX543") # 'T2P',
      normal.cells <- c('T2P','TT1','IMR90')
      rep.cells <- c("A549","T2Prep","H1650rep")
      list.files(paste0('/Volumes/proj-tracerx-lung/tctProjects/mcgranahanLab/mdietzen/Repli-Seq/TRACERx/Release_w1000_20210305/'))
      PDCs <- c("LTX1000", "LTX543")
    }else if(cancer_type == 'BRCA'){ 
      reference <- c('HMEC')
      cancer.cells <- c('MDA453', 'SK-BR3','T47D','MCF-7') # (ENCODE)
      cells.encode <- c('T47D', 'MCF-7')
      all.cells <- c('MCF10A', 'MDA453', 'SK-BR3','T47D', 'MCF-7') # 'HMEC', 
      normal.cells <- c('HMEC','MCF10A')
    }
    out <- list(reference, cancer.cells, all.cells, normal.cells)
    names(out) <- c('reference', 'cancer.cells', 'all.cells', 'normal.cells')
    return(out)
  })
  names(cell.line_list) <- cancer.types
  cell.line_list[[1]]
}

# load expressed genes in TCGA:
TCGA.tumour_expressed.genes <- readRDS(paste0(output_dir, 'TCGA.tumour_expressed.genes.RDS'))

# load ARTextreme.genes:
ARTgenes.extreme_list <- readRDS(paste0(data_dir,'resultsARTgenes.extreme_list.RDS'))

# shared ART:
ARTgenes.consist_type_list <- lapply(cancer.types, function(ca){
  print(ca)
  data <- fread(paste0(data_dir, ca,"_shared", "ARTgenes.txt"))
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

## ----------------------------------------------------
## function: 
# Bootstrap test: compare mean.log2FC of earlier vs. (randomly selected genes in late of normal)
# bootstrapping: Function --> 
# --> whether log2FC of earlier_genes (late in normal, early in cancer) > of any late genes (in normal)
# --> whether log2FC of later_genes (early in normal, late in cancer) < of any early genes (in normal)

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
  # # test: # new note:
  # genes_test.ART = sub.ARTgenes
  # ref.genes.df = sub.Refgenes
  # y.axis = c('log2.FC')[1]
  # column.ART = c('timing')[1]
  # column.ref = c('refTiming')[1]
  # niter = iteration
  # timing.ART = c('earlier','later')
  # timing.ref = c('early','late')
  # replace_or.not = c(TRUE, FALSE)[1]
  
  # process data:
  # column.log2FC <- paste0(DEG.package.to.use, '_log2FoldChange')
  genes_test.ART[1:2,]
  ref.genes.df[1:2,]
  colnames(genes_test.ART)[match(column.ART, colnames(genes_test.ART))] <- 'column.ART'
  colnames(ref.genes.df)[match(column.ref, colnames(ref.genes.df))] <- 'column.ref'
  
  colnames(genes_test.ART)[match(y.axis, colnames(genes_test.ART))] <- 'y.axis'
  colnames(ref.genes.df)[match(y.axis, colnames(ref.genes.df))] <- 'y.axis'
  
  genes_test.ART[1:2,]
  ref.genes.df[1:2,]
  
  # bootstrapping: ===> Use cancer log2FC!
  ### eg: 1-earlier genes - Up: exclude those in late normal
  # earlier genes: mean.logFC_a; randomly selected late in normal: mean.logFC_b
  # --> calculate & plot 95% CI of mean.logFC_b, then plot & vs. mean.logFC_a
  genes_test.ART$column.ART %>% unique()
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
      # freq.random <- nrow(bg.sub[bg.sub$DE_change=='Up',])/nrow(bg.sub)
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
    out.list <- list(out, freq.random)
    names(out.list) <- c('boot.df', 'freq.random')
    return(out.list)
  })
  
  # combine <- Reduce(rbind, boot.per.cancer_list)
  names(boot.per.cancer_list) <- timing.ART
  combine <- boot.per.cancer_list
  return(combine)
}


## simple fisher: =======
# fisher:
fisher_expr.ARTgenes <- function(ARTgenes = ARTgenes, #cancer_genes, # nonCancer_genes, 
                                 use.gene.list = NULL, # or: test.genes
                                 ARTcolumn = c('timing'), 
                                 DEGcolumn = 'change_to.use',
                                 column.ref = c('refTiming')[1],
                                 DEG.change = c("Up","Down")[2],
                                 DEG.reference = "Stable"){
  # test:
  # ARTcolumn = c('timing')
  # DEGcolumn = 'change_to.use'
  # DEG.change = c("Up","Down")[2]
  # DEG.reference = "Stable"
  
  ARTgenes[1:2,]
  colnames(ARTgenes)[match(ARTcolumn,colnames(ARTgenes))] <- 'ARTcolumn'
  colnames(ARTgenes)[match(DEGcolumn,colnames(ARTgenes))] <- 'DEGcolumn'
  colnames(ARTgenes)[match(column.ref,colnames(ARTgenes))] <- 'column.ref'
  
  #classify genes are cancer or non_cancer
  if(!is.null(use.gene.list)){
    ARTgenes <- ARTgenes%>%filter(.$gene_name %in% use.gene.list)
  }
  
  ##run tests with normal early genes as reference
  ARTgenes[1:2,]
  early_ARTgenes_raw <- ARTgenes[ARTgenes$column.ref == 'early' & 
                                   ARTgenes$DEGcolumn%in%c(DEG.change, DEG.reference),] %>%
    group_by(ARTcolumn, DEGcolumn) %>%
    summarise(count = n(), .groups='keep')
  
  if(all(early_ARTgenes_raw$count[early_ARTgenes_raw$ARTcolumn%in%c('early', 'late')] != 0)){
    not_altered_counts <- early_ARTgenes_raw[early_ARTgenes_raw$ARTcolumn %in% c('early', 'late'),]
    early_ARTgenes     <- early_ARTgenes_raw[! early_ARTgenes_raw$ARTcolumn %in% c('early', 'late'),]
    not_altered_counts
    early_ARTgenes
    timings <- early_ARTgenes$ARTcolumn[duplicated(early_ARTgenes$ARTcolumn)]
    timings <- timings[!grepl('_noSwitch', timings)]
    if(length(timings) != 0){
      early_tests <- lapply(timings, function(x){
        # ma <- matrix(c(early_ARTgenes$count[early_ARTgenes$ARTcolumn == x & early_ARTgenes$DEGcolumn == DEG.change],
        #                early_ARTgenes$count[early_ARTgenes$ARTcolumn == x & early_ARTgenes$DEGcolumn == DEG.reference],
        #                not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.change],
        #                not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.reference]),
        #              ncol = 2, byrow = F,
        #              dimnames =  list(c(DEG.change, DEG.reference), c('altered', 'not_altered')))
        num11 <- early_ARTgenes$count[early_ARTgenes$ARTcolumn == x & early_ARTgenes$DEGcolumn == DEG.change]
        num22 <- early_ARTgenes$count[early_ARTgenes$ARTcolumn == x & early_ARTgenes$DEGcolumn == DEG.reference]
        num33 <- not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.change]
        num44 <- not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.reference]
        if(is_empty(num11)){num1 <- 0 }else{ num1 <- num11}
        if(is_empty(num22)){num2 <- 0 }else{ num2 <- num22}
        if(is_empty(num33)){num3 <- 0 }else{ num3 <- num33}
        if(is_empty(num44)){num4 <- 0 }else{ num4 <- num44}
        ma <- matrix(c(num1, num2, num3, num4), ncol = 2, byrow = F,
                     dimnames =  list(c(DEG.change, DEG.reference), c('altered', 'not_altered')))
        test <- fisher.test(ma)
        data.frame(timing = x, reference = 'early', oddsRatio = test$estimate, lowerCI = test$conf.int[1], upperCI = test$conf.int[2], pvalue = test$p.value,
                   altered_DEG = ma[1,1], nonAltered_DEG = ma[1,2],
                   altered_nonDEG = ma[2,1], nonAltered_nonDEG = ma[2,2])
        # OR <- (ma[1,1]/ma[1,2])/(ma[2,1]/ma[2,2])
        # data.frame(timing = x, reference = 'early', oddsRatio = OR #test$estimate, #lowerCI = test$conf.int[1], upperCI = test$conf.int[2], pvalue = test$p.value,
        #            # ,altered_DEG = ma[1,1], nonAltered_DEG = ma[1,2], 
        #            # altered_nonDEG = ma[2,1], nonAltered_nonDEG = ma[2,2]
        #            )
      })
      early_tests <- Reduce(rbind, early_tests)
    } else {
      early_tests <- c()
    }
  }
  early_tests
  
  #run tests with normal late genes as reference
  late_ARTgenes_raw <- ARTgenes[ARTgenes$column.ref == 'late' & 
                                  ARTgenes$DEGcolumn%in%c(DEG.change, DEG.reference),] %>%
    group_by(ARTcolumn, DEGcolumn) %>%
    summarise(count = n(), .groups='keep')
  late_ARTgenes_raw
  
  if(all(late_ARTgenes_raw$count[late_ARTgenes_raw$ARTcolumn%in%c('early', 'late')] != 0)){
    not_altered_counts <- late_ARTgenes_raw[late_ARTgenes_raw$ARTcolumn%in%c('early', 'late'),]
    late_ARTgenes     <- late_ARTgenes_raw[! late_ARTgenes_raw$ARTcolumn%in%c('early', 'late'),]
    timings <- late_ARTgenes$ARTcolumn[duplicated(late_ARTgenes$ARTcolumn)]
    timings <- timings[!grepl('_noSwitch', timings)]
    if(length(timings) != 0){
      late_tests <- lapply(timings, function(x){
        # x=timings[1]
        # ma <- matrix(c(late_ARTgenes$count[late_ARTgenes$ARTcolumn == x & late_ARTgenes$DEGcolumn == DEG.change],
        #                late_ARTgenes$count[late_ARTgenes$ARTcolumn == x & late_ARTgenes$DEGcolumn == DEG.reference],
        #                not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.change],
        #                not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.reference]),
        #              ncol = 2, byrow = F,
        #              dimnames =  list(c(DEG.change, DEG.reference), c('altered', 'not_altered')))
        num11.L <- late_ARTgenes$count[late_ARTgenes$ARTcolumn == x & late_ARTgenes$DEGcolumn == DEG.change]
        num22.L <- late_ARTgenes$count[late_ARTgenes$ARTcolumn == x & late_ARTgenes$DEGcolumn == DEG.reference]
        num33.L <- not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.change]
        num44.L <- not_altered_counts$count[not_altered_counts$DEGcolumn == DEG.reference]
        if(is_empty(num11.L)){num1.L <- 0 }else{ num1.L <- num11.L}
        if(is_empty(num22.L)){num2.L <- 0 }else{ num2.L <- num22.L}
        if(is_empty(num33.L)){num3.L <- 0 }else{ num3.L <- num33.L}
        if(is_empty(num44.L)){num4.L <- 0 }else{ num4.L <- num44.L}
        ma <- matrix(c(num1.L, num2.L, num3.L, num4.L), ncol = 2, byrow = F,
                     dimnames =  list(c(DEG.change, DEG.reference), c('altered', 'not_altered')))
        test <- fisher.test(ma)
        data.frame(timing = x, reference = 'late', oddsRatio = test$estimate, lowerCI = test$conf.int[1], upperCI = test$conf.int[2], pvalue = test$p.value,
                   altered_DEG = ma[1,1], nonAltered_DEG = ma[1,2],
                   altered_nonDEG = ma[2,1], nonAltered_nonDEG = ma[2,2])
      })
      late_tests <- Reduce(rbind, late_tests)
    } else {
      late_tests <- c()
    }
  }
  late_tests
  #combine tests
  df <- rbind(early_tests, late_tests)
  if(!is.null(df)){
    df <- df %>% dplyr::mutate(DEG_comparison = paste(DEG.change, DEG.reference, sep = '-'), .after='reference')
    rownames(df) <- NULL
    return(df)
  }else{
    return(NULL)
  }
  
}

## ----------------------------------------------------
# boot of logFC
ART.column.used <- c('timing')
data.sets <- c('TCGA')
DE.package <- c('DESeq2')
iteration <- 100000
y.axis_used <- c('log2.FC')
FC_cutoff <- 2 #(not log2.FC_cutoff) 
P.value_cutoff <- 0.05

ARTgenes.DEG_comb.list <- lapply(data.sets, function(data.set){
  print(data.set)
  to.comb_DEG.df <- DEG.TCGA_filter
  DEG.out_filter.list <- list(BRCA = DEG.outputs.TCGA_list[['BRCA']][['DEG_output.comb']],
                              LUAD = DEG.outputs.TCGA_list[['LUAD']][['DEG_output.comb']])
  cancer.types <- names(DEG.out_filter.list)
  
  #### Process DEG results with ART: (define up/down/stable) --------
  DEG.out_filter_sub <- lapply(cancer.types, function(x){
    print(paste0('DEG:', x))
    if(data.set=="TRACERx"){
      sub.deg <- DEG.out_filter.list[[x]]%>%filter(coef=="Tumour - Normal")
    }else{
      sub.deg <- DEG.out_filter.list[[x]]%>%
        dplyr::select(c('gene_id', grep(DE.package, colnames(.), value = T, ignore.case = T) ))
    }
    sub.deg <- sub.deg%>%
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
    return(ARTgenes)
  })
  names(comb.ARTgenes.DEG_list) <- cancer.types
  
  #final output:
  return(comb.ARTgenes.DEG_list)
})
names(ARTgenes.DEG_comb.list) <- data.sets

## run boot: 
ref.used <- c('vs.allRT_')[1] 
remove.ART.in.ref <- c(TRUE)[1]

boot.logFC_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  cancers <- names(ARTgenes.DEG_comb.list[[data.set]])
  
  to.comb_DEG.df <- DEG.TCGA_filter
  DEG.out_filter.list <- list(BRCA = DEG.outputs.TCGA_list[['BRCA']][['DEG_output.comb']],
                              LUAD = DEG.outputs.TCGA_list[["LUAD"]][['DEG_output.comb']])
  
  #### Process DEG results with ART: (define up/down/stable) --------
  DEG.out_filter_sub <- lapply(cancers, function(x){
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
  names(DEG.out_filter_sub) <- cancers
  
  # run boot of logFC:
  boot.list <- lapply(cancers, function(ca){
    print(ca)
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[ca]]
    ARTgenes <- ARTgenes.DEG_comb.list[[data.set]][[ca]]
    ARTgenes[1:2,]
    
    if(ref.used=='vs.consistRT_'){
      ref.genes.df <- ARTgenes
    }else if(ref.used=="vs.allRT_"){
      # if all.genes in normal as ref: --> use DEG for all genes:
      ARTgenes.extreme_list[[ca]][1:2,]
      DEG.out_filter_sub[[ca]][1:2,]
      ref.genes.df <- ARTgenes.extreme_list[[ca]] %>% 
        dplyr::select(gene_name,chr,start,stop,normal_l2r) %>% distinct() %>%
        mutate(refTiming = if_else(normal_l2r>0, 'early', 'late')) %>%
        left_join(DEG.out_filter_sub[[ca]], by=c('gene_name'='gene')) %>%
        filter(gene_name %in% sub.expressed.TCGA) %>%
        filter(!is.na(log2.FC))
      ref.genes.df[1:2,]
    }
    
    ref.genes.df$refTiming %>% unique()
    ref.genes.df$change_to.use %>% unique()
    
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
                                       remove.ART.in.ref = remove.ART.in.ref)
    return(boot.out)
  })
  # out:
  # out <- Reduce(r,bind, boot.list)
  names(boot.list) <- cancers
  out <- boot.list
  return(out)
})
names(boot.logFC_ART_list) <- data.sets

saveRDS(boot.logFC_ART_list, file = paste0(output_dir,'boot.logFC_rm.ART_matched.',use.matched,
                                           '_',ref.used,data.sets,'_', iteration/1000, 'k.RDS'))

## 
file.used <- "boot.logFC_rm.ART_matched.TRUE_vs.allRT_TCGA_100k.RDS"
boot.logFC_ART_list <- readRDS(paste0(output_dir,file.used))

# count % DEG among earlier/later ========
frac.DEG_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  cancers <- names(ARTgenes.DEG_comb.list[[data.set]])
  # frac of DEG per cancer:
  frac.list <- lapply(cancers, function(ca){
    print(ca)
    ARTgenes <- ARTgenes.DEG_comb.list[[data.set]][[ca]]
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

# simple.fisher test: fisher: earlier vs. late & Up/Down vs. Stable========
fisher.DEG_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  cancers <- names(ARTgenes.DEG_comb.list[[data.set]])
  # fisher: earlier vs. late & Up/Down vs. Stable -->
  fisher.list <- lapply(cancers, function(ca){
    print(ca)
    ARTgenes <- ARTgenes.DEG_comb.list[[data.set]][[ca]]
    ARTgenes[1:2,]
    ARTgenes$DEGcolumn %>% unique()
    ARTgenes$timing %>% unique()
    fisher <- lapply(c("Up", "Down"), function(d){
      df <- fisher_expr.ARTgenes(ARTgenes = ARTgenes, 
                                 use.gene.list = NULL, 
                                 ARTcolumn = c('timing'), 
                                 DEGcolumn = 'change_to.use',
                                 column.ref = c('refTiming')[1],
                                 DEG.change = d,
                                 DEG.reference = "Stable")
    })
    fisher <- Reduce(rbind, fisher) %>% mutate(cancer.type=ca)
    fisher
    return(fisher)
  })
  # out:
  out <- Reduce(rbind, fisher.list)
  out
  return(out)
})
names(fisher.DEG_ART_list) <- data.sets

# load boot.logFC data: for plotting: =======
# annot.color:
{
  display.brewer.all()
  display.brewer.pal(8,'Set2')
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

# calculate p value for bootstrapping analyses:
data.set <- 'TCGA'
plot.boot.logFC <- rbind(boot.logFC_ART_list[[data.set]][['BRCA']][['earlier']][[1]],
                         boot.logFC_ART_list[[data.set]][['BRCA']][['later']][[1]],
                         boot.logFC_ART_list[[data.set]][['LUAD']][['earlier']][[1]],
                         boot.logFC_ART_list[[data.set]][['LUAD']][['later']][[1]])
boot.p.value_dt <- cbind(boot.logFC_ART_list[[data.set]][['BRCA']][['earlier']][[2]],
                         boot.logFC_ART_list[[data.set]][['BRCA']][['later']][[2]],
                         boot.logFC_ART_list[[data.set]][['LUAD']][['earlier']][[2]],
                         boot.logFC_ART_list[[data.set]][['LUAD']][['later']][[2]])
colnames(boot.p.value_dt) <- c('BRCA_earlier.vs.late','BRCA_later.vs.early',
                               'LUAD_earlier.vs.late','LUAD_later.vs.early')
p.value.df <- data.frame()
for (i in colnames(boot.p.value_dt)) {
  # i=colnames(boot.p.value_dt)[1]
  cancer1 <- unlist(strsplit(i,split='_'))[1]
  timing1 <- unlist(strsplit(unlist(strsplit(i,split='_'))[2],split='.vs.'))[1]
  mean1 <- plot.boot.logFC %>% filter(cancer.type==cancer1&timing==timing1)%>%pull(mean_original)
  freq.boot <- boot.p.value_dt[,i]
  if (timing1 == 'earlier') {
    p1 <- length(freq.boot[mean1<=freq.boot])/iteration
  }else{
    p1 <- p1 <- length(freq.boot[mean1>=freq.boot])/iteration
  }
  dt <- data.frame(cancer.type = cancer1, timing = timing1,
                   mean_original = mean1, p.boot = p1)
  p.value.df <- rbind(p.value.df, dt)
}

# make plots:
plot.frac.DEG <- frac.DEG_ART_list[[data.set]] %>%
  mutate(change_to.use=factor(change_to.use,levels = c('Up','Down','Stable'))) %>%
  mutate(timing=factor(timing,levels = names(color.key)))%>%
  mutate(cancer.type=factor(cancer.type,levels = c(cancer.types)))

plot.frac.DEG <- plot.frac.DEG[order(plot.frac.DEG$cancer.type,
                                     plot.frac.DEG$timing,
                                     plot.frac.DEG$change_to.use),] %>%
  group_by(cancer.type, timing) %>%
  mutate(cumulative = cumsum(frac), midpoint = cumulative-frac/2)

plot.fisher.DEG <- fisher.DEG_ART_list[[data.set]] %>% rename(comparison=DEG_comparison) %>%
  mutate(comparison=sub('-Stable','',comparison))

# plots: facets by earlier/later:
plot.boot.logFC <- plot.boot.logFC %>% 
  pivot_longer(cols=c('mean_original', 'mean_iter'), names_to='data.type', values_to='mean.value') %>%
  as.data.frame() %>%
  mutate(x.label = if_else(data.type=='mean_original', timing, if_else(data.type=='mean_iter', refTiming, 'other'))) %>%
  mutate(y.observed = if_else(data.type=='mean_original', mean.value, NULL),
         y.random = if_else(data.type=='mean_iter', mean.value, NULL)) %>%
  mutate(y.lowerCI = if_else(data.type=='mean_iter', lowerCI, NULL),
         y.upperCI = if_else(data.type=='mean_iter', upperCI, NULL),
         numb.label = if_else(data.type=='mean_original', 'N.Gene', 'N.rdGene') )
plot.boot.logFC[1:4,]

p_boot_by.cancer2 <- plot.boot.logFC %>% 
  ggplot(aes(x=factor(x.label, levels = c("earlier","late","later","early")),# factor(cancer.type, levels = cancer.types), 
             color=x.label)) +
  geom_errorbar(aes(ymin=y.lowerCI, ymax=y.upperCI, colour=factor(x.label, levels = names(color.key))), width = 0.1) +
  geom_point(aes(y = y.observed, colour = factor(x.label, levels = names(color.key)),
                 shape = paste0('Observed(mean)_ART')), size=7) +
  geom_point(aes(y = y.random, colour=factor(refTiming,levels=names(color.key)),
                 shape = paste0('Random(mean)_Reference') ), size=5) +
  geom_text(aes(y=y.observed, label=paste0('(N=', n.test.genes,')')), 
            nudge_y=-0.03, stat = 'identity',check_overlap=T) +
  scale_colour_manual(name = paste0('rep.Timing'), values = c(color.key)) +
  scale_shape_manual(name =paste0('Data'), values = c(18,16)) +
  facet_grid(~factor(cancer.type, levels = cancer.types),scales="free") +
  xlab('') + ylab('Log2.Fold.Change\nin tumour vs. normal') +
  theme_bw() + theme_hz +
  theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
        panel.grid = element_blank(),
        strip.text = element_text(size=12))
p_boot_by.cancer2 

pdf(paste0(output_dir, 'p.boot_',file.used,'.pdf'), 
    width = 9, height = 6)
print(p_boot_by.cancer2)
dev.off()

write_csv(p.value.df, paste0(output_dir,'p.value_boot.log2FC_',file.used,'.csv'))

# pie: earlier/later per.cancer
plot.frac.DEG[1:4,]
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
print(p.pie_list[[1]])

pdf(paste0(output_dir, 'p.pie_DEG.in.ART_', data.set,'.pdf'), 
    width = 8, height = 5)
for (i in names(p.pie_list)) {
  print(p.pie_list[[i]])
}
dev.off()

rm(data.dir)

# ====================================================================================================
# Figure 6B: Per.cancer: 
# mean.CN/Ploidy (vs. random in reference): similar to what was initially done in logFC of gene expression
## load CN/Ploidy data:
data_dir
gene.CN_TCGA.list <- lapply(cancer.types, function(x){
  if(x=='BRCA'){
    path <- paste0(data_dir, x, 'subset_mean_cnDiff_perGene.rds')
  }else{
    path <- paste0(data_dir, x, '_mean_cnDiff_perGene.rds')
  }
  out <- readRDS(path)
  return(out)
})
names(gene.CN_TCGA.list) <- cancer.types

# load expressed genes in TCGA:
TCGA.tumour_expressed.genes <- readRDS(paste0(output_dir, 'TCGA.tumour_expressed.genes.RDS'))

# load ARTextreme.genes:
ARTgenes.extreme_list <- readRDS(paste0(data_dir,'resultsARTgenes.extreme_list.RDS'))

# shared ART:
ARTgenes.consist_type_list <- lapply(cancer.types, function(ca){
  print(ca)
  data <- fread(paste0(data_dir, ca,"_shared", "ARTgenes.txt"))
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

# boot of logFC ===========
ART.column.used <- c('timing')
data.sets <- c('TCGA')
iteration <- 100000
y.axis_used <- c('mean_diff')

ARTgenes.CNA_comb.list <- lapply(data.sets, function(data.set){
  print(data.set)
  gene.CN_used.list <- gene.CN_TCGA.list
  cancer.types <- names(gene.CN_used.list)
  
  #### Process CNA results with ART: --------
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

## run boot: ---------
ref.used <- c('vs.allRT_')[1] 
genes.used <- c('all', 'expressed')[2]
remove.ART.in.ref <- c(TRUE, FALSE)[1]

boot.CN.ploidy_ART_list <- lapply(data.sets, function(data.set){
  print(data.set)
  gene.CN_used.list <- gene.CN_TCGA.list
  cancers <- names(gene.CN_used.list)
  
  #### Process CNA results with ART: --------
  # run boot of CN.ploidy:
  boot.list <- lapply(cancers, function(ca){
    # ca=cancers[1]
    print(ca)
    sub.expressed.TCGA <- TCGA.tumour_expressed.genes[[ca]]
    ARTgenes <- ARTgenes.CNA_comb.list[[data.set]][[ca]] %>% filter(!is.na(mean_diff))
    ARTgenes[1:2,]
    # is.na(ARTgenes$mean_diff) %>% table()
    sub.CNA <- gene.CN_used.list[[ca]]
    sub.CNA[1:2,]
    
    if(ref.used=='vs.consistRT_'){
      ref.genes.df <- ARTgenes
    }else if(ref.used=="vs.allRT_"){
      # if all.genes in normal as ref: --> use CNA for all genes:
      ARTgenes.extreme_list[[ca]][1:2,]
      ref.genes.df <- ARTgenes.extreme_list[[ca]] %>% 
        dplyr::select(gene_name,chr,start,stop,normal_l2r) %>% distinct() %>%
        mutate(refTiming = if_else(normal_l2r>0, 'early', 'late')) %>%
        left_join(sub.CNA, by=c('gene_name', 'chr', 'start', 'stop'))
      ref.genes.df[1:2,]
    }
    ref.genes.df <- ref.genes.df %>% filter(!is.na(mean_diff))
    
    if(genes.used=="expressed"){
      ref.genes.df <- ref.genes.df %>% filter(gene_name %in% sub.expressed.TCGA)
      ARTgenes <- ARTgenes %>% filter(gene_name %in% sub.expressed.TCGA)
    }
    
    # run boot:
    ARTgenes[1:2,]
    ref.genes.df[1:2,]
    boot.out <- boot.log2FC_per.cancer(genes_test.ART=ARTgenes,
                                       ref.genes.df = ref.genes.df,
                                       y.axis = c('mean_diff')[1],
                                       column.ART = c('timing')[1],
                                       column.ref = c('refTiming')[1],
                                       niter = iteration,
                                       timing.ART = c('earlier','later'),
                                       timing.ref = c('early','late'),
                                       replace_or.not = c(TRUE, FALSE)[1],
                                       remove.ART.in.ref = remove.ART.in.ref)
    return(boot.out)
  })
  names(boot.list) <- cancers
  # out:
  return(boot.list)
})
names(boot.CN.ploidy_ART_list) <- data.sets

saveRDS(boot.CN.ploidy_ART_list, file = paste0(output_dir,'boot.CN.ploidy-1_',genes.used, '_',data.sets,'_',
                                               ref.used,'removeARTinRef.',remove.ART.in.ref,'_', iteration/1000, 'k_BRCAsubset.RDS'))

# load boot.CN.ploidy data: for plotting: =======
data.set <- c("TCGA")[1] 

plot.boot.CN.ploidy <- rbind(boot.CN.ploidy_ART_list[[data.set]][['BRCA']][['earlier']][[1]],
                             boot.CN.ploidy_ART_list[[data.set]][['BRCA']][['later']][[1]],
                             boot.CN.ploidy_ART_list[[data.set]][['LUAD']][['earlier']][[1]],
                             boot.CN.ploidy_ART_list[[data.set]][['LUAD']][['later']][[1]])
plot.boot.CN.ploidy[1:2,]

# calculate p value for bootstrapping analyses:
data.set <- 'TCGA'
boot.p.value_dt <- cbind(boot.CN.ploidy_ART_list[[data.set]][['BRCA']][['earlier']][[2]],
                         boot.CN.ploidy_ART_list[[data.set]][['BRCA']][['later']][[2]],
                         boot.CN.ploidy_ART_list[[data.set]][['LUAD']][['earlier']][[2]],
                         boot.CN.ploidy_ART_list[[data.set]][['LUAD']][['later']][[2]])
colnames(boot.p.value_dt) <- c('BRCA_earlier.vs.late','BRCA_later.vs.early',
                               'LUAD_earlier.vs.late','LUAD_later.vs.early')
p.value.df <- data.frame()
for (i in colnames(boot.p.value_dt)) {
  # i=colnames(boot.p.value_dt)[1]
  cancer1 <- unlist(strsplit(i,split='_'))[1]
  timing1 <- unlist(strsplit(unlist(strsplit(i,split='_'))[2],split='.vs.'))[1]
  mean1 <- plot.boot.CN.ploidy %>% filter(cancer.type==cancer1&timing==timing1)%>%pull(mean_original)
  freq.boot <- boot.p.value_dt[,i]
  if (timing1 == 'earlier') {
    p1 <- length(freq.boot[mean1<=freq.boot])/iteration
  }else{
    p1 <- p1 <- length(freq.boot[mean1>=freq.boot])/iteration
  }
  dt <- data.frame(cancer.type = cancer1, timing = timing1,
                   mean_original = mean1, p.boot = p1)
  p.value.df <- rbind(p.value.df, dt)
}

# annot.color:
{
  display.brewer.all()
  display.brewer.pal(8,'Set2')
  color.key <- c(brewer.pal(9,'PiYG')[c(1,2, 8,9)])
  names(color.key) <- c("early", "earlier","later","late")
  color.key
  
  col.CNA <- c(rev(brewer.pal(8,'Set2')[c(1:2)]), 'darkgrey')
  names(col.CNA) <- c('Up','Down','Stable')
  
  theme_hz <- theme(axis.text.x = element_text(angle=90, hjust = 0.5,size = 12, color = 'black'),
                    axis.text.y = element_text(size = 12, color = 'black'),
                    axis.title = element_text(size = 14, color = 'black'),
                    legend.text = element_text(size = 10, color = 'black'),
                    legend.title = element_text(size = 12, color = 'black', face = 'bold'),
                    plot.title = element_text(size = 15, color = 'black', face = 'bold'))
  
}

# facets by earlier/later:
{
  plot.boot.CN.ploidy <- plot.boot.CN.ploidy %>% 
    pivot_longer(cols=c('mean_original', 'mean_iter'), names_to='data.type', values_to='mean.value') %>%
    as.data.frame() %>%
    mutate(x.label = if_else(data.type=='mean_original', timing, if_else(data.type=='mean_iter', refTiming, 'other'))) %>%
    mutate(y.observed = if_else(data.type=='mean_original', mean.value, NULL),
           y.random = if_else(data.type=='mean_iter', mean.value, NULL)) %>%
    mutate(y.lowerCI = if_else(data.type=='mean_iter', lowerCI, NULL),
           y.upperCI = if_else(data.type=='mean_iter', upperCI, NULL),
           numb.label = if_else(data.type=='mean_original', 'N.Gene', 'N.rdGene') )
  plot.boot.CN.ploidy[1:2,]
  
  p_boot_by.cancer2 <- plot.boot.CN.ploidy %>% 
    ggplot(aes(x=factor(x.label, levels = c("earlier","late","later","early")),# factor(cancer.type, levels = cancer.types), 
               color=x.label)) +
    geom_errorbar(aes(ymin=y.lowerCI, ymax=y.upperCI, colour=factor(x.label, levels = names(color.key))), width = 0.1) +
    geom_point(aes(y = y.observed, colour = factor(x.label, levels = names(color.key)),
                   shape = paste0('Observed(mean)_ART')), size=7) +
    geom_point(aes(y = y.random, colour=factor(refTiming,levels=names(color.key)),
                   shape = paste0('Random(mean)_Reference') ), size=5) +
    geom_text(aes(y=y.observed, label=paste0('(N=', n.test.genes,')')), 
              nudge_y=-0.005, stat = 'identity',check_overlap=T) +
    scale_colour_manual(name = paste0('rep.Timing'), values = c(color.key)) +
    scale_shape_manual(name =paste0('Data'), values = c(18,16)) +
    facet_grid(~factor(cancer.type, levels = cancer.types),scales="free") +
    xlab('') + ylab('Mean CN/Ploidy in tumour') +
    theme_bw() + theme_hz +
    theme(axis.text.x = element_text(angle=45,hjust=1,vjust=1),
          panel.grid = element_blank(),
          strip.text = element_text(size=12))
  p_boot_by.cancer2
  
}
p.title <- paste0('p.boot.wide_CN.ploidy-1_', genes.used) 
pdf(paste0(output_dir, p.title,'_',ref.used, data.set,'_BRCAsubtype.pdf'), 
    width = 9, height = 5)
print(p_boot_by.cancer2)
dev.off()

write_csv(p.value.df, paste0(output_dir,'p.value_',file.used,'.csv'))


###############################################################################################################################
# Figure 6C:
###############################################################################################################################
cancer.types <- c("BRCA", "LUAD")

# Load overlaping.ART genes:
overlap.ARTgenes_list <- lapply(cancer.types, function(x){
  print(x)
  data <- readRDS(paste0(data_dir,x, '_overlapping_ARTgenes.rds')) %>% 
    mutate(cancerType = x)
  return(data)
})
names(overlap.ARTgenes_list) <- cancer.types
overlap.ARTgenes_list[[1]][1:2,]

overlap.ART_cancer.genes_list <- lapply(cancer.types, function(cancer){
  overlap.ART <- overlap.ARTgenes_list[[cancer]] %>% filter(class %in% c('shared', 'recurrent'))
  # "onco","tsg","cancerType_cancerGene"
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
{
  display.brewer.all()
  display.brewer.pal(11, 'PiYG')
  annot.col <- list(RT = c(brewer.pal(11, 'PiYG')[c(3,9)]))
  names(annot.col$RT) <- c("earlier", "later")
  
  theme_hz <- theme(axis.text.x = element_text(angle=90, hjust = 0.5,size = 12, color = 'black'),
                    axis.text.y = element_text(size = 12, color = 'black'),
                    axis.title = element_text(size = 14, color = 'black'),
                    legend.text = element_text(size = 10, color = 'black'),
                    legend.title = element_text(size = 12, color = 'black', face = 'bold'),
                    plot.title = element_text(size = 15, color = 'black', face = 'bold'))
}


# geneType:
plot_data.cancer <- plot_data %>% filter(geneType=="Cancer genes")
plot_data.cancer[1:2,]
p.cancer_genes <- plot_data.cancer %>%
  ggplot(aes(x = factor(cancerType,levels = c(cancer.types)), #factor(cancerType,levels = cancer.types), 
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

pdf(paste0(output_dir, 'bar_frac_cancerGenes.pdf'), width = 7, height = 4, useDingbats = F)
print(p.cancer_genes)
dev.off()

#####################################################################
### Figure 6C: Cancer genes ===========
# Load overlaping.ART genes:
cancer.types <- c('BRCA', 'LUAD')
overlap.ARTgenes_list <- lapply(cancer.types, function(x){
  print(x)
  normal <- normal_cellLines[[x]]
  data <- readRDS(paste0(data_dir, x, '_overlapping_ARTgenes.rds')) %>% # genes with recurrent/shared/unique ART or consistent RT (identified in Figure 2E)
    mutate(cancerType = x)
  return(data)
})
names(overlap.ARTgenes_list) <- cancer.types

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

# Figure 6C:
plot_data.cancer <- plot_data %>% filter(geneType=="Cancer genes")
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

pdf(paste0(output_dir, 'bar_frac_cancerGenes.pdf'), width = 7, height = 4, useDingbats = F)
print(p.cancer_genes)
dev.off()

#####################################################################
#####################################################################
### Supplementary Figure 11A: Proportions of genes with consistent RT or shared/recurrent ART (in BRCA, LUAD) ======
all.cellLines <- data.frame(cellLine = c(unlist(cancer_cellLines), normal_cellLines)) %>%
  rownames_to_column(var = 'cancer.type') %>% 
  mutate(cancer.type = substr(cancer.type, 1,4)) %>%
  mutate(tissue.type = if_else(cellLine %in% normal_cellLines, 'normal', 'cancer'))

# load ARTextreme.genes:
# ARTgenes.extreme_list <- readRDS(paste0(data_dir, 'ARTgenes.per.cell_list.RDS'))
ARTgenes.extreme_list <- readRDS(paste0(data_dir, 'resultsARTgenes.extreme_list.RDS'))

# shared ARTgenes:
ARTgenes.consist_type_list <- lapply(cancer.types, function(x){
  print(x)
  # genes with shared/recurrent ART or consistent RT in BRCA, LUAD
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
  
  for (i in cancer.types[2:length(cancer.types)]) {
    df <- ARTgenes.consist_type_list[[i]] %>% dplyr::select('gene_name','timing')
    colnames(df)[2] <- i
    ARTgenes.consist_df <- ARTgenes.consist_df %>% 
      full_join(df, by=c('gene_name'))
  }
  ARTgenes.consist_df[1:2,]
  count_ART <- apply(ARTgenes.consist_df[, 5:(4+length(cancer.types)) ], 1, function(x) length(x[x%in%c('earlier','later')]))
  
  timings <- apply(ARTgenes.consist_df[,5:(4+length(cancer.types))], 1, function(x){
    if( length(unique(x))==1 ){
      y <- unique(x)
    }else{
      y <- 'discordant'
    }
    return(y)
  })
  ARTgenes.consist_df$timing <- timings
  ARTgenes.consist_df$count_ART <- count_ART
  
  keep <- apply(ARTgenes.consist_df[,5:(4+length(cancer.types))], 1, function(x) !anyNA(x) )
  genes.consistCRT <- ARTgenes.consist_df[keep, ]%>%
    filter(count_ART==0 & timing!='discordant')
  
}
ARTgenes.consist_df[1:2,]
genes.sharedART <- ARTgenes.consist_df%>%filter(count_ART==length(cancer.types))

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

pdf(paste0(output_dir,'figS11_Frac.sharedART.genes_per.cancer','.pdf'), 
    width = 7, height = 6)
print(p.frac)
dev.off()

## Supplementary Figure 11B: Oncogenes and Tumour suppressor genes =====
cancer.types <- c("BRCA","LUAD")

all.cellLines <- data.frame(cellLine = c(unlist(cancer_cellLines), normal_cellLines)) %>%
  rownames_to_column(var = 'cancer.type') %>% 
  mutate(cancer.type = substr(cancer.type, 1,4)) %>%
  mutate(tissue.type = if_else(cellLine %in% normal_cellLines, 'normal', 'cancer'))
all.cellLines[1:4,]

# load shared ARTgenes:
ARTgenes.consist_type_list <- lapply(cancer.types, function(x){
  print(x)
  data <- fread( paste0(data_dir, x,"_sharedARTgenes",".txt") ) %>%
    unique
  return(data)
})
names(ARTgenes.consist_type_list) <- cancer.types

# shared art.genes in all cancer:
{
  ARTgenes.consist_df <- ARTgenes.consist_type_list[[1]] %>% 
    dplyr::select('gene_name','chr','start','stop','timing') 
  ARTgenes.consist_df[1:2,]
  for (i in cancer.types[2:length(names(ARTgenes.consist_type_list))]) {
    df <- ARTgenes.consist_type_list[[i]] %>% dplyr::select('gene_name','timing')
    colnames(df)[2] <- i
    ARTgenes.consist_df <- ARTgenes.consist_df %>% 
      full_join(df, by=c('gene_name'))
  }
  ARTgenes.consist_df[1:2,]
  colnames(ARTgenes.consist_df)[5:(4+length(names(ARTgenes.consist_type_list)))] <- cancer.types
  count_ART <- apply(ARTgenes.consist_df[,5:(4+length(names(ARTgenes.consist_type_list))) ], 1, function(x) length(x[x%in%c('earlier','later')]))
  
  timings <- apply(ARTgenes.consist_df[,5:(4+length(names(ARTgenes.consist_type_list)))], 1, function(x){
    if( length(unique(x))==1 ){
      y <- unique(x)
    }else{
      y <- 'discordant'
    }
    return(y)
  })
  ARTgenes.consist_df[1:2,]
  ARTgenes.consist_df$timing <- timings
  ARTgenes.consist_df$count_ART <- count_ART
  
  keep <- apply(ARTgenes.consist_df[,5:(4+length(names(ARTgenes.consist_type_list)))], 1, function(x) !anyNA(x) )
  genes.consistCRT <- ARTgenes.consist_df[keep, ]%>%
    filter(count_ART==0 & timing!='discordant')
}
genes.sharedART <- ARTgenes.consist_df%>%filter(count_ART>=2)

## Analyses: ------------
{
  essen.breast <- read_xlsx(paste0(data_dir,'20220203_Essential_breast.xlsx')) %>%
    as.data.frame() %>% pull(`Gene/Compound`) %>% unique()
  essen.lung <- read_xlsx(paste0(data_dir,'20220203_Essential_lung.xlsx')) %>%
    as.data.frame() %>% pull(`Gene/Compound`) %>% unique()
  
  driver.lung <- fread(paste0(data_dir,'20220203_lung_drivers.csv')) %>% 
    pull(Gene_Symbol) %>% unique()
  driver.breast <- fread(paste0(data_dir,'20220208_breast.driver.csv')) %>%
    pull(gene) %>% unique()
  driver.pancan <- fread(paste0(data_dir,'20220203_pancan_drivers.csv')) %>%
    pull(gene) %>% unique()
  
  cancer.og.tsg <- fread(paste0(data_dir,'20220203_pan.driver_og.tsg.csv'))
  cancer.og <- cancer.og.tsg%>%filter(oncogene==T) %>% pull(Gene) %>% unique()
  cancer.tsg <- cancer.og.tsg%>%filter(tumour_suppressor==T) %>% pull(Gene) %>% unique()
  
  to.comb_all.geneSet <- data.frame(gene = unique(c(cancer.og, cancer.tsg, driver.lung, driver.breast, driver.pancan,
                                                    essen.breast, essen.lung)) ) %>%
    mutate(OG.TSG = if_else(gene%in%cancer.og, 'OG', if_else(gene%in%cancer.tsg, 'TSG', NULL)),
           driver.gene = if_else(gene%in% intersect(driver.lung, driver.breast), 'lung&breast',
                                 if_else(gene%in%driver.lung, 'lung', if_else(gene%in%driver.breast, 'breast', NULL))),
           driver.pancan = if_else(gene%in%driver.pancan, 'pancancer',NULL),
           essential.type = if_else(gene%in% intersect(essen.lung, essen.breast), 'lung&breast',
                                    if_else(gene%in%essen.lung, 'lung', if_else(gene%in%essen.breast, 'breast', NULL)))
    ) 
}

shared.genes_CAG_list <- lapply(cancer.types, function(x){
  print(x)
  tissue <- if_else(x=='BRCA','breast','lung')
  sub.shared.all <- ARTgenes.consist_type_list[[x]] %>%
    filter(timing %in% c('earlier','later')) %>%
    left_join(to.comb_all.geneSet,by=c('gene_name'='gene')) 
  sub.shared.all[1:2,]
  freq.all <- sub.shared.all %>% group_by(timing) %>%
    summarise(count_all = n())
  
  sub.shared.test <- sub.shared.all %>%
    filter(!is.na(OG.TSG) | grepl(tissue, driver.gene)) %>%
    dplyr::select(c('gene_name','timing','OG.TSG','driver.gene'))
  sub.shared.test$is.OT.TSG <- if_else(!is.na(sub.shared.test$OG.TSG),TRUE,NULL)
  sub.shared.test$driver.gene <- if_else(grepl(tissue, sub.shared.test$driver.gene),tissue,NULL)
  sub.shared.test <- sub.shared.test %>%
    pivot_wider(names_from = 'OG.TSG', values_from='is.OT.TSG') %>%
    as.data.frame()
  sub.shared.test <- sub.shared.test[,c(1:(ncol(sub.shared.test)-1))]
  sub.shared.test[1:2,]
  
  tests <- colnames(sub.shared.test)[3:ncol(sub.shared.test)]
  tmp.list <- lapply(tests, function(t){
    # t=tests[1]
    sub <- sub.shared.test[,c('gene_name','timing', t)]
    colnames(sub)[3] <- 'test'
    sub[1:4,]
    freq.sub <- sub %>% filter(!is.na(test)) %>% 
      group_by(timing) %>%
      summarise(count = n(), 
                gene=paste(gene_name[order(gene_name,decreasing = F)],collapse = '\n')) %>% 
      mutate(test = t)
    freq.sub
  })
  freq.tab <- Reduce(rbind, tmp.list) %>%
    left_join(freq.all) %>% mutate(cancer.type = x)
  freq.tab[1:2,]
  return(freq.tab)
})

shared.genes_CAG <- Reduce(rbind, shared.genes_CAG_list) %>%
  as.data.frame() %>%
  mutate(freq = count/count_all) %>%
  mutate(text = paste(count, count_all, sep='/'))

# make plots: Supplenmentary Figure 11B:
annot.color <- c("OG"=brewer.pal(9,'Set1')[1],
                 "TSG"=brewer.pal(9,'Set1')[2],
                 "driver.gene"=brewer.pal(9,'Set1')[3])
names(annot.color)

plot.data.all <- shared.genes_CAG

# earlier:
plot.data_earlier <- plot.data.all %>% filter(grepl('earlier', timing, ignore.case = T)) 
p_bar_earlier <- plot.data_earlier %>%
  ggplot(aes(x = factor(test, levels = names(annot.color)), y = freq, fill = test)) + 
  geom_bar(stat = 'identity', position = 'dodge', colour = 'black') + 
  geom_text(aes(label = gene), position = position_dodge(width = 0.9), vjust = 1.2, size = 2, colour='white') +
  facet_grid(.~ cancer.type, scales = 'free_x', space = 'free_x') + 
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  scale_fill_manual(name = 'Genes', values =annot.color )+
  xlab('') + ylab('% genes with Earlier.timing') + 
  facet_grid(.~ factor(cancer.type,levels = cancer.types), scales = 'free_x', space = 'free_x') + 
  theme_bw() + theme_hz +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_rect(fill = '#fde0ef50'), 
        panel.grid = element_blank(),
        strip.text = element_text(size=12))
p_bar_earlier

#later:
plot.data_later <- plot.data.all %>% filter(grepl('later', timing, ignore.case = T)) 
p_bar_later <- plot.data_later %>%
  ggplot(aes(x = factor(test, levels = names(annot.color)), y = freq, fill = test)) + 
  geom_bar(stat = 'identity', position = 'dodge', colour = 'black') + 
  geom_text(aes(label = gene), position = position_dodge(width = 0.9), vjust = 0, size = 2, colour='white') +
  scale_y_reverse(expand = c(0.1,0,0,0)) + #!
  scale_fill_manual(name = 'Genes', values =annot.color )+
  xlab('') + ylab('% genes with Later.timing') + 
  facet_grid(.~ factor(cancer.type,levels = cancer.types), scales = 'free_x', space = 'free_x') + 
  theme_bw() + theme_hz +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = '#e6f5d050'), 
        panel.grid = element_blank(),
        strip.text = element_text(size=12))
p_bar_later

# align legends
## function:
{
  ### helper functions to align legends with ggplot ###
  
  library(gtable)
  
  #Extract the legend from a ggplot object with this function:
  g_legend <- function(a.gplot){
    tmp    <- ggplot_gtable(ggplot_build(a.gplot))
    leg    <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if(length(leg) == 0){
      return(NULL)
    }
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  #Then use the following function on a list of legends to return a single grob that will contain all the legends aligned:
  align.legends.fun  <- function(legend.list){
    aligned.leg <- legend.list[[1]]
    # Haoran added this if:
    if(length(legend.list) > 1){
      for(i in 2:length(legend.list)){
        leg1        <- legend.list[[i]]$grobs[[1]]
        leg_a       <- gtable_add_rows(aligned.leg, pos = nrow(aligned.leg) - 1, heights = sum(leg1$heights))
        leg_final   <- gtable_add_grob(leg_a, leg1, t = nrow(leg_a) - 1, l = 3)
        aligned.leg <- leg_final
      }
    }
    
    return(aligned.leg)
  }
  
}

# align plots:
{
  p_bar_earlier
  p_bar_later
  
  plot.list <- list(p_bar_earlier+ theme(legend.position = 'none'),
                    p_bar_later)
  
  # 1:length(plot.list)
  legends.list <- lapply(c(2), function(x){
    sub.legend <- g_legend(plot.list[[x]])
  })
  length(legends.list)
  
  library(gtable) #align.legends.fun: >=2 legend.lists
  aligend.legends1 <- align.legends.fun(legends.list[1:length(legends.list)])
  
  tmp.plot.list <- list(tmp1 = plot.list[[1]]+theme(legend.position = 'none')+theme(plot.margin = unit(c(0.5, 0.5, -0.2, 0.2), "cm")),
                        tmp2 = plot.list[[2]]+theme(legend.position = 'none')+theme(plot.margin = unit(c(-0.2, 0.5, 1, 0.2), "cm")))
  
  tmp.plot.list <- lapply(c(1:length(tmp.plot.list)), function(x){
    tmp.plot2 <- ggplotGrob(tmp.plot.list[[x]])
    return(tmp.plot2)
  })
  
  width.list <- lapply(c(1:length(tmp.plot.list)), function(x){
    sub.width <- tmp.plot.list[[x]]$widths
  })
  width.common <- sapply(width.list, grid::unit.pmax)
  
  for (plot.each in tmp.plot.list) {
    plot.each$widths <- width.common
  }
  
  g <- gtable_rbind(tmp.plot.list[[1]],tmp.plot.list[[2]])
  id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
  
  g$heights[id_panels_h] <- grid::unit(c(4,4), "null")
  
  grid.arrange(g, aligend.legends1, ncol = 2, widths = c(9, 1.2))
  
}

# output: Supplenmentary Figure 11B:
pdf(paste0(output_dir,'count.CAG_per.cancer','.pdf'), 
    width = 12, height = 8.5)
grid.arrange(g, aligend.legends1, ncol = 2, widths = c(9, 1.2))
dev.off()





