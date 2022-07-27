############################################################################################################
#############               Hierarchical clustering of Log2(E/L) of 30 cell-lines              ############# 
############################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create Figure 1 of the manuscript "Replication timing alterations shape the genomic and transcriptomic landscape during breast and lung cancer evolution"
# Data accessibility statement can be found in the manuscript.

#library and options
options(stringsAsFactors = F)
library(vroom)
library(ggplot2) 
library(fpc)
library(dplyr)
library(dendextend)
library(ComplexHeatmap)
library(ggpubr)
library(grid)
library(gridExtra)
library(RColorBrewer)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved



############################
#######     Main     #######
############################

#------- Figure 1 A -------#

#load RPKM counts and log2ratio in 50kb windows for T2P
early_RPKM <- as.data.frame(vroom(paste0(data_dir ,'/T2P_Early.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))
late_RPKM  <- as.data.frame(vroom(paste0(data_dir ,'/T2P_Late.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))
log2ratio  <- as.data.frame(vroom(paste0(data_dir, '/T2P.loess300000.bedGraph'), col_names = c('chr', 'start', 'stop', 'score')))

#plot per chromosome (chromsome 3 shown in Figure 1 A)
pdf(paste0(output_dir, '/T2P.RPKMcounts.log2ratio.pdf'), width = 10, height = 3)
for(chr in paste0('chr', c(1:22, 'X'))) {
  sub        <- early_RPKM[early_RPKM[,1]==chr,]
  x.range    <- seq(min(sub$start), max(sub$start), by = as.numeric(WINDOWSIZE))
  full.early <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Early')
  full.early$score[full.early$start %in% sub$start] <- sub$score
  full.early$start <- full.early$start / 1000000
  
  sub        <- late_RPKM[late_RPKM[,1]==chr,]
  full.late <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Late')
  full.late$score[full.late$start %in% sub$start] <- -1*sub$score
  full.late$start <- full.late$start / 1000000
  
  count.data        <- rbind(full.early, full.late)
  count.data$Timing <- factor(count.data$Timing, levels = c('Early', 'Late'))
  
  log2ratio.chr       <- log2ratio[log2ratio$chr == chr,]
  log2ratio.chr$start <- log2ratio.chr$start / 1000000
  
  min.y <- max(min(c(min(log2ratio.chr$score), min(count.data$score))), -10)
  min.y <- ifelse((floor(min.y) %% 2) == 0, floor(min.y), floor(min.y) - 1)
  max.y <- min(max(c(max(log2ratio.chr$score), max(count.data$score))), 10)
  max.y <- ifelse((ceiling(max.y) %% 2) == 0, ceiling(max.y), ceiling(max.y) + 1)
  y.axis.breaks <- c(rev(seq.int(0, min.y, by = -2)), seq.int(0, max.y, by = 2))
  y.axis.breaks <- unique(y.axis.breaks)
  y.axis.labels <- sub('-', '', y.axis.breaks)
  
  p <- ggplot() + 
    geom_area(data = count.data, mapping = aes(x = start, y = score, fill = Timing, group = Timing), position = 'identity') +
    guides(fill=guide_legend(title="RPKM counts")) +
    geom_line(data = log2ratio.chr, mapping = aes(x=log2ratio.chr$start, y = log2ratio.chr$score, colour='black'), size = 0.2) +
    scale_colour_manual(name = '', values = c('black'='black'), labels = c('log2-ratio')) +
    scale_y_continuous(breaks = y.axis.breaks, labels = y.axis.labels, limits = c(min.y, max.y)) +
    theme_bw() + 
    ylab('') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()


#------- Figure 1 B -------#

#load RPKM counts and log2ratio in 50kb windows for T47D
early_RPKM <- as.data.frame(vroom(paste0(data_dir ,'/T47D_Early.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))
late_RPKM  <- as.data.frame(vroom(paste0(data_dir ,'/T47D_Late.filteredRPKM.bedGraph'),  col_names =  c('chr', 'start', 'stop', 'score')))
log2ratio  <- as.data.frame(vroom(paste0(data_dir, '/T47D.loess300000.bedGraph'), col_names = c('chr', 'start', 'stop', 'score')))

#plot per chromosome (chromsome 3 shown in Figure 1 A)
pdf(paste0(output_dir, '/T47D.RPKMcounts.log2ratio.pdf'), width = 10, height = 3)
for(chr in paste0('chr', c(1:22, 'X'))) {
  sub        <- early_RPKM[early_RPKM[,1]==chr,]
  x.range    <- seq(min(sub$start), max(sub$start), by = as.numeric(WINDOWSIZE))
  full.early <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Early')
  full.early$score[full.early$start %in% sub$start] <- sub$score
  full.early$start <- full.early$start / 1000000
  
  sub        <- late_RPKM[late_RPKM[,1]==chr,]
  full.late <- data.frame(chr = chr, start = x.range, score = 0, Timing = 'Late')
  full.late$score[full.late$start %in% sub$start] <- -1*sub$score
  full.late$start <- full.late$start / 1000000
  
  count.data        <- rbind(full.early, full.late)
  count.data$Timing <- factor(count.data$Timing, levels = c('Early', 'Late'))
  
  log2ratio.chr       <- log2ratio[log2ratio$chr == chr,]
  log2ratio.chr$start <- log2ratio.chr$start / 1000000
  
  min.y <- max(min(c(min(log2ratio.chr$score), min(count.data$score))), -10)
  min.y <- ifelse((floor(min.y) %% 2) == 0, floor(min.y), floor(min.y) - 1)
  max.y <- min(max(c(max(log2ratio.chr$score), max(count.data$score))), 10)
  max.y <- ifelse((ceiling(max.y) %% 2) == 0, ceiling(max.y), ceiling(max.y) + 1)
  y.axis.breaks <- c(rev(seq.int(0, min.y, by = -2)), seq.int(0, max.y, by = 2))
  y.axis.breaks <- unique(y.axis.breaks)
  y.axis.labels <- sub('-', '', y.axis.breaks)
  
  p <- ggplot() + 
    geom_area(data = count.data, mapping = aes(x = start, y = score, fill = Timing, group = Timing), position = 'identity') +
    guides(fill=guide_legend(title="RPKM counts")) +
    geom_line(data = log2ratio.chr, mapping = aes(x=log2ratio.chr$start, y = log2ratio.chr$score, colour='black'), size = 0.2) +
    scale_colour_manual(name = '', values = c('black'='black'), labels = c('log2-ratio')) +
    scale_y_continuous(breaks = y.axis.breaks, labels = y.axis.labels, limits = c(min.y, max.y)) +
    theme_bw() + 
    ylab('') + xlab('Chromosome Position (mb)') + 
    ggtitle(paste0('Chromosome ', sub('chr', '', chr)))
  plot(p)
}
dev.off()



#------- Figure 1 C -------#

#read in info to cell-lines
tissue_info <- read.table(paste0(data_dir, '/tissueInfo_cellLines_20210309.tsv'), header = T, sep = '\t')
tissue_info$cellLine <- sub('A549$', 'A549encode', tissue_info$cellLine)
tissue_info$cellLine <- sub('A549rep', 'A549', tissue_info$cellLine)

#read in cell-line mutations
mutTable <- readRDS(paste0(data_dir, '/mutTable_DepMap_20210802.rds'))

#read in log2-ratios for the 30 cell-lines
log2ratio_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))

#heatmap with  euclidean distance and ward.2 criterion
cluster_method  <- 'ward.D2'
dist_matrix     <- dist(t(log2ratio_df[,-1*c(1:3)]))
hc              <- hclust(dist_matrix, method = cluster_method)
mat             <- as.matrix(dist_matrix)
o               <- rownames(mat)
mat             <- mat[hc$order, hc$order]
mat[lower.tri(mat)] <- NA
mat             <- mat[o, o]

#tissue annotations
annotation_col <- data.frame(cellLine = colnames(mat))
annotation_col <- annotation_col %>% left_join(tissue_info[,c('cellLine', 'tissueType', 'tissue', 'Morphology', 'RepliSeq_dataset')])
annotation_col$RepliSeq_dataset <- sub('Haoran', 'IN-HOUSE', annotation_col$RepliSeq_dataset)
colnames(annotation_col) <-c('cellLine', 'TissueType', 'Tissue', 'Morphology', 'Source')
annotation_col$BreastSubtype <- 'not_applicable'
annotation_col$BreastSubtype[annotation_col$cellLine == 'MCF-7']  <- 'ER+'
annotation_col$BreastSubtype[annotation_col$cellLine == 'T47D']   <- 'ER+PR+'
annotation_col$BreastSubtype[annotation_col$cellLine == 'SK-BR3'] <- 'HER2+'
annotation_col$BreastSubtype[annotation_col$cellLine == 'MDA453'] <- 'TNBC'
rownames(annotation_col) <- annotation_col$cellLine
annotation_col <- annotation_col[,-1]
annotation_col <- annotation_col[,rev(colnames(annotation_col))]

tissue_col <- brewer.pal(n = 10, 'Paired')[c(1:6, 9:10, 7:8)]
names(tissue_col) <- sort(unique(annotation_col$Tissue))

morphology_col <- brewer.pal(n = 9, 'Blues')[c(2,4,6,8,9)]
names(morphology_col) <- sort(unique(annotation_col$Morphology))

#driver annotations
driver_genes <- c('TP53','KRAS', 'EGFR', 'PTEN', 'PIK3CA', 'SKT11', 'KEAP1')

driver_matrix <- matrix('not_applicable', ncol = length(driver_genes), nrow = nrow(mat), dimnames = list(rownames(mat), driver_genes))
driver_matrix[rownames(driver_matrix) %in% unique(mutTable$cellLine),] <- 'noMut'
for( x in driver_genes){
  sub <- mutTable[mutTable$Hugo_Symbol %in% x,]
  driver_matrix[rownames(driver_matrix) %in% unique(sub$cellLine), x] <- 'Mut'
}
driver_matrix[rownames(driver_matrix) == 'A549encode',] <- driver_matrix[rownames(driver_matrix) == 'A549',]

#annotation row
annotation_row <- data.frame(cellLine = colnames(mat), AnalysisGroup = NA)
annotation_row$AnalysisGroup[annotation_row$cellLine %in% c('T2P', 'H1650', 'H1792', 'H2009', 'A549')] <- 'LUAD'
annotation_row$AnalysisGroup[annotation_row$cellLine %in% c('HBEC3', 'H2170', 'H520', 'SW900')] <- 'LUSC'
annotation_row$AnalysisGroup[annotation_row$cellLine %in% c('HMEC', 'MDA453', 'SK-BR3', 'MCF-7', 'T47D')] <- 'BRCA'
rownames(annotation_row) <- annotation_row$cellLine
annotation_row <- annotation_row[,-1, drop = F]

annotation_row$colour <- 'black'
annotation_row$colour[ annotation_row$AnalysisGroup == 'LUAD'] <- '#b2182b'
annotation_row$colour[ annotation_row$AnalysisGroup == 'LUSC'] <- '#2166ac'
annotation_row$colour[ annotation_row$AnalysisGroup == 'BRCA'] <- '#fd8d3c'


#annotation colours
ann_colors = list(
  TissueType = c('Normal' = '#80cdc1', 'Cancer' = '#01665e'),
  Tissue = tissue_col,
  Morphology = morphology_col,
  Source = c('IN-HOUSE' = '#ae017e', 'ENCODE' = '#fcc5c0'),
  BreastSubtype = c('ER+' = '#d9f0a3', 'ER+PR+' = '#78c679', 'HER2+' = '#238443', 'TNBC' = '#004529', 'not_applicable' = '#f0f0f0'),
  DriverGenes = c('noMut' = '#bdbdbd', 'Mut' = '#67001f', 'not_applicable' = '#f0f0f0')
)

#plot
p <- ComplexHeatmap::Heatmap(mat, name = 'Euclidian Distance',
                             col = c(rep('#a50026', 50), colorRampPalette(rev(brewer.pal(n = 7, name = "YlOrRd")))(100)),
                             clustering_distance_rows = dist_matrix, clustering_method_rows = cluster_method, row_dend_reorder = F,
                             clustering_distance_columns = dist_matrix,  clustering_method_columns = cluster_method, column_dend_reorder = F,
                             row_split = 2, column_split = 2, row_gap = unit(2, "mm"), column_gap = unit(2, "mm"),
                             show_row_dend = F, show_column_names = T, row_title = ' ', column_title = ' ',
                             rect_gp = gpar(col = "white", lwd = 1), na_col = "white", show_parent_dend_line = FALSE,
                             top_annotation = HeatmapAnnotation(TissueType = annotation_col$TissueType,
                                                                Tissue = annotation_col$Tissue,
                                                                Morphology = annotation_col$Morphology,
                                                                Source = annotation_col$Source,
                                                                BreastSubtype = annotation_col$BreastSubtype,
                                                                DriverGenes = driver_matrix,
                                                                col = ann_colors, na_col = "white",
                                                                gp = gpar(col = "white"), border = T,
                                                                gap = unit(c(0.5,0.5,0.5,0.5,2,2), "mm"),
                                                                simple_anno_size = unit(0.4, "cm")),
                             row_names_gp = gpar(col = annotation_row$colour),
                             column_names_gp = gpar(col = annotation_row$colour))

pdf(paste0(output_dir, 'heatmap_dist_', cluster_method, '_all.pdf'), width = 10, height = 10)
draw(p, merge_legends = T)
dev.off()


#calculate stability of clustering
clusterStability <- clusterboot(dist_matrix, B = 100, bootmethod = 'boot', clustermethod = disthclustCBI, method = cluster_method, k = 2, seed = 123)
write.table(clusterStability$bootmean, file = paste0(output_dir, 'bootmean_jaccard_', cluster_method, '_all.txt'), quote = F, row.names = F, col.names = F)









