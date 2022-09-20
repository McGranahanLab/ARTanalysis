############################################################################################################
#############               Hierarchical clustering of RT signal of 30 cell-lines              ############# 
############################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create Figure 1 of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution".
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
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved



#load BRCA mutation data
load(paste0(data_dir, '/560Breast_subset_mutTable.RData'))

#load BRCA mutation data
load(paste0(data_dir, '/560Breast_subset_cnTable.RData'))

#read in info to cell-lines
tissue_info <- read.table(paste0(data_dir, '/tissueInfo_cellLines_20210309.tsv'), header = T, sep = '\t')
tissue_info$cellLine <- sub('A549$', 'A549encode', tissue_info$cellLine)
tissue_info$cellLine <- sub('A549rep', 'A549', tissue_info$cellLine)

#read in cell-line mutations
mutTable <- readRDS(paste0(data_dir, '/mutTable_DepMap_20210802.rds'))

#read in replication timing signal (log2-ratio) for the 30 cell-lines
repTiming_hg19_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))



#################################
#######     Functions     #######
#################################

#function to calculate mutation frequency in gained, lost and no event regions per sample
mutLoad_cnEvents <- function(sample_mutTable, sample_cnTable, gain_threshold = log2(2.5/2), loss_threshold = log2(1.5/2)){
  
  #calculate log2 value
  sample_cnTable$logCN <- log2(sample_cnTable$nTotal / sample_cnTable$ploidy)
  
  #classify gains and losses based on thresholds
  sample_cnTable$event <- 'no_event'
  sample_cnTable$event[sample_cnTable$logCN >= gain_threshold] <- 'gain'
  sample_cnTable$event[sample_cnTable$logCN <= loss_threshold] <- 'loss'
  
  #overlap mutations with different events
  mut_gr     <- GRanges(seqnames = sample_mutTable$Chrom, IRanges(start = sample_mutTable$Pos, end =  sample_mutTable$Pos))
  
  if(sum(sample_cnTable$event == 'gain') > 0){
    gain_gr    <- makeGRangesFromDataFrame(sample_cnTable[sample_cnTable$event == 'gain',])
    overlap <- findOverlaps(gain_gr, mut_gr)
    gain_mutDensity <- length(unique(subjectHits(overlap))) / sum(width(gain_gr))
  } else {
    gain_mutDensity <- NA
  }
  
  if(sum(sample_cnTable$event == 'loss') > 0){
    loss_gr    <- makeGRangesFromDataFrame(sample_cnTable[sample_cnTable$event == 'loss',])
    overlap <- findOverlaps(loss_gr, mut_gr)
    loss_mutDensity <- length(unique(subjectHits(overlap))) / sum(width(loss_gr))
  } else{
    loss_mutDensity <- NA
  }
  
  if(sum(sample_cnTable$event == 'no_event') > 0){
    noEvent_gr <- makeGRangesFromDataFrame(sample_cnTable[sample_cnTable$event == 'no_event',])
    overlap <- findOverlaps(noEvent_gr, mut_gr)
    noEvent_mutDensity <- length(unique(subjectHits(overlap))) / sum(width(noEvent_gr))
  } else {
    noEvent_mutDensity <- NA
  }
  
  #create output
  output_df <- data.frame(sample = sample_mutTable$Sample[1], gain_mutDensity, loss_mutDensity, noEvent_mutDensity)
  return(output_df)
}


#function to calculate copy number adjusted mutation load in 50kb windows
cnMut_cnTotal.mutationLoad.bin <- function(sample.mutTable, 
                                           chromosomes = paste0('chr', c(1:22)), 
                                           binSize = 50000, 
                                           genome = BSgenome.Hsapiens.UCSC.hg19,
                                           lspan = 500000){
  
  chr.length    <- seqlengths(genome)
  mutLoad.table <- c()
  for(chr in chromosomes){
    chr.size <- chr.length[chr]
    
    #divide chromosome into binSize windows
    seq.bin <- c(seq(from = 1, to = chr.size, by = binSize), chr.size)
    bins_gr <- GRanges(seqnames = chr, IRanges(start = seq.bin[1:(length(seq.bin)-1)], end = seq.bin[2:length(seq.bin)]-1))
    
    #subset mutation table
    sub.mutTable    <- sample.mutTable[sample.mutTable$Chrom == chr,]
    sub.mutTable_gr <- GRanges(seqnames = sub.mutTable$Chrom, IRanges(start = sub.mutTable$Pos, end = sub.mutTable$Pos))
    
    #count mutations in bins
    overlap        <- findOverlaps(bins_gr, sub.mutTable_gr)
    if(length(overlap) == 0){
      count.bins          <- as.data.frame(bins_gr)
      count.bins$mutLoad  <- NA
      count.bins$mutCount <- NA
    } else{
      counts.values   <- lapply(unique(queryHits(overlap)), function(x){
        mut_bin <- sub.mutTable[subjectHits(overlap)[queryHits(overlap) == x],]
        mut_cn_count <- round(mut_bin[,'cnMut']) / mut_bin[,'cnTotal']
        mut_cn_count[mut_bin[,'cnTotal'] == 0] <- 0
        data.frame(mutLoad = sum(mut_cn_count, na.rm = T), mutCount = nrow(mut_bin))
      })
      counts.values <- Reduce(rbind, counts.values)
      rownames(counts.values) <- unique(queryHits(overlap))
      count.bins              <- as.data.frame(bins_gr)
      count.bins$mutLoad <- count.bins$mutCount <- NA
      count.bins$mutCount[as.numeric(rownames(counts.values))] <- as.numeric(counts.values$mutCount)
      count.bins$mutLoad[as.numeric(rownames(counts.values))] <- as.numeric(counts.values$mutLoad)
      count.bins$mutLoad <- (count.bins$mutLoad / count.bins$width) * binSize
    }
    colnames(count.bins)[1] <- 'chr'
    
    
    
    if(!is.na(lspan)){
      #smooth
      chromlspan <- lspan/(max(count.bins[,3])-min(count.bins[,2]))
      
      count.bins$smooth_mutLoad <- predict(loess(mutLoad~start, data = count.bins, span=chromlspan), count.bins$start)
      count.bins$smooth_mutLoad[count.bins$smooth_mutLoad < 0] <- 0
      count.bins$smooth_mutLoad[is.na(count.bins$mutLoad)] <- NA
    }
    
    mutLoad.table <- rbind(mutLoad.table, count.bins)
  }
  
  mutLoad.table$start <- mutLoad.table$start - 1
  return(mutLoad.table)
  
}



############################
#######     Main     #######
############################

#------- Figure 1 A -------#
#--> the same code was used to create Figure 1 B-C for LUAD and LUSC

#calculate mutation frequency in gained, lost and neutral genomic regions per sample
cn_mutLoad_perSample <- lapply(unique(copyNumber_df$sample), function(x){
  print(x)
  sample_cnTable  <- copyNumber_df[copyNumber_df$sample == x,]
  sample_mutTable <- mutTable[grep(x, mutTable$Sample),]
  
  mutDensity <- data.frame(mutLoad_cnEvents(sample_mutTable, sample_cnTable), nMuts = nrow(sample_mutTable))
  return(mutDensity)
})
cn_mutLoad_perSample <- Reduce(rbind, cn_mutLoad_perSample)

#boxplot and paired t-test
plot_data <- reshape2::melt(cn_mutLoad_perSample, id.vars = c('sample', 'nMuts'))
plot_data$variable <- sub('_mutDensity', '', plot_data$variable)
plot_data$variable <- sub('noEvent', 'neutral', plot_data$variable)
plot_data$variable <- factor(plot_data$variable, levels = c('loss', 'neutral', 'gain'))

pdf(paste0(output_dir, 'boxplot_mutDensity_cnEvents.pdf'), width = 4, height = 4, useDingbats = F)
ggplot(plot_data, aes(x = variable, y = value)) +
  geom_line(aes(group = sample), colour = '#e0e0e0') +
  ggbeeswarm::geom_quasirandom(aes(colour = variable), width = 0.3, alpha = 0.8) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  stat_compare_means(paired = T, comparisons = list(c('loss', 'neutral'), c('neutral', 'gain')), method.args = list(alternative = "less")) +
  scale_colour_manual(values = c('loss' = '#2166ac', 'neutral' = '#878787', 'gain' = '#b2182b')) +
  xlab('') + ylab('Mutation Density') +
  labs(title = "Mutation Densities in Regions\nwith Copy Number Events",
       subtitle = '(Paired t-test with alternative less)') +
  scale_y_continuous(trans = 'log10') +
  theme_bw() +
  theme(legend.position = 'none', panel.grid = element_blank())
dev.off()


#------- Figure 1 D -------#

#barplot with number of tumours
plot_data <- data.frame(cancerType = c('BRCA', 'LUAD', 'LUSC'), 
                        value = c(482, 470, 319))

pdf(paste0(output_dir, 'bar_nTumours_cancerType.pdf'), width = 4, height = 3, useDingbats = F)
ggplot(plot_data, aes(x = cancerType, y = value, fill = cancerType)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = c('LUAD' = '#4292c6', 'LUSC' = '#41ab5d', 'BRCA' = '#fd8d3c'), guide = F) +
  geom_text(aes(label = value), colour = 'white', position = position_stack(vjust = 0.5)) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  xlab('') + ylab(' Number of tumours') +
  theme_classic()
dev.off()


#barplot with mutation load in 5Mb bins across the genome
#--> the same code was used to plot the mutation load across the genome for LUAD and LUSC
BRCA_mutLoad <- cnMut_cnTotal.mutationLoad.bin(mutTable, binSize = 5000000, lspan = NA)
colnames(BRCA_mutLoad)[3] <- 'stop'

plot_data     <- rbind(data.frame(BRCA_mutLoad, cancerType = 'BRCA'))
plot_data$bin <- paste(plot_data$chr, plot_data$start, plot_data$stop, sep = ':')
bins_to_use   <- plot_data$bin #--> overlap bins from LUAD, LUSC and BRCA and only use those. The bins vary slightly because of the genome build hg19 used for BRCA and hg38 used for LUAD and LUSC

plot_data$chr <- factor(plot_data$chr, levels = paste0('chr', c(1:22)))
plot_data$bin <- factor(paste(plot_data$chr, plot_data$start, plot_data$stop, sep = ':'), levels = bins_to_use)
plot_data$col <- 'straight'
plot_data$col[plot_data$chr %in% c(paste0('chr', c(1,3,5,7,9,11,13,15,17,19,21)))] <- 'odd'

#get postions for lines and text
temp <- plot_data[plot_data$cancerType == 'BRCA',]
pos.lines <- cumsum(table(temp$chr))
index.tick <- round(as.numeric(pos.lines - (table(temp$chr) / 2)))

#create plot
pdf(paste0(output_dir, 'mutLoad_alongGenome_5Mb_BRCA.pdf'), width = 10, height = 4)
ggplot(plot_data, aes(x = bin, y = mutLoad, fill = col)) + 
  geom_bar(stat = 'identity') + 
  geom_vline(xintercept = pos.lines + 0.5, linetype = 'dashed', size = 0.2) +
  facet_grid(cancerType ~ ., scale = "free") + 
  scale_fill_manual(values = c('straight' = '#4d4d4d', 'odd' = '#bababa'), guide = F) +
  scale_x_discrete(breaks = levels(plot_data$bin)[index.tick], labels = paste0('chr', c(1:22))) + 
  scale_y_continuous(expand  = c(0,0)) +
  xlab('5Mb Bins') + ylab('Mutation Load') +
  ggtitle('') +
  theme_bw() + 
  theme(legend.position = 'none', panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#------- Figure 1 E -------#
#identify consereved replication timing regions only using ENCODE data (with NA values allowed)
normal_cellLines <- c("HBEC3", "T2P", "TT1", "HMEC", "MCF10A", "BG02", "BJ", "HUVEC", "IMR90", "keratinocyte")
ENCODE_cellLines <- c("GM06990", "GM12801", "GM12812", "GM12813", "GM12878", "K562", "SK-N-MC", "SK-N-SH",  "MCF-7", "T47D", "HeLa-S3", "BG02", "HUVEC",       
                      "Caki2", "G401", "HepG2", "A549encode", "H460", "IMR90",  "LNCAP", "BJ", "keratinocyte")
ENCODE_cellLines <- ENCODE_cellLines[ENCODE_cellLines %in% colnames(repTiming_hg19_df)]
normal_cellLines <- normal_cellLines[normal_cellLines %in% ENCODE_cellLines]

log2ratio_ENCODE <- repTiming_hg19_df[,c('chr', 'start', 'stop', ENCODE_cellLines)]
log2ratio_ENCODE$mean <- rowMeans(log2ratio_ENCODE[,ENCODE_cellLines], na.rm = T)

normal_log2ratio_ENCODE <- log2ratio_ENCODE[,c('chr', 'start', 'stop', normal_cellLines)]
normal_log2ratio_ENCODE$repTiming           <- 'unknown'
ww_nonConserved                         <- apply(normal_log2ratio_ENCODE[,normal_cellLines], 1, function(x) sum(sign(x[!is.na(x)]) > 0) > 0 & sum(sign(x[!is.na(x)]) < 0) > 0)
normal_log2ratio_ENCODE$repTiming[ww_nonConserved] <- 'non conserved'
ww_early                                <- apply(normal_log2ratio_ENCODE[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) > 0))
normal_log2ratio_ENCODE$repTiming[ww_early] <- 'early'
ww_late                                 <- apply(normal_log2ratio_ENCODE[,normal_cellLines], 1, function(x) sum(is.na(x)) == 0 && all(sign(x) < 0))
normal_log2ratio_ENCODE$repTiming[ww_late]  <- 'late'

log2ratio_ENCODE <- log2ratio_ENCODE %>% left_join(normal_log2ratio_ENCODE[,c('chr', 'start', 'stop', 'repTiming')])

#plot pie chart
plot_data <- normal_log2ratio_ENCODE %>%
  filter(repTiming != 'unknown') %>%
  group_by(repTiming) %>%
  summarise(count = n()) %>%
  mutate(fraction = count / sum(count) * 100) %>%
  mutate(cumulative = cumsum(fraction),
         midpoint = cumulative - fraction / 2,
         label = paste0(round(fraction,1), '%'))


pdf(paste0(output_dir, 'pie_conservedRT_normalCells_ENCODE.pdf'), width = 5, height = 4)
ggplot(plot_data, aes(x = 1, weights = fraction, fill = repTiming)) +
  geom_bar(width = 1, position = position_stack(reverse = T), colour = 'white') +
  coord_polar(theta = "y") +
  scale_fill_manual(name = 'Replication Timing', values = c('early' = '#c51b7d', 'late' = '#4d9221', 'non conserved' = '#bababa')) +
  geom_text(aes(x = rev(c(0.9, 1, 1.1)), y = midpoint, label = label), size = 5, colour = c('white', 'white', 'black')) + 
  xlab('') + ylab('') +
  ggtitle(paste0('Conserved Replication Timing\nin Normal Cells')) +
  theme_void() + 
  theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = 'top')
dev.off()



#------- Figure 1 F -------#

#calculate mutation load in 50kb bins
#--> the same code was used for LUAD and LUSC
BRCA_mutLoad <- cnMut_cnTotal.mutationLoad.bin(mutTable, binSize = 50000, lspan = NA)

#explained variance of model between mutation load and consereved or non consereved regions
#--> the same code was used to calculate and plot the variance explained in LUAD and LUSC
corr_BRCA <- log2ratio_ENCODE %>%
  select(chr, start, stop, mean, repTiming) %>%
  dplyr::filter(repTiming %in% c('early', 'late', 'non conserved')) %>% 
  left_join(BRCA_mutLoad, by = c('chr', 'start', 'stop' = 'end'))

corr_conserved_BRCA    <- summary(lm(mutLoad ~ mean, data = corr_BRCA[corr_BRCA$repTiming != 'non conserved',]))$r.squared
corr_nonconserved_BRCA <- summary(lm(mutLoad ~ mean, data = corr_BRCA[corr_BRCA$repTiming == 'non conserved',]))$r.squared

#plot
plot_data <- rbind(data.frame(cancerType = 'BRCA', type = c('conserved', 'non conserved'), value = c(corr_conserved_BRCA, corr_nonconserved_BRCA)))

pdf(paste0(output_dir, 'bar_r2_mutLoad_conservedRT_ENCODEonly.pdf'), width = 5, height = 4)
ggplot(plot_data, aes(x = cancerType, y = value, fill = type)) + 
  geom_bar(stat = 'identity', position = 'dodge', colour = 'white') + 
  scale_y_continuous(expand = c(0,0, 0.01, 0)) +
  scale_fill_manual(name = 'RT in normal cells', values = c('conserved' = 'black', 'non conserved' = 'gray')) +
  xlab('') + ylab('variance in mutation load explained by log2(E/L)') +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()



#------- Figure 1 G -------#

#heatmap with  euclidean distance and ward.2 criterion
log2ratio_df    <- repTiming_hg19_df
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
annotation_row$colour[ annotation_row$AnalysisGroup == 'LUAD'] <- '#4292c6'
annotation_row$colour[ annotation_row$AnalysisGroup == 'LUSC'] <- '#41ab5d'
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









