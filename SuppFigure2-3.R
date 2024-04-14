############################################################################################
#############               Find thresholds for ART indetification             ############# 
############################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create SuppFigure 2 - 3 of the manuscript 
# "Replication timing alterations are associated with mutation acquisition 
# during tumour evolution in breast and lung cancer".
# Data accessibility statement can be found in the manuscript.

#library and options
options(stringsAsFactors = F)
library(dplyr)
library(ggplot2)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(vroom)

#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

#paired replicates
replicates <- c('T2P' = 'T2Prep', 'H1650' = 'H1650rep', 'A549' = 'A549encode')

#load altered replication timing regions
ARTregions_list <- readRDS(paste0(data_dir, '/ARTregions_full.rds'))

#read in info to cell-lines
tissue_info <- read.table(paste0(data_dir, '/tissueInfo_cellLines_20210309.tsv'), header = T, sep = '\t')
tissue_info$cellLine <- sub('A549$', 'A549encode', tissue_info$cellLine)
tissue_info$cellLine <- sub('A549rep', 'A549', tissue_info$cellLine)

#load RT values in 50kb windows
repTiming_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))
normal_cellLines <- colnames(repTiming_df)[colnames(repTiming_df) %in% tissue_info$cellLine[tissue_info$tissueType == 'Normal']]




##########################
###     Functions      ###    
##########################

#function to calculate correlation between 2 replicates and plot 2d density plot
fun_corr_replicates <- function(rep1_file, rep2_file, title = NULL){
  rep1 <- data.frame(vroom(rep1_file, col_names = c('chr', 'start', 'stop', 'l2r.value')))
  rep2 <- data.frame(vroom(rep2_file, col_names = c('chr', 'start', 'stop', 'l2r.value')))
  
  plot_data <- rep1 %>% left_join(rep2, by = c('chr', 'start', 'stop'))
  
  ggplot(plot_data, aes(x = l2r.value.x, y = l2r.value.y)) + 
    geom_bin2d(show.legend = F) + 
    stat_cor() +
    geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
    scale_fill_gradient(low = 'white', high = '#67001f') +
    xlab('Replicate 1') + ylab('Replicate 2') +
    ggtitle(title) + 
    theme_bw() + theme(panel.grid = element_blank())
}

#####################
###     Main      ###    
#####################

#--------- SuppFigure 2 A-C---------#
p_raw <- lapply(1:length(replicates), function(i){
  rep1 <- names(replicates[i])
  rep2 <- as.character(replicates[i])
  
  rep1_file <- paste0(data_dir, rep1, '.l2r.bedGraph')
  rep2_file <- paste0(data_dir, rep2, '.l2r.bedGraph')
  
  title <- paste0(rep1, ' loess smoothed Log2-Ratios')
  fun_corr_replicates(rep1_file, rep2_file, title)
})

p_qnorm <- lapply(1:length(replicates), function(i){
  rep1 <- names(replicates[i])
  rep2 <- as.character(replicates[i])
  
  rep1_file <- paste0(data_dir, rep1, '.qnorm_l2r.bedGraph')
  rep2_file <- paste0(data_dir, rep2, '.qnorm_l2r.bedGraph')
  
  title <- paste0(rep1, ' loess smoothed Log2-Ratios')
  fun_corr_replicates(rep1_file, rep2_file, title)
})

p_loess <- lapply(1:length(replicates), function(i){
  rep1 <- names(replicates[i])
  rep2 <- as.character(replicates[i])
  
  rep1_file <- paste0(data_dir, rep1, '.loess300000.bedGraph')
  rep2_file <- paste0(data_dir, rep2, '.loess300000.bedGraph')
  
  title <- paste0(rep1, ' loess smoothed Log2-Ratios')
  fun_corr_replicates(rep1_file, rep2_file, title)
})


#plot
pdf(paste0(output_dir, 'corr_replicates_50kb.pdf'), width = 4, height = 4)
p_raw
p_qnorm
p_loess
dev.off()



#--------- SuppFigure 3 A---------#
log2ratio_df <- repTiming_df
fract_RT <- lapply(colnames(log2ratio_df)[4:ncol(log2ratio_df)], function(x){
  sub <- log2ratio_df[,x] 
  sub <- sub[!is.na(sub)]
  data.frame(cellLine = x, early = round(sum(sub > 0) / length(sub), 2), late = round(sum(sub < 0) / length(sub), 2))
})
fract_RT <- Reduce(rbind, fract_RT)

fract_RT$cancer <- 'Cancer'
fract_RT$cancer[fract_RT$cellLine %in% normal_cellLines] <- 'Normal'

fract_RT <- fract_RT[order(fract_RT$late),]
fract_RT <- fract_RT[order(fract_RT$cancer),]

plot_data <- reshape2::melt(fract_RT, id.vars = c('cellLine', 'cancer'))
plot_data$cellLine <- factor(plot_data$cellLine, levels = fract_RT$cellLine)

colour_cellLines <- setNames(rep('black', length(plot_data$cellLine)), plot_data$cellLine)
colour_cellLines[names(colour_cellLines) %in% tissue_info$cellLine[tissue_info$RepliSeq_dataset == 'ENCODE']] <- '#636363'

pdf(paste0(output_dir, 'bar_fractionRT_allCellLines.pdf'), width = 10, height = 4)
ggplot(plot_data, aes(x = cellLine, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', colour = 'black') + 
  geom_text(aes(label = value), colour = 'white', position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(name = 'Replication Timing', values = c('early' = '#c51b7d', 'late' = '#4d9221')) +
  facet_grid(.~cancer, scales = 'free_x', space = 'free_x') + 
  ylab('Proportion of genome') +
  xlab('') +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = colour_cellLines), legend.position = 'top')
dev.off()





#--------- SuppFifure 3 B---------#
#read in log2-ratios
for(i in 1:length(replicates)){
  rep1 <- names(replicates[i])
  rep2 <- as.character(replicates[i])
  
  rep1_file <- paste0(data_dir, rep1, '.loess300000.bedGraph')
  rep2_file <- paste0(data_dir, rep2, '.loess300000.bedGraph')
  
  rep1 <- data.frame(vroom(rep1_file, col_names = c('chr', 'start', 'stop', rep1)))
  rep2 <- data.frame(vroom(rep2_file, col_names = c('chr', 'start', 'stop', rep2)))
  
  combined_df <- rep1 %>% full_join(rep2, by = c('chr', 'start', 'stop'))
  
  if(i == 1){
    log2ratio_df <- combined_df
  } else {
    log2ratio_df <- log2ratio_df %>% full_join(combined_df, by = c('chr', 'start', 'stop'))
  }
}

#calculate difference between replicates
diff_df <- lapply(1:length(replicates), function(i){
  data      <- log2ratio_df[,c('chr', 'start', 'stop', names(replicates)[i], replicates[i])]
  data$diff <- data[,names(replicates)[i]] - data[,replicates[i]]
  data      <- data[!is.na(data$diff),]
  data$timing <- 'earlier'
  data$timing[data$diff > 0] <- 'later'
  data$cellLine <-  names(replicates)[i]
  data[,-1*c(4:5)]
})
diff_df <- Reduce(rbind, diff_df)


#calculate 99%-CI (using qunatiles)
quantiles_df <- diff_df %>%
  group_by(timing) %>%
  summarise(lowCI = quantile(diff, 0.01),
            highCI = quantile(diff, 0.99))
#--> use |2| threhsold?

min_max_df <- diff_df %>%
  group_by(timing) %>%
  summarise(lowCI = min(diff),
            highCI = max(diff))


#plot densities
plot_data <- diff_df
plot_data$cellLine <- 'all'
plot_data <- rbind(plot_data, diff_df)
plot_data1 <- plot_data
plot_data1$timing <- 'all'
plot_data <- rbind(plot_data, plot_data1)
plot_data$cellLine <- factor(plot_data$cellLine, levels = c('all', 'T2P', 'H1650', 'A549'))

pdf(paste0(output_dir, 'density_diffReplicates.pdf'), width = 8, height = 5)
ggplot(plot_data, aes(x = diff, fill = cellLine, colour = cellLine)) + 
  geom_density(alpha = 0.3) + 
  facet_grid(timing ~.) + 
  geom_vline(xintercept = c(-2, 2), linetype = 'dotted', colour = '#b2182b') +
  xlab('diff(log2ratio)') + ggtitle('Difference in Log2-ratio between Replicates') +
  theme_bw()
dev.off()




#--------- SuppFifure 3 C---------#
ARTregions_df <- Reduce(rbind, ARTregions_list)
ARTregions_df$cancerType <- 'LUAD'
ARTregions_df$cancerType[ARTregions_df$normal_cancer %in% paste0('HBEC3-', c('H520', 'H2170', 'SW900'))] <- 'LUSC'
ARTregions_df$cancerType[ARTregions_df$normal_cancer %in% paste0('HMEC-', c('SK-BR3', 'MCF-7', 'T47D', 'MDA453'))] <- 'BRCA'

#create polygons to plot density and colour ARTclass
dens_list <- lapply(c('LUAD', 'LUSC', 'BRCA'), function(x){
  dens <- density(ARTregions_df$diff[ARTregions_df$cancerType == x])
  dens <- data.frame(x = dens$x, y = dens$y)
  dens$ARTclass <- 'not_altered'
  dens$ARTclass[dens$x < -2] <- 'earlier'
  dens$ARTclass[dens$x > 2]   <- 'later'
  
  dens <- do.call("rbind", lapply(split(dens, dens$ARTclass), function(df) {
    df <- rbind(df[1,], df, df[nrow(df),])
    df$y[c(1, nrow(df))] <- 0
    df
  }))
  
  dens$cancerType <- x
  return(dens)
})
dens_df <- Reduce(rbind, dens_list)

#plot
pdf(paste0(output_dir, 'density_diffL2R_cancerTypes.pdf'), width = 6, height = 6)
ggplot(dens_df, aes(x, y)) + 
  geom_polygon(aes(fill = factor(ARTclass), color = factor(ARTclass))) +
  scale_fill_manual(name = "ARTclass", values = c('earlier' = '#de77ae', 'later' = '#7fbc41', 'not_altered' = '#bababa')) +
  scale_colour_manual(name = "ARTclass", values = c( 'earlier' = '#de77ae', 'later' = '#7fbc41', 'not_altered' = '#bababa'), guide = guide_none()) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = c(-2,2), linetype = 'dotted') +
  facet_grid(cancerType~.) +
  ylab('density') + 
  xlab('difference log2-ratio between normal and cancer') +
  ggtitle('Density of Log2-Ratio Differences') +
  theme_bw()
dev.off() 

