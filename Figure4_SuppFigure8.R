########################################################################################################################
#############               Mutation distribution in altered replication timing (ART) regions              ############# 
########################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create Figure 3 and SuppFigure 7  of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution".
# --> BRCA data used as an example (LUAD and LUSC data Data from the 100,000 Genomes Project are held in a secure research environment and are available to registered users. 
#     See https://www.genomicsengland.co.uk/research/academic for further information.)
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
library(ggpubr)
library(ComplexUpset)
library(signal)
library(ggnewscale)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

#load BRCA mutation data
load(paste0(data_dir, '/560Breast_subset_mutTable.RData'))

#download Supplementary Table 1 from  https://www.nature.com/articles/nature17676
clinical_data <- read.table(paste0(data_dir, '/Supplementary Table 1 CLINICAL.PATHOLOGY.DATA.FREEZE.ANALYSIS.v4.032015.csv', sep = ',',  header = T))
colnames(clinical_data) <- clinical_data[1,]
clinical_data           <- clinical_data[-1,]
clinical_data           <- clinical_data[clinical_data$sample_name %in% mutTable$patient,]

#load shared altered replication timing
overlap_repTiming <- read.table(paste0(data_dir, '/BRCA_sharedARTregions.txt'), header = T)
recurrentART      <- read.table(paste0(data_dir, '/BRCA_consistentARTregions.txt'), header = T)

#load log2(E/L) values in 50kb windows
repTiming_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))



#################################
#######     Functions     #######
#################################

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



#function to estimate mean mutLoad distribution in repTiming regions
meanMutload_distribution_sampling <- function(mutLoad_df, altered_repTiming, iter = 10000, nMuts = NULL, zTrans = F){
  
  set.seed(123)
  
  #add mutation load
  if(zTrans){
    data <- altered_repTiming %>% 
      left_join(mutLoad_df, by = c('chr', 'start', 'stop' = 'end')) %>%
      dplyr::filter(!is.na(mutLoad)) %>%
      mutate(mutLoad = (mutLoad - mean(mutLoad)) / sd(mutLoad))
  } else {
    if(is.null(nMuts)){
      data <- altered_repTiming %>% 
        left_join(mutLoad_df, by = c('chr', 'start', 'stop' = 'end')) %>%
        dplyr::filter(!is.na(mutLoad)) %>%
        mutate(mutLoad = mutLoad / sum(mutLoad))
    } else {
      data <- altered_repTiming %>% 
        left_join(mutLoad_df, by = c('chr', 'start', 'stop' = 'end')) %>%
        dplyr::filter(!is.na(mutLoad)) %>%
        mutate(mutLoad = mutLoad / nMuts)
    }
  }
  
  #use minumum number of bins for sampling
  nSample <- min(table(data$timing))
  
  #iterations 
  meanMutLoad_iterations <- sapply(unique(data$timing), function(t){
    mutLoad <- data$mutLoad[data$timing == t]
    iterations <- sapply(1:iter, function(i){
      mean(sample(mutLoad, nSample, replace = T))
    })
    df <- data.frame(iterations)
    return(df)
  })
  meanMutLoad_iterations <- Reduce(cbind, meanMutLoad_iterations)
  meanMutLoad_iterations <- data.frame(meanMutLoad_iterations)
  colnames(meanMutLoad_iterations) <- unique(data$timing)
  
  return(meanMutLoad_iterations)
  
}


#function to identify replication timing events
#--> group adjacent bins with same replication timing change together
identify_repTimingEvents <- function(repTiming_df, maxGapBin = 1){
  
  width <- repTiming_df$stop[1] - repTiming_df$start[1]
  
  #combine bins
  ARTevents <- lapply(unique(repTiming_df$timing), function(x){
    sub       <- repTiming_df[repTiming_df$timing == x,]
    
    if(!is.null(maxGapBin)){
      ww        <- which((sub$start[2:nrow(sub)] - sub$stop[1:(nrow(sub)-1)]) == (maxGapBin * width))
      ww        <- ww[sub$chr[ww] == sub$chr[ww+1]]
      sub$stop[ww] <- sub$start[ww+1] 
    }
    
    sub_gr       <- reduce(makeGRangesFromDataFrame(sub))
    
    #return
    df <- data.frame(sub_gr)
    df$timing <- x
    return(df)
  })
  ARTevents <- Reduce(rbind, ARTevents)
  
  return(ARTevents)
}


#function to calculate average mutLoad and log2-ratio around certain area around a center point
averageBehaviour_genomicArea <- function(centers_df, repTiming_normal_cancer, mutLoad_df, windowLength = 500000, binsize = 50000,
                                         title = 'Average behaviour around center', flip_mutLoad = FALSE){
  
  #create vector of positions based on windowLength and binSize
  seq_bins <- seq(-1 * windowLength, windowLength, binsize)
  
  #assign values in bin around center position
  cancer_log2ratio <- normal_log2ratio <- mutLoad_matrix <- matrix(NA, nrow = nrow(centers_df), ncol = length(seq_bins),
                                                                   dimnames = list(paste0(centers_df[,1], ':', centers_df[,2], '-', centers_df[,3]), seq_bins))
  
  for(i in 1:nrow(centers_df)){
    bin <- paste0(centers_df[i,1], ':', centers_df[i,2], '-', centers_df[i,3])
    df  <- data.frame(index = seq_bins, chr = centers_df[i,1], start = seq_bins + centers_df[i,2], end = seq_bins + centers_df[i,2] + binsize)
    df  <- df %>% left_join(repTiming_normal_cancer, by = c('chr' = 'chr', 'start' = 'start', 'end' = 'stop'))
    
    cancer_log2ratio[bin,] <- df$cancer
    normal_log2ratio[bin,] <- df$normal
    
    df <-  df %>% left_join(mutLoad_df[,c('chr', 'start', 'end', 'mutLoad')], by = c('chr', 'start', 'end'))
    mutLoad_matrix[bin,] <- df$mutLoad
  }
  
  #calculate average across postions
  average_df <- data.frame(bin_index = seq_bins, 
                           cancer_log2ratio = colMeans(cancer_log2ratio, na.rm = T),
                           normal_log2ratio = colMeans(normal_log2ratio, na.rm = T),
                           mutLoad = colMeans(mutLoad_matrix, na.rm = T))
  rownames(average_df) <- NULL
  output <- average_df
  
  #calculate standard error across positions
  se_df <- data.frame(bin_index = seq_bins, 
                      cancer_log2ratio = apply(cancer_log2ratio, 2, function(x) sd(x, na.rm = T) / sqrt(sum(!is.na(x)))),
                      normal_log2ratio = apply(normal_log2ratio, 2, function(x) sd(x, na.rm = T) / sqrt(sum(!is.na(x)))),
                      mutLoad = apply(mutLoad_matrix, 2, function(x) sd(x, na.rm = T) / sqrt(sum(!is.na(x)))))
  rownames(se_df) <- NULL
  output_se <- se_df
  
  if(flip_mutLoad){
    average_df$mutLoad <- -1 * average_df$mutLoad
  }
  
  average_df[,2:ncol(average_df)] <- apply(average_df[,2:ncol(average_df)], 2, function(x) (x - mean(x, na.rm = T)) / sd(x, na.rm = T))
  average_df[,2:ncol(average_df)] <- apply(average_df[,2:ncol(average_df)], 2, function(x) sgolayfilt(x, p = 3, n = 7, m = 0))
  
  se_df[,2:ncol(se_df)] <- apply(se_df[,2:ncol(se_df)], 2, function(x) sgolayfilt(x, p = 3, n = 7, m = 0))
  
  
  #plot
  plot_data <- reshape2::melt(average_df, id.vars = 'bin_index') %>%
    left_join(reshape2::melt(se_df, id.vars = 'bin_index'), by = c('bin_index', 'variable'))
  colnames(plot_data) <- c('bin_index', 'variable', 'mean', 'se')
  
  plot_data$ymin <- plot_data$mean - plot_data$se
  plot_data$ymax <- plot_data$mean + plot_data$se
  
  plot_data_area <- rbind(data.frame(start = seq_bins[1], end = 0 - 10000, type = 'pre'),
                          data.frame(start = 0 - 10000, end = 0 + 10000, type = 'event'),
                          data.frame(start =0 + 10000, end = seq_bins[length(seq_bins)], type = 'post')) 
  
  p <- ggplot() + 
    geom_rect(data = plot_data_area, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = type), alpha = 0.3) +
    scale_fill_manual(name = '', values = c('pre' = '#fdae61', 'event' = '#878787', 'post' = '#fee090'), guide = 'none') +
    new_scale('fill') +
    geom_ribbon(data = plot_data, aes(x = bin_index, ymin = ymin, ymax = ymax, fill = variable), alpha = 0.2) +
    scale_fill_manual(name = '', values = c('mutLoad' = 'black', 'cancer_log2ratio' = '#b2182b', 'normal_log2ratio' = '#2166ac'),
                      labels = c('mutLoad' = 'mutation load', 'cancer_log2ratio' = 'cancer log2-ratio', 'normal_log2ratio' = 'normal log2-ratio')) +
    geom_line(data = plot_data, aes(x = bin_index, y = mean, colour = variable, linetype = variable)) + 
    scale_x_continuous(expand = c(0,0), breaks = c(seq_bins[1], 0, seq_bins[length(seq_bins)]), 
                       labels = c(paste0(-1*(windowLength/1000), 'kb'), 'center', paste0((windowLength/1000), 'kb'))) + 
    scale_colour_manual(name = '', values = c('mutLoad' = 'black', 'cancer_log2ratio' = '#b2182b', 'normal_log2ratio' = '#2166ac'),
                        labels = c('mutLoad' = 'mutation load', 'cancer_log2ratio' = 'cancer log2-ratio', 'normal_log2ratio' = 'normal log2-ratio')) +
    scale_linetype_manual(name = '', values = c('mutLoad' = 'dashed', 'cancer_log2ratio' = 'solid', 'normal_log2ratio' = 'solid'),
                          labels = c('mutLoad' = 'mutation load', 'cancer_log2ratio' = 'cancer log2-ratio', 'normal_log2ratio' = 'normal log2-ratio')) +
    
    xlab('') + ylab('') +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) 
  
  #output
  return(list(mean = output, se = output_se, plot = p))
  
}



############################
#######     Main     #######
############################

#--------- Figure 7 A ---------#
BRCA_mutLoad <- cnMut_cnTotal.mutationLoad.bin(mutTable, lspan = NA)

plot_data <- BRCA_mutLoad %>% 
  left_join(overlap_repTiming, by = c('chr', 'start', 'end' = 'stop')) %>%
  dplyr::filter(!is.na(mutLoad)) %>%
  dplyr::filter(!is.na(timing)) %>%
  mutate(timing = factor(timing, levels = c('early', 'earlier', 'later', 'late')))
  
pdf(paste0(output_dir, 'violin_mutLoad_BRCA.pdf'), width = 7, height = 5)
ggplot(plot_data, aes(x = timing, y = mutLoad, fill = timing)) + 
  geom_violin() + 
  geom_boxplot(outlier.shape = NA, width = 0.2, alpha = 0.5, fill = 'white') + 
  stat_compare_means(method = 'wilcox', comparisons = list(c('early', 'late'), c('earlier', 'later'), c('early', 'later'), c('late', 'earlier')),
                     label = 'p.signif') +
  scale_fill_manual(name = 'Replication Timing', values = brewer.pal(n = 11, 'PiYG')[c(1,3,9,11)]) +
  scale_y_continuous(trans = 'log10') +
  xlab('') + ylab('mutation load in 50kb bins') +
  theme_bw() + theme(legend.position = 'none', axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#--------- Figure 7 B ---------#
consistent_repTiming_df <- recurrentART
cellLines               <- c('HMEC', 'MCF-7', 'MDA453', 'SK-BR3', 'T47D')
repTiming_df            <- repTiming_df[,c('chr', 'start', 'stop', cellLines)]
mutLoad                 <- BRCA_mutLoad

BRCA_repTimingDomains <- identify_repTimingEvents(consistent_repTiming_df)
repTimingDomains      <- BRCA_repTimingDomains
earlier_domains  <- repTimingDomains[repTimingDomains$timing %in% 'earlier',]
later_domains    <- repTimingDomains[repTimingDomains$timing %in% 'later',]
early_domains    <- repTimingDomains[repTimingDomains$timing %in% 'early',]
late_domains     <- repTimingDomains[repTimingDomains$timing %in% 'late',]

repTiming_normal_cancer        <- repTiming_df
repTiming_normal_cancer$normal <- repTiming_normal_cancer[,cellLines[1]]
repTiming_normal_cancer$cancer <- rowMeans(repTiming_normal_cancer[,cellLines[-1]], na.rm = T)

average_earlierDomains <- averageBehaviour_genomicArea(earlier_domains, repTiming_normal_cancer, mutLoad, title = 'Average behaviour at recurrent \nlate to early ART regions')
average_laterDomains   <- averageBehaviour_genomicArea(later_domains, repTiming_normal_cancer, mutLoad, title = 'Average behaviour at recurrent \nearly to late ART regions')
average_earlyDomains   <- averageBehaviour_genomicArea(early_domains, repTiming_normal_cancer, mutLoad, title = 'Average behaviour at consistent \nearly RT regions')
average_lateDomains    <- averageBehaviour_genomicArea(late_domains, repTiming_normal_cancer, mutLoad, title = 'Average behaviour at consistent \nlate ART regions')

pdf(paste0(output_dir, 'averageBehvaiour_repTimingDomains_BRCA.pdf'), width = 5, height = 4)
average_earlierDomains$plot
average_laterDomains$plot
average_earlyDomains$plot
average_lateDomains$plot
dev.off()



#--------- Figure 7 C ---------#
#calculate mutation load
clonal_mutTable   <- mutTable[!is.na(mutTable$absolute.ccf),]
clonal_mutTable   <- clonal_mutTable[clonal_mutTable$absolute.ccf > 0.95,]
mutLoad           <- cnMut_cnTotal.mutationLoad.bin(clonal_mutTable, lspan = NA)

# sample 1000 bins from each category 10000 times (boostrapping) #
meanMutLoad_df <- meanMutload_distribution_sampling(mutLoad, overlap_repTiming, zTrans = T)

#mirrored densityplot
colours <- brewer.pal(n = 11, 'PiYG')[c(2,3,9,10)]

pdf(paste0(output_dir, 'BRCA_mirrorDist_meanMutload.pdf'), width = 5, height = 5)
ggplot(meanMutLoad_df) + 
  geom_density(aes(x = early, y = ..density..), fill = colours[1], colour = colours[1], alpha = 0.8) +
  geom_density(aes(x = late, y = ..density..), fill = 'white', colour = colours[4], linetype = 'dashed') +
  geom_density(aes(x = earlier, y = ..density..), fill = colours[2], colour = colours[2], alpha = 0.8) +
  geom_density(aes(x = early, y = -..density..), fill = 'white', colour = colours[1], linetype= 'dashed') +
  geom_density(aes(x = late, y = -..density..), fill = colours[4], colour = colours[4], alpha = 0.8) +
  geom_density(aes(x = later, y = -..density..), fill = colours[3], colour = colours[3], alpha = 0.8) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(expand = c(0.05,0,0.05,0), breaks = c(-20, 0, 20), labels = c(20, 0, 20)) +
  xlab('z-transformed mean(mutLoad)') + 
  ggtitle("Bootstrapped mean mutation load") +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()





#--------- SuppFigure 7 ---------#
#calculate mutation load per subtype
# SK-BR3 = HER2+ #
samples <- clinical_data$sample_name[which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'positive')]
mutLoad_HER2pos  <- cnMut_cnTotal.mutationLoad.bin(mutTable[mutTable$patient %in% samples,])

# MDA453 = TNBC #
samples <- clinical_data$sample_name[which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'negative')]
mutLoad_TNBC  <- cnMut_cnTotal.mutationLoad.bin(mutTable[mutTable$patient %in% samples,])

# HER2 neg #
#--> MCF-7 = ER+  and T47D = ER+ and PR+
samples <- clinical_data$sample_name[unique(c(which(clinical_data$final.ER == 'positive' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'negative'),
                                              which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'positive' & clinical_data$final.HER2 == 'negative'),
                                              which(clinical_data$final.ER == 'positive' & clinical_data$final.PR == 'positive' & clinical_data$final.HER2 == 'negative')))]
mutLoad_HER2neg       <- cnMut_cnTotal.mutationLoad.bin(mutTable[mutTable$patient %in% samples,])

mutLoad_subtype_list <- list('SK-BR3' = mutLoad_HER2pos, 'MDA453' = mutLoad_TNBC, 'MCF-7' = mutLoad_HER2neg, 'T47D' = mutLoad_HER2neg)


#sample 1000 bins from each category 10000 times
cellLines <- c('HMEC', 'MCF-7', 'MDA453', 'SK-BR3', 'T47D')
sample_meanMutLoad_switch_subtype <- lapply(cellLines[-1], function(j){
  
  set.seed(123)
  
  #classify not_altered regions in early or late
  x <- ARTregions_list[[j]]
  x$timing <- x$ARTclass
  x$timing[x$timing == 'not_altered' & x$normal_l2r < 0] <- 'late'
  x$timing[x$timing == 'not_altered' & x$normal_l2r > 0] <- 'early'
  
  #add mutation load
  data <- x %>% 
    left_join(mutLoad_subtype_list[[j]], by = c('chr', 'start', 'stop' = 'end')) %>%
    dplyr::filter(!is.na(mutLoad)) %>%
    mutate(mutLoad = (mutLoad - mean(mutLoad)) / sd(mutLoad))
  
  #use minumum number of bins for sampling
  nSample <- min(table(data$timing))
  
  #iterations 
  meanMutLoad_iterations <- sapply(unique(data$timing), function(t){
    mutLoad <- data$mutLoad[data$timing == t]
    iterations <- sapply(1:10000, function(i){
      mean(sample(mutLoad, nSample, replace = T))
    })
    df <- data.frame(iterations)
    return(df)
  })
  meanMutLoad_iterations <- Reduce(cbind, meanMutLoad_iterations)
  meanMutLoad_iterations <- data.frame(meanMutLoad_iterations)
  colnames(meanMutLoad_iterations) <- unique(data$timing)
  meanMutLoad_iterations$cellLine <- j
  
  return(meanMutLoad_iterations)
  
})
sample_meanMutLoad_switch_subtype <- Reduce(rbind, sample_meanMutLoad_switch_subtype)


#mirrored plot
colours <- brewer.pal(n = 11, 'PiYG')[c(2,3,9,10)]

pdf(paste0(output_dir, 'mirrorDist_meanMutload_subtypes.pdf'), width = 5, height = 5)
ggplot(sample_meanMutLoad_switch_subtype) + 
  geom_density(aes(x = early, y = ..density..), fill = colours[1], colour = colours[1], alpha = 0.8) +
  geom_density(aes(x = late, y = ..density..), fill = 'white', colour = colours[4], linetype = 'dashed') +
  geom_density(aes(x = earlier, y = ..density..), fill = colours[2], colour = colours[2], alpha = 0.8) +
  geom_density(aes(x = early, y = -..density..), fill = 'white', colour = colours[1], linetype= 'dashed') +
  geom_density(aes(x = late, y = -..density..), fill = colours[4], colour = colours[4], alpha = 0.8) +
  geom_density(aes(x = later, y = -..density..), fill = colours[3], colour = colours[3], alpha = 0.8) +
  geom_hline(yintercept = 0) +
  facet_wrap(cellLine ~ .) + 
  scale_y_continuous(expand = c(0.05,0,0.05,0), breaks = c(-20, 0, 20), labels = c(20, 0, 20)) +
  xlab('z-transformed mean(mutLoad)') + 
  ggtitle("Bootstrapped mean mutation load (subtype)") +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()


#number of patients per subtye 
HER2pos_samples <- clinical_data$sample_name[which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'positive')]
HER2neg_samples <- clinical_data$sample_name[unique(c(which(clinical_data$final.ER == 'positive' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'negative'),
                                                      which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'positive' & clinical_data$final.HER2 == 'negative'),
                                                      which(clinical_data$final.ER == 'positive' & clinical_data$final.PR == 'positive' & clinical_data$final.HER2 == 'negative')))]
TNBC_samples <- clinical_data$sample_name[which(clinical_data$final.ER == 'negative' & clinical_data$final.PR == 'negative' & clinical_data$final.HER2 == 'negative')]


nPatients_subtypes <- data.frame(subtype = c('HER2+', 'HER2-', 'TNBC'),
                                 n = c(length(HER2pos_samples), length(HER2neg_samples), length(TNBC_samples)))
nPatients_subtypes <- nPatients_subtypes[order(nPatients_subtypes$n, decreasing = T),] 
nPatients_subtypes$subtype <- factor(nPatients_subtypes$subtype, levels = nPatients_subtypes$subtype)

pdf(paste0(output_dir, 'bar_nPatients_subtypes.pdf'), width = 3, height = 3)
ggplot(nPatients_subtypes, aes(x = subtype, y = n, fill = subtype)) +
  geom_bar(stat = 'identity', colour = 'black') + 
  scale_fill_manual(values = brewer.pal(n = 3, name = 'Dark2')) +
  scale_y_continuous(expand = c(0,0,0.01,0)) +
  xlab('') + ylab('# patients') +
  theme_bw() + theme(legend.position = 'none', panel.grid = element_blank())
dev.off()


#number of mutations per patient
nMuts_perPatient <- lapply(unique(mutTable$patient), function(x) {
  data.frame(patient = x, nMuts = sum(mutTable$patient == x))
})
nMuts_perPatient <- Reduce(rbind, nMuts_perPatient)

nMuts_perPatient$subtype <- 'unknown'
nMuts_perPatient$subtype[nMuts_perPatient$patient %in% HER2pos_samples] <- 'HER2+'
nMuts_perPatient$subtype[nMuts_perPatient$patient %in% HER2neg_samples] <- 'HER2-'
nMuts_perPatient$subtype[nMuts_perPatient$patient %in% TNBC_samples] <- 'TNBC'
nMuts_perPatient <- nMuts_perPatient %>% filter(subtype != 'unknown') 
nMuts_perPatient$subtype <- factor(nMuts_perPatient$subtype, levels = levels(nPatients_subtypes$subtype))

pdf(paste0(output_dir, 'box_nMuts_perPatient_subtypes.pdf'), width = 3, height = 3, useDingbats = F)
ggplot(nMuts_perPatient, aes(x = subtype, y = nMuts)) + 
  geom_quasirandom(aes(colour = subtype), alpha = 0.7) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  scale_colour_manual(values = brewer.pal(n = 3, name = 'Dark2')) +
  scale_y_continuous(trans = 'log10') +
  xlab('') + ylab('# mutations per patient') +
  theme_bw() + theme(legend.position = 'none', panel.grid = element_blank())
dev.off()




