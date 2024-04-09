############################################################################################################################################## 
#############                        Simulations for timing of ART relative to mutation accumulation in MRCA                     ############# 
############################################################################################################################################## 
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create SuppFigure 8  of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution in breast and lung cancer".
# --> BRCA data used as an example (LUAD Data from the 100,000 Genomes Project are held in a secure research environment and are available to registered users. 
#     See https://www.genomicsengland.co.uk/research/academic for further information.)
# Data accessibility statement can be found in the manuscript.

#options and libraries
options(stringsAsFactors = F)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggridges)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

#load BRCA mutation data (file was too big, so had to be split in 3 to be uploaded)
mutTable_1 <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_1.rds'))
mutTable_2 <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_2.rds'))
mutTable_3 <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_3.rds'))
mutTable   <- rbind(mutTable_1, mutTable_2, mutTable_3)
mutTable   <- mutTable[!is.na(mutTable$absolute.ccf),]
mutTable   <- mutTable[mutTable$absolute.ccf > 0.95,]

#load shared altered replication timing
overlap_repTiming <- read.table('/camp/lab/swantonc/working/dietzem/RepliSeq/Data/cohortTables/BRCA_sharedARTregions.txt', header = T)



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


#function to simulate timing of replication timing alterations based on mutation accumulation
simulation_ARTtiming_relativeMuts <- function(nBins = 1000, percLate = 0.66, ARTfraction = 0.15, mutRate = 0.1,
                                              mutRate_factor_late = 2, iterations = 10000, EarlierTiming = 0.1, LaterTiming = 0.1, 
                                              timing = 'punctuated', gap = 100){
  
  set.seed(123)
  
  ### early ###
  #initiate list with replication timing and mutation informations
  nEarly     <- round(nBins * (1-percLate))
  early_list <- lapply(1:nEarly, function(i){ data.frame(repTiming = 'early', mutRate = mutRate, nMuts = 0)})
  
  #create mutations until beginning of ART event
  Later_iteration   <- iterations * LaterTiming
  early_list <- lapply(early_list, function(x){
    x$nMuts  <- sum(sample(x = c(0,1), size = Later_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
    return(x)
  })
  
  #simulate ART event
  nLater     <- round(nBins * ARTfraction / 2)
  ww_later   <- 1:nLater
  
  if(timing == 'punctuated'){
    #update ART entries
    early_list[ww_later]   <- lapply(early_list[ww_later], function(x) {
      x$repTiming <- 'later'
      x$mutRate   <- mutRate * mutRate_factor_late
      return(x)
    })
    #create mutations after ART event
    early_list <- lapply(early_list, function(x){
      x$nMuts <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Later_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
  } else if(timing == 'linear'){
    #mutation accumulation for ART events with gap-times iterations between next ART event
    early_list[ww_later]   <- lapply(1:length(ww_later), function(i) {
      x           <- early_list[[ww_later[i]]]
      x$nMuts     <- x$nMuts + sum(sample(x = c(0,1), size = (i-1)*gap, prob = c(1-x$mutRate, x$mutRate), replace = T))
      x$repTiming <- 'later'
      x$mutRate   <- mutRate * mutRate_factor_late
      x$nMuts     <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Later_iteration - (i-1)*gap, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
    #mutation accumulation for ART events 
    ww_non_later <- c(1:length(early_list))[!c(1:length(early_list)) %in% ww_later]
    early_list[ww_non_later]   <- lapply(early_list[ww_non_later], function(x) {
      x$nMuts <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Later_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
  }
  
  
  ### late ###
  #initiate list with replication timing and mutation informations
  nLate     <- nBins - nEarly
  late_list <- lapply(1:nLate, function(i){ data.frame(repTiming = 'late', mutRate = mutRate * mutRate_factor_late, nMuts = 0)})
  
  #create mutations until ART event
  Earlier_iteration   <- iterations * EarlierTiming
  late_list <- lapply(late_list, function(x){
    x$nMuts <- sum(sample(x = c(0,1), size = Earlier_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
    return(x)
  })
  
  #simulate ART event
  nEarlier   <- round(nBins * ARTfraction / 2)
  ww_earlier   <- 1:nEarlier
  
  if(timing == 'punctuated'){
    #update ART entries
    late_list[ww_earlier]   <- lapply(late_list[ww_earlier], function(x) {
      x$repTiming <- 'earlier'
      x$mutRate   <- mutRate 
      return(x)
    })
    #create mutations after ART event
    late_list <- lapply(late_list, function(x){
      x$nMuts <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Earlier_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
  } else if(timing == 'linear'){
    #mutation accumulation for ART events with gap-times iterations between next ART event
    late_list[ww_earlier]   <- lapply(1:length(ww_earlier), function(i) {
      x           <- late_list[[ww_earlier[i]]]
      x$nMuts     <- x$nMuts + sum(sample(x = c(0,1), size = (i-1)*gap, prob = c(1-x$mutRate, x$mutRate), replace = T))
      x$repTiming <- 'earlier'
      x$mutRate   <- mutRate
      x$nMuts     <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Earlier_iteration - (i-1)*gap, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
    #mutation accumulation for ART events 
    ww_non_earlier <- c(1:length(late_list))[!c(1:length(late_list)) %in% ww_earlier]
    late_list[ww_non_earlier]   <- lapply(late_list[ww_non_earlier], function(x) {
      x$nMuts <- x$nMuts + sum(sample(x = c(0,1), size = iterations - Earlier_iteration, prob = c(1-x$mutRate, x$mutRate), replace = T))
      return(x)
    })
  }
  
  
  #output
  iteration_list <- c(early_list, late_list)
  iteration_df   <- Reduce(rbind, iteration_list)
  return(iteration_df)
  
}



#function to estimate mean mutLoad distribution in repTiming regions
meanMutload_distribution_sampling <- function(simulated_Muts, iter = 10000){
  
  set.seed(123)
  
  #z-tranform mutation load
  data <- simulated_Muts %>% 
    mutate(mutLoad = (nMuts - mean(nMuts)) / sd(nMuts))
  
  #use minumum number of bins for sampling
  nSample <- min(table(data$repTiming))
  
  #iterations 
  meanMutLoad_iterations <- sapply(unique(data$repTiming), function(t){
    mutLoad <- data$mutLoad[data$repTiming == t]
    iterations <- sapply(1:iter, function(i){
      mean(sample(mutLoad, nSample, replace = T))
    })
    df <- data.frame(iterations)
    return(df)
  })
  meanMutLoad_iterations <- Reduce(cbind, meanMutLoad_iterations)
  meanMutLoad_iterations <- data.frame(meanMutLoad_iterations)
  colnames(meanMutLoad_iterations) <- unique(data$repTiming)
  
  return(meanMutLoad_iterations)
  
}



#function to test if observed and simulated distributions overlap
test_oberserved_simulated_mutDist <- function(observed_data, simulated_data){
  
  #calculate difference between early and earlier and late and later
  diff_observed  <- data.frame(earlier = abs(observed_data$early - observed_data$earlier) / abs(observed_data$early - observed_data$late),
                               later = abs(observed_data$late - observed_data$later) / abs(observed_data$early - observed_data$late))
  diff_simulated <- data.frame(earlier = abs(simulated_data$early - simulated_data$earlier) / abs(simulated_data$early - simulated_data$late),
                               later = abs(simulated_data$late - simulated_data$later) / abs(simulated_data$early - simulated_data$late))
  
  #calculate 95%-CI for simulated data (ground truth)
  lowCI  <- apply(diff_simulated, 2, function(i) quantile(i, 0.025))
  highCI <- apply(diff_simulated, 2, function(i) quantile(i, 0.975))
  
  #count how many observed values don't fall within the CI and divide by number of iterations
  pvalue_earlier <- sum(diff_observed$earlier < lowCI['earlier'] | diff_observed$earlier > highCI['earlier']) / nrow(diff_observed)
  pvalue_later   <- sum(diff_observed$later < lowCI['later'] | diff_observed$later > highCI['later']) / nrow(diff_observed)
  
  #output
  df <- rbind(data.frame(repTiming = 'earlier', lowCI = lowCI['earlier'], highCI = highCI['earlier'], 
                         count = sum(diff_observed$earlier < lowCI['earlier'] | diff_observed$earlier > highCI['earlier']),
                         niter = nrow(diff_observed), pvalue = pvalue_earlier),
              data.frame(repTiming = 'later', lowCI = lowCI['later'], highCI = highCI['later'], 
                         count = sum(diff_observed$later < lowCI['later'] | diff_observed$later > highCI['later']),
                         niter = nrow(diff_observed), pvalue = pvalue_later))
  return(df)
}

############################
######       Main     ###### 
############################

#---------- observed data ----------#
#calculate mutation load
mutLoad       <- cnMut_cnTotal.mutationLoad.bin(mutTable)

#assign replication timing
mutLoad <- mutLoad %>% left_join(overlap_repTiming, by = c('chr' = 'chr', 'start' = 'start', 'end' = 'stop'))
mutLoad <- mutLoad[!is.na(mutLoad$timing) & !is.na(mutLoad$mutLoad),]

#calculate mean mutation load distributions in different replication timing regions
mutLoad$nMuts     <- mutLoad$mutLoad
mutLoad$repTiming <- mutLoad$timing
observed_meanMutLoad_df <- meanMutload_distribution_sampling(mutLoad)

#calculate differences between early and earlier and late and later
diff_observed_meanMutLoad_df         <- observed_meanMutLoad_df
diff_observed_meanMutLoad_df$earlier <- abs(diff_observed_meanMutLoad_df$early - diff_observed_meanMutLoad_df$earlier)
diff_observed_meanMutLoad_df$later   <- abs(diff_observed_meanMutLoad_df$late - diff_observed_meanMutLoad_df$later)
diff_observed_meanMutLoad_df         <- diff_observed_meanMutLoad_df[,c('earlier', 'later')]


#---------- mutation rate in consistent RT regions ----------#
#assign replication timing
repTiming_df <- overlap_repTiming
mutLoad <- mutLoad %>% left_join(repTiming_df, by = c('chr' = 'chr', 'start' = 'start', 'end' = 'stop'))
mutLoad <- mutLoad[!is.na(mutLoad$timing) & !is.na(mutLoad$mutLoad),]

#calculate average mutation load per Mb for each timing
mutLoad_perMB_perRepTiming <- mutLoad %>%
  group_by(timing) %>%
  summarise(mutCount = sum(mutLoad),
            regionSize_perMb = sum(width) / 1000000) %>%
  mutate(mutLoad_perMb = mutCount / regionSize_perMb)

#plot pie of fraction of mutations in early and late replicated regions adjusted for size
#--> SuppFigure 8A
plot_data <- data.frame(timing = c('early', 'late'),
                        fraction = c(mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'early'] / (mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'late'] + mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'early']),
                                     mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'late'] / (mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'late'] + mutLoad_perMB_perRepTiming$mutLoad_perMb[mutLoad_perMB_perRepTiming$timing == 'early'])))

plot_data <- plot_data %>%
  mutate(prop = fraction * 100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop ) %>%
  mutate(timing = factor(timing, levels = rev(c('early', 'late'))))

pdf(paste0(output_dir, 'pie_frac_consistentRT_mutLoad_BRCA.pdf'), width =4, height = 4)
ggplot(plot_data, aes(x="", y=prop, fill=timing)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(x = 1, y = ypos, label = paste0(round(prop), '%')), color = 'white', size = 5) +
  scale_fill_manual(name = 'Replication Timing', values = c('early' = '#c51b7d', 'late' = '#4d9221')) +
  ggtitle('Proportion of mutations corrected for region size') +
  theme_void() + theme(plot.title = element_text(hjust = 0.5), legend.position = 'top')
dev.off()



#---------- simulated data ----------#
#those parameters changed for LUAD and LUSC --> see Methods of the manuscript
mutRate              <- 0.1
mutRate_factor_late  <- 1.3
ARTfraction          <- 0.15

#simulations
simulated_meanMutLoad_list <- lapply(seq(0, 1, 0.01), function(s){
  meanMutload_distribution_sampling(simulation_ARTtiming_relativeMuts(EarlierTiming = s, LaterTiming = s))
})
names(simulated_meanMutLoad_list) <- paste0('ARTtiming_',seq(0, 1, 0.01))

#plot simulated distributions
#--> SuppFigure 8C
plot_data <- lapply(paste0('ARTtiming_',seq(0, 1, 0.1)), function(x){
  meanMutLoad <- simulated_meanMutLoad_list[[x]]
  meanMutLoad <- reshape2::melt(meanMutLoad)
  meanMutLoad$ARTtiming <- sub('ARTtiming_', '', x)
  return(meanMutLoad)
})
plot_data            <- Reduce(rbind, plot_data)
colnames(plot_data)  <- c('repTiming', 'meanMutLoad', 'ARTtiming')
plot_data$repTiming  <- factor(plot_data$repTiming, levels = c('early', 'earlier', 'later', 'late'))

pdf(paste0(output_dir, 'simulations_ART_mutDist_BRCA.pdf'), width = 6, height = 10)
ggplot(plot_data, aes(x = meanMutLoad, fill = repTiming)) +
  geom_histogram(position = 'identity', alpha = 1, bins = 500) +
  scale_fill_manual(name = 'Replication Timing', values = brewer.pal(n = 11, 'PiYG')[c(1,3,9,11)]) +
  facet_grid(ARTtiming ~ .) +
  xlab('zTrans(mean(mutLoad))') +
  ggtitle(paste0("Different frations of mutations accumulated before ART")) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#plot ridge of differences
#--> SuppFigure 8C
plot_data <- lapply(paste0('ARTtiming_',seq(0, 1, 0.1)), function(x){
  meanMutLoad <- simulated_meanMutLoad_list[[x]]
  meanMutLoad$earlier <- abs(meanMutLoad$earlier - meanMutLoad$early) / abs(meanMutLoad$late - meanMutLoad$early)
  meanMutLoad$later   <- abs(meanMutLoad$later - meanMutLoad$late) / abs(meanMutLoad$late - meanMutLoad$early)
  meanMutLoad <- reshape2::melt(meanMutLoad[,c('earlier', 'later')])
  meanMutLoad$ARTtiming <- sub('ARTtiming_', '', x)
  return(meanMutLoad)
})
plot_data            <- Reduce(rbind, plot_data)
colnames(plot_data)  <- c('repTiming', 'meanMutLoad', 'ARTtiming')
plot_data$repTiming  <- factor(plot_data$repTiming, levels = c('earlier', 'later'))
plot_data$ARTtiming  <- factor(plot_data$ARTtiming, levels = rev(seq(0, 1, 0.1)))

pdf(paste0(output_dir, 'simulations_ART_mutDist_ridge_diffMeanMtLoad_BRCA.pdf'), width = 6, height = 5)
ggplot(plot_data, aes(x = meanMutLoad, y = ARTtiming, fill = repTiming)) + 
  geom_density_ridges(alpha = 0.8) +
  scale_fill_manual(name = 'RepTiming Change', values = c('later' = '#7fbc41', 'earlier' = '#de77ae'),
                    labels = c('later' = 'Early-to-Late minus Late', 'earlier' = 'Late-to-Early minus Early')) +
  facet_grid(.~repTiming) +
  xlab('Difference in bootstrapped mean mutaion load') + ylab('fractions of iterations before ART') +
  theme_bw() + theme(legend.position = 'none')
dev.off()


#test which simulated distribution is significantly the same as the observed data
test_similarity_distribution <- lapply(names(simulated_meanMutLoad_list), function(x){
  data <- test_oberserved_simulated_mutDist(observed_meanMutLoad_df, simulated_meanMutLoad_list[[x]])
  data$ARTtiming <- sub('ARTtiming_', '', x)
  return(data)
})
test_similarity_distribution <- Reduce(rbind, test_similarity_distribution)
rownames(test_similarity_distribution) <- NULL


#plot fraction overlaps
#--> SuppFigure 8F
plot_data <- test_similarity_distribution
plot_data$frac_verlap <- (1 - plot_data$pvalue)

pdf(paste0(output_dir, 'ART_punctuating_homogenous_overlapDistributions.pdf'), width = 8, height = 4)
ggplot(plot_data, aes(x = factor(ARTtiming), y = frac_verlap, fill = repTiming)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  geom_hline(yintercept = 0.95, linetype = 'dashed') +
  scale_fill_manual(name = 'ART timing', values = c('later' = '#7fbc41', 'earlier' = '#de77ae')) +
  scale_x_discrete(breaks = seq(0, 1, 0.1)) +
  xlab("Fraction of Mutations accumulated before ART") +
  ylab("Fraction overlap between distributions") +
  ggtitle('Overlap between observed abd simulated distribution of differences \nof mean mutation load between altered and non-altered regions') +
  theme_bw() + theme(legend.position = 'top', plot.title = element_text(hjust = 0.5))
dev.off()

