################################################################################################################
#############               Association chromatin compartments and Replication Timing              ############# 
################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and Chris Barrington
# and run in R version 3.5.1

# Description:
# Script to create Figure 4 and SuppFigure 7 F  of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution in breast and lung cancer".

#libraries and options
options(stringsAsFactors = F)
options(useDingbats = F)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(data.table)
library(RColorBrewer)
library(grid)
library(dplyr)
library(ggalluvial)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(ggnewscale)
library(ggridges)
library(readr)

#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

#load BRCA mutation data
load(paste0(data_dir, '/560Breast_subset_mutTable.RData'))

#subset for clonal mutations
mutTable   <- mutTable[!is.na(mutTable$absolute.ccf.0.95),]
mutTable   <- mutTable[mutTable$absolute.ccf.0.95 >= 1,]

#load Hi-C data
HiC_list <- readRDS(paste0(data_dir, '/HiC_BRCA.rds'))
HMEC_insitu <- HiC_list[["HMEC_insitu"]]
T47D_insitu <- HiC_list[["T47D_insitu"]]
HMEC_intact <- HiC_list[["HMEC_intact"]]
MCF7_intact <- HiC_list[["MCF7_intact"]]

#load shared altered replication timing
overlap_repTiming <- read.table(paste0(data_dir, '/BRCA_sharedARTregions.txt'), header = T)

#load log2-ratios
log2ratio_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))

#load cell-line ART 
ART_list <- readRDS(paste0(data_dir, '/ARTregions.rds'))
MCF7_ART <- ART_list[["MCF-7"]]
T47D_ART <- ART_list[["T47D"]]

#load BRCA intersect between Hi-C and ART
BRCA_intersect_CC_ART <- read_tsv(paste0(data_dir, 'overlap_HiC_ART_BRCA.tsv'))
LUAD_intersect_CC_ART <- read_tsv(paste0(data_dir, 'overlap_HiC_ART_LUAD.tsv'))



####################################
#######      Functions       #######
####################################

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



###############################
#######      Main       #######
###############################

#--------- Figure 4 A ---------#
plot_data_normal <- BRCA_intersect_CC_ART %>%
  dplyr::filter(hic_cell_type == "Normal") %>%
  dplyr::select(chrom, start, end, normal_cancer, timing_class, hic_cell_type, compartment) %>%
  tidyr::unite(col=id, chrom,start,end)

plot_data_cancer <- BRCA_intersect_CC_ART %>%
  filter(hic_cell_type == "Tumour") %>%
  dplyr::select(chrom, start, end, normal_cancer, timing_class, hic_cell_type, compartment) %>%
  tidyr::unite(col=id, chrom,start,end)
  
plot_data <- plot_data_normal %>%
  left_join(plot_data_cancer, by = c("normal_cancer", "id", "timing_class")) %>%
  dplyr::select(id, normal_cancer, timing_class, compartment.x, compartment.y) %>%
  setNames(c("id", "normal_cancer", "timing_class", "compartment_normal", "compartment_cancer")) %>%
  mutate(normal_cancer = sub("HMEC-", "", normal_cancer)) %>%
  group_by(normal_cancer, timing_class, compartment_normal, compartment_cancer) %>%
  summarise(count = n())

plot_data$normal_cancer      <- factor(plot_data$normal_cancer, levels = c("MCF-7", "T47D"))
plot_data$timing_class       <- factor(plot_data$timing_class, levels = c("earlier", "early", "late", "later"))
plot_data$compartment_normal <- factor(plot_data$compartment_normal, levels = c("A", "B"))
plot_data$compartment_cancer <- factor(plot_data$compartment_cancer, levels = c("A", "B"))

pdf(paste0(output_dir, 'alluvial_CC_in_ART_perCellLine_BRCA.pdf'), width = 8, height = 6)
ggplot(plot_data, aes(y = count, axis1 = compartment_normal, axis2 = compartment_cancer)) +
  geom_alluvium(aes(fill = compartment_normal), width = 1/4) +
  geom_stratum(width = 1/4, aes(fill = compartment_normal), color = "black") +
  scale_x_discrete(limits = c("compartment_normal", "compartment_cancer"), expand = c(.2, .2), labels = c("normal", "tumour")) +
  scale_fill_manual(name = '', values = c('A' = '#e08214', 'B' = '#35978f')) +
  facet_wrap(normal_cancer~timing_class, scales='free_y', nrow = 2) +
  guides(colour='none', fill=guide_legend(title.position='top')) +
  theme_bw() +
  theme(legend.position='bottom',
        axis.title.x=element_blank(),
        panel.border=element_rect(colour='black', fill=NA))
dev.off()


#--------- Figure 4 B-D ---------#

#identify ACC regions
ACC_MCF7 <- MCF7_intact %>%
  left_join(HMEC_intact, by = c("chr", "start", "stop")) %>%
  setNames(c("chr", "start", "stop", "tumour_score", "normal_score")) %>%
  filter(!is.na(normal_score)) %>%
  filter(normal_score != 0 & tumour_score != 0) %>%
  mutate(diff_score = normal_score - tumour_score,
         tumour_CC = ifelse(tumour_score > 0, 'A', 'B'),
         normal_CC = ifelse(normal_score > 0, 'A', 'B'),
         tumour_normal = "MCF7_HMEC")

ACC_T47D <- T47D_insitu %>%
  left_join(HMEC_insitu, by = c("chr", "start", "stop")) %>%
  setNames(c("chr", "start", "stop", "tumour_score", "normal_score")) %>%
  filter(!is.na(normal_score)) %>%
  filter(normal_score != 0 & tumour_score != 0) %>%
  mutate(diff_score = normal_score - tumour_score,
         tumour_CC = ifelse(tumour_score > 0, 'A', 'B'),
         normal_CC = ifelse(normal_score > 0, 'A', 'B'),
         tumour_normal = "T47D_HMEC")

#density of score differences
rbind(ACC_MCF7, ACC_T47D) %>%
  ggplot(. , aes(x = diff_score)) +
  geom_density(fill = "gray") +
  geom_vline(xintercept = c(-0.03, 0.03), linetype = 'dashed') +
  facet_grid(tumour_normal ~ .) +
  ggtitle("Hi-C score normal - cancer") +
  theme_bw()

#classify regions with abs(diff) > 0.03 as ACC
ACC_MCF7 <- ACC_MCF7 %>%
  mutate(ACC = normal_CC,
         ACC = ifelse(normal_CC == "A" & tumour_CC == "B", "A>B", ACC),
         ACC = ifelse(normal_CC == "B" & tumour_CC == "A", "B>A", ACC)) %>%
  filter((ACC %in% c("A>B", "B>A") & abs(diff_score) > 0.03) | (ACC %in% c("A", "B") & abs(diff_score) < 0.03))

ACC_T47D <- ACC_T47D %>%
  mutate(ACC = normal_CC,
         ACC = ifelse(normal_CC == "A" & tumour_CC == "B", "A>B", ACC),
         ACC = ifelse(normal_CC == "B" & tumour_CC == "A", "B>A", ACC)) %>%
  filter((ACC %in% c("A>B", "B>A") & abs(diff_score) > 0.03) | (ACC %in% c("A", "B") & abs(diff_score) < 0.03))


#plot bar plot with fraction of ACC
plot_data <- rbind(ACC_MCF7, ACC_T47D) %>%
  group_by(tumour_normal) %>%
  count(ACC) %>%
  mutate(total = sum(n),
         precent = (n / total) * 100) %>%
  filter(ACC %in% c("A>B", "B>A")) %>%
  mutate(cell_line = sub("_HMEC", "", tumour_normal))

plot_data$ACC <- factor(plot_data$ACC, levels = c("B>A", "A>B"))

pdf(paste0(output_dir, 'bar_ACC_perCellLine_BRCA.pdf'), width = 4, height = 4)
ggplot(plot_data, aes(x = cell_line, y = precent, fill = ACC)) +
  geom_bar(stat = 'identity', colour = 'black') +
  scale_fill_manual(name = "ACC", values = c("B>A" = "#fdb863", "A>B" = "#80cdc1")) +
  geom_text(aes(label = paste(round(precent, 1), "%")), position = position_stack(vjust = 0.5), colour = "white") +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  ggtitle("BRCA") +
  xlab("") +
  ylab("% genome with ACC") +
  theme_bw() +
  theme(panel.grid = element_blank())
dev.off()

#replication timing in ACC regions
ACC_MCF7_RT <- ACC_MCF7 %>%
  left_join(log2ratio_df[, c("chr", "start", "stop", "HMEC", "MCF-7")], by = c("chr", "start", "stop")) %>%
  filter(!is.na(HMEC)) %>%
  filter(!is.na(`MCF-7`)) %>% 
  mutate(normal_RT = ifelse(HMEC > 0, 'early', 'late'),
         tumour_RT = ifelse(`MCF-7` > 0, 'early', 'late')) %>%
  dplyr::select(-HMEC, -`MCF-7`)

ACC_T47D_RT <- ACC_T47D %>%
  left_join(log2ratio_df[, c("chr", "start", "stop", "HMEC", "T47D")], by = c("chr", "start", "stop")) %>%
  filter(!is.na(HMEC)) %>%
  filter(!is.na(T47D)) %>% 
  mutate(normal_RT = ifelse(HMEC > 0, 'early', 'late'),
         tumour_RT = ifelse(T47D > 0, 'early', 'late')) %>%
  dplyr::select(-HMEC, -T47D)

plot_data <- rbind(ACC_MCF7_RT, ACC_T47D_RT) %>%
  mutate(cell_line = sub("_HMEC", "", tumour_normal)) %>%
  dplyr::select(cell_line, ACC, normal_RT, tumour_RT) %>%
  group_by(cell_line, ACC, normal_RT, tumour_RT) %>%
  summarise(count = n())

plot_data$cell_line <- factor(plot_data$cell_line, levels = c("MCF7", "T47D"))
plot_data$ACC <- factor(plot_data$ACC, levels = c("B>A", "A", "B", "A>B"))
plot_data$normal_RT <- factor(plot_data$normal_RT, levels = c("early", "late"))
plot_data$tumour_RT <- factor(plot_data$tumour_RT, levels = c("early", "late"))

pdf(paste0(output_dir, 'alluvial_RT_in_ACC_perCellLine_BRCA.pdf'), width = 8, height = 6)
ggplot(plot_data, aes(y = count, axis1 = normal_RT, axis2 = tumour_RT)) +
  geom_alluvium(aes(fill = normal_RT), width = 1/4) +
  geom_stratum(width = 1/4, aes(fill = normal_RT), color = "black") +
  scale_x_discrete(limits = c("normal_RT", "tumour_RT"), expand = c(.2, .2), labels = c("normal", "tumour")) +
  scale_fill_manual(name = '', values = c('early' = '#c51b7d', 'late' = '#4d9221')) +
  facet_wrap(cell_line~ACC, scales='free_y', nrow = 2) +
  guides(colour='none', fill=guide_legend(title.position='top')) +
  theme_bw() +
  theme(legend.position='bottom',
        axis.title.x=element_blank(),
        panel.border=element_rect(colour='black', fill=NA))
dev.off()

#calculate mutLoad
mutLoad <- cnMut_cnTotal.mutationLoad.bin(mutTable, lspan = NA)

# mutation distributions in ACC regions
meanMutLoad_MCF7_df <- meanMutload_distribution_sampling(mutLoad, 
                                                         ACC_MCF7 %>% 
                                                           dplyr::select(chr, start, stop, ACC) %>% 
                                                           setNames(c("chr", "start", "stop", "timing")), 
                                                         zTrans = T)
meanMutLoad_T47D_df <- meanMutload_distribution_sampling(mutLoad, 
                                                         ACC_T47D %>% 
                                                           dplyr::select(chr, start, stop, ACC) %>% 
                                                           setNames(c("chr", "start", "stop", "timing")), 
                                                         zTrans = T)

meanMutLoad_df <- rbind(meanMutLoad_MCF7_df %>% mutate(cellLine = 'MCF-7',
                                                       type = "chromatin compartment"),
                        meanMutLoad_T47D_df %>% mutate(cellLine = 'T47D',
                                                       type = "chromatin compartment"))


# mutation distributions in ART regions
ART_MCF7 <- MCF7_ART %>%
  mutate(timing = ifelse(normal_l2r > 0, "early", "late"),
         timing = ifelse(ARTclass != "not_altered", ARTclass, timing)) %>%
  dplyr::select(chr, start, stop, timing)

meanMutLoad_MCF7_ART_df <- meanMutload_distribution_sampling(mutLoad, ART_MCF7, zTrans = T)

ART_T47D <- T47D_ART %>%
  mutate(timing = ifelse(normal_l2r > 0, "early", "late"),
         timing = ifelse(ARTclass != "not_altered", ARTclass, timing)) %>%
  dplyr::select(chr, start, stop, timing)

meanMutLoad_T47D_ART_df <- meanMutload_distribution_sampling(mutLoad, ART_T47D, zTrans = T)

meanMutLoad_ART_df <- rbind(meanMutLoad_MCF7_ART_df %>% mutate(cellLine = 'MCF-7',
                                                               type = 'replication timing'),
                            meanMutLoad_T47D_ART_df %>% mutate(cellLine = 'T47D',
                                                               type = "replication timing"))


#plot mirrored plot
colours <- c("#e08214",  "#fdb863", "#80cdc1", "#35978f")
p_chromatin <- ggplot(meanMutLoad_df) + 
  geom_density(aes(x = A, y = ..density..), fill = colours[1], colour = colours[1], alpha = 0.8) +
  geom_density(aes(x = B, y = ..density..), fill = 'white', colour = colours[4], linetype = 'dashed') +
  geom_density(aes(x = `B>A`, y = ..density..), fill = colours[2], colour = colours[2], alpha = 0.8) +
  geom_density(aes(x = A, y = -..density..), fill = 'white', colour = colours[1], linetype= 'dashed') +
  geom_density(aes(x = B, y = -..density..), fill = colours[4], colour = colours[4], alpha = 0.8) +
  geom_density(aes(x = `A>B`, y = -..density..), fill = colours[3], colour = colours[3], alpha = 0.8) +
  geom_hline(yintercept = 0) +
  facet_grid(type ~ cellLine) + 
  scale_y_continuous(expand = c(0.05,0,0.05,0), breaks = c(-20, 0, 20), labels = c(20, 0, 20)) +
  xlab('z-transformed mean(mutLoad)') + 
  ggtitle("BRCA") +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank())

colours <- brewer.pal(n = 11, 'PiYG')[c(2,3,9,10)]
p_repTiming <- ggplot(meanMutLoad_ART_df) + 
  geom_density(aes(x = early, y = ..density..), fill = colours[1], colour = colours[1], alpha = 0.8) +
  geom_density(aes(x = late, y = ..density..), fill = 'white', colour = colours[4], linetype = 'dashed') +
  geom_density(aes(x = earlier, y = ..density..), fill = colours[2], colour = colours[2], alpha = 0.8) +
  geom_density(aes(x = early, y = -..density..), fill = 'white', colour = colours[1], linetype= 'dashed') +
  geom_density(aes(x = late, y = -..density..), fill = colours[4], colour = colours[4], alpha = 0.8) +
  geom_density(aes(x = later, y = -..density..), fill = colours[3], colour = colours[3], alpha = 0.8) +
  geom_hline(yintercept = 0) +
  facet_grid(type ~ cellLine) + 
  scale_y_continuous(expand = c(0.05,0,0.05,0), breaks = c(-20, 0, 20), labels = c(20, 0, 20)) +
  xlab('z-transformed mean(mutLoad)') + 
  theme_bw() + theme(panel.grid = element_blank())

pdf(paste0(mutDist_output_dir, 'mirrorDist_meanMutload_ACC_ART_BRCA.pdf'), width = 5, height = 5)
cowplot::plot_grid(p_chromatin, p_repTiming, ncol = 1)
dev.off()


#--------- Figure 4 E and SuppFigure 7 F ---------#
# univariate models #
# MCF-F (CC and ACC)
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(ACC_MCF7, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ACC))

normal_CC_MCF7 <- summary(lm(mutLoad ~ normal_score, data = input_data))$r.squared
cancer_CC_MCF7 <- summary(lm(mutLoad ~ tumour_score, data = input_data))$r.squared

# T47D (CC and ACC)
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(ACC_T47D, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ACC))

normal_CC_T47D <- summary(lm(mutLoad ~ normal_score, data = input_data))$r.squared
cancer_CC_T47D <- summary(lm(mutLoad ~ tumour_score, data = input_data))$r.squared

# MCF-F (RT and ART)
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(MCF7_ART, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ARTclass))

normal_RT_MCF7 <- summary(lm(mutLoad ~ normal_l2r, data = input_data))$r.squared
cancer_RT_MCF7 <- summary(lm(mutLoad ~ cancer_l2r, data = input_data))$r.squared

# T47D (RT and ART)
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(T47D_ART, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ARTclass))

normal_RT_T47D <- summary(lm(mutLoad ~ normal_l2r, data = input_data))$r.squared
cancer_RT_T47D <- summary(lm(mutLoad ~ cancer_l2r, data = input_data))$r.squared

#plot univariat models
plot_data <- rbind(data.frame(cell_line = "MCF-7", 
                              tissue = c("normal", "cancer"),
                              class = "CC", 
                              rsquare = c(normal_CC_MCF7, cancer_CC_MCF7)),
                   data.frame(cell_line = "T47D", 
                              tissue = c("normal", "cancer"),
                              class = "CC", 
                              rsquare = c(normal_CC_T47D, cancer_CC_T47D)),
                   data.frame(cell_line = "MCF-7", 
                              tissue = c("normal", "cancer"),
                              class = "RT", 
                              rsquare = c(normal_RT_MCF7, cancer_RT_MCF7)),
                   data.frame(cell_line = "T47D", 
                              tissue = c("normal", "cancer"),
                              class = "RT", 
                              rsquare = c(normal_RT_T47D, cancer_RT_T47D)))

pdf(paste0(mutDist_output_dir, 'bar_varianceExplained_univariateModels_BRCA.pdf'), width = 5, height = 4)
ggplot(plot_data %>% filter(class %in% c("CC", "RT")), 
                      aes(x = class, y = rsquare, fill = tissue)) +
  geom_bar(stat = 'identity', position = "dodge", colour = 'black') +
  scale_fill_manual(name = "", values = c("cancer" = "#525252", "normal" = "#bdbdbd")) +
  facet_grid(.~cell_line) +
  xlab('') + ylab('variance in mutation load explained by') +
  scale_y_continuous(expand = c(0,0, 0.05, 0)) +
  ggtitle("Unaltered and altered regions") +
  theme_bw()
dev.off()


# multivariate models #
#MCF-7
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(ACC_MCF7, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ACC)) %>%
  left_join(MCF7_ART, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ARTclass)) %>%
  mutate(tumour_score = (tumour_score - mean(tumour_score)) / sd(tumour_score),
         normal_score = (normal_score - mean(normal_score)) / sd(normal_score),
         cancer_l2r = (cancer_l2r - mean(cancer_l2r)) / sd(cancer_l2r),
         normal_l2r = (normal_l2r - mean(normal_l2r)) / sd(normal_l2r))

total_MCF7  <- summary(lm(mutLoad ~ tumour_score + cancer_l2r + normal_score + normal_l2r, data = input_data))

#T47D
input_data <- mutLoad %>%
  filter(!is.na(mutLoad)) %>%
  left_join(ACC_T47D, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ACC)) %>%
  left_join(T47D_ART, by = c("chr", "start", "end" = "stop")) %>%
  filter(!is.na(ARTclass)) %>%
  mutate(tumour_score = (tumour_score - mean(tumour_score)) / sd(tumour_score),
         normal_score = (normal_score - mean(normal_score)) / sd(normal_score),
         cancer_l2r = (cancer_l2r - mean(cancer_l2r)) / sd(cancer_l2r),
         normal_l2r = (normal_l2r - mean(normal_l2r)) / sd(normal_l2r))

total_T47D  <- summary(lm(mutLoad ~ tumour_score + cancer_l2r + normal_score + normal_l2r, data = input_data))

#plot
plot_data <- rbind(data.frame(coef = rownames(total_MCF7$coefficients),
                              total_MCF7$coefficients,
                              cell_line = "MCF7"),
                   data.frame(coef = rownames(total_T47D$coefficients),
                              total_T47D$coefficients,
                              cell_line = "T47D"))
colnames(plot_data) <- c("coef", "estimate", "stderr", "t", "p", "cell_line")
rownames(plot_data) <- NULL
plot_data           <- plot_data[plot_data$coef != "(Intercept)",]
plot_data           <- plot_data %>%
  mutate(p_signif = ifelse(p < 0.05, "*", "n.s."),
         p_signif = ifelse(p < 0.01, "**", p_signif),
         p_signif = ifelse(p < 0.001, "***", p_signif))

pdf(paste0(mutDist_output_dir, 'bar_coef_multivariateReg_meanMutload_BRCA.pdf'), width = 5, height = 4)
ggplot(plot_data, aes(x = coef, y = abs(estimate))) +
  geom_errorbar(aes(x=coef, ymin=abs(estimate), ymax=abs(estimate)+stderr), width = 0.5) +
  geom_bar(stat = "identity", colour = 'black', fill = "darkgray") +
  facet_grid(.~cell_line) +
  scale_y_continuous(expand = c(0,0, 0.05, 0)) +
  scale_x_discrete(label = c("cancer_l2r" = "cancer RT", "normal_l2r" = "normal RT",
                             "normal_score" = "normal CC", "tumour_score" = "cancer CC")) +
  geom_text(aes(label = p_signif), position = position_stack(vjust = 0.5)) +
  xlab("") +
  ylab("abs(estimate) + stderr") +
  ggtitle("Multivariate linear regression") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#----------- Figure 4 F-G -----------#
plot_data <- LUAD_intersect_CC_ART %>%
  group_by(timing_class) %>%
  count(compartment) %>%
  mutate(percent = (n / sum(n)) * 100)

pdf(paste0(mutDist_output_dir, 'bar_compartments_ART_LUAD.pdf'), width = 5, height = 4)
ggplot(plot_data, aes(x = timing_class, y = n, fill = compartment)) +
  geom_bar(stat = "identity", position = "dodge", colour = 'black') +
  scale_fill_manual(name = 'Compartment', values = c('A' = '#e08214', 'B' = '#35978f')) +
  ylab('Number of 50kb bins') +
  theme_bw() +
  theme(legend.position='bottom',
        axis.title.x=element_blank(),
        panel.border=element_rect(colour='black', fill=NA))

ggplot(plot_data, aes(x = timing_class, y = percent, fill = compartment)) +
  geom_bar(stat = "identity", colour = 'black') +
  scale_fill_manual(name = 'Compartment', values = c('A' = '#e08214', 'B' = '#35978f')) +
  ylab('Proportion (%)') +
  theme_bw() +
  theme(legend.position='bottom',
        axis.title.x=element_blank(),
        panel.border=element_rect(colour='black', fill=NA))
dev.off()


