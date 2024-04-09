###############################################################################################################################
#############                   Mutational signatures and APOBEC omikli events in RT and ART regions              ############# 
###############################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create Figure 5 and SuppFigure 9-10 of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution in breast and lung cancer".
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
library(Biostrings)
library(cowplot)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved

#load BRCA mutation data (file was too big, so had to be split in 3 to be uploaded)
mutTable_1  <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_1.rds'))
mutTable_2  <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_2.rds'))
mutTable_3  <- readRDS(paste0(data_dir, '/560Breast_subset_mutTable_3.rds'))
mutTable    <- rbind(mutTable_1, mutTable_2, mutTable_3)
mutTable_gr <- GRanges(seqnames = mutTable$Chrom, IRanges(start = mutTable$Pos, end = mutTable$Pos))

#load shared altered replication timing
overlap_repTiming <- read.table(paste0(data_dir, '/BRCA_sharedARTregions.txt'), header = T)

#SBS signatures
exposures <- readRDS(paste0(data_dir, 'exposures_SBS_BRCA.rds'))
signatures <- readRDS(paste0(data_dir, 'signatures_SBS_BRCA.rds'))

#hyperClust results
hyperClust_table <- readRDS(paste0(data_dir, 'hyperClust_table_BRCA.rds'))
hyperClust_table <- hyperClust_table[hyperClust_table$sample %in% unique(mutTable$patient),]

#load cancer genes
pancan_genes <- read.table(paste0(data_dir, '20220203_pan.driver_og.tsg.csv'), sep = ',', header = T)




##################################
#####        Functions        #####   
##################################

#function to count trinucleotide context motifs in regions of interest relative to whole genome ###
trinuc_count_fun <- function(regionInterest_df, 
                             genome_name = "BSgenome.Hsapiens.UCSC.hg19", 
                             chromosomes = paste0('chr', 1:22)){
  
  #trinucleotide patterns
  sigs <- deconstructSigs::signatures.genome.cosmic.v3.may2019
  trinuc_pattern <- data.frame(trinuc = colnames(sigs),
                               first_base = substr(colnames(sigs), start = 1, stop = 1),
                               second_base = substr(colnames(sigs), start = 3, stop = 3),
                               third_base = substr(colnames(sigs), start = 7, stop = 7))
  trinuc_pattern$trinuc_pattern <- paste0(trinuc_pattern$first_base, trinuc_pattern$second_base, trinuc_pattern$third_base)
  trinuc_context <- DNAStringSet(unique(unique(trinuc_pattern$trinuc_pattern)))
  reverseComp_trinuc_context <- reverseComplement(trinuc_context)
  
  #count patterns
  context_counts        <- rep(0, length(unique(trinuc_pattern$trinuc_pattern)))
  names(context_counts) <- unique(trinuc_pattern$trinuc_pattern)
  for(chr in chromosomes){
    print(chr)
    chr_genome   <- getBSgenome(genome_name)[[chr]]
    chr_region   <- regionInterest_df[regionInterest_df[,1] == chr,,drop = F]
    if(nrow(chr_region) == 0){ next }
    seq_segments <- Views(chr_genome, start = chr_region[,2], end = chr_region[,3])
    counts_trinuc   <- sapply(data.frame(trinuc_context, stringsAsFactors = F)[,1], function(x) vcountPattern(x, seq_segments))
    counts_reverse  <- sapply(data.frame(reverseComp_trinuc_context, stringsAsFactors = F)[,1], function(x) vcountPattern(x, seq_segments))
    
    if(class(counts_trinuc) == 'matrix'){
      counts <- colSums(counts_trinuc) + colSums(counts_reverse)
    } else {
      counts <- counts_trinuc + counts_reverse
    }
    context_counts <- context_counts + counts[match(names(context_counts), names(counts))]
  }
  
  return(context_counts)
  
}


#function to plot SBS profiles
plot_SBS <- function(signatures_df, title = ''){
  
  #add groups
  colnames(signatures_df) <- c('channel', 'value')
  signatures_df$group <- substr(as.character(signatures_df$channel), start = 3, stop = 5)
  signatures_df <- signatures_df %>%
    mutate(channel = factor(channel, levels = signatures_df$channel),
           group = factor(group, levels = unique(signatures_df$group)),
           value = value * 100)
  
  #set colours
  colours   <- setNames(c('#1ebff0', '#050708', '#e62725', '#cbcacb', '#a1cf64', '#edc8c5'), unique(signatures_df$group))
  strip_name_colours <- c('black','white','white','black','black','black')
  xlabels   <- paste0(substr(as.character(signatures_df$channel), start = 1, stop = 1),
                      substr(as.character(signatures_df$channel), start = 3, stop = 3),
                      substr(as.character(signatures_df$channel), start = 7, stop = 7))
  
  #plot
  p <- ggplot(signatures_df, aes(x = channel, y = value, fill = group)) +
    geom_bar(stat = 'identity') + 
    facet_grid(.~ group, space = 'free_x', scales = 'free_x') +
    scale_fill_manual(name = '', values = colours, guide = 'none') +
    xlab('') +
    ylab('% SBS') +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    scale_x_discrete(labels = xlabels) +
    ggtitle(title) +
    theme_classic() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                            strip.text = element_text(face = 'bold'), plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
  
  #change colours and labels of strips
  g       <- ggplot_gtable(ggplot_build(p))
  striprt <- which(grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name))
  k <- 1
  for (i in striprt) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colours[k]
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col  <- colours[k]
    
    t <- which(grepl('text', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[t]]$children[[1]]$gp$col <- strip_name_colours[k]
    
    k <- k+1
  }
  
  return(g)
}


#function to calssify APOBEC mutations
call.APOBECmut = function(mutTable, flanking.size = 5, flanking.needed = FALSE, genome = BSgenome.Hsapiens.UCSC.hg19){
  
  #add flanking sequence if needed
  if(flanking.needed){
    context  <- Biostrings::getSeq(genome, name = mutTable$Chrom, start = (mutTable$Pos-flanking.size), end = (mutTable$Pos+flanking.size))
    flanking <- paste(toupper(substr(as.character(context), 1, flanking.size)), 
                      tolower(substr(as.character(context), (flanking.size + 1), (flanking.size + 1))), 
                      toupper(substr(as.character(context), flanking.size + 2, flanking.size + 2 + flanking.size)), sep = "")
    mutTable$flanking <- flanking
  } 
  mutTable_flanking <- mutTable
  mutTable_flanking$APOBECmut <- 'notAPOBEC'
  
  #identify c > T or G mutations
  index.C <- which(mutTable_flanking$Ref=='C'&mutTable_flanking$Alt%in%c('T','G')| mutTable_flanking$Ref=='G'&mutTable_flanking$Alt%in%c('A','C'))
  
  #APOBEC mutations fullfill stringent TCW motif (W either A or T)
  index.TpCpW         <- which(substr(mutTable_flanking$flanking, flanking.size,flanking.size)=='T'
                               &substr(mutTable_flanking$flanking, flanking.size+2,flanking.size+2)%in%c('A','T')
                               &mutTable_flanking$Ref=='C'
                               | substr(mutTable_flanking$flanking, flanking.size,flanking.size)%in%c('A','T')
                               &substr(mutTable_flanking$flanking, flanking.size+2,flanking.size+2)=='A'
                               &mutTable_flanking$Ref=='G')
  
  #not sure if APOBEC or not when W is C or G, so considered as none of both groups 
  index.unknown      <- which(substr(mutTable_flanking$flanking, flanking.size,flanking.size)=='T'
                              &substr(mutTable_flanking$flanking, flanking.size+2,flanking.size+2)%in%c('C','G')
                              &mutTable_flanking$Ref=='C'
                              | substr(mutTable_flanking$flanking, flanking.size,flanking.size)%in%c('C','G')
                              &substr(mutTable_flanking$flanking, flanking.size+2,flanking.size+2)=='A'
                              &mutTable_flanking$Ref=='G')
  
  mutTable_flanking$APOBECmut[intersect(index.C, index.TpCpW)]   <- 'APOBEC'
  mutTable_flanking$APOBECmut[intersect(index.C, index.unknown)] <- 'Unknown'
  
  #no flanking information --> NA
  mutTable_flanking$APOBECmut[is.na(mutTable_flanking$flanking)] <- NA
  
  return(mutTable_flanking)
}



####################################
#####         Main             #####   
####################################

#--------------- create input files for HDP_sigExtraction pipeline ---------------#
#assign repTiming to mutations
overlap_repTiming_gr     <- makeGRangesFromDataFrame(overlap_repTiming_df)
overlap                  <- findOverlaps(overlap_repTiming_gr, mutTable_gr)
mutation_df              <- mutTable[subjectHits(overlap),]
mutation_df$timing       <- overlap_repTiming_df$timing[queryHits(overlap)]

#calculate 96-trinuc counts per timing and adjust for background
input_matrix <- c()
for(time in c('early', 'earlier', 'later', 'late')){
  print(time)
  sub_mutations <- mutation_df[mutation_df$timing %in% time,]
  input         <- deconstructSigs::mut.to.sigs.input(sub_mutations, sample.id = 'patient', chr = 'Chrom', pos = 'Pos', ref = 'Ref', alt = 'Alt')
  
  trinuc_genome  <- readRDS(paste0(output_dir, 'trinucContext_genome.rds'))
  trinuc_region  <- trinuc_count_fun(regionInterest_df = overlap_repTiming_df[overlap_repTiming_df$timing %in% time,])
  trinuc_genome$region_count <- trinuc_region[match(trinuc_genome$trinuc_context, names(trinuc_region))]
  trinuc_genome$ratio        <- trinuc_genome$genome_count / trinuc_genome$region_count
  norm_ratio                 <- trinuc_genome[,'ratio',drop = F]
  rownames(norm_ratio)       <- trinuc_genome$trinuc_context
  
  norm.it <- function(col, trimer.ratio){
    trimer  <- paste(substr(colnames(col), 1, 1), substr(colnames(col), 3, 3), substr(colnames(col), 7, 7), sep = "")
    new.col <- col*trimer.ratio[trimer,]
    return(new.col)
  }
  
  norm_input           <- sapply(colnames(input), function(x) {norm.it(input[,x,drop=F], trimer.ratio = norm_ratio)})
  norm_input           <- data.frame(norm_input, row.names = paste(rownames(input), time, sep = '_'))
  colnames(norm_input) <- colnames(input)
  norm_input           <- norm_input / rowSums(norm_input)
  norm_input           <- round(norm_input * rowSums(input))
  
  input_matrix <- rbind(input_matrix, norm_input)
}

#save matrix
saveRDS(input_matrix, file = paste0(output_dir, 'input_96matrix_patient_sharedRepTiming.rds'))

#treeLayer
treeLayer_df <- data.frame(sample = rownames(input_matrix),
                           repTiming = matrix(unlist(strsplit(rownames(input_matrix), '_')), ncol = 2, byrow = T)[,2])
treeLayer_df <- treeLayer_df[order(treeLayer_df$repTiming),]
rownames(treeLayer_df) <- NULL
saveRDS(treeLayer_df, file = paste0(output_dir, 'treeLayers_patient_sharedRepTiming.rds'))

####### priors file ########
# genome breast signatures
sigs           <- read.table(paste0(data_dir, 'COSMIC_v3.2_SBS_GRCh37.txt'), header = T)
rownames(sigs) <- sigs$Type
sigs           <- sigs[,-1]
prior_sigs     <- as.matrix(sigs[,paste0('SBS', c(1,2,3,5,6,8,13))])
deconstructSigs_sigs <- t(deconstructSigs::signatures.genome.cosmic.v3.may2019)
prior_sigs           <- prior_sigs[match(rownames(deconstructSigs_sigs), rownames(prior_sigs)),]

saveRDS(prior_sigs, file = paste0(output_dir, 'priorSBS_breast_genome.rds'))

#--> run HDP_sigExtraction (https://github.com/McGranahanLab/HDP_sigExtraction) pipeline to extract mutation signatures in RT and ART regions


#---------------- Figure 5 A ----------------#
#--> Figure 5 B was created in the same way using the LUAD data from the 100,000 Genomes Project
sigActivity_threshold <- 0.1
exposures[exposures < sigActivity_threshold] <- 0
order_signatures <- colnames(exposures)[order(as.numeric(gsub('SBS|b', '', colnames(exposures))))]
exposures        <- exposures[,order_signatures]

#find paired exposures in different RT and ART regions
plot_weights <- exposures
plot_weights$sample    <- rownames(plot_weights)
plot_weights$patient   <- matrix(unlist(strsplit(rownames(plot_weights), '_')), ncol = 2, byrow = T)[,1]
plot_weights$repTiming <- matrix(unlist(strsplit(rownames(plot_weights), '_')), ncol = 2, byrow = T)[,2]
plot_weights <- reshape2::melt(plot_weights, id.vars = c('patient', 'repTiming', 'sample'))
plot_weights$value[plot_weights$value == 0] <- NA
plot_weights$repTiming <- factor(plot_weights$repTiming, levels = c('early', 'earlier', 'later', 'late'))
plot_weights$variable  <- factor(plot_weights$variable, levels = order_signatures)

overlap_patients <- intersect(intersect(plot_weights$patient[plot_weights$repTiming == 'early'], plot_weights$patient[plot_weights$repTiming == 'earlier']),
                              intersect(plot_weights$patient[plot_weights$repTiming == 'late'], plot_weights$patient[plot_weights$repTiming == 'later']))
overlap_patients <- unique(overlap_patients)

paired_exposures <- plot_weights[plot_weights$sample %in% paste0(overlap_patients, '_early'),] %>%
  dplyr::select(patient, variable, value)
for(t in c('earlier', 'later', 'late')){
  data <- plot_weights[plot_weights$sample %in% paste0(overlap_patients, '_', t),] %>%
    dplyr::select(patient, variable, value)
  paired_exposures <- paired_exposures %>%
    left_join(data, by = c('patient', 'variable'))
}
colnames(paired_exposures)[grep('value', colnames(paired_exposures))] <- c('early', 'earlier', 'later', 'late')

ww <- apply(paired_exposures[,c('early', 'earlier', 'later', 'late')], 1, function(x) all(is.na(x)))
paired_exposures <- paired_exposures[!ww,]

plot_data <- paired_exposures %>%
  mutate(early = ifelse(is.na(early), 0, early),
         earlier = ifelse(is.na(earlier), 0, earlier),
         later = ifelse(is.na(later), 0, later),
         late = ifelse(is.na(late), 0, late),
         early_late = ifelse(late - early == 0, NA, late - early),
         earlier_later = ifelse(later - earlier == 0, NA, later - earlier)) %>%
  group_by(variable) %>%
  summarise(early_late = median(early_late, na.rm = T),
            earlier_later = median(earlier_later, na.rm = T))

ww <- apply(paired_exposures[,c('early', 'earlier', 'later', 'late')], 1, function(x) all(is.na(x)))
plot_data_size <- paired_exposures[!ww,] %>%
  group_by(variable) %>%
  summarise(active = n() / length(overlap_patients) * 100)

plot_data <- plot_data %>% left_join(plot_data_size, by = 'variable')

sig_exclude <- as.character(plot_data_size$variable[plot_data_size$active < 1])
if(length(sig_exclude) > 0){
  plot_data <- plot_data[!plot_data$variable %in% sig_exclude,]
}

pdf(paste0(output_dir, 'scatter_medianExposure_ART.pdf'), width = 4, height = 3.5, useDingbats = F)
max_value <- max(abs(c(plot_data$early_late, plot_data$earlier_later)))
ggplot(plot_data, aes(x = early_late, y = earlier_later)) + 
  geom_rect(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = '#7fbc41', alpha = 0.01) +
  geom_rect(xmin = 0, xmax = -Inf, ymin = 0, ymax = -Inf, fill = '#de77ae', alpha = 0.01) +
  geom_point(aes(size = active, fill = early_late), shape = 21) + 
  scale_fill_gradient2(low = '#c51b7d', mid = '#b2abd2', high = '#4d9221', midpoint = 0, guide = 'none') +
  scale_size_continuous(name = '% patients active', breaks = seq(0,100,1), guide = 'none') +
  ggrepel::geom_text_repel(aes(label = variable))+
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(c(-1*max_value, max_value)) + 
  ylim(c(-1*max_value, max_value)) +
  xlab('median exposure difference RT') +
  ylab('median exposure difference ART') +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()



#---------------- SuppFigure 9 ----------------#
pdf(paste0(output_dir, "SBSprofiles_BRCA.pdf"), width = 10, height = 5)
lapply(colnames(signatures), function(x){
  plot_data <- data.frame(channel = rownames(signatures), value = signatures[,x])
  grid.draw(plot_SBS(plot_data, title = x))
  grid.newpage()
})
dev.off()



#---------------- SuppFigure 10 A-B ----------------#
input_matrix <- readRDS(paste0(data_dir, 'input_96matrix_patient_sharedRepTiming.rds'))

plot_weights <- exposures
plot_weights$unknown   <- 1-rowSums(plot_weights)
plot_weights$sample    <- rownames(plot_weights)
plot_weights$patient   <- matrix(unlist(strsplit(rownames(plot_weights), '_')), ncol = 2, byrow = T)[,1]
plot_weights$repTiming <- matrix(unlist(strsplit(rownames(plot_weights), '_')), ncol = 2, byrow = T)[,2]
plot_weights <- reshape2::melt(plot_weights, id.vars = c('patient', 'repTiming', 'sample'))
plot_weights$type <- 'weights'

plot_counts <- data.frame(sample = rownames(input_matrix),
                          patient = matrix(unlist(strsplit(rownames(input_matrix), '_')), ncol = 2, byrow = T)[,1],
                          repTiming = matrix(unlist(strsplit(rownames(input_matrix), '_')), ncol = 2, byrow = T)[,2],
                          variable = 'counts',
                          value = rowSums(input_matrix),
                          type = 'counts')
plot_counts <- plot_counts[rownames(input_matrix) %in% rownames(exposures),]
plot_counts <- plot_counts %>%
  group_by(repTiming) %>%
  arrange(value, .by_group = T) %>%
  data.frame()

plot_data <- rbind(plot_weights, plot_counts)
plot_data$sample   <- factor(plot_data$sample, levels = plot_counts$sample)
plot_data$variable <- factor(plot_data$variable, levels = c(order_signatures,'unknown','counts'))

pdf(paste0(output_dir, 'bar_exposures.pdf'), width = 10, height = 4, useDingbats = F)
ggplot(plot_data, aes(x = sample, y = value, fill = variable)) + 
  geom_bar(stat = 'identity', position = position_stack(reverse = T)) + 
  scale_fill_manual(name = 'Signatures', values = c(brewer.pal(n = 12, 'Paired'), '#d8b365', '#e0e0e0', '#878787')) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(type ~ repTiming, scales = 'free', space = 'free_x') + 
  xlab('') + ylab('') +
  theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
dev.off()




#---------------- Figure 5 C ----------------#
#--> run hyperClust (https://github.com/davidmasp/hyperclust) and use results here

#annotate mutation cluster
mutTable <- mutTable %>% left_join(hyperClust_table, by = c('patient' = 'sample', 'Chrom' = 'seqnames', 'Pos' = 'start'))
mutTable$event_type[is.na(mutTable$event_type)] <- 'unclustered'

#annotate A3B mutations
mutTable <- call.APOBECmut(mutTable, flanking.needed = T)

#annotate timing
mutTable_gr             <- GRanges(seqnames = mutTable$Chrom, IRanges(start = mutTable$Pos, end = mutTable$Pos))
consistent_repTiming_gr <- makeGRangesFromDataFrame(overlap_repTiming)
overlap                 <- findOverlaps(mutTable_gr, consistent_repTiming_gr)
mutTable$repTiming      <- NA
mutTable$repTiming[queryHits(overlap)] <- overlap_repTiming$timing[subjectHits(overlap)]

#exclude mutations with no timing information
mutTable_all <- mutTable
mutTable     <- mutTable[!is.na(mutTable$repTiming),]

repTiming_size <- overlap_repTiming %>%
  group_by(timing) %>%
  mutate(width = stop - start) %>%
  summarise(size = sum(width) / 1000000)

sub_mutTable <- mutTable[mutTable$Ref %in% c('C', 'G'),]
sub_mutTable$mutType <- 'other'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'G') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'C')] <- 'C>G'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'T') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'A')] <- 'C>T'

clusterMuts_frequencies <- sub_mutTable %>%
  dplyr::filter(mutType != 'other') %>%
  group_by(event_type, APOBECmut, repTiming) %>%
  summarise(n = n())

clusterMuts_frequencies <- clusterMuts_frequencies %>%
  left_join(repTiming_size, by = c('repTiming' = 'timing'))

clusterMuts_frequencies$frequency <- clusterMuts_frequencies$n / clusterMuts_frequencies$size

omikli <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('omikli')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC')) %>%
  mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

unclustered <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('unclustered')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC')) %>%
  mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

plot_data <- rbind(omikli[c('event_type', 'repTiming', 'frequency')],
                   unclustered[c('event_type', 'repTiming', 'frequency')])
plot_data$repTiming <- factor(plot_data$repTiming, levels = c('early', 'earlier', 'later', 'late'))

pdf(paste0(output_dir, 'bar_APOBEComikli_muts_perMB_BRCA.pdf'), width = 5, height = 5)
ggplot(plot_data, aes(x = repTiming, y = frequency, fill = repTiming)) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'white') + 
  scale_fill_manual(name = 'RepTiming', values = c('early' = '#c51b7d', 'earlier' = '#de77ae', 'later' = '#7fbc41', 'late' = '#4d9221')) +
  facet_grid(event_type~., scales = 'free_y') +
  xlab('') + ylab('APOBEC mutations per Mb') +
  ggtitle('BRCA') +
  scale_y_continuous(expand = c(0,0,0.01, 0)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



#---------------- Figure 5 D ----------------#
mutTable$cancerGene <- ifelse(mutTable$Gene %in% pancan_genes$Gene, 'cancerGene', 'non_cancerGene')

repTiming_size <- overlap_repTiming %>%
  group_by(timing) %>%
  mutate(width = stop - start) %>%
  summarise(size = sum(width) / 1000000)

sub_mutTable <- mutTable[mutTable$Ref %in% c('C', 'G'),]
sub_mutTable$mutType <- 'other'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'G') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'C')] <- 'C>G'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'T') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'A')] <- 'C>T'

clusterMuts_frequencies <- sub_mutTable %>%
  dplyr::filter(!CDS %in% '-') %>%
  dplyr::filter(cancerGene %in% 'cancerGene') %>%
  dplyr::filter(mutType != 'other') %>%
  dplyr::filter(event_type != 'kataegis') %>%
  group_by(event_type, APOBECmut, repTiming) %>%
  summarise(n = n())

clusterMuts_frequencies <- clusterMuts_frequencies %>%
  left_join(repTiming_size, by = c('repTiming' = 'timing'))

clusterMuts_frequencies$frequency <- clusterMuts_frequencies$n / clusterMuts_frequencies$size

#plot omikli
plot_data <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('omikli')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC')) %>%
  mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

p_omikli <- ggplot(plot_data, aes(x = repTiming, y = n, fill = repTiming)) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'white') + 
  scale_fill_manual(name = 'RepTiming', values = c('early' = '#c51b7d', 'earlier' = '#de77ae', 'later' = '#7fbc41', 'late' = '#4d9221')) +
  scale_y_continuous(expand = c(0,0,0.1,0)) +
  xlab('') + ylab('# A3 omikli mutations \nin cancer genes') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA))

#plot unclustered
plot_data <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('unclustered')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC')) %>%
  mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

p_unclustered <- ggplot(plot_data, aes(x = repTiming, y = n, fill = repTiming)) +
  geom_bar(stat = 'identity', position = position_dodge(), colour = 'white') + 
  scale_fill_manual(name = 'RepTiming', values = c('early' = '#c51b7d', 'earlier' = '#de77ae', 'later' = '#7fbc41', 'late' = '#4d9221')) +
  scale_y_reverse(expand = c(0.1,0,0,0)) +
  scale_x_discrete(position = "top") +
  xlab('') + ylab('A3 unclustered mutations\n in cancer genes') +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = 'none', axis.text.x = element_text(angle = 90, vjust = -0.5, hjust = 0.5, size = 12),
                     panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA))


#cancer gene enrichment
count <- mutTable %>%
  dplyr::filter(!CDS %in% '-') %>%
  dplyr::filter(cancerGene %in% 'cancerGene') %>%
  dplyr::filter(!event_type %in% 'kataegis') %>%
  dplyr::filter(!is.na(repTiming)) %>%
  group_by(APOBECmut, event_type, repTiming) %>%
  summarise(count = n()) %>%
  tidyr::spread(event_type, count, fill = 0)

fisher_results <- lapply(unique(count$repTiming), function(x){
  m <- matrix(c(count$omikli[count$repTiming == x & count$APOBECmut == 'APOBEC'],
                count$unclustered[count$repTiming == x & count$APOBECmut == 'APOBEC'],
                count$omikli[count$repTiming == x & count$APOBECmut == 'notAPOBEC'],
                count$unclustered[count$repTiming == x & count$APOBECmut == 'notAPOBEC']), 
              ncol = 2, byrow = F, dimnames = list(c('omikli', 'unclustered'), c('APOBEC', 'nonAPOBEC')))
  test <- fisher.test(m)
  data.frame(repTiming = x, p_value = test$p.value, oddRatio = test$estimate, lowCI = test$conf.int[1], highCI = test$conf.int[2])
})
fisher_results <- Reduce(rbind, fisher_results)
rownames(fisher_results) <- NULL
fisher_results$repTiming <- factor(fisher_results$repTiming, levels = c('early', 'earlier', 'later', 'late'))

p_fisher <- ggplot(fisher_results) + 
  geom_pointrange(aes(x = repTiming, y = oddRatio, ymin = lowCI, ymax = highCI, colour = repTiming)) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  scale_colour_manual(name = 'RepTiming', values = c('early' = '#c51b7d', 'earlier' = '#de77ae', 'later' = '#7fbc41', 'late' = '#4d9221')) +
  scale_y_continuous(trans = 'log10') +
  xlab('') + ylab('Odds-Ratio ') + 
  theme_bw() + theme(axis.text.x = element_blank(), legend.position = 'none', panel.background = element_rect(fill = "transparent",colour = NA),
                     plot.background = element_rect(fill = "transparent",colour = NA))

#combine plots
p_omikli      <- ggplotGrob(p_omikli + theme(plot.margin = unit(c(0, 0.5, 0, 0.1), "cm")))
p_fisher      <- ggplotGrob(p_fisher + theme(plot.margin = unit(c(-0.2, 0.5, -0.5, 0.1), "cm")))
p_unclustered <- ggplotGrob(p_unclustered + theme(plot.margin = unit(c(0, 0.5, 0, 0.1), "cm")))

plot_list <- list(p_omikli, p_fisher, p_unclustered)
all_widths <- lapply(plot_list, function(x) {x$widths})
plot_list_alignedWidths <- lapply(plot_list, function(x){
  x$widths <- do.call(unit.pmax, all_widths)
  return(x)
})

#plot
full_list <- rev(plot_list_alignedWidths)
heights    <- c(0.3,0.2, 0.3)
ypos      <- c(0, cumsum(heights)) + 0.05
full_plot <- ggdraw()
for(x in 1:length(full_list)){
  full_plot <- full_plot + draw_plot(full_list[[x]], x = 0, y = ypos[x], width = 1, height = heights[x])
}
full_plot <- full_plot + draw_plot_label(c('BRCA'), x = 0, y = 0.9, size = 20)

pdf(paste0(output_dir, 'combined_oddratioCancerGene_A3clusterPerMb_inCancerGenes.pdf'), useDingbats = F, width = 3, height = 7)
plot(full_plot)
dev.off()


### cancer genes that harbour APOBEC omikli mutations and are earlier replicated ###
genes <- mutTable %>%
  dplyr::filter(!CDS %in% '-') %>%
  dplyr::filter(cancerGene %in% 'cancerGene') %>%
  dplyr::filter(event_type %in% 'omikli') %>%
  dplyr::filter(APOBECmut %in% 'APOBEC') %>%
  dplyr::filter(!is.na(repTiming))

patients_perGene <- genes %>%
  group_by(Gene) %>%
  summarise(nPatient = length(unique(Sample))) %>%
  arrange(desc(nPatient)) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)))



#---------------- SuppFigure 10 C-D ----------------#
repTiming_size <- overlap_repTiming %>%
  group_by(timing) %>%
  mutate(width = stop - start) %>%
  summarise(size = sum(width) / 1000000)

sub_mutTable <- mutTable[mutTable$Ref %in% c('C', 'G'),]
sub_mutTable$mutType <- 'other'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'G') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'C')] <- 'C>G'
sub_mutTable$mutType[(sub_mutTable$Ref %in% 'C' & sub_mutTable$Alt %in% 'T') | (sub_mutTable$Ref %in% 'G' & sub_mutTable$Alt %in% 'A')] <- 'C>T'

clusterMuts_frequencies <- sub_mutTable %>%
  dplyr::filter(mutType != 'other') %>%
  dplyr::group_by(event_type, APOBECmut, repTiming, mutType) %>%
  dplyr::summarise(n = n())

clusterMuts_frequencies <- clusterMuts_frequencies %>%
  left_join(repTiming_size, by = c('repTiming' = 'timing'))

clusterMuts_frequencies$frequency <- clusterMuts_frequencies$n / clusterMuts_frequencies$size

#Omikli
plot_data <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('omikli', 'unclustered')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC', 'notAPOBEC')) %>%
  dplyr::filter(mutType %in% c('C>G', 'C>T')) %>%
  dplyr::mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

plot_data <- plot_data %>%
  dplyr::select(event_type, APOBECmut, mutType, repTiming, frequency) %>%
  tidyr::spread(repTiming, frequency)

plot_data[,c('early', 'earlier', 'later', 'late')] <- plot_data[,c('early', 'earlier', 'later', 'late')] / plot_data$late
plot_data <- reshape2::melt(plot_data, is.vars = c('event_type', 'APOBECmut', 'mutType'))
plot_data$variable <- factor(plot_data$variable, levels = rev(c('early', 'earlier', 'later', 'late')))
plot_data$value    <- log2(plot_data$value)

p_omikli <- ggplot(plot_data) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(aes(x = variable, y = value, colour = event_type, shape = mutType)) +
  geom_line(aes(x = variable, y = value, colour = event_type, group = interaction(mutType, event_type))) +
  scale_colour_manual(name = 'eventType', values = c('omikli' = '#1f78b4', 'unclustered' = '#b2df8a')) +
  facet_grid(.~APOBECmut) +
  xlab('Replication Timing') + ylab('log2(mutations per Mb / late mutations per Mb)') +
  ggtitle('Omikli') +
  theme_bw()


#Kataegis
plot_data <- clusterMuts_frequencies %>%
  dplyr::filter(event_type %in% c('kataegis', 'unclustered')) %>%
  dplyr::filter(APOBECmut %in% c('APOBEC', 'notAPOBEC')) %>%
  dplyr::filter(mutType %in% c('C>G', 'C>T')) %>%
  dplyr::mutate(repTiming = factor(repTiming, levels = c('early', 'earlier', 'later', 'late')))

plot_data <- plot_data %>%
  dplyr::select(event_type, APOBECmut, mutType, repTiming, frequency) %>%
  tidyr::spread(repTiming, frequency)

plot_data[,c('early', 'earlier', 'later', 'late')] <- plot_data[,c('early', 'earlier', 'later', 'late')] / plot_data$late
plot_data <- reshape2::melt(plot_data, is.vars = c('event_type', 'APOBECmut', 'mutType'))
plot_data$variable <- factor(plot_data$variable, levels = rev(c('early', 'earlier', 'later', 'late')))
plot_data$value    <- log2(plot_data$value)

p_kataegis <- ggplot(plot_data) + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_point(aes(x = variable, y = value, colour = event_type, shape = mutType)) +
  geom_line(aes(x = variable, y = value, colour = event_type, group = interaction(mutType, event_type))) +
  scale_colour_manual(name = 'eventType', values = c('kataegis' = '#a6cee3', 'unclustered' = '#b2df8a')) +
  facet_grid(.~APOBECmut) +
  xlab('Replication Timing') + ylab('log2(mutations per Mb / late mutations per Mb)') +
  ggtitle('Kataegis') +
  theme_bw()

#plot
pdf(paste0(output_dir, 'points_relativeMutCounts_clusterEvent_BRCA.pdf'), width = 8, height = 4, useDingbats = F)
p_omikli
p_kataegis
dev.off()










