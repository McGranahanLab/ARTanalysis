#################################################################################################################################
#############               Analysis of altered replication timing (ART) regions in BRCA, LUAD and LUSC             ############# 
#################################################################################################################################
# written by Michelle Dietzen (m.dietzen@ucl.ac.uk) and run in R version 3.5.1

# Description:
# Script to create Figure 2 and SuppFigure 4 - 6 of the manuscript "Replication timing alterations impact mutation acquisition during tumour evolution".
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
library(ComplexUpset)


#parameters
data_dir   <- '.' #set full path to the directory where the data for this analysis has been saved
output_dir <- '.' #set full path to the directory where the results for this analysis should be saved



#load altered replication timing regions
ARTregions_list <- readRDS(paste0(data_dir, '/ARTregions.rds'))

#load ART regions identfied with wrong normals
TT1_artRegions_list   <- readRDS(paste0(data_dir, 'TT1_ARTregions.rds'))
IMR90_artRegions_list <- readRDS(paste0(data_dir, 'IMR90_ARTregions.rds'))

#load RT values in 50kb windows
repTiming_df <- readRDS(paste0(data_dir, '/cohort_50kb_l2r.rds'))




#################################
#######     Functions     #######
#################################

### function to plot altered replication timing regions per chromosome ###
plot_ARTregions_heatmap <- function(repTiming_df, ref, test, chromosomes = paste0('chr', c(1:22))) {
  
  #heatmap
  plot_data            <- repTiming_df[!repTiming_df$ARTclass %in% c('not_altered'),]
  plot_data$ARTclass   <- factor(plot_data$ARTclass, levels = c('earlier', 'later'))
  plot_data            <- plot_data[plot_data$chr %in% chromosomes,]
  plot_data$chr        <- factor(plot_data$chr, levels = chromosomes)
  index_names          <- setNames(seq(0, length(chromosomes) - 1), chromosomes)
  plot_data$ystart     <- index_names[plot_data$chr]
  plot_data$yend       <- index_names[plot_data$chr] + 1
  
  #get start and stop points of chromosomes
  chromosome.length <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)[chromosomes] 
  hline.data        <- data.frame(start = 0, end = chromosome.length, y = seq(1,length(chromosomes)))
  if(length(chromosomes) > 22){
    hline.data[22, 2] <- chromosome.length[23]
  }
  hline.data[21, 2] <- chromosome.length[22]
  hline.data[19, 2] <- chromosome.length[20]
  vline.data        <- data.frame(x = chromosome.length, start = seq(0,(length(chromosomes) - 1)), end = seq(1,length(chromosomes)))
  
  #colour timeChanges differently
  colours        <- brewer.pal(n = 11, 'PiYG')[c(3,9)]
  names(colours) <- c('earlier', 'later')
  labels         <- c('late normal --> early tumour',  'early normal --> late tumour')
  names(labels) <- c('earlier', 'later')
  
  #heatmap plot
  heatmap <- ggplot() +
    geom_rect(data = plot_data, aes(xmin = as.numeric(start) / 1000000, xmax = as.numeric(stop) / 1000000, ymin = ystart, ymax = yend, fill = ARTclass)) + 
    scale_fill_manual(name = 'Altered RepTiming', values = colours, labels = labels) +
    scale_x_continuous(expand = c(0, 0)) + 
    scale_y_continuous(breaks = seq(0.5, (length(chromosomes) - 0.5)), labels = chromosomes, expand = c(0, 0)) + 
    geom_segment(data = hline.data, mapping = aes(x = start, xend = end / 1000000, y = y, yend = y), colour = 'black', size = 0.3) +
    geom_segment(data = vline.data, mapping = aes(x = x / 1000000, xend = x / 1000000, y = start, yend = end), colour = 'black', size = 0.3) +
    xlab('Genomic Position (mb)') + ylab('') +
    theme_classic() +
    theme(strip.text.y = element_blank()) +
    theme(axis.ticks.y = element_blank()) +
    theme(panel.grid = element_blank()) +
    theme(axis.line = element_line(size = 0.3)) +
    theme(legend.text = element_text(size = 8))
  
  
  #barplot with fraction of chromosomes that are different
  count_chr <- repTiming_df %>%
    group_by(chr) %>%
    summarise(total = n())
  
  plot_data_bar <- plot_data %>%
    group_by(chr, ARTclass) %>%
    summarise(count = n()) %>%
    left_join(count_chr) %>%
    mutate(fraction = (count / total) * 100,
           ARTclass = factor(ARTclass, levels = levels(plot_data$ARTclass)),
           chr = factor(chr, levels = chromosomes))
  
  bar <- ggplot(plot_data_bar, aes(x = chr, y = fraction, fill = ARTclass)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(name = 'repTiming Change', values = colours, labels = labels) +
    scale_y_reverse(expand = c(0, 0)) +
    geom_hline(yintercept = 0.001, color = "black", size = 0.2) +
    ylab('% Chromosome') + xlab("") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = 'left') + 
    theme(panel.grid = element_blank()) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(axis.line.x = element_line(size = 0.3), axis.line.y = element_blank(), legend.position = 'none')
  
  #combine both plots
  heatmap <- ggplotGrob(heatmap + theme(plot.margin = unit(c(1, 0.5, 0.5, 0), "cm")))
  bar     <- ggplotGrob(bar + theme(plot.margin = unit(c(1, 0.1, 0.5, 0.5), "cm")))
  
  heatmap$heights <- grid::unit.pmax(heatmap$heights, bar$heights)
  bar$heights     <- grid::unit.pmax(heatmap$heights, bar$heights)
  
  combined.list        <- list(bar, heatmap)
  combined.plot.layout <- matrix(c(rep(1,3), rep(2,12)), nrow = 1)
  combined.plot        <- arrangeGrob(grobs = combined.list, ncol = dim(combined.plot.layout)[2], nrow = 1, layout_matrix = combined.plot.layout, top = paste0('Altered Replication Timing \n', ref, ' - ', test))
  
  return(combined.plot)
}


### function to plot pie of altered replication timing regions ###
plot_ARTregions_pie <- function(repTiming_df, ref, test) {
  
  #count regions
  plot_data       <- repTiming_df[!repTiming_df$ARTclass %in% c('not_altered'),]
  plot_data_count <- plot_data %>%
    group_by(ARTclass) %>%
    summarise(count = n()) %>%
    mutate(total = nrow(repTiming_df),
           fraction = (count / total) * 100)
  plot_data_count <- rbind(plot_data_count, data.frame(ARTclass = 'not_altered', count = nrow(repTiming_df) - sum(plot_data_count$count),
                                                       total = sum(plot_data_count$count), fraction =  (nrow(repTiming_df) - sum(plot_data_count$count)) * 100 / nrow(repTiming_df)))
  plot_data_count$ARTclass <- factor(plot_data_count$ARTclass, levels = c('earlier', 'later','not_altered'))
  plot_data_count <- plot_data_count[match(levels(plot_data_count$ARTclass), plot_data_count$ARTclass),]
  
  plot_data_count <- plot_data_count %>% 
    mutate(cumulative = cumsum(fraction),
           midpoint = cumulative - fraction / 2,
           label = paste0(round(fraction,1), '%'))
  
  
  #set colour for changes
  colours        <- c(brewer.pal(n = 11, 'PiYG')[c(3,9)], '#bababa')
  names(colours) <- c('earlier',  'later', 'not_altered')
  labels         <- c('late normal --> early tumour',  'early normal --> late tumour', 'not altered')
  names(labels)  <- c('earlier', 'later', 'not_altered')
  
  #plot pie
  p <- ggplot(plot_data_count, aes(x = 1, weights = fraction, fill = ARTclass)) +
    geom_bar(width = 1, position = position_stack(reverse = T), colour = 'white') +
    coord_polar(theta = "y") +
    scale_fill_manual(name = 'Altered RepTiming', values = colours, labels = labels) +
    geom_label(aes(x = rev(c(1, 1.2, 1)), y = midpoint, label = label), size = 4) + 
    xlab('') + ylab('') +
    ggtitle(paste0('Altered Replication Timing \n', ref, ' - ', test)) +
    theme_void() + 
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5))
  
  #return plot and fractions
  repTiming_fractions <- plot_data_count[,1:4]
  colnames(repTiming_fractions) <- c('ARTclass', 'nBins', 'totalBins', 'fraction')
  
  return(list(fractions = repTiming_fractions, pie = p))
  
}


#function copied from https://github.com/cran/spatialEco/blob/master/R/impute.loess.R
impute.loess <- function(y, s = 0.2, smooth = FALSE) {
  x.length = length(y)
  if(length(y[!is.na(y)]) < 6) {
    warning("Fewer than 6 real-value observations, assigning NA")
    y <- rep(NA, x.length)
  } else {
    x <- 1:x.length
    p <- suppressWarnings(stats::loess(y ~ x, span = s, 
                                       data.frame(x = x, y = y)))
    if (smooth == TRUE) {
      y <- stats::predict(p, x)
    } else {
      na.idx <- which(is.na(y))
      if(length(na.idx) >= 1) {
        y[na.idx] <- stats::predict(p, data.frame(x = na.idx))
      }
    }
  }
  return(y)
}



#function to check if the observed overlap is significantly bigger than expected by chance
check_random_overlaps <- function(ARTregions_list, niter = 1000){
  
  set.seed(123)
  
  iterations <- lapply(1:niter, function(iter){
    print(iter)
    
    #randomly sample earlier and later regions from normal later / early regions
    random_ARTregions_list <- lapply(ARTregions_list, function(x){
      
      #reference early 
      ARTregions <- x %>%
        filter(ARTclass != 'not_altered') %>%
        filter(normal_l2r > 0)
      ARTregions_gr <- makeGRangesFromDataFrame(ARTregions)
      ARTregions_gr <- GenomicRanges::reduce(ARTregions_gr)
      
      potential_ARTregions_gr <- makeGRangesFromDataFrame(x[x$normal_l2r > 0,])
      potential_ARTregions_combined_gr <- GenomicRanges::reduce(potential_ARTregions_gr)
      
      ramdon_index <- sample(1:length(potential_ARTregions_gr), length(ARTregions_gr), replace = F)
      random_gr <- GRanges(seqnames = seqnames(potential_ARTregions_gr)[ramdon_index], IRanges(start = start(potential_ARTregions_gr)[ramdon_index], width = width(ARTregions_gr)))
      overlap   <- findOverlaps(random_gr, potential_ARTregions_combined_gr, type = 'within')
      index     <- c(1:length(ARTregions_gr))[!1:length(ARTregions_gr) %in% queryHits(overlap)]
      
      random_refEarly_df <- as.data.frame(random_gr[queryHits(overlap)])
      
      while(length(index) != 0){
        potential_ARTregions_gr <- setdiff(potential_ARTregions_gr, random_gr[queryHits(overlap)])
        potential_ARTregions_combined_gr  <- GenomicRanges::reduce(potential_ARTregions_gr)
        ARTregions_gr           <- ARTregions_gr[index]
        
        ramdon_index <- sample(1:length(potential_ARTregions_gr), length(ARTregions_gr), replace = F)
        random_gr <- GRanges(seqnames = seqnames(potential_ARTregions_gr)[ramdon_index], IRanges(start = start(potential_ARTregions_gr)[ramdon_index], width = width(ARTregions_gr)))
        overlap   <- findOverlaps(random_gr, potential_ARTregions_combined_gr, type = 'within')
        index     <- c(1:length(ARTregions_gr))[!1:length(ARTregions_gr) %in% queryHits(overlap)]
        
        random_refEarly_df <- rbind(random_refEarly_df, as.data.frame(random_gr[queryHits(overlap)]))
      }
      
      #reference late 
      ARTregions <- x %>%
        filter(ARTclass != 'not_altered') %>%
        filter(normal_l2r < 0)
      ARTregions_gr <- makeGRangesFromDataFrame(ARTregions)
      ARTregions_gr <- GenomicRanges::reduce(ARTregions_gr)
      
      potential_ARTregions_gr <- makeGRangesFromDataFrame(x[x$normal_l2r < 0,])
      potential_ARTregions_combined_gr <- GenomicRanges::reduce(potential_ARTregions_gr)
      
      ramdon_index <- sample(1:length(potential_ARTregions_gr), length(ARTregions_gr), replace = F)
      random_gr <- GRanges(seqnames = seqnames(potential_ARTregions_gr)[ramdon_index], IRanges(start = start(potential_ARTregions_gr)[ramdon_index], width = width(ARTregions_gr)))
      overlap   <- findOverlaps(random_gr, potential_ARTregions_combined_gr, type = 'within')
      index     <- c(1:length(ARTregions_gr))[!1:length(ARTregions_gr) %in% queryHits(overlap)]
      
      random_refLate_df <- as.data.frame(random_gr[queryHits(overlap)])
      
      while(length(index) != 0){
        potential_ARTregions_gr <- setdiff(potential_ARTregions_gr, random_gr[queryHits(overlap)])
        potential_ARTregions_combined_gr  <- GenomicRanges::reduce(potential_ARTregions_gr)
        ARTregions_gr           <- ARTregions_gr[index]
        
        ramdon_index <- sample(1:length(potential_ARTregions_gr), length(ARTregions_gr), replace = F)
        random_gr <- GRanges(seqnames = seqnames(potential_ARTregions_gr)[ramdon_index], IRanges(start = start(potential_ARTregions_gr)[ramdon_index], width = width(ARTregions_gr)))
        overlap   <- findOverlaps(random_gr, potential_ARTregions_combined_gr, type = 'within')
        index     <- c(1:length(ARTregions_gr))[!1:length(ARTregions_gr) %in% queryHits(overlap)]
        
        random_refLate_df <- rbind(random_refLate_df, as.data.frame(random_gr[queryHits(overlap)]))
      }
      
      random_ARTevents <- rbind(data.frame(random_refEarly_df, ref_timing = 'early'),
                                data.frame(random_refLate_df,  ref_timing = 'late'))
      
      overlap        <- findOverlaps(makeGRangesFromDataFrame(random_ARTevents), makeGRangesFromDataFrame(x), minoverlap = 10)
      random_ARTbins <- x[subjectHits(overlap), 1:3]
      random_ARTbins$ref_timing <- random_ARTevents$ref_timing[queryHits(overlap)]
      
      return(random_ARTbins)
    })
    names(random_ARTregions_list) <- names(ARTregions_list)
    
    #find overlapping events
    overlap_random_ARTregions <- random_ARTregions_list[[1]]
    for(i in 2:length(random_ARTregions_list)){
      ARTregions <- random_ARTregions_list[[i]]
      overlap_random_ARTregions <- overlap_random_ARTregions %>% full_join(ARTregions, by = c('chr','start','stop'))
    }
    
    #number of bins overlapping
    n_cellLines <- length(random_ARTregions_list)
    overlap_random_ARTregions$class <- apply(overlap_random_ARTregions[4:ncol(overlap_random_ARTregions)], 1, function(i){
      if(sum(!is.na(i)) == 1){
        return('unique')
      } else if(sum(!is.na(i)) == n_cellLines){
        return('recurrent')
      } else {
        return('shared')
      }
    })
    
    #output
    output <- data.frame(recurrent = sum(overlap_random_ARTregions$class == 'recurrent'),
                         shared = sum(overlap_random_ARTregions$class == 'shared'),
                         unique = sum(overlap_random_ARTregions$class == 'unique'))
    return(output)
  })
  iterations <- Reduce(rbind, iterations)
  return(iterations)
}




############################
#######     Main     #######
############################

#-------- Figure 2 A and SuppFigure 4 --------#
ARTregions_df <- Reduce(rbind, ARTregions_list)

ARTclass_fractions <-c()
for(x in unique(ARTregions_df$normal_cancer)){
  repTiming_df <- ARTregions_df[ARTregions_df$normal_cancer == x,]
  
  #normal and cancer cellLines
  normal <- unlist(strsplit(repTiming_df$normal_cancer[1], '-'))[1]
  if(length(unlist(strsplit(repTiming_df$normal_cancer[1], '-'))) > 2){
    cancer <- paste(unlist(strsplit(repTiming_df$normal_cancer[1], '-'))[2:length(unlist(strsplit(repTiming_df$normal_cancer[1], '-')))], collapse = '-')
  } else {
    cancer <- unlist(strsplit(repTiming_df$normal_cancer[1], '-'))[2]
  }
  print(cancer)
  
  #plot heatmap
  p_heatmap <- plot_ARTregions_heatmap(repTiming_df, ref = normal, test = cancer, chromosomes = paste0('chr', c(1:22)))
  pdf(paste0(output_dir, normal, '/', cancer, '/ARTregions_heatmap_perChr.pdf'), width = 9, height = 6)
  grid.draw(p_heatmap)
  dev.off()
  
  #plot pie
  pie_fractions <- plot_ARTregions_pie(repTiming_df, ref = normal, test = cancer)
  pdf(paste0(output_dir, normal, '/', cancer, '/pie_fraction_ARTregions_exclNoSwitch.pdf'), width = 7, height = 5)
  print(pie_fractions[[2]])
  dev.off()
  
  #combine fractions for all cell-lines
  fractions_df        <- pie_fractions[[1]]
  fractions_df$normal <- normal
  fractions_df$cancer <- cancer
  ARTclass_fractions  <- rbind(ARTclass_fractions, fractions_df)
}



#-------- Figure 2 B --------#
ARTclass_fractions$cancerType <- 'BRCA'
ARTclass_fractions$cancerType[ARTclass_fractions$normal == 'T2P']   <- 'LUAD'
ARTclass_fractions$cancerType[ARTclass_fractions$normal == 'HBEC3'] <- 'LUSC'
ARTclass_fractions$cancerType <- factor(ARTclass_fractions$cancerType, levels = c('LUAD', 'LUSC', 'BRCA'))

plot_data          <- ARTclass_fractions[ARTclass_fractions$ARTclass != 'not_altered',]
plot_data$ARTclass <- factor(plot_data$ARTclass, levels = c('earlier', 'later'))
plot_data$label    <- paste0(round(plot_data$fraction, 1), '%')

colours        <- brewer.pal(n = 11, 'PiYG')[c(3,9)]
names(colours) <- c('earlier', 'later')
labels         <- c('late normal --> early tumour',  'early normal --> late tumour')
names(labels) <- c('earlier', 'later')

pdf(paste0(output_dir, '/bar_fraction_ARTregions.pdf'), width = 8, height = 4)
ggplot(plot_data, aes(x = cancer, y = fraction, fill = ARTclass)) + 
  geom_bar(stat = 'identity', colour = 'black') + 
  scale_fill_manual(name = 'Altered RepTiming', values = colours, labels = labels) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), colour = 'white', size = 3) +
  facet_grid( ~ cancerType, scale = 'free_x', space = 'free_x') +
  scale_y_continuous(expand = c(0,0,0,0.2)) +
  xlab('') + ylab('% genome altered repTiming') + 
  labs(title = 'Altered Replication Timing Regions',
       subtitle = '(References: LUAD = T2P, LUSC = HBEC3, BRCA = HMEC)') +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
dev.off()


#-------- Figure 2 C --------#

T2P_artRegions   <- ARTregions_list$H1650
TT1_artRegions   <- TT1_artRegions_list$H1650
IMR90_artRegions <- IMR90_artRegions_list$H1650

#classify concordant and discordant timing
combined_df <- T2P_artRegions[c('chr',  'start',   'stop', 'ARTclass')] %>% 
  dplyr::full_join(TT1_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop')) %>%
  dplyr::full_join(IMR90_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop'))

combined_df <- combined_df[,c(1:3, 4:ncol(combined_df))]
colnames(combined_df)[4:6] <- c('T2P', 'TT1', 'IMR90')

ww <- apply(combined_df[4:6], 1, function(x) any(is.na(x)))
combined_df <- combined_df[!ww,]

ww <- apply(combined_df[4:6], 1, function(x) all(x[!is.na(x)] == 'not_altered'))
combined_df <- combined_df[!ww,]

combined_df$timing <- 'discordant'
combined_df$timing <- sapply(1:nrow(combined_df), function(i) ifelse(all(combined_df[i,4:6] %in% 'earlier'), 'late-to-early', combined_df$timing[i])) 
combined_df$timing <- sapply(1:nrow(combined_df), function(i) ifelse(all(combined_df[i,4:6] %in% 'later'), 'early-to-late', combined_df$timing[i])) 

#IMR90
IMR90_unique <- combined_df %>%
  filter(IMR90 != 'not_altered') %>%
  filter(T2P == 'not_altered') %>%
  filter(TT1 == 'not_altered')

domains <- IMR90_unique
domains <- makeGRangesFromDataFrame(domains)
domains <- reduce(domains)
domains <- domains[order(width(domains), decreasing = T)]
domains <- as.data.frame(domains)

i <- 5

plot_data <- repTiming_df[,c('chr', 'start', 'stop', 'T2P', 'TT1', 'IMR90', 'H1650')] %>%
  filter(chr == domains$seqnames[i]) %>%
  filter(start > domains$start[i] - 1000000) %>%
  filter(stop < domains$end[i] + 1000000) %>%
  data.frame()

plot_data[,c( 'T2P', 'TT1', 'IMR90', 'H1650')] <- apply(plot_data[,c( 'T2P', 'TT1', 'IMR90', 'H1650')] , 2, function(x) impute.loess(x))

plot_data_lines <- plot_data %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop')) %>%
  mutate(type = ifelse(variable == 'H1650', 'cancer', 'normal'),
         variable = sub('_l2r', '', variable),
         start = start / 1000000)


pdf(paste0(output_dir, '/H1650example_IMR90discordant.pdf'), width = 6, height = 3)
ggplot(plot_data_lines, aes(x = start, y = value)) + 
  geom_rect(xmin = domains$start[i] / 1000000, xmax = domains$end[i] / 1000000, ymin = -Inf, ymax = Inf, fill = '#d9d9d9') +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = variable, linetype = variable), width = 2) + 
  scale_colour_manual(name = 'Cell-line', values = c('IMR90' = '#000000','TT1' = '#000000', 'T2P' = '#a50f15', 'H1650' = '#08519c')) +
  scale_linetype_manual(name = 'Cell-line', values = c('IMR90' = 'dotted', 'TT1' = 'dashed', 'T2P' = 'solid', 'H1650' = 'solid')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab('Genomic Position (mb)') + ylab('RT signal') + 
  ggtitle(paste0('IMR90 discordant region on Chromsome ', sub('chr', '',  domains$seqnames[1]))) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()   


#TT1
TT1_unique <- combined_df %>%
  filter(IMR90 == 'not_altered') %>%
  filter(T2P == 'not_altered') %>%
  filter(TT1 != 'not_altered')

domains <- TT1_unique
domains <- makeGRangesFromDataFrame(domains)
domains <- reduce(domains)
domains <- domains[order(width(domains), decreasing = T)]
domains <- as.data.frame(domains)

plot_data <- repTiming_df[,c('chr', 'start', 'stop', 'T2P', 'TT1', 'IMR90', 'H1650')] %>%
  filter(chr == domains$seqnames[2]) %>%
  filter(start > domains$start[2] - 1000000) %>%
  filter(stop < domains$end[2] + 1000000) %>%
  data.frame()

plot_data[,c( 'T2P', 'TT1', 'IMR90', 'H1650')] <- apply(plot_data[,c( 'T2P', 'TT1', 'IMR90', 'H1650')] , 2, function(x) impute.loess(x))

plot_data_lines <- plot_data %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop')) %>%
  mutate(type = ifelse(variable == 'H1650', 'cancer', 'normal'),
         variable = sub('_l2r', '', variable),
         start = start / 1000000)

missedART <- plot_data %>%
  filter(T2P > 0)


pdf(paste0(output_dir, '/H1650example_TT1discordant.pdf'), width = 6, height = 3)
ggplot(plot_data_lines, aes(x = start, y = value)) + 
  geom_rect(xmin = domains$start[2] / 1000000, xmax = domains$end[2] / 1000000, ymin = -Inf, ymax = Inf, fill = '#d9d9d9') +
  geom_rect(xmin = missedART$start[1] / 1000000, xmax = missedART$stop[nrow(missedART)] / 1000000, ymin = -Inf, ymax = Inf, fill = '#fee8c8') +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = variable, linetype = variable), width = 2) + 
  scale_colour_manual(name = 'Cell-line', values = c('IMR90' = '#000000','TT1' = '#000000', 'T2P' = '#a50f15', 'H1650' = '#08519c')) +
  scale_linetype_manual(name = 'Cell-line', values = c('IMR90' = 'dotted', 'TT1' = 'dashed', 'T2P' = 'solid', 'H1650' = 'solid')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab('Genomic Position (mb)') + ylab('RT signal') + 
  ggtitle(paste0('TT1 discordant region on Chromsome ', sub('chr', '',  domains$seqnames[1]))) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()   



#-------- Figure 2 D and SuppFigure 5 --------#
# Upset plots
upset_plots <- lapply(c('H1650', 'H1792', 'H2009', 'A549'), function(x){
  
  #read in artRegions with different cell-of-origins
  T2P_artRegions   <- ARTregions_list[[x]]
  TT1_artRegions   <- TT1_artRegions_list[[x]]
  IMR90_artRegions <- IMR90_artRegions_list[[x]]
  
  #create binary input matrix
  combined_df <- T2P_artRegions[c('chr',  'start',   'stop', 'ARTclass')] %>% 
    dplyr::full_join(TT1_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop')) %>%
    dplyr::full_join(IMR90_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop'))
  
  combined_df <- combined_df[,c(1:3, 4:ncol(combined_df))]
  colnames(combined_df)[4:6] <- c('T2P', 'TT1', 'IMR90')
  
  ww <- apply(combined_df[4:6], 1, function(x) any(is.na(x)))
  combined_df <- combined_df[!ww,]
  
  ww <- apply(combined_df[4:6], 1, function(x) all(x[!is.na(x)] == 'not_altered'))
  combined_df <- combined_df[!ww,]
  
  timing <- sapply(1:nrow(combined_df), function(i){
    data <- as.character(combined_df[i,c('T2P', 'TT1', 'IMR90')])
    data <- data[!is.na(data)]
    data <- data[data != "not_altered"]
    if(length(unique(data)) == 1){
      return(unique(data))
    } else {
      return('discordant')
    }
  })
  combined_df$timing <- timing
  
  input       <- combined_df
  input$T2P   <- ifelse(input$T2P == 'not_altered', 0, 1)
  input$TT1   <- ifelse(input$TT1 == 'not_altered', 0, 1)
  input$IMR90 <- ifelse(input$IMR90 == 'not_altered', 0, 1)
  
  #plot
  p <- upset(input, c('T2P', 'TT1', 'IMR90'), name = '',
             base_annotations=list('# 50kb bins' = intersection_size(counts=FALSE, mapping=aes(fill=timing)) 
                                   + scale_fill_manual(name = 'Replication Timing', values = c('earlier' = '#de77ae', 'later' = '#7fbc41'))
                                   + ggtitle(x)))
  return(p)
  
})

pdf(paste0(output_dir, '/upset_artRegions_wrongNormals_LUAD.pdf'), width = 8, height = 5, useDingbats = F)
upset_plots
dev.off()


# Pie Charts 
pie_plots <- lapply(c('H1650', 'H1792', 'H2009', 'A549'), function(x){
  
  #read in artRegions with different cell-of-origins
  T2P_artRegions   <- ARTregions_list[[x]]
  TT1_artRegions   <- TT1_artRegions_list[[x]]
  IMR90_artRegions <- IMR90_artRegions_list[[x]]
  
  #create binary input matrix
  combined_df <- T2P_artRegions[c('chr',  'start',   'stop', 'ARTclass')] %>% 
    dplyr::full_join(TT1_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop')) %>%
    dplyr::full_join(IMR90_artRegions[,c('chr',  'start',   'stop', 'ARTclass')], by = c('chr', 'start', 'stop'))
  
  combined_df <- combined_df[,c(1:3, 4:ncol(combined_df))]
  colnames(combined_df)[4:6] <- c('T2P', 'TT1', 'IMR90')
  
  ww <- apply(combined_df[4:6], 1, function(x) any(is.na(x)))
  combined_df <- combined_df[!ww,]
  
  ww <- apply(combined_df[4:6], 1, function(x) all(x[!is.na(x)] == 'not_altered'))
  combined_df <- combined_df[!ww,]
  
  input       <- combined_df
  input$T2P   <- ifelse(input$T2P == 'not_altered', 0, 1)
  input$TT1   <- ifelse(input$TT1 == 'not_altered', 0, 1)
  input$IMR90 <- ifelse(input$IMR90 == 'not_altered', 0, 1)
  
  input$sum <- rowSums(input[,c('T2P', 'TT1', 'IMR90')])
  input$class <- 'shared'
  input$class[input$sum == 1] <- 'unique'
  input$class[input$sum == 3] <- 'recurrent'
  
  plot_data <- input %>%
    group_by(class) %>%
    summarise(count = n()) %>%
    mutate(fraction = count / sum(count) * 100) %>%
    mutate(cumulative = cumsum(fraction),
           midpoint = cumulative - fraction / 2,
           label = paste0(round(fraction,1), '%'))
  
  #plot pie
  p <- ggplot(plot_data, aes(x = 1, weights = fraction, fill = class)) +
    geom_bar(width = 1, position = position_stack(reverse = T), colour = 'white') +
    coord_polar(theta = "y") +
    geom_label(aes(x = c(1.2, 1.1, 1), y = midpoint, label = label, colour = class), size = 5) + 
    scale_fill_manual(name = 'Altered RepTiming', values = c('recurrent' = '#762a83', 'shared' = '#c2a5cf', 'unique' = '#bababa')) +
    scale_colour_manual(name = 'Altered RepTiming', values = c('recurrent' = 'white', 'shared' = 'black', 'unique' = 'black')) +
    xlab('') + ylab('') +
    ggtitle(x) +
    theme_void() + 
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5))
  
  return(p)
  
})

pdf(paste0(output_dir, '/pies_artRegions_wrongNormals_LUAD.pdf'), width = 4, height = 3, useDingbats = F)
pie_plots
dev.off()



#-------- Figure 2 E and SuppFigure 6 --------#
overlappingART_list <- readRDS(paste0(data_dir, 'overlappingART.rds'))

#pie charts
pies <- lapply(names(overlappingART_list), function(x){
  overlappingART <- overlappingART_list[[x]]
  plot_data <- overlappingART %>%
    filter(!class %in% c('unknown', 'not_altered')) %>%
    group_by(class) %>%
    summarise(count = n()) %>%
    group_by() %>%
    mutate(total = sum(count),
           fraction = (count / total) *100,
           cumulative = cumsum(fraction),
           midpoint = cumulative - fraction / 2,
           label = paste0(round(fraction), '%'))
  
  ggplot(plot_data, aes(x = 1, weights = fraction, fill = class)) +
    geom_bar(width = 1, position = position_stack(reverse = T), colour = 'white') +
    coord_polar(theta = "y") +
    geom_label(aes(x = c(1.2, 1.1, 1), y = midpoint, label = label, colour = class) , size = 5) + 
    scale_fill_manual(name = '', values = c('recurrent' = '#762a83', 'shared' = '#c2a5cf', 'unique' = '#bababa')) +
    scale_colour_manual(name = '', values = c('recurrent' = 'white', 'shared' = 'black', 'unique' = 'black')) +
    xlab('') + ylab('') +
    ggtitle(x) +
    theme_void() + 
    theme(axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5))
})

pdf(paste0(output_dir,'/pies_overlappingART.pdf'), width = 4, height = 3, useDingbats = F)
pies
dev.off()





# upset-plots
upsetPlots <- lapply(names(overlappingART_list), function(x){
  
  if(x == 'LUAD'){
    cancer <- c('H1650', 'H1792', 'H2009', 'A549')
  } else if (x == 'LUSC'){
    cancer <- c('H520', 'H2170', 'SW900')
  } else if (x == 'BRCA'){
    cancer <- c('MCF-7', 'MDA453', 'SK-BR3', 'T47D')
  }
  
  overlappingART <- overlappingART_list[[x]]
  input <- overlappingART %>%
    filter(!class %in% c('unknown', 'not_altered')) %>%
    select(chr, start, stop, paste0(cancer, '_timing'), 'timing')
  input[,paste0(cancer, '_timing')] <- apply(input[,paste0(cancer, '_timing')], 2, function(x) ifelse(x %in% c('earlier', 'later'), 1, 0))
  colnames(input) <- sub('_timing', '', colnames(input))
  
  upset(input, cancer, name = '',
        base_annotations=list('# 50kb bins' = intersection_size(counts=FALSE, mapping=aes(fill=timing)) 
                              + scale_fill_manual(name = 'Replication Timing', values = c('earlier' = '#de77ae', 'later' = '#7fbc41')) 
                              + ggtitle(x)))
  
})

pdf(paste0(output_dir, '/upset_ARTregions.pdf'), width = 7, height = 4, useDingbats = F)
upsetPlots
dev.off()


#check random overlapps
LUAD_random_overlaps   <- check_random_overlaps(ARTregions_list[c('H1650', 'H1792', 'H2009', 'A549')])
LUSC_random_overlaps   <- check_random_overlaps(ARTregions_list[c('H520', 'H2170', 'SW900')])
BRCA_random_overlaps   <- check_random_overlaps(ARTregions_list[c('MCF-7', 'MDA453', 'SK-BR3', 'T47D')])

LUAD_overlappingART <- overlappingART_list$LUAD
LUSC_overlappingART <- overlappingART_list$LUSC
BRCA_overlappingART <- overlappingART_list$BRCA

LUAD_plot_data <- LUAD_random_overlaps %>%
  mutate(overlap = recurrent + shared) %>%
  summarise(mean = mean(overlap), 
            lowCI = quantile(overlap, 0.025),
            highCI = quantile(overlap, 0.975)) %>%
  mutate(cancerType = 'LUAD') %>%
  mutate(observed = sum(LUAD_overlappingART$class %in% c('recurrent', 'shared')))

LUSC_plot_data <- LUSC_random_overlaps %>%
  mutate(overlap = recurrent + shared) %>%
  summarise(mean = mean(overlap), 
            lowCI = quantile(overlap, 0.025),
            highCI = quantile(overlap, 0.975)) %>%
  mutate(cancerType = 'LUSC') %>%
  mutate(observed = sum(LUSC_overlappingART$class %in% c('recurrent', 'shared')))

BRCA_plot_data <- BRCA_random_overlaps %>%
  mutate(overlap = recurrent + shared) %>%
  summarise(mean = mean(overlap), 
            lowCI = quantile(overlap, 0.025),
            highCI = quantile(overlap, 0.975)) %>%
  mutate(cancerType = 'BRCA') %>%
  mutate(observed = sum(BRCA_overlappingART$class %in% c('recurrent', 'shared')))

plot_data <- rbind(LUAD_plot_data, LUSC_plot_data, BRCA_plot_data)

pdf(paste0(output_dir,'/random_ART_overlap_pointrange.pdf'), width = 5, height = 3, useDingbats = F)
ggplot(plot_data) + 
  geom_errorbar(aes(x = cancerType, ymin = lowCI, ymax = highCI), width = 0.2) + 
  geom_point(aes(x = cancerType, y = mean, colour = 'random (mean)', shape = 'random (mean)'), fill = 'black') + 
  geom_point(aes(x = cancerType, y = observed, colour = 'observed', shape = 'observed'), size = 3) + 
  scale_colour_manual(name = '', values = c('random (mean)'= 'black',  'observed' = '#b2182b')) + 
  scale_shape_manual(name = '', values = c('random (mean)' = 21,  'observed' = 8)) + 
  coord_flip() +
  xlab('') + ylab('# overlapping ART 50kb bins') +
  theme_bw() + theme(legend.position = 'top')
dev.off()


#examples for recurrent, shared and unique ART
# recurrent BRCA #
BRCA_overlappingART <- overlappingART_list$BRCA
BRCA_overlappingART <- BRCA_overlappingART %>%
  group_by(chr) %>%
  arrange(desc(start))

recurrent_domains <- BRCA_overlappingART[BRCA_overlappingART$class == 'recurrent',]
recurrent_domains <- makeGRangesFromDataFrame(recurrent_domains)
recurrent_domains <- reduce(recurrent_domains)
recurrent_domains <- recurrent_domains[order(width(recurrent_domains), decreasing = T)]
recurrent_domains <- as.data.frame(recurrent_domains)

i <- 5

plot_data <- BRCA_overlappingART %>%
  filter(chr == recurrent_domains$seqnames[i]) %>%
  filter(start > recurrent_domains$start[i] - 1000000) %>%
  filter(stop < recurrent_domains$end[i] + 1000000) %>%
  data.frame()

plot_data[,grep('l2r', colnames(plot_data))] <- apply(plot_data[,grep('l2r', colnames(plot_data))] , 2, function(x) impute.loess(x))

plot_data_lines <- plot_data[,c('chr', 'start', 'stop', grep('l2r', colnames(plot_data), value = T))] %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop')) %>%
  mutate(type = ifelse(variable == 'HHMEC_l2r', 'normal', 'cancer'),
         variable = sub('_l2r', '', variable),
         start = start / 1000000)


pdf(paste0(output_dir, '/BRCAexample_recurrentART.pdf'), width = 6, height = 3)
ggplot(plot_data_lines, aes(x = start, y = value)) + 
  geom_rect(xmin = recurrent_domains$start[i] / 1000000, xmax = recurrent_domains$end[i] / 1000000, ymin = -Inf, ymax = Inf, fill = '#762a83') +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = variable, linetype = variable), width = 2) + 
  scale_colour_manual(name = 'Cell-line', values = c('HMEC' = 'black', 'MCF.7' = '#a63603', 'MDA453' = '#d94801', 'SK.BR3' = '#f16913', 'T47D' = '#fd8d3c'),
                      labels = c('HMEC' = 'HMEC', 'MCF.7' = 'MCF-7', 'MDA453' = 'MDA453', 'SK.BR3' = 'SK-BR3', 'T47D' = 'T47D')) +
  scale_linetype_manual(name = 'Cell-line', values = c('HMEC' = 'dashed', 'MCF.7' = 'solid', 'MDA453' = 'solid', 'SK.BR3' = 'solid', 'T47D' = 'solid'),
                        labels = c('HMEC' = 'HMEC', 'MCF.7' = 'MCF-7', 'MDA453' = 'MDA453', 'SK.BR3' = 'SK-BR3', 'T47D' = 'T47D')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab('Genomic Position (mb)') + ylab('RT signal') + 
  ggtitle(paste0('BRCA - Recurrent ART on Chromsome ', sub('chr', '',  recurrent_domains$seqnames[1]))) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()                   


# shared LUAD #
LUAD_overlappingART <- overlappingART_list$LUAD
LUAD_overlappingART <- LUAD_overlappingART %>%
  group_by(chr) %>%
  arrange(desc(start))

shared_domains <- LUAD_overlappingART[LUAD_overlappingART$class == 'shared',]
shared_domains <- makeGRangesFromDataFrame(shared_domains)
shared_domains <- reduce(shared_domains)
shared_domains <- shared_domains[order(width(shared_domains), decreasing = T)]
shared_domains <- as.data.frame(shared_domains)

i <- 2

plot_data <- LUAD_overlappingART %>%
  filter(chr == shared_domains$seqnames[i]) %>%
  filter(start > shared_domains$start[i] - 1000000) %>%
  filter(stop < shared_domains$end[i] + 1000000) %>%
  data.frame()

plot_data[,grep('l2r', colnames(plot_data))] <- apply(plot_data[,grep('l2r', colnames(plot_data))] , 2, function(x) impute.loess(x))

plot_data_lines <- plot_data[,c('chr', 'start', 'stop', grep('l2r', colnames(plot_data), value = T))] %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop')) %>%
  mutate(type = ifelse(variable == 'T2P_l2r', 'normal', 'cancer'),
         variable = sub('_l2r', '', variable),
         start = start / 1000000)


pdf(paste0(output_dir, '/LUADexample_sharedART.pdf'), width = 6, height = 3)
ggplot(plot_data_lines, aes(x = start, y = value)) + 
  geom_rect(xmin = shared_domains$start[i] / 1000000, xmax = shared_domains$end[i] / 1000000, ymin = -Inf, ymax = Inf, fill = '#c2a5cf') +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = variable, linetype = variable), width = 2) + 
  scale_colour_manual(name = 'Cell-line', values = c('T2P' = 'black', 'H1650' = '#08519c', 'H1792' = '#2171b5', 'H2009' = '#4292c6', 'A549' = '#6baed6')) +
  scale_linetype_manual(name = 'Cell-line', values = c('T2P' = 'dashed', 'H1650' = 'solid', 'H1792' = 'solid', 'H2009' = 'solid', 'A549' = 'solid')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab('Genomic Position (mb)') + ylab('RT signal') + 
  ggtitle(paste0('LUAD - Shared ART on Chromsome ', sub('chr', '',  shared_domains$seqnames[i]))) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()                   


# unique LUSC #
LUSC_overlappingART <- overlappingART_list$LUSC
LUSC_overlappingART <- LUSC_overlappingART %>%
  group_by(chr) %>%
  arrange(desc(start))

unique_domains <- LUSC_overlappingART[LUSC_overlappingART$class == 'unique',]
unique_domains <- makeGRangesFromDataFrame(unique_domains)
unique_domains <- reduce(unique_domains)
unique_domains <- unique_domains[order(width(unique_domains), decreasing = T)]
unique_domains <- as.data.frame(unique_domains)

i <- 6

plot_data <- LUSC_overlappingART %>%
  filter(chr == unique_domains$seqnames[i]) %>%
  filter(start > unique_domains$start[i] - 1000000) %>%
  filter(stop < unique_domains$end[i] + 1000000) %>%
  data.frame()

plot_data[,grep('l2r', colnames(plot_data))] <- apply(plot_data[,grep('l2r', colnames(plot_data))] , 2, function(x) impute.loess(x))

plot_data_lines <- plot_data[,c('chr', 'start', 'stop', grep('l2r', colnames(plot_data), value = T))] %>%
  reshape2::melt(id.vars = c('chr', 'start', 'stop')) %>%
  mutate(type = ifelse(variable == 'HBEC_l2r', 'normal', 'cancer'),
         variable = sub('_l2r', '', variable),
         start = start / 1000000)


pdf(paste0(output_dir, '/LUSCexample_uniqueART.pdf'), width = 6, height = 3)
ggplot(plot_data_lines, aes(x = start, y = value)) + 
  geom_rect(xmin = unique_domains$start[i] / 1000000, xmax = unique_domains$end[i] / 1000000, ymin = -Inf, ymax = Inf, fill = '#bababa') +
  geom_hline(yintercept = 0) +
  geom_line(aes(colour = variable, linetype = variable), width = 2) + 
  scale_colour_manual(name = 'Cell-line', values = c('HBEC3' = 'black', 'H520' = '#006d2c', 'H2170' = '#238b45', 'SW900' = '#41ab5d')) +
  scale_linetype_manual(name = 'Cell-line', values = c('HBEC3' = 'dashed', 'H520' = 'solid', 'H2170' = 'solid', 'SW900' = 'solid')) +
  scale_x_continuous(expand = c(0,0)) +
  xlab('Genomic Position (mb)') + ylab('RT signal') + 
  ggtitle(paste0('LUSC - Unique ART on Chromsome ', sub('chr', '',  unique_domains$seqnames[1]))) +
  theme_bw() + theme(panel.grid = element_blank())
dev.off()  


