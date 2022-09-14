
#####################################################################
## Supplementary Figure 3E/F/G: Conserved RT and ART
#####################################################################
# written by Haoran Zhai (haoran.zhai.17@ucl.ac.uk) and run in R version 4.0.2

# Description:
# Scripts to create Supplementary Figure 3E/F/G in the manuscript named "Replication timing alterations impact mutation acquisition during tumour evolution".
# Data accessibility statement can be found in the manuscript.

#libraries
options(stringsAsFactors = F)
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

### Supplementary Figure 3E: conserved RT (in all 10 normal) ========
tissue_info <- read.table(paste0(data_dir, 'tissueInfo_cellLines_20210309.tsv'), header = T, sep = '\t')
tissue_info$cellLine <- sub('A549$', 'A549encode', tissue_info$cellLine)
tissue_info$cellLine <- sub('A549rep', 'A549', tissue_info$cellLine)

#load log2-ratios
log2ratio_df        <- readRDS(paste0(data_dir,'cohort_50kb_l2r.rds'))
normal_cellLines    <- colnames(log2ratio_df)[colnames(log2ratio_df) %in% tissue_info$cellLine[tissue_info$tissueType == 'Normal']]
normal_log2ratio_df <- log2ratio_df[,c('chr', 'start', 'stop', normal_cellLines)]

#load log2-ratios per gene
log2ratio_genes        <- readRDS(paste0(data_dir, 'cohort_meanL2R_genes.rds') )
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
### Supplementary Figure 3F-G: ART per cancer type being conserved RT or not ========
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

# Supplementary Figure 3F:
pdf(paste0(output_dir, 'bar_frac.conserRT_ART.pdf'), width=6.5, height = 5)
print(p_in.ART)
dev.off()

# Supplementary Figure 3G:
pdf(paste0(output_dir, 'bar_frac.ART_conserved.RT.pdf'), width=7, height = 5)
print(p_in.conserved.RT)
dev.off()

