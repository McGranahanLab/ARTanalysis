# Data file descriptions

For more detailed information about the data have a look at the methods section of the paper: "Replication timing alterations are associated with mutation acquisition during tumour evolution in breast and lung cancer" by Dietzen and Zhai et al.


* 20220203_pan.driver_og.tsg.csv: List of pan-cancer driver genes including annotation of tumour suppressor genes and oncogenes
             
* 20220203_pancan_drivers.csv: List of driver genes for different cancer types

* 20220203_lung_drivers.csv: List of genes identified as driver genes for lung cancer  
                        
* 20220208_breast.driver.csv: List of genes identified as driver genes for breast cancer 

* tissueInfo_cellLines_20210309.tsv: Table with tissue information of cell lines used in this study.

* A549encode.l2r.bedGraph: Raw replication timing signal from the A549 cell line sequenced by ENCODE

* A549.l2r.bedGraph: Raw replication timing signal from the A549 cell line sequenced by us (IN_STUDY)

* H1650.l2r.bedGraph: Raw replication timing signal from the H1650 cell line sequenced by us (IN_STUDY)

* H1650rep.l2r.bedGraph: Raw replication timing signal from the H1650 cell line replicate sequenced by us (IN_STUDY)

* T2P.l2r.bedGraph: Raw replication timing signal from the T2P cell line sequenced by us (IN_STUDY)

* T2Prep.l2r.bedGraph: Raw replication timing signal from the T2P cell line replicate sequenced by us (IN_STUDY)

* A549encode.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the A549 cell line sequenced by ENCODE

* A549.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the A549 cell line sequenced by us (IN_STUDY)

* H1650.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the H1650 cell line sequenced by us (IN_STUDY)

* H1650rep.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the H1650 cell line replicate sequenced by us (IN_STUDY)

* T2P.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the T2P cell line sequenced by us (IN_STUDY)

* T2Prep.qnorm_l2r.bedGraph: Quantile normalised replication timing signal from the T2P cell line replicate sequenced by us (IN_STUDY)

* A549encode.loess300000.bedGraph: Loess smoothed replication timing signal from the A549 cell line sequenced by ENCODE

* A549.loess300000.bedGraph: Loess smoothed replication timing signal from the A549 cell line sequenced by us (IN_STUDY)

* H1650.loess300000.bedGraph: Loess smoothed replication timing signal from the H1650 cell line sequenced by us (IN_STUDY)

* H1650rep.loess300000.bedGraph: Loess smoothed replication timing signal from the H1650 cell line replicate sequenced by us (IN_STUDY)

* T2P.loess300000.bedGraph: Loess smoothed replication timing signal from the T2P cell line sequenced by us (IN_STUDY)

* T2Prep.loess300000.bedGraph: Loess smoothed replication timing signal from the T2P cell line replicate sequenced by us (IN_STUDY)

* ARTregions_full.rds: ART regions including no switch regions for 4 BRCA and 4 LUAD cell lines. No switch regions are genomic regions that pass the threshold in replication timing signal difference between normal and cancer, but those regions were classified as either early in both normal anc caner or late.

* ARTregions.rds: ART regions excluding no switch regions for 4 BRCA and 4 LUAD cell lines.

* overlappingART.rds: ART regions classified as unique, shared and recurrent for the LUAD and BRCA cell lines

* BRCA_consistentARTregions.txt: ART regions that were classified as early, late, early-to-late or late-to-early consistently across all 4 BRCA cell lines

* LUAD_consistentARTregions.txt: ART regions that were classified as early, late, early-to-late or late-to-early consistently across all 4 LUAD cell lines

* BRCA_overlapping_ARTgenes.rds: ART genes classified as unique, shared and recurrent for the BRCA cell lines

* LUAD_overlapping_ARTgenes.rds: ART genes classified as unique, shared and recurrent for the LUAD cell lines

* resultsARTgenes.extreme_list.RDS: ART and consistent RT per cancer cell line compared to their matched tissue-of-origin (T2P as reference for LUAD, HMEC as reference for BRCA)

* BRCA_sharedARTregions.txt: ART regions that were classified as ART within at least 2 out of the 4 BRCA cell lines

* LUAD_sharedARTregions.txt: ART regions that were classified as ART within at least 2 out of the 4 LUAD cell lines

* BRCA_sharedARTgenes.txt: Genes classified as ART presenting in at least 2 out of the 4 BRCA cell lines

* LUAD_sharedARTgenes.txt: Genes classified as ART presenting in at least 2 out of the 4 LUAD cell lines

* IMR90_ARTregions.rds: ART regions for the 4 LUAD cell lines when using IMR90 as normal reference

* TT1_ARTregions.rds: ART regions for the 4 LUAD cell lines when using TT1 as normal reference

* MCF10A_ARTregions.rds: ART regions for the 4 BRCA cell lines when using MCF10A as normal reference

* cohort_50kb_l2r.rds: Replication timing signal (log2-ratio) for 31 cell-lines (16 ENCODE and 15 IN-STUDY)

* cohort_meanL2R_genes.rds: Mean replication timing signal (log2-ratio) per gene for 31 cell-lines (16 ENCODE and 15 IN-STUDY)

* COSMIC_v3.2_SBS_GRCh37.txt: Mutational signature profiles downloaded from COSMIC

* 560Breast_subset_cnTable.RData: Copy number information from 482 lobular and ductal breast cancer tumours from the breast cancer study: "Landscape of somatic mutations in 560 breast cancer whole-genome sequences" by Nik-Zainal et al. (2016)

* 560Breast_subset_mutTable_[1-3].rds: Mutation information from 482 lobular and ductal breast cancer tumours from the breast cancer study: "Landscape of somatic mutations in 560 breast cancer whole-genome sequences" by Nik-Zainal et al. (2016). The data combined data was too large to upload on github which is why it had been split into 3 files.

* exposures_SBS_BRCA.rds: Mutation signature exposures for the 482 BRCA tumours across unaltered RT and ART regions.

* signatures_SBS_BRCA.rds: Mutation signature profiles extracted from the 482 BRCA tumours across unaltered RT and ART regions. 

* input_96matrix_patient_sharedRepTiming.rds: 96 trinucleotide counts for the 482 BRCA tumours in unaltered RT and ART regions as input to the HDP signature extraction pipeline

* hyperClust_table_BRCA.rds: HyperClust results to identify kataegis and omikli events for the 482 BRCA tumours

* mutTable_DepMap_20210802.rds: Mutation information for most cancer cell lines (ENCODE and IN-STUDY) downloaded from DepMap

* HiC_BRCA.rds: Hi-C data for 2 BRCA cell lines and the normal (tissue-of-origin) cell line HMEC

* overlap_HiC_ART_BRCA.tsv: ART regions from T47D and MCF-7 (relative to HMEC) classified as chromatin compartent A and B using Hi-C data. 

* overlap_HiC_ART_LUAD.tsv: ART regions from A549 (relative to T2P) classified as chromatin compartent A and B using Hi-C data. 

* PDC_ART.rds: ART regions for the 2 TRACERx patient derived cell lines

* PDC_mutTable.rds: Mutation information from the WGS data of the 2 TRACERx patient derived cell lines

* mmc2_Histologic.BRCA.inTCGA.xlsx: Data S1 from (https://www.cell.com/cell-genomics/fulltext/S2666-979X(21)00083-5#supplementaryMaterial) containing histological annotations of BRCA tumours which has been used to filter for ductal and lobular subtypes in the TCGA analysis (https://www.cell.com/cell-genomics/fulltext/S2666-979X(21)00083-5#supplementaryMaterial)

* LUAD_mean_cnDiff_perGene.rds: The mean copy number value per gene was calculated by overlapping the copy number segments of TCGA tumours with a certain gene and averaging the copy number relative to the overlapping size of the segments using WES data of LUAD tumours in TCGA

* BRCAsubset_mean_cnDiff_perGene.rds: The mean copy number value per gene which was calculated by overlapping the copy number segments of TCGA tumours with a certain gene and averaging the copy number relative to the overlapping size of the segments using WES data of ductal and lobular BRCA tumours in TCGA

* CN.cohort_sample.list.RDS: The list of samples related to LUAD and BRCA in TCGA with exome-sequencing data available and passing the quality control





















