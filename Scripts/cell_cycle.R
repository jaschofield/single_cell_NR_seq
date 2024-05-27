## Cell cycle calculations
library(tidyverse)
library(pheatmap)
library(purrr)

## bring in yeast coactivator class dataset from Donczew et al. 2020
donczew_data <- read.csv(file = "//donczew_2020_data.csv")
colnames(donczew_data) <- c("UAS_name", "gene", "TATA", "class", "tail_dep", "endo_expression", "Spt3",
                        "Spt7", "Taf1", "Taf13", "Bdf1Bdf2", "Esa1", "Med1", "Med5", "Med10",
                        "Med14", "Med15", "TID", "Mot1", "Bur6", "NCB2", "Toa1", "Ssl2")

## bring in cell cycle data from Spellman et al. 1998
spellman_DF <- read.table(file = "orf_data.txt")
colnames(spellman_DF) <- c("gene_ID", "Gene_old", "Score", "plotX", "plotY", "ATAN", "cycle")
spellman_DF <- left_join(spellman_DF, donczew_data, by = "gene_ID") %>% filter(!is.na(gene_ID))

## Select highly cell cycle reglated yeast genes from each cycle
spellman_filt_G1 <- spellman_DF %>% filter(cycle %in% c("G1") & Score > 2) %>% dplyr::select(gene_ID) %>% data.frame()
spellman_filt_G2 <- spellman_DF %>% filter(cycle %in% c("G2") & Score > 2) %>% dplyr::select(gene_ID) %>% data.frame()
spellman_filt_S <- spellman_DF %>% filter(cycle == "S" & Score > 2) %>% dplyr::select(gene_ID) %>% data.frame()
spellman_filt_M_G1 <- spellman_DF %>% filter(cycle %in% c("M/G1") & Score > 2) %>% dplyr::select(gene_ID) %>% data.frame()
spellman_filt_M <- spellman_DF %>% filter(cycle %in% c("M") & Score > 2) %>% dplyr::select(gene_ID) %>% data.frame()


name_joiner <- donczew_data %>% dplyr::select(gene_ID, gene)
colnames(name_joiner) <- c("UAS_name", "EF")

## bring in filtered dataset (10 minute timepoint used in manuscript)
Nascent_filt <- read.csv(file = "Nascent_filt.csv")

## Combine nascent RNA count dataframe with naming column
Nascent_filt_form <- left_join(Nascent_filt, name_joiner)

## Generate scores for each cell based on highly cell cycle regulated genes
G1_nascent_scores <- Nascent_filt_form %>% filter(UAS_name %in% spellman_filt_G1$gene_ID) %>% dplyr::select(-UAS_name, -EF) %>% colSums()
G2_nascent_scores <- Nascent_filt_form %>% filter(UAS_name %in% spellman_filt_G2$gene_ID) %>% dplyr::select(-UAS_name, -EF) %>% colSums()
S_nascent_scores <- Nascent_filt_form %>% filter(UAS_name %in% spellman_filt_S$gene_ID) %>% dplyr::select(-UAS_name, -EF) %>% colSums()
M_nascent_scores <- Nascent_filt_form %>% filter(UAS_name %in% spellman_filt_M$gene_ID) %>% dplyr::select(-UAS_name, -EF) %>% colSums()
M_G1_nascent_scores <- Nascent_filt_form %>% filter(UAS_name %in% spellman_filt_M_G1$gene_ID) %>% dplyr::select(-UAS_name, -EF) %>% colSums()


scored_nascent_cells <- data.frame(G1_nascent_scores, G2_nascent_scores, S_nascent_scores, M_nascent_scores, M_G1_nascent_scores)

## normalizing cell scores to the max count in each cycle and generating a fraction of expression for each cycle
## Hard coded values are max scores from 10 minute sample
scored_nascent_cells <- scored_nascent_cells %>%
  mutate(G1_norm = G1_nascent_scores / 48,
         G2_norm = G2_nascent_scores / 59,
         S_norm = S_nascent_scores / 115,
         M_norm = M__nascent_scores / 36,
         M_G1_norm = M_G1_nascent_scores / 39) %>%
  mutate(G1_frac = G1_norm / (G1_norm + G2_norm + S_norm + M_norm + M_G1_norm),
         G2_frac = G2_norm / (G1_norm + G2_norm + S_norm + M_norm + M_G1_norm),
         S_frac = S_norm / (G1_norm + G2_norm + S_norm + M_norm + M_G1_norm),
         M_frac = M_norm / (G1_norm + G2_norm + S_norm + M_norm + M_G1_norm),
         M_G1_frac = M_G1_norm / (G1_norm + G2_norm + S_norm + M_norm + M_G1_norm))

## clustering cells based on the fraction of expression in each cycle
cycle_cor <- scored_nascent_10min_cells %>% dplyr::select(contains("frac"))
cycle_cor <- na.omit(cycle_cor)
clust_out <- pheatmap::pheatmap(cycle_cor, kmeans_k = 5)
cycle_assignment_clusters <- clust_out$kmeans$cluster %>% data.frame()
colnames(cycle_assignment_clusters) <- "cluster"
cycle_assignment_clusters <- cycle_assignment_clusters %>%
  mutate(cycle_assignment = ifelse(cluster == 1, "S",
                                   ifelse(cluster == 2, "G1",
                                          ifelse(cluster == 3, "M",
                                                 ifelse(cluster == 4, "G2",
                                                        ifelse(cluster == 5, "M_G1", NA))))))

## Making lists of cells assigned to each cycle
G1_cells <- cycle_assignment_clusters %>% filter(cycle_assignment == "G1") %>% rownames() %>% gsub("\\.", "-", .)
S_cells <- cycle_assignment_clusters %>% filter(cycle_assignment == "S") %>% rownames() %>% gsub("\\.", "-", .)
G2_cells <- cycle_assignment_clusters %>% filter(cycle_assignment == "G2") %>% rownames() %>% gsub("\\.", "-", .)
M_cells <- cycle_assignment_clusters %>% filter(cycle_assignment == "M") %>% rownames() %>% gsub("\\.", "-", .)
M_G1_cells <- cycle_assignment_clusters %>% filter(cycle_assignment == "M_G1") %>% rownames() %>% gsub("\\.", "-", .)

Nascent_G1 <- Nascent_filt %>% dplyr::select(EF, G1_cells)
Nascent_S <- Nascent_filt %>% dplyr::select(EF, S_cells)
Nascent_G2 <- Nascent_filt %>% dplyr::select(EF, G2_cells)
Nascent_M_G1 <- Nascent_filt %>% dplyr::select(EF, M_G1_cells)
Nascent_M <- Nascent_filt %>% dplyr::select(EF, M_cells)

## Calculating Fon for cells assigned to each cycle, using same ten-fold above background cutoff as full dataset
G1_ratios <- T10_sums %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>% filter(cell_code %in% G1_cells) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(G1_ratios) <- c("gene", "G1_cells")
G2_ratios <- T10_sums %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>% filter(cell_code %in% G2_cells) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(G2_ratios) <- c("gene", "G2_cells")
S_ratios <- T10_sums %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>% filter(cell_code %in% S_cells) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(S_ratios) <- c("gene", "S_cells")
M_ratios <- T10_sums %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>% filter(cell_code %in% M_cells) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(M_ratios) <- c("gene", "M_cells")
M_G1_ratios <- T10_sums %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>% filter(cell_code %in% M_G1_cells) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(M_G1_ratios) <- c("gene", "M_G1_cells")

cell_cycle_Fon <- list(G1_ratios, S_ratios, G2_ratios, M_ratios, M_G1_ratios) %>%
  purrr::reduce(full_join, by = "gene")

## Fon calculation, hard coded values are number of cells generated in clustering for the 10 minute sample
cell_cycle_Fon <- cell_cycle_Fon %>% replace(is.na(.), 0) %>%
  mutate(G1_Fon = G1_cells/2941,
         G2_Fon = G2_cells/3403,
         M_Fon = M_cells/3683,
         M_G1_Fon = M_G1_cells/3317,
         S_Fon = S_cells/2789)

## Calculating Fano for cells assigned to each cycle
G1_var <- Nascent_G1 %>% mutate(row_wise_var = rowVars(as.matrix(Nascent_G1[,c(2:ncol(Nascent_G1))]))) %>% dplyr::select(EF, row_wise_var)
colnames(G1_var) <- c("gene", "G1_var")
G1_mean <- Nascent_G1 %>% dplyr::select(-EF) %>% rowMeans()
G1_mean <- data.frame(Nascent_G1$EF, G1_mean)
colnames(G1_mean) <- c("gene", "G1_mean")
G1_fano_DF <- full_join(G1_var, G1_mean) %>% mutate(G1_fano = G1_var/G1_mean)

G2_var <- Nascent_G2 %>% mutate(row_wise_var = rowVars(as.matrix(Nascent_G2[,c(2:ncol(Nascent_G2))]))) %>% dplyr::select(EF, row_wise_var)
colnames(G2_var) <- c("gene", "G2_var")
G2_mean <- Nascent_G2 %>% dplyr::select(-EF) %>% rowMeans()
G2_mean <- data.frame(Nascent_G2$EF, G2_mean)
colnames(G2_mean) <- c("gene", "G2_mean")
G2_fano_DF <- full_join(G2_var, G2_mean) %>% mutate(G2_fano = G2_var/G2_mean)

M_var <- Nascent_M %>% mutate(row_wise_var = rowVars(as.matrix(Nascent_M[,c(2:ncol(Nascent_M))]))) %>% dplyr::select(EF, row_wise_var)
colnames(M_var) <- c("gene", "M_var")
M_mean <- Nascent_M %>% dplyr::select(-EF) %>% rowMeans()
M_mean <- data.frame(Nascent_M$EF, M_mean)
colnames(M_mean) <- c("gene", "M_mean")
M_fano_DF <- full_join(M_var, M_mean) %>% mutate(M_fano = M_var/M_mean)

S_var <- Nascent_S %>% mutate(row_wise_var = rowVars(as.matrix(Nascent_S[,c(2:ncol(Nascent_S))]))) %>% dplyr::select(EF, row_wise_var)
colnames(S_var) <- c("gene", "S_var")
S_mean <- Nascent_S %>% dplyr::select(-EF) %>% rowMeans()
S_mean <- data.frame(Nascent_S$EF, S_mean)
colnames(S_mean) <- c("gene", "S_mean")
S_fano_DF <- full_join(S_var, S_mean) %>% mutate(S_fano = S_var/S_mean)

M_G1_var <- Nascent_M_G1 %>% mutate(row_wise_var = rowVars(as.matrix(Nascent_M_G1[,c(2:ncol(Nascent_M_G1))]))) %>% dplyr::select(EF, row_wise_var)
colnames(M_G1_var) <- c("gene", "M_G1_var")
M_G1_mean <- Nascent_M_G1 %>% dplyr::select(-EF) %>% rowMeans()
M_G1_mean <- data.frame(Nascent_M_G1$EF, M_G1_mean)
colnames(M_G1_mean) <- c("gene", "M_G1_mean")
M_G1_fano_DF <- full_join(M_G1_var, M_G1_mean) %>% mutate(M_G1_fano = M_G1_var/M_G1_mean)

cell_cycle_fano_df <- full_join(G1_fano_DF, G2_fano_DF)
cell_cycle_fano_df <- full_join(cell_cycle_fano_df, M_fano_DF)
cell_cycle_fano_df <- full_join(cell_cycle_fano_df, S_fano_DF)
cell_cycle_fano_df <- full_join(cell_cycle_fano_df, M_G1_fano_DF)
