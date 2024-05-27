## Joining cellranger BAM with TL csv
library(tidyverse)
library(Rsamtools)
library(purrr)

## bring in TimeLapse-seq pipeline counts
TL_counts <- read.csv(file = "/path_to_TimeLapse/TimeLapse.dir/sample_counts.csv")

## bring in BAM file from CellRanger
CR_bam <- scanBam(param = ScanBamParam(what = c("qname"), tag = c("CB", "UB")), file = "/path_to_Cellranger/outs/possorted_genome_bam.bam")

## Check formatting
names(CR_bam[[1]])

## Bring in list of cell-associated barcodes from CellRanger
CR_barcodes <- read.csv(file = "/path_to_Cellranger/outs/filtered_feature_bc_matrix/barcodes.tsv", header = FALSE)

## unlist function
.unlist <- function (x){
  x1 <- x[[1L]]
  if (is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

#store names of BAM fields
bam_field <- names(CR_bam[[1]])

#go through each BAM field and unlist
list <- lapply(bam_field, function(y) .unlist(lapply(CR_bam, "[[", y)))

#store as data frame
bam_df <- do.call("DataFrame", list)
names(bam_df) <- bam_field
bam_df <- data.frame(bam_df)
bam_df_filt <- bam_df %>% filter(tag %in% CR_barcodes$V1)
colnames(bam_df_filt) <- c("qname", "cell_code", "mol_code")

## Join reads together by their read name and filter for gene-associated reads
joined_df <- left_join(TL_counts, bam_df_filt, by = "qname")

TL_filt <- joined_df %>% filter(!is.na(EF))

write.csv(TL_filt, file = "TL_filt.csv")

## mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## Select for the most commonly appearing observation and collapse
TL_collapsed <- TL_filt %>% dplyr::select(cell_code, mol_code, nT, TC, EF) %>%
  group_by(cell_code, mol_code, EF) %>% filter(nT == getmode(nT)) %>% distinct(.) %>% ungroup()

import_prob_table <- read.csv(file = "/path/Binom_prob_table.csv")
TL_collapsed <- left_join(TL_collapsed, import_prob_table, by = "nT") %>% filter(!is.na(X1mut_prob))

## execute simulation to sample out potential background mutations and collapse into nascent counts
set.seed(1)
TL_mut <- TL_collapsed %>% rowwise() %>% mutate(nascent = ifelse(TC == 0, 0,
                                                                           ifelse(TC == 1, sample(c(0, 1), size = 1, replace = TRUE, prob = c(1-`X1mut_prob`, `X1mut_prob`)),
                                                                                  ifelse(TC == 2, sample(c(0, 1), size = 1, replace = TRUE, prob = c(1-`X2mut_prob`, `X2mut_prob`)),
                                                                                         ifelse(TC == 3, sample(c(0, 1), size = 1, replace = TRUE, prob = c(1-`X3mut_prob`, `X3mut_prob`)), 1
))))) %>%
dplyr::select(cell_code, nascent, EF) %>% pivot_wider(names_from = cell_code, values_from = nascent, values_fn = sum)
head(TL_mut)


write.csv(TL_mut, file = "Sample_nascent_counts.csv")

## collapse into total counts
TL_collapsed_tot <- TL_filt %>% dplyr::select(cell_code, mol_code, nT, TC, EF) %>%
  group_by(cell_code, mol_code, EF) %>% filter(nT == getmode(nT)) %>% distinct(.) %>% ungroup() %>%
  mutate(nascent = ifelse(TC >-1, 1, 0)) %>% dplyr::select(cell_code, nascent, EF) %>%
  pivot_wider(names_from = cell_code, values_from = nascent, values_fn = sum)
head(TL_collapsed_tot)
write.csv(TL_collapsed_tot, file = "Sample_total_counts.csv")

## create cB for half life analysis
TL_cB <- TL_filt %>% dplyr::select(cell_code, mol_code, nT, TC, EF) %>%
  group_by(cell_code, mol_code, EF) %>% filter(nT == getmode(nT)) %>% distinct(.) %>% ungroup() %>% 
  mutate(sample = "M15_plus")

TL_cB <- TL_cB %>% group_by(sample, nT, EF, TC) %>% summarize(n = n())

write.csv(TL_cB, file = "Sample_cB.csv")
