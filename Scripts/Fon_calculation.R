## Fon calculation on full dataset
library(tidyverse)

## Bring in pre-processed TL_filt.csv
TL_filt <- read.csv(file = "TL_filt.csv")

## Bring in cell IDs from doublet filtering analysis
cell_IDs <- read.csv(file = "Total_singlet_cell_IDs.csv")

## mode function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

## collapse into Tsum DF for Fon calculations
TL_collapsed_Tsum <- TL_filt %>% filter(EF %in% cell_IDs) %>% dplyr::select(cell_code, mol_code, nT, TC, EF) %>%
  group_by(cell_code, mol_code, EF) %>% filter(nT == getmode(nT)) %>% distinct(.) %>% ungroup() %>%
  group_by(cell_code, EF) %>% summarise(totN = sum(nT), totMut = sum(TC))
head(TL_collapsed_Tsum)
write.csv(TL_collapsed_Tsum, file = "Sample_Tsums.csv")

TL_collapsed_ratios <- TL_collapsed_Tsum %>% filter(!grepl("ambig", EF)) %>% filter(!grepl("feature", EF)) %>%
  mutate(mut_ratio = totMut/totN) %>% filter(mut_ratio > 0.006636853) %>% group_by(EF) %>% summarise(n())
colnames(TL_collapsed_ratios) <- c("gene", "On_cells")
