### Doublet removal from single cell datasets
library(tidyverse)
library(Seurat)
library(DoubletFinder)

Sample_total <- read.csv(file = "/path_to/sample_total.csv")
Sample_nascent <- read.csv(file = "/path_to/sample_nascent.csv")

cellct_total <- Sample_total %>% filter(!grepl('ambig', EF)) %>% filter(!grepl("feature", EF)) %>% dplyr::select(-EF, -X, -NA.) %>% colSums(na.rm = TRUE)
cellct_total <- data.frame(cellct_total)
cellct_total$cellID <- rownames(cellct_total)


Sample_total_filt <- Sample_total %>% filter(!grepl('ambig', EF)) %>% dplyr::select(cellct_total$cellID) %>% replace(is.na(.), 0)
seu_total <- CreateSeuratObject(Sample_total_filt)
seu_total <- NormalizeData(seu_total)
seu_total <- FindVariableFeatures(seu_total, selection.method = "vst", nfeatures = 2000)
seu_total <- ScaleData(seu_total)
seu_total <- RunPCA(seu_total)
seu_total <- RunUMAP(seu_total, dims = 1:10)
sweep.res.list_total <- paramSweep_v3(seu_total, PCs = 1:10, sct = FALSE)
sweep.stats_total <- summarizeSweep(sweep.res.list_total, GT = FALSE)
bcmvn_total <- find.pK(sweep.stats_total)

### performing using defaults
annotations <- seu_total@meta.data$ClusteringResults
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(seu_total@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

seu_total <- doubletFinder_v3(seu_total, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
seu_res_total_DF <- data.frame(seu_total$DF.classifications)
seu_res_total_DF_DR <- seu_res_total_DF %>% filter(seu_total.DF.classifications == "Singlet")
cellfilt_total <- rownames(seu_res_total_DF_DR)

write.csv(cellfilt_total, file = "Total_singlet_cell_IDs.csv")

Total_doublet_removed <- Sample_total %>% dplyr::select(EF, cellfilt_total) %>% replace(is.na(.), 0)
Nascent_doublet_removed <- Sample_nacent %>% dplyr::select(EF, cellfilt_total) %>% replace(is.na(.), 0)
