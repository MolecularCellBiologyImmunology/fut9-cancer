
# Install Packages if needed and load them (R version 3.6.2)
PACKAGES_TOTAL <- c("readr", "dplyr", "tidyr", "pheatmap", "tibble", "readxl", "tibble", "xlsx", "ggplot2", "gplots", "rstudioapi", "Seurat", "matrixStats")
installation_needed <- unlist(lapply(PACKAGES_TOTAL, require, character.only = TRUE, quietly = TRUE))
installation_needed <- PACKAGES_TOTAL[installation_needed == FALSE]
if(length(installation_needed) > 0) {
  
  ## From CRAN
  chooseCRANmirror(ind = 8)
  for (package in installation_needed) {if (package %in% PACKAGES_TOTAL) {try({install.packages(pkgs = package, repos='http://cran.us.r-project.org', dependencies = TRUE)})}}
  installation_needed <- unlist(lapply(PACKAGES_TOTAL, require, character.only = TRUE, quietly = TRUE))
  installation_needed <- PACKAGES_TOTAL[installation_needed == FALSE]
  
  ## From Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
  for (package in installation_needed) {if (package %in% PACKAGES_TOTAL) {try({BiocManager::install(package, update = FALSE)})}}
  installation_needed <- unlist(lapply(PACKAGES_TOTAL, require, character.only = TRUE, quietly = TRUE))
  installation_needed <- PACKAGES_TOTAL[installation_needed == FALSE]
  if(length(installation_needed) > 0) {stop("Not All Packages could be installed and loaded. Check what version of R you are using or contact the corresponding author.")}
}
rm(PACKAGES_TOTAL, installation_needed)



# Read Data and add cell type (Normal / Tumor) to cell names
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
normal <- read_delim(file = "GSE81861_CRC_NM_all_cells_COUNT.csv", delim = ",", col_names = TRUE)
colnames(normal) <- gsub(colnames(normal), pattern = "__", replacement = "_Normal_")
tumor <- read_delim(file = "GSE81861_CRC_tumor_all_cells_COUNT.csv", delim = ",", col_names = TRUE)
colnames(tumor) <- gsub(colnames(tumor), pattern = "__", replacement = "_Tumor_")
data <- full_join(x = normal, y = tumor, by = "X1")
rm(normal, tumor)



# Clean Up Gene Names
genes <- data$X1
genes_clean <- list()
for (gene in genes) {genes_clean <- c(genes_clean, strsplit(x = gene, split = "_")[[1]][2])}
genes_clean <- make.unique(names = unlist(genes_clean))
data$X1 <- genes_clean
data <- column_to_rownames(data, "X1")
rm(gene, genes, genes_clean)



# Clean up Cell names
cells <- colnames(data)
celltypes <- list()
for (cell in cells) {split <- strsplit(x = cell, split = "_")[[1]]; celltypes <- c(celltypes, paste(split[2], split[3]))}
celltypes <- unlist(celltypes)
celltypes <- as.factor(celltypes)
rm(cell, cells, split)



# Even though it seems this data is after their filtering steps,
# Reproduce Filtering Steps from the paper just to be certain
data_logical <- data > 0
FODG <- colSums(data_logical)
names(FODG) <- colnames(data_logical)
sort(FODG)
FODG_to_filter <- FODG < 1000
data_filtered <- data[,names(FODG_to_filter[FODG_to_filter == FALSE])]
TOTAL <- colSums(data_filtered)
sort(TOTAL)
TOTAL_to_filter <- TOTAL < 100000
data_filtered <- data_filtered[,names(TOTAL_to_filter[TOTAL_to_filter == FALSE])]
rm(data_logical, FODG, FODG_to_filter, TOTAL, TOTAL_to_filter)



# Identify FUT9 Expressing Cells
FUT9 <- as.data.frame(t(data_filtered["FUT9",]))
FUT9$Celltype <- celltypes
FUT9 <- FUT9[order(-FUT9$FUT9),]
FUT9 <- FUT9[FUT9$FUT9 > 0,] 
celltypestable <- subset(as.data.frame(table(FUT9)), Freq > 0)
write.xlsx(x = celltypestable, file = "FUT9 Expression per Celltype.xlsx", col.names = TRUE, row.names = FALSE)
FUT9 <- FUT9[(FUT9$Celltype == "Normal Epithelial") | (FUT9$Celltype == "Tumor Epithelial"),]
FUT9$Celltype <- droplevels(FUT9$Celltype)
barplot(height = FUT9$FUT9, col = FUT9$Celltype)
legend("topright", legend = levels(FUT9$Celltype), col=1:nlevels(FUT9$Celltype), pch = 16)
p <- ggplot(FUT9, aes(factor(Celltype), FUT9)) 
p <- p + geom_boxplot() + facet_wrap(~Celltype, scale="free")
p + geom_jitter(shape=16, position=position_jitter(0.2))       



# Group Cells by Normal/Tumor and FUT9+/-, average all expression for the heatmap
FUT9NEGATIVE <- as.data.frame(t(data_filtered["FUT9",]))
FUT9NEGATIVE$Celltype <- celltypes
FUT9NEGATIVE <- FUT9NEGATIVE[FUT9NEGATIVE$FUT9 == 0,]
NFP <- rownames(FUT9[FUT9$Celltype == "Normal Epithelial",])
TFP <- rownames(FUT9[FUT9$Celltype == "Tumor Epithelial",])
NFN <- rownames(FUT9NEGATIVE[FUT9NEGATIVE$Celltype == "Normal Epithelial",])
TFN <- rownames(FUT9NEGATIVE[FUT9NEGATIVE$Celltype == "Tumor Epithelial",])
NFP <- rowMeans(data_filtered[,NFP])
TFP <- rowMeans(data_filtered[,TFP])
NFN <- rowMeans(data_filtered[,NFN])
TFN <- rowMeans(data_filtered[,TFN])
AverageSampleExpressions <- data.frame(NFN, NFP, TFN, TFP)
rm(NFN, NFP, TFN, TFP, FUT9, FUT9NEGATIVE)



# Cluster Data
clusterdata <- data_filtered[,grep(x = celltypes, pattern = "Epithelial")]
metadata <- data.frame(row.names = colnames(clusterdata), 
                       "Cell Type" = celltypes[grep(x = celltypes, pattern = "Epithelial")], 
                       "FUT9" = factor(as.logical(clusterdata["FUT9",] > 0), labels = c('Negative', 'Positive') ))
metadata <- tidyr::unite(metadata, col = "metadata", sep = " - ", remove = FALSE)
SeuratObject <- CreateSeuratObject(counts = as.matrix(clusterdata), project = "Courtois et al., FUT9", meta.data = metadata)
SeuratObject <- NormalizeData(object = SeuratObject)
SeuratObject <- FindVariableFeatures(object = SeuratObject)
SeuratObject <- ScaleData(object = SeuratObject)
SeuratObject <- RunPCA(object = SeuratObject)
SeuratObject <- FindNeighbors(object = SeuratObject)
SeuratObject <- FindClusters(object = SeuratObject)
SeuratObject <- RunTSNE(object = SeuratObject)
DimPlot(object = SeuratObject, reduction = "tsne", pt.size = 1, cols = rainbow(9))
DimPlot(object = SeuratObject, reduction = "tsne", pt.size = 3, group.by = "metadata", shape.by = "FUT9", cols = c("lightblue","blue","pink","red"))
rm(clusterdata, metadata)



# Group Genes by relevant pathways, average all expression for the heatmap
setwd("../Genesets")
Genesets <- list()
Genesets_Miranda <- readxl::read_excel(path = "Miranda et al (2019) supplemental tables.xlsx", sheet = "S1A Curated stemness signature")
Genesets[["Stemcells_Curated"]] <- as.vector(na.omit(Genesets_Miranda$`Curated  (withouth immune and proliferative genes, annotated in TCGA)`))
Genesets_Files <- list.files(pattern = ".txt")
for (Geneset_File in Genesets_Files) {Genesets[[Geneset_File]] <- as.character(unlist(read.table(file = Geneset_File, skip = 2)))}
AverageSampleGeneExpressions <- data.frame()
DeviationSampleGeneExpressions <- data.frame()
for (Geneset in Genesets) {
  Data_Geneset <- AverageSampleExpressions[rownames(AverageSampleExpressions) %in% Geneset,]
  AverageSampleGeneExpressions <- rbind(AverageSampleGeneExpressions, t(colMeans(Data_Geneset)))
  DeviationSampleGeneExpressions <- rbind(DeviationSampleGeneExpressions, t(colSds(as.matrix(Data_Geneset))))
}
rownames(AverageSampleGeneExpressions) <- unlist(lapply(X = sub(pattern = ".txt", x = names(Genesets), replacement = ""), FUN = tolower))
rownames(DeviationSampleGeneExpressions) <- unlist(lapply(X = sub(pattern = ".txt", x = names(Genesets), replacement = ""), FUN = tolower))
rm(Genesets_Miranda, Data_Geneset, Geneset, Genesets, Geneset_File, Genesets_Files)



# Make Final Heatmap
colnames(AverageSampleGeneExpressions) <- c("Healthy FUT9 Negative (N = 198)", "Healthy FUT9 Positive (N = 5)", "Tumor FUT9 Negative (N = 261)", "Tumor FUT9 Positive (N = 11)")
colnames(DeviationSampleGeneExpressions) <- c("Healthy FUT9 Negative (N = 198)", "Healthy FUT9 Positive (N = 5)", "Tumor FUT9 Negative (N = 261)", "Tumor FUT9 Positive (N = 11)")
AverageSampleGeneExpressions <- AverageSampleGeneExpressions[c("stemcells_curated","affar_yy1_targets_up", "pid_myc_pathway", "hsf1_01", "benporath_sox2_targets", "kegg_hedgehog_signaling_pathway", "go_notch_signaling_pathway", "go_canonical_wnt_signaling_pathway"),]
DeviationSampleGeneExpressions <- DeviationSampleGeneExpressions[c("stemcells_curated","affar_yy1_targets_up", "pid_myc_pathway", "hsf1_01", "benporath_sox2_targets", "kegg_hedgehog_signaling_pathway", "go_notch_signaling_pathway", "go_canonical_wnt_signaling_pathway"),]
rownames(AverageSampleGeneExpressions) <- c("stemcells_curated","affar_yy1_targets_up", "pid_myc_pathway", "hsf1_01", "benporath_sox2_targets", "kegg_hedgehog_signaling_pathway", "go_notch_signaling_pathway", "go_canonical_wnt_signaling_pathway")
rownames(DeviationSampleGeneExpressions) <- c("stemcells_curated","affar_yy1_targets_up", "pid_myc_pathway", "hsf1_01", "benporath_sox2_targets", "kegg_hedgehog_signaling_pathway", "go_notch_signaling_pathway", "go_canonical_wnt_signaling_pathway")
pheatmap(mat = AverageSampleGeneExpressions, main = "Average FUT9 Expression", kmeans_k = NA, color = bluered(100), scale = "row", angle_col = 315, fontsize_col = 11, cluster_cols = FALSE, cluster_rows = FALSE)
pheatmap(mat = DeviationSampleGeneExpressions, main = "Standard Deviation of FUT9 Expression", kmeans_k = NA, color = bluered(100), scale = "row", angle_col = 315, fontsize_col = 11, cluster_cols = FALSE, cluster_rows = FALSE)
