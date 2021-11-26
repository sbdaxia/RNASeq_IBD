# create folder for revision pdf figures
dir.create("PDF_for_revision", showWarnings = F)

## Figure 2a, heatmap across genes from two models
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(mosaic)
library(org.Mm.eg.db)
library(AnnotationDbi)

tct_inf <- read_csv("../TCT/Significant genes_I_V vs NI_V.csv")
tnf_inf <- read_csv("../TNF/Significant genes_I_V vs NI_V.csv")

overlap <- intersect(tct_inf$X1, tnf_inf$X1)

tct_inf_overlap <- tct_inf %>% filter(X1 %in% overlap)
tnf_inf_overlap <- tnf_inf %>% filter(X1 %in% overlap)


overlap_matrix <- cbind(tct_inf_overlap[8:18], tnf_inf_overlap[8:17])
colnames(overlap_matrix) <- c(paste0("TCT_", colnames(overlap_matrix)[1:11]), paste0("TNF_", colnames(overlap_matrix)[12:21]))
my_palette <- colorRampPalette(c("blue", "white", "red"))(128)
heatmap_matrix <- as.matrix(overlap_matrix)
rownames(heatmap_matrix) <- tct_inf_overlap$X1

pdf('PDF_for_revision/Inflammation_overlap_heatmap.pdf',
    width = 22,
    height = 22)
out <- heatmap.2(heatmap_matrix,
          main = "Overlapped inflammation genes\npadj < 0.1 LFC > 0.25",
          density.info = "none",
          key = TRUE,
          lwid = c(2,7),
          lhei = c(1,8),
          col=my_palette,
          margins = c(8,10),
          symbreaks = TRUE,
          trace = "none",
          dendrogram = "row",
          labRow = rownames(heatmap_matrix),
          cexCol = 1,
          Colv = "NA")
dev.off()

# Overlap order
overlap_roworder <- rownames(heatmap_matrix)[out$rowInd]


merged <- unique(c(tct_inf$X1, tnf_inf$X1))
tct_vst <- read_csv("../TCT/VST_I_V vs NI_V.csv")
tnf_vst <- read_csv("../TNF/VST_I_V vs NI_V.csv")
tct_vst_merged <- tct_vst %>% filter(X1 %in% merged)
tnf_vst_merged <- tnf_vst %>% filter(X1 %in% merged)
merged_overlap <- intersect(tct_vst_merged$X1, tnf_vst_merged$X1)


# TCT specific DE gene order
tct_inf_unique <- tct_inf %>% filter(! X1 %in% overlap) %>% filter(X1 %in% merged_overlap)
heatmap_matrix <- as.matrix(tct_inf_unique[,8:18])
rownames(heatmap_matrix) <- tct_inf_unique$X1
pdf('PDF_for_revision/Inflammation_TCT_unique_heatmap.pdf',
    width = 22,
    height = 22)
out <- heatmap.2(heatmap_matrix,
                 main = "TCT unique inflammation genes\npadj < 0.1 LFC > 0.25",
                 density.info = "none",
                 key = TRUE,
                 lwid = c(2,7),
                 lhei = c(1,8),
                 col=my_palette,
                 margins = c(8,10),
                 symbreaks = TRUE,
                 trace = "none",
                 dendrogram = "row",
                 labRow = rownames(heatmap_matrix),
                 cexCol = 1,
                 Colv = "NA")
dev.off()

tct_unique_roworder <- rownames(heatmap_matrix)[out$rowInd]


# TNF specific DE gene order
tnf_inf_unique <- tnf_inf %>% filter(! X1 %in% overlap) %>% filter(X1 %in% merged_overlap)
heatmap_matrix <- as.matrix(tnf_inf_unique[,8:17])
rownames(heatmap_matrix) <- tnf_inf_unique$X1
pdf('PDF_for_revision/Inflammation_TNF_unique_heatmap.pdf',
    width = 22,
    height = 22)
out <- heatmap.2(heatmap_matrix,
                 main = "TNF unique inflammation genes\npadj < 0.1 LFC > 0.25",
                 density.info = "none",
                 key = TRUE,
                 lwid = c(2,7),
                 lhei = c(1,8),
                 col=my_palette,
                 margins = c(8,10),
                 symbreaks = TRUE,
                 trace = "none",
                 dendrogram = "row",
                 labRow = rownames(heatmap_matrix),
                 cexCol = 1,
                 Colv = "NA")
dev.off()

tnf_unique_roworder <- rownames(heatmap_matrix)[out$rowInd]


all_order <- c(overlap_roworder, tct_unique_roworder, tnf_unique_roworder)


tct_vst_merged <- tct_vst_merged %>% filter(X1 %in% merged_overlap) %>% as.data.frame()
tnf_vst_merged <- tnf_vst_merged %>% filter(X1 %in% merged_overlap) %>% as.data.frame()

for (i in 1:dim(tct_vst_merged)[1]) {
  tct_vst_merged[i,2:12] <- zscore(as.numeric(tct_vst_merged[i,2:12]))
}

for (i in 1:dim(tnf_vst_merged)[1]) {
  tnf_vst_merged[i,2:11] <- zscore(as.numeric(tnf_vst_merged[i,2:11]))
}

overlap_matrix <- cbind(tct_vst_merged[2:12], tnf_vst_merged[2:11])
colnames(overlap_matrix) <- c(paste0("TCT_", colnames(overlap_matrix)[1:11]), paste0("TNF_", colnames(overlap_matrix)[12:21]))
overlap_matrix$ID <- tnf_vst_merged$X1
all_order <- as.data.frame(all_order)
overlap_matrix_ordered <- left_join(all_order, overlap_matrix, by = c("all_order" = "ID"))
overlap_gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                           keys = all_order$all_order,
                                           columns = c("SYMBOL"),
                                           keytype = "ENSEMBL")

non_duplicates_idx <- which(duplicated(overlap_gene_name$ENSEMBL) == FALSE)
overlap_gene_name <- overlap_gene_name[non_duplicates_idx, ]

heatmap_matrix <- as.matrix(overlap_matrix_ordered[,-1])
rownames(heatmap_matrix) <- overlap_gene_name$SYMBOL

pdf('PDF_for_revision/Inflammation_merged_heatmap.pdf',
    width = 75,
    height = 150)
heatmap.2(heatmap_matrix,
          main = "Merged inflammation genes\npadj < 0.1 LFC > 0.25",
          density.info = "none",
          key = TRUE,
          lwid = c(2,7),
          lhei = c(1,8),
          col=my_palette,
          margins = c(8,10),
          symbreaks = TRUE,
          trace = "none",
          dendrogram = "none",
          labRow = rownames(heatmap_matrix),
          cexCol = 1,
          Colv = "NA",
          Rowv = FALSE)
dev.off()



# Generate DE table for supplement with I_V vs NI_V comparison
de_tct <- read_csv("../TCT/Differential Analysis_I_V vs NI_V.csv")
norm_tct <- read_csv("../TCT/normalized_raw_gene_counts_I_V vs NI_V.csv")
tb_tct <- left_join(de_tct, norm_tct, by = c("X1" = "X1"))
gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                           keys = tb_tct$X1,
                                           columns = c("SYMBOL"),
                                           keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]
tb_tct <- left_join(tb_tct, gene_name, by = c("X1" = "ENSEMBL"))
write_csv(tb_tct, "PDF_for_revision/TCT_IV_NIV.csv")

de_tnf <- read_csv("../TNF/Differential Analysis_I_V vs NI_V.csv")
norm_tnf <- read_csv("../TNF/normalized_raw_gene_counts_I_V vs NI_V.csv")
tb_tnf <- left_join(de_tnf, norm_tnf, by = c("X1" = "X1"))
gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = tb_tnf$X1,
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]
tb_tnf <- left_join(tb_tnf, gene_name, by = c("X1" = "ENSEMBL"))
write_csv(tb_tnf, "PDF_for_revision/TNF_IV_NIV.csv")

# Generate DE table for supplement with I_Mk2i vs I_V comparison
de_tct <- read_csv("../TCT/Differential Analysis_I_MKI vs I_V.csv")
norm_tct <- read_csv("../TCT/normalized_raw_gene_counts_I_MKI vs I_V.csv")
tb_tct <- left_join(de_tct, norm_tct, by = c("X1" = "X1"))
gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = tb_tct$X1,
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]
tb_tct <- left_join(tb_tct, gene_name, by = c("X1" = "ENSEMBL"))
write_csv(tb_tct, "PDF_for_revision/TCT_IMKI_IV.csv")

de_tnf <- read_csv("../TNF/Differential Analysis_I_MKI vs I_V.csv")
norm_tnf <- read_csv("../TNF/normalized_raw_gene_counts_I_MKI vs I_V.csv")
tb_tnf <- left_join(de_tnf, norm_tnf, by = c("X1" = "X1"))
gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = tb_tnf$X1,
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]
tb_tnf <- left_join(tb_tnf, gene_name, by = c("X1" = "ENSEMBL"))
write_csv(tb_tnf, "PDF_for_revision/TNF_IMKI_IV.csv")


# NI_V vs I_VS vs I_Mki heatmap
## TCT
# rearrange the TCT master count matrix
rawCountTable_transform_detected_TCT <- read_csv("../TCT/VST_TCT.csv")
row_names <- rawCountTable_transform_detected_TCT$X1
rawCountTable_transform_detected_TCT$X1 <- NULL
rawCountTable_transform_detected_reshape <- cbind(rawCountTable_transform_detected_TCT[,19:23], rawCountTable_transform_detected_TCT[,7:12], rawCountTable_transform_detected_TCT[,1:6])
rownames(rawCountTable_transform_detected_reshape) <- row_names

# plot the heatmap
dif_analysis <- read_csv("../TCT/Differential Analysis_I_MKI vs I_V.csv")
sig_dif <- subset(dif_analysis, dif_analysis$padj < 0.1 & abs(dif_analysis$log2FoldChange) > 0.25)
sig_count <- rawCountTable_transform_detected_reshape[rownames(rawCountTable_transform_detected_reshape) %in% sig_dif$X1,]
sig_dif <- cbind(sig_dif, sig_count)
sig_dif$X1 <- NULL
for (i in 1:dim(sig_dif)[1]) {
  sig_dif[i,7:23] <- zscore(as.numeric(sig_dif[i,7:23]))
}

gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = rownames(sig_dif),
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]

heatmap_matrix <- as.matrix(sig_dif[,7:23])
rownames(heatmap_matrix) <- gene_name$SYMBOL
my_palette <- colorRampPalette(c("blue", "white", "red"))(128)

pdf('PDF_for_revision/TCT MK2i targets in TCT NI_V vs I_V vs I_MK2i RNASeq.pdf',
    width = 9,
    height = 11)
heatmap.2(heatmap_matrix,
          main = "TCT MK2i targets\nTCT NI_V vs I_V vs I_MK2i\npadj < 0.1 LFC > 0.25",
          density.info = "none",
          key = TRUE,
          lhei = c(1,7),
          lwid = c(1,5),
          col=my_palette,
          cexCol = 1,
          margins = c(8,8),
          trace = "none",
          dendrogram = "row",
          keysize = 2,
          labRow = rownames(heatmap_matrix),
          Colv = "NA")
dev.off()

## TNF
# rearrange the TCT master count matrix
rawCountTable_transform_detected_TNF <- read_csv("../TNF/VST_TNF.csv")
row_names <- rawCountTable_transform_detected_TNF$X1
rawCountTable_transform_detected_TNF$X1 <- NULL
rawCountTable_transform_detected_reshape <- cbind(rawCountTable_transform_detected_TNF[,17:21], rawCountTable_transform_detected_TNF[,7:11], rawCountTable_transform_detected_TNF[,1:6])
rownames(rawCountTable_transform_detected_reshape) <- row_names

# plot the heatmap
dif_analysis <- read_csv("../TNF/Differential Analysis_I_MKI vs I_V.csv")
sig_dif <- subset(dif_analysis, dif_analysis$padj < 0.1 & abs(dif_analysis$log2FoldChange) > 0.25)
sig_count <- rawCountTable_transform_detected_reshape[rownames(rawCountTable_transform_detected_reshape) %in% sig_dif$X1,]
sig_dif <- cbind(sig_dif, sig_count)
sig_dif$X1 <- NULL
for (i in 1:dim(sig_dif)[1]) {
  sig_dif[i,7:22] <- zscore(as.numeric(sig_dif[i,7:22]))
}

gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = rownames(sig_dif),
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]

heatmap_matrix <- as.matrix(sig_dif[,7:22])
rownames(heatmap_matrix) <- gene_name$SYMBOL

pdf('PDF_for_revision/TNF MK2i targets in TNF NI_V vs I_V vs I_MK2i RNASeq.pdf',
    width = 9,
    height = 11)
heatmap.2(heatmap_matrix,
          main = "TNF MK2i targets\nTNF NI_V vs I_V vs I_MK2i\npadj < 0.1 LFC > 0.25",
          density.info = "none",
          key = TRUE,
          lhei = c(1,7),
          lwid = c(1,5),
          col=my_palette,
          cexCol = 1,
          margins = c(8,8),
          trace = "none",
          dendrogram = "row",
          keysize = 2,
          labRow = rownames(heatmap_matrix),
          Colv = "NA")
dev.off()


# Generate csv file ccontaining unified DEGs from TCT and TNF inflammation
de_tct <- read_csv("../TCT/Differential Analysis_I_V vs NI_V.csv")
de_tnf <- read_csv("../TNF/Differential Analysis_I_V vs NI_V.csv")
sig_tct <- de_tct %>% filter(padj < 0.05) %>% dplyr::select(X1) %>% unlist()
sig_tnf <- de_tnf %>% filter(padj < 0.05) %>% dplyr::select(X1) %>% unlist()
sig_uni <- unique(c(sig_tct, sig_tnf))
uni_tct <- de_tct %>% filter(X1 %in% sig_uni)
uni_tnf <- de_tnf %>% filter(X1 %in% sig_uni)
colnames(uni_tct) <- paste0("TCT_", colnames(uni_tct))
colnames(uni_tnf) <- paste0("TNF_", colnames(uni_tnf))

uni_tbl <- inner_join(uni_tct, uni_tnf, by = c("TCT_X1" = "TNF_X1"))

gene_name <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = uni_tbl$TCT_X1,
                                   columns = c("SYMBOL"),
                                   keytype = "ENSEMBL")
non_duplicates_idx <- which(duplicated(gene_name$ENSEMBL) == FALSE)
gene_name <- gene_name[non_duplicates_idx, ]

uni_tbl <- left_join(uni_tbl, gene_name, by = c("TCT_X1"="ENSEMBL"))

uni_tbl <- as.data.frame(uni_tbl)
rownames(uni_tbl) <- uni_tbl$TCT_X1
uni_tbl <- uni_tbl[,-1]

uni_tbl$SameDirection <- FALSE

uni_tbl$SameDirection[(uni_tbl$TCT_log2FoldChange * uni_tbl$TNF_log2FoldChange) > 0] <- TRUE


write.csv(uni_tbl, "PDF_for_revision/Uni_DEGs.csv")

sum(uni_tbl$SameDirection)
