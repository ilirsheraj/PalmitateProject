library(ensembldb)
library(EnsDb.Mmusculus.v79)
library(tximport)
library(readr)
library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(tidyverse)

my_dir <- getwd()

files <- list.files(path=paste0(my_dir, "/Salmon"), 
                    pattern="quant.sf$", recursive=TRUE, full.names = TRUE, include.dirs = TRUE)
files

# Create the datatable
mydata <- data.frame(sample = c("S46551_CTR1", 
                                "S46552_CTR2", 
                                "S46553_CTR3",
                                "S46554_PALM1",
                                "S46555_PALM2",
                                "S46556_PALM3"),
                     condition = as.factor(c("Control",
                                             "Control",
                                             "Control",
                                             "Palmitate",
                                             "Palmitate",
                                             "Palmitate")),
                     short_let = as.factor(rep(c("CTR", "PA"), each=3)))

mydata <- cbind(mydata, files = files)
mydata

# Create Output Directories
# Output directories
out_dir   <- file.path(getwd(), "DiffExp")
plot_dir  <- file.path(getwd(), "Plots")
dir.create(out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
#############################################################################
# TXImport
edb <- EnsDb.Mmusculus.v79
tx2gene <- transcripts(edb, return.type="DataFrame", columns=c("tx_id", "gene_id", "gene_name", "entrezid"))

txi <- tximport(files,type ="salmon", tx2gene = tx2gene, ignoreTxVersion=TRUE)

dds <- DESeqDataSetFromTximport(txi, colData = mydata, design = ~condition)

dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 3
dds <- dds[filter,]
vsd <- rlog(dds, blind = TRUE)


# Plot PCA
pcaData <- plotPCA(vsd, intgroup= "condition", returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf(file = file.path(plot_dir, "pca_RNASeq.pdf"), width = 3, height = 3)
ggplot(pcaData, aes(PC1, PC2, color=short_let)) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() + 
  theme(legend.position = c(0.12, 0.9),
        legend.title=element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_point(size=4) +
  # geom_text(aes(label=mydata$sample), vjust=-1, size=3) +
  # xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  # ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  theme(panel.background = element_rect(colour = "black", size=0.8, fill=NA))
dev.off()

# Get the raw data
raw_counts <- counts(dds)
colnames(raw_counts) <- mydata$sample
write.csv(raw_counts, file = file.path(out_dir, "rnaseq_raw_counts.csv"))

################################################################################
# Differential Expression
dds <- DESeq(dds)
diff_exp <- results(dds, contrast=c("condition", "Palmitate", "Control"))
summary(diff_exp, alpha = 0.05)
# Palmitate vs Control

# Apply fold change shrinkage for better visualization
library(apeglm)
res <- lfcShrink(dds, coef="condition_Palmitate_vs_Control", type="ashr", res = diff_exp)
summary(res, alpha= 0.05)

# Annotate the genes
# Convert ensemble IDs into symbols
ensembl_genenames = rownames(diff_exp)
mouse_genenames <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ensembl_genenames, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
diff_exp$symbol = cbind(ensembl_genenames, mouse_genenames[match(ensembl_genenames,mouse_genenames$GENEID),1])[,2]
write.csv(diff_exp, file = file.path(out_dir, "palmitate_vs_control_diff_expression_annotated.csv"))

# Do the same for log-shrinked data
ensembl_genenames = rownames(res)
mouse_genenames <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ensembl_genenames, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
res$symbol = cbind(ensembl_genenames, mouse_genenames[match(ensembl_genenames,mouse_genenames$GENEID),1])[,2]
write.csv(res, file = file.path(out_dir, "palmitate_vs_control_lfc_diff_expression_annotated.csv"))

##################################################################################
# GO Term Enrichment
diff_exp <- data.frame(diff_exp)
rnaseq_up <- diff_exp %>% dplyr::filter(log2FoldChange > 0.5 & padj < 0.05)
rnaseq_down <- diff_exp %>% dplyr::filter(log2FoldChange < -0.5 & padj < 0.05)

# Call "conventient_functions.R" in R
go_terms <- go_enrich_func(diff_exp, adjp=0.05, lfc=0.5, "deseq_lfc05", trim_val=0.7)

# Plot
pdf(file = file.path(plot_dir, "Deseq_GO_Terms.pdf"))
for (i in 1:2) {
  print(go_terms[[i]])
}
dev.off()

# KEGG Pathway
deseq_kegg <- kegg_enrich_func(diff_exp, adjp=0.05, lfc=0.5, "trimgalore_deseq_lfc1")

pdf(file = file.path(plot_dir, "Deseq_KEGG.pdf"))
for (i in 1:2) {
  print(deseq_kegg[[i]])
}
dev.off()

# Reactome Pathways
deseq_reactome <- reactome_enrich_func(diff_exp, adjp=0.05, lfc=0.5, "trimgalore_deseq_lfc1")

pdf(file = file.path(plot_dir, "Deseq_reactome.pdf"))
for (i in 1:2) {
  print(deseq_reactome[[i]])
}
dev.off()

################################################################################
# Select Genes for Volcano Plot
goterm_genes <- read.csv("Data/genes_selected_from_goterms_PA.csv")
head(goterm_genes)

# Convert them into Esnemble (Works better than others)
library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Gene Symbols
gene_symbols <- goterm_genes$gene_name

# Retrieve Ensembl IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "mgi_symbol",
  values = gene_symbols,
  mart = mart
)

head(results)

merged_data <- merge(goterm_genes, results, by.x = "gene_name", by.y = "mgi_symbol", all.x = TRUE)
# write.csv(merged_data, file = "Data/annotated_genes_from_go.csv")

# After some modification, read it back again and use for Annotation
goterm_genes <- read.csv("Data/annotated_genes_from_go.csv")
head(goterm_genes)

# Unify the symbols with the rest of the data
ensembl_genenames = goterm_genes$ensembl_gene_id
mouse_genenames <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ensembl_genenames, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
goterm_genes$symbol = cbind(ensembl_genenames, mouse_genenames[match(ensembl_genenames,mouse_genenames$GENEID),1])[,2]
write.csv(goterm_genes, file = "Data/go_term_genes_volcano_annotated.csv")

# Convert logFC shrank into data frame
diff_exp_annotated <- data.frame(res)
goterm_genes <- goterm_genes[,c(4,2)]
diff_exp_annotated <- left_join(diff_exp_annotated, goterm_genes, by = "symbol")
diff_exp_annotated <- diff_exp_annotated[order(diff_exp_annotated[,"GO_Term"],na.last=F),]
diff_exp_annotated <- diff_exp_annotated %>% dplyr::filter(abs(log2FoldChange) < 6)

library(ggnewscale)
bs <- diff_exp_annotated
bs <- bs %>% mutate(GO_Term = ifelse(is.na(GO_Term), 0, GO_Term))
bs$GO_Term[bs$GO_Term == 0 & bs$padj < 0.05 & abs(bs$log2FoldChange) > 0.5] <- "SIGNIFICANT"
bs$GO_Term[bs$GO_Term == 0] <- "NONSIGNIFICANT"

ggplot(bs, aes(x=log2FoldChange, y=-log10(padj), size = ifelse(abs(log2FoldChange) > 0.5, abs(log2FoldChange), 0.1))) +
  geom_point(data = subset(bs, GO_Term %in% c("SIGNIFICANT", "NONSIGNIFICANT")), aes(color = GO_Term)) +
  scale_colour_manual(values = c("lightblue", "gray90"),
                      limits = c("SIGNIFICANT", "NONSIGNIFICANT")) +
  guides(color=guide_legend(title = "SIGNIFICANCE", override.aes = list(size = 6, face = "bold"))) +
  new_scale_color() +
  geom_point(data = subset(bs, GO_Term %in% c("APOPTOSIS", "CARBOXYLIC ACID", "CELL CYCLE", "LIPID CATABOLISM", 
                                              "MITOCHONDRIA", "NEGATIVE REGULATION OF T CELL ACTIVATION")), aes(color = GO_Term)) +
  scale_color_manual(values = c("#990099", "#303030", "#006000", "#7e0021", "#ffd320", "blue"),
                     limits = c("APOPTOSIS", "CARBOXYLIC ACID", "CELL CYCLE", "LIPID CATABOLISM", 
                                "MITOCHONDRIA", "NEGATIVE REGULATION OF T CELL ACTIVATION")) +
  annotate("text", x=2, y=25, label= "UP\n 849", size=5) +
  annotate("text", x=-2, y=25, label= "Down\n 468", size = 5) +
  theme_light() +
  ggtitle("Palmitate vs Control") +
  xlab("log2 fold change") + 
  ylab("-log10(padj)") +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(1.5), hjust = 0.5, face = "bold"),
        axis.title = element_text(size = rel(1.25))) +
  guides(color=guide_legend(title = "FUNCTION", override.aes = list(size = 6, face = "bold"))) +
  guides(size = FALSE) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 10, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_vline(xintercept = 0.5, color = "grey", size=0.7, linetype = "dashed") +
  geom_vline(xintercept = -0.5, color = "grey", size=0.7, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "grey", size = 0.7, linetype = "dashed") +
  xlim(c(-3, 3))

ggsave(file.path(plot_dir, "RNASeq_Volcano_final.pdf"), width = 10, height = 8, dpi = 600)

##################################################################################
# Heatmaps: First from GO Terms
rld <- data.frame(assay(rlog(dds, blind=FALSE)))
mydata$name = as.factor(c("CTR1", "CTR2", "CTR3", "PALM1", "PALM2", "PALM3"))
colnames(rld) <- mydata$name

go_all <- read.csv("Data/go_term_genes_volcano_annotated.csv")

go_all <- go_all |> dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

all_genes_rlog <- rld %>% dplyr::filter(rownames(.) %in% go_all$ensembl_gene_id)

library(pheatmap)
my_gene_col <- go_all %>% dplyr::filter(ensembl_gene_id %in% rownames(all_genes_rlog)) %>% dplyr::select(ensembl_gene_id, GO_Term)
rownames(my_gene_col) <- my_gene_col$ensembl_gene_id
colnames(my_gene_col)[2] <- "GO"
my_gene_col <- my_gene_col %>% dplyr::select(GO)

# my_gene_col$GO <- gsub(
#   "NEGATIVE REGULATION OF T CELL ACTIVATION",
#   "NEGATIVE REGULATION\nOF T CELL ACTIVATION",
#   my_gene_col$GO
# )

annoCol<-list(GO=c(APOPTOSIS="#990099", `LIPID CATABOLISM`="#7e0021", 
                   `CARBOXYLIC ACID`="#303030", `CELL CYCLE` = "#006000",
                   MITOCHONDRIA =  "#ffd320",
                   `NEGATIVE REGULATION OF T CELL ACTIVATION`="blue"))

breaks <- seq(-1.5, 1.5, length.out = 100)
pdf(file = file.path(plot_dir, "heatmap_for_goterm.pdf"), width = 8, height = 6)
a <- pheatmap::pheatmap(all_genes_rlog, scale="row", 
                        clustering_distance_rows="correlation", 
                        clustering_distance_cols = "correlation", 
                        show_rownames =F, 
                        show_colnames =F, 
                        cluster_cols=T, 
                        fontsize_row=4, 
                        annotation_row = my_gene_col,
                        angle_col = 45, 
                        breaks = breaks,
                        annotation_colors = annoCol, 
                        color = colorRampPalette(c("blue", "white", "red"))(100))

print(a)
dev.off()
#######################################
# For a number of selected genes
library(ComplexHeatmap)
genes_heatmap <- read.csv("Data/gene_selected_goterms_heatmap.csv")
head(genes_heatmap)
gene_symbols <- genes_heatmap$gene_name

# Retrieve Ensembl IDs
results <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "mgi_symbol",
  values = gene_symbols,
  mart = mart
)

merged_data <- merge(genes_heatmap, results, by.x = "gene_name", by.y = "mgi_symbol", all.x = TRUE)
# write.csv(merged_data, file = "Data/heatmap_genes_annotated.csv")
# Some that were not annotated were manually checked and annotated accordingly

merged_data <- read.csv("Data/heatmap_genes_annotated.csv")
merged_data <- merged_data %>% dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)
head(merged_data)
# Unify the annotation
ensembl_genenames = merged_data$ensembl_gene_id
mouse_genenames <- ensembldb::select(EnsDb.Mmusculus.v79, keys=ensembl_genenames, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
merged_data$symbol = cbind(ensembl_genenames, mouse_genenames[match(ensembl_genenames,mouse_genenames$GENEID),1])[,2]
write.csv(merged_data, file = "Data/heatmap_genes_annotated_new.csv")

go_all <- merged_data

all_genes_rlog <- rld %>% dplyr::filter(rownames(.) %in% go_all$ensembl_gene_id)

rlog_new <- all_genes_rlog %>%  rownames_to_column("GENEID")
colnames(go_all)[3] <- "GENEID"

go_all_rlog <- inner_join(go_all, rlog_new, by = "GENEID")

gene_rlog <- go_all_rlog[,-c(1,2)]
rownames(gene_rlog) <- gene_rlog$GENEID
gene_rlog <- gene_rlog[,-1]

row_labels <- structure(go_all_rlog$gene_name, names = go_all_rlog$GENEID)
split_labels <- structure(go_all_rlog$GO_Term, names = go_all_rlog$GENEID)

genes_mat <- scale(t(all_genes_rlog))
all_genes_rlog <- data.frame(t(genes_mat))


pdf(file = file.path(plot_dir, "smaller_heatmap.pdf"), width = 3.5, height = 8)

Heatmap(as.matrix(all_genes_rlog), row_labels = row_labels[rownames(all_genes_rlog)], 
        cluster_rows = T, cluster_columns = T, row_split = split_labels[rownames(all_genes_rlog)],
        row_names_gp = gpar(col = c("#990099", "#303030", "#006000", "#7e0021", "#ffd200", "blue"), fontsize=9),
        row_title = NULL, row_gap = unit(1, "mm"), name = "z-score", #column_names_rot = 45,
        column_names_gp = gpar(fontsize = 12, col = "black"), cluster_row_slices = F,
        col = circlize::colorRamp2(c(-1.5, 0, 1.5), c("darkblue", "white", "red")),
        heatmap_legend_param = list(
          color_bar = "continuous"
        ))

dev.off()

# Do it with pheatmap too
go_all <- merged_data
my_gene_col <- go_all %>% dplyr::filter(ensembl_gene_id %in% rownames(all_genes_rlog)) %>% dplyr::select(ensembl_gene_id, GO_Term)
rownames(my_gene_col) <- my_gene_col$ensembl_gene_id
colnames(my_gene_col)[2] <- "GO"
my_gene_col <- my_gene_col %>% dplyr::select(GO)

# my_gene_col$GO <- gsub(
#   "NEGATIVE REGULATION OF T CELL ACTIVATION",
#   "NEGATIVE REGULATION\nOF T CELL ACTIVATION",
#   my_gene_col$GO
# )

annoCol<-list(GO=c(APOPTOSIS="#990099", `LIPID CATABOLISM`="#7e0021", 
                   `CARBOXYLIC ACID`="#303030", `CELL CYCLE` = "#006000",
                   MITOCHONDRIA =  "#ffd320",
                   `NEGATIVE REGULATION OF T CELL ACTIVATION`="blue"))

breaks <- seq(-1.5, 1.5, length.out = 100)
pdf(file = file.path(plot_dir, "smaller_heatmap_for_goterm.pdf"), width = 8, height = 6)
a <- pheatmap::pheatmap(all_genes_rlog, #scale="row", 
                        clustering_distance_rows="correlation", 
                        clustering_distance_cols = "correlation", 
                        show_rownames =F, 
                        show_colnames =F, 
                        cluster_cols=T, 
                        fontsize_row=4, 
                        annotation_row = my_gene_col,
                        angle_col = 45, 
                        breaks = breaks,
                        annotation_colors = annoCol, 
                        color = colorRampPalette(c("blue", "white", "red"))(100))

print(a)
dev.off()
####################################################################################
# Dot Plots for Commonly Up and Downregulted Genes in ATACSeq and RNASeq
# ATACSeq CSAW Differential Accessibility
all_peaks <- read.csv("ATACSeq/csaw_macs2_chipseeker_annotated_batch_peaks.csv")

atac_up <- all_peaks %>% dplyr::filter(logFC > 0.5 & FDR < 0.1)
atac_up_uq <- atac_up %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)

atac_down <- all_peaks %>% dplyr::filter(logFC < -0.5 & FDR < 0.1)
atac_down_uq <- atac_down %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>% na.omit(ENSEMBL)

library(ggVennDiagram)

# List of items
x <- list(RNASeq_Down = rownames(rnaseq_down), 
          ATAC_Down = atac_down_uq$ENSEMBL)

# 2D Venn diagram
pdf(file = file.path(plot_dir, "rsaseq_atacseq_down_venn.pdf"), width = 3, height = 2.5)
a <- ggVennDiagram(x, label = "count",  set_color = c("red","blue"),
                   category.names = c("RNASeq", "ATACSeq")) + 
  scale_color_manual(values = c("red","blue")) + 
  coord_flip() +
  theme(legend.position = "none")

print(a)
dev.off()

# Genes Upregulated
xx <- list(RNASeq_Up = rownames(rnaseq_up), 
           ATAC_Up = atac_up_uq$ENSEMBL)

pdf(file = file.path(plot_dir, "rsaseq_atacseq_up_venn.pdf"), width = 3, height = 2.5)
a <- ggVennDiagram(xx, label = "count",  set_color = c("red","blue"), lwd = 0.8, lty = 1, 
                   category.names = c("RNASeq", "ATACSeq")) + 
  scale_color_manual(values = c("red","blue")) + 
  coord_flip() +
  theme(legend.position = "none") # + scale_fill_gradient(low = "deepskyblue", high = "pink")

print(a)
dev.off()

# Pick Promoters only to find the intersection of genes and then pick all peaks for those genes
atacseq_promoters <- all_peaks %>% dplyr::filter(annotation == "Promoter")

atac_up_prom <- atacseq_promoters %>% dplyr::filter(logFC > 0.5 & FDR < 0.1)
atac_up_prom <- atac_up_prom %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)

atac_down_prom <- atacseq_promoters %>% dplyr::filter(logFC < -0.5 & FDR < 0.1)
atac_down_prom <- atac_down_prom %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE) %>% na.omit(ENSEMBL)

# List of items
x <- list(RNASeq_Down = rownames(rnaseq_down), 
          ATAC_Down = atac_down_prom$ENSEMBL)

# 2D Venn diagram
pdf(file = file.path(plot_dir, "rsaseq_atacseq_down_promoters_venn.pdf"), width = 3, height = 2.5)
a <- ggVennDiagram(x, label = "count",  set_color = c("red","blue"), lwd = 0.8, lty = 1, 
                   category.names = c("RNASeq", "ATACSeq")) + 
  scale_color_manual(values = c("red","blue")) + 
  coord_flip() +
  theme(legend.position = "none") #+ scale_fill_gradient(low = "deepskyblue", high = "pink")

print(a)
dev.off()

##############################################################################
common_genes <- intersect(rownames(rnaseq_down), atac_down_prom$ENSEMBL)

common_go <- enrichGO(gene = common_genes,
                      OrgDb         = org.Mm.eg.db,
                      # universe      = rownames(diff_exp),
                      keyType       = 'ENSEMBL',
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      # qvalueCutoff  = 0.05,
                      readable      = TRUE)

common_trim <- simplify(common_go, cutoff=0.7, by="p.adjust", select_fun=min)

write.csv(common_trim, file = file.path(out_dir, "atacseq_rnaseq_promoters_common_goterms_cc.csv"))

pdf(file = file.path(plot_dir, "Common_ATAC_RNASeq_Promoters_GO.pdf"), width = 4.5, height = 4)
dotplot(down_trim, showCategory=20, font.size=8, title = "Control Upregulated GO")
dev.off()
###############################################################################

# Genes Upregulated
xx <- list(RNASeq_Up = rownames(rnaseq_up), 
           ATAC_Up = atac_up_prom$ENSEMBL)

pdf(file = file.path(plot_dir, "rsaseq_atacseq_up_promoters_venn.pdf"), width = 3, height = 2.5)
a <- ggVennDiagram(xx, label = "count",  set_color = c("red","blue"), lwd = 0.8, lty = 1, 
                   category.names = c("RNASeq", "ATACSeq")) + 
  scale_color_manual(values = c("red","blue")) + 
  coord_flip() +
  theme(legend.position = "none") #+ scale_fill_gradient(low = "deepskyblue", high = "pink")

print(a)
dev.off()

################################################################################
# Use ATACSeq Promoters only
goterm_genes <- read.csv("Data/go_term_genes_volcano_annotated.csv")
go_all <- goterm_genes %>% dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

promoters_annotated <- left_join(atacseq_promoters, go_all, by = c("SYMBOL" = "gene_name"))
# Keep the most significant ones
promoters_annotated <- promoters_annotated %>% dplyr::arrange(FDR) %>% dplyr::distinct(ENSEMBL, .keep_all = TRUE)
promoters_annotated <- promoters_annotated[order(promoters_annotated[,"GO_Term"],na.last=F),]


# Modify the above plot
promoters_annotated <- promoters_annotated %>% mutate(lab = ifelse(is.na(GO_Term), 0, SYMBOL))
promoters_annotated <- promoters_annotated %>% mutate(GO_Term = ifelse(is.na(GO_Term), 0, GO_Term))
promoters_annotated$GO_Term[promoters_annotated$GO_Term == 0 & promoters_annotated$FDR < 0.1 & abs(promoters_annotated$logFC) > 0.5] <- "SIGNIFICANT"
promoters_annotated$GO_Term[promoters_annotated$GO_Term == 0] <- "NONSIGNIFICANT"

library(ggrepel)
options(ggrepel.max.overlaps = Inf)
ggplot(promoters_annotated, aes(x=logFC, y=-log10(FDR), size = ifelse(abs(logFC) > 0.5, abs(logFC), 0.1))) +
  geom_point(data = subset(promoters_annotated, GO_Term %in% c("SIGNIFICANT", "NONSIGNIFICANT")), aes(color = GO_Term)) +
  scale_colour_manual(values = c("lightblue", "gray90"),
                      limits = c("SIGNIFICANT", "NONSIGNIFICANT")) +
  guides(color=guide_legend(title = "SIGNIFICANCE", override.aes = list(size = 6, face = "bold"))) +
  new_scale_color() +
  geom_point(data = subset(promoters_annotated, GO_Term %in% c("APOPTOSIS", "CARBOXYLIC ACID", "CELL CYCLE", "LIPID CATABOLISM", 
                                                               "MITOCHONDRIA", "NEGATIVE REGULATION OF T CELL ACTIVATION")), aes(color = GO_Term)) +
  scale_color_manual(values = c("#990099", "#303030", "#006000", "#7e0021", "#ffd320", "blue"),
                     limits = c("APOPTOSIS", "CARBOXYLIC ACID", "CELL CYCLE", "LIPID CATABOLISM", 
                                "MITOCHONDRIA", "NEGATIVE REGULATION OF T CELL ACTIVATION")) +
  geom_vline(xintercept = 0.5, color = "grey", size=0.7, linetype = "dashed") +
  geom_vline(xintercept = -0.5, color = "grey", size=0.7, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), color = "grey", size = 0.7, linetype = "dashed") +
  geom_label_repel(aes(x=logFC, y=-log10(FDR), label = ifelse(abs(logFC) > 0.5 & FDR < 0.1 & lab != 0, lab, "")),
                   size = 5) +
  theme_light() +
  ggtitle("") +
  xlab("log2 fold change") + 
  ylab("-log10(padj)") +
  theme(legend.position = "",
        plot.title = element_text(size = rel(1.5), hjust = 0.5, face = "bold"),
        axis.title = element_text(size = rel(1.25))) +
  guides(color=guide_legend(title = "FUNCTION", override.aes = list(size = 6, face = "bold"))) +
  guides(size = FALSE) +
  theme(legend.title = element_text(face = "bold")) +
  theme(legend.text = element_text(size = 10, face = "bold")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylim(c(0,8)) + xlim(c(-2.5, 1.5))

ggsave(file.path(plot_dir,"ATAC_Volcano_final.pdf"), width = 5, height = 7, dpi = 600)

##################################################################################
# Capture the common genes with all peaks
## RNASEQ and ATACSeq upregulated
# Use RNASeq annotation as reference, otherwise it sucks
rnaseq_annot <- diff_exp |> rownames_to_column("ENSEMBL")
# rnaseq_annot <- rnaseq_annot |> dplyr::select(ENSEMBL, symbol)

rnaseq_up_atac_up <- atac_up_prom %>% dplyr::filter(ENSEMBL %in% rownames(rnaseq_up))
rnaseq_up_atac_up <- left_join(rnaseq_up_atac_up, rnaseq_annot, by = "ENSEMBL")
rnaseq_up_atac_up_clean <- rnaseq_up_atac_up %>% dplyr::select(symbol, logFC, annotation)

## RNASEQ and ATACSeq downregulated
rnaseq_down_atac_down <- atac_down_prom %>% dplyr::filter(ENSEMBL %in% rownames(rnaseq_down))
rnaseq_down_atac_down <- left_join(rnaseq_down_atac_down, rnaseq_annot, by = "ENSEMBL")
rnaseq_down_atac_down_clean <- rnaseq_down_atac_down %>% dplyr::select(symbol, logFC, annotation)

# Continue from up
rnaseq_down_atac_down_all <- atac_down %>% dplyr::filter(ENSEMBL %in% rnaseq_down_atac_down$ENSEMBL)
rnaseq_down_atac_down_all <- left_join(rnaseq_down_atac_down_all, rnaseq_annot, by = "ENSEMBL")
rnaseq_down_atac_down_all_clean <- rnaseq_down_atac_down_all %>% dplyr::select(symbol, logFC, annotation)

rnaseq_atacseq_intersection <- rbind(rnaseq_up_atac_up_clean, rnaseq_down_atac_down_all_clean)

my_cols <- rnaseq_atacseq_intersection$annotation
my_cols <- sub("\\(.*", "", my_cols)

rnaseq_atacseq_intersection$annotation <- my_cols

rnaseq_atacseq_intersection <- rnaseq_atacseq_intersection %>% dplyr::arrange(logFC)
rnaseq_atacseq_intersection$symbol <- factor(rnaseq_atacseq_intersection$symbol, levels = unique(rnaseq_atacseq_intersection$symbol))

d <- rnaseq_atacseq_intersection %>% ggplot(aes(symbol, logFC, color = annotation)) +
  geom_point(size = 4) + coord_flip() + theme_bw() + 
  # scale_color_manual() +
  geom_hline(yintercept = 0.0, color = "red", size=0.7, linetype = "dashed") +
  scale_color_manual(values = c("blue", "red", "yellow", "cyan"),
                     limits = c("Distal Intergenic", "Promoter", "Intron ", "Exon ")) + 
  guides(color=guide_legend(title = "Annotation", override.aes = list(size = 5, face = "bold"))) +
  theme(axis.text.y=element_blank()) + xlab("") + ylab("ATACSeq Log2 FC")

# Now the other
## Set the order of the genes according to ATACSeq
atac_up_rnaseq_up <- rnaseq_up_atac_up %>% dplyr::select(symbol, log2FoldChange)
atac_down_rnaseq_down <- rnaseq_down_atac_down %>% dplyr::select(symbol, log2FoldChange)
atac_intersect_rnaseq <- rbind(atac_up_rnaseq_up, atac_down_rnaseq_down)
atac_intersect_rnaseq$symbol <- factor(atac_intersect_rnaseq$symbol, levels = unique(rnaseq_atacseq_intersection$symbol))

e <- atac_intersect_rnaseq %>% ggplot(aes(symbol, log2FoldChange)) +
  geom_point(size = 4, color = "red") + coord_flip() + theme_bw() +
  geom_hline(yintercept = 0.0, color = "red", size=0.7, linetype = "dashed") +
  theme(axis.text = element_text(face="bold")) + xlab("") + ylab("RNASeq Log2 FC")

library(gridExtra)
pdf(file = file.path(plot_dir, "dotplots_rnaseq_atacseq_intersection.pdf"), width = 7, height = 7)
f <- grid.arrange(e,d, ncol = 2)
dev.off()

# Modify more
d_no_legend <- d + theme(legend.position = "none")
# Extract the legend
legend <- cowplot::get_legend(
  d + theme(legend.position = "right", legend.justification = "center"))

# Align plots and legend
b <- cowplot::plot_grid(
  cowplot::plot_grid(e, d_no_legend, ncol = 2, align = "h", axis = "tb"),
  legend,
  ncol = 2,
  rel_widths = c(7, 1.5) # Adjust relative widths for plot and legend
)

pdf(file = file.path(plot_dir, "dotplots_rnaseq_atacseq_intersection2.pdf"), width = 7, height = 7)
print(b)
dev.off()
################################################################################
# Boxplot
promoters <- rnaseq_atacseq_intersection %>% dplyr::arrange(logFC) %>% 
  dplyr::filter(annotation == "Promoter" & logFC < 0) %>% dplyr::distinct(symbol)
promoters$Group <- rep("Common", nrow(promoters))

promoters_boxplot <- left_join(promoters_annotated, promoters, by = c("SYMBOL" = "symbol"))
promoters_boxplot <- promoters_boxplot %>% mutate(Group = ifelse(is.na(Group), 0, Group))
promoters_boxplot$lab <- ifelse(promoters_boxplot$Group == 0, "Rest", "Down")

library(ggsignif)
promoters_boxplot %>% ggplot(aes(lab, logFC, fill = lab)) + geom_boxplot() + 
  theme_classic() + scale_fill_manual(values=c("dodgerblue", "royalblue4")) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(face = "bold",colour = "black", size = 12, hjust=0.5),
        axis.title.y = element_text(size = 12, face = "bold", color = "black")) +
  theme(legend.position = "none") +
  ylab("logFold Change PA/CTR") + 
  geom_signif(comparisons = list(c("Down", "Rest")), map_signif_level=TRUE)

ggsave(file.path(plot_dir,"boxplot.pdf"), width = 3, height = 4, dpi = 600)
######################################################################################
## RADAR Plots for GO Terms
library(ggradar)
library(scales)
# Upregulated GO Terms
palm_up_goterms <- read.csv("SpiderPlots/palmitate_upregulated_goterms_spider.csv")
palm_up_goterms <- palm_up_goterms %>% dplyr::mutate(GR = Count/792, logqval = round(-log10(qvalue),1))
palm_up_goterms <- palm_up_goterms %>% arrange(-Count) %>% dplyr::mutate(GR = round(GR, 3))
palm_up_table <- palm_up_goterms %>% dplyr::select(Description, GR, logqval)
palm_up_table$Palmitate <- rep("Upregulated", nrow(palm_up_table))

# Downregulated GO Terms
palm_down_goterms <- read.csv("SpiderPlots/palmitate_downregulated_goterms_spider.csv")
palm_down_goterms <- palm_down_goterms %>% dplyr::mutate(GR = Count/442, logqval = round(-log10(qvalue),1))
palm_down_goterms <- palm_down_goterms %>% arrange(-Count) %>% dplyr::mutate(GR = round(GR, 3))
palm_down_table <- palm_down_goterms %>% dplyr::select(Description, GR, logqval)
palm_down_table$Palmitate <- rep("Downregulated", nrow(palm_down_table))

# Join the tables
palm_data <- rbind(palm_up_table, palm_down_table)

# # Capitalize
# palm_data$Description <- title_case_bio(palm_data$Description)

# Cut long text to shorter pieces
your_data <- palm_data %>%
  mutate(
    Description_wrapped = sapply(Description, wrap_label_4line, width = 15, USE.NAMES = FALSE),
    x = factor(seq_along(Description), levels = seq_along(Description))
  )

# Arrange for plotting
group_ranges <- your_data %>%
  group_by(Palmitate) %>%
  summarise(
    start = min(as.numeric(x)),
    end = max(as.numeric(x)),
    .groups = "drop"
  )


# Define the outer radius
circle_breaks <- c(0.02, 0.04, 0.06)
r_max   <- max(circle_breaks)
# where the labels will sit (inside)
label_r <- r_max


p <- ggplot(your_data, aes(x = factor(x), y = GR, fill = logqval)) +
  # Background shading for groups
  geom_rect(
    data = group_ranges,
    aes(xmin = start - 0.5, xmax = end + 0.5, ymin = -Inf, ymax = Inf),
    fill = ifelse(group_ranges$Palmitate == "Upregulated", "darkturquoise", "pink"),
    alpha = 0.45,
    inherit.aes = FALSE
  ) +
  # Bars
  geom_col(width = 0.7, color = "black") +
  
  # Color scale
  scale_fill_gradient(low = "lightblue", high = "darkblue", name = "logqval") +
  
  # Make the panel a bit smaller and keep labels visible
  scale_y_continuous(limits = c(0, r_max + 0.01), breaks = circle_breaks, expand = c(0, 0)) +
  coord_polar(start = 0, clip = "off") +  # << important: don't clip labels
  
  # Remove default x labels (we'll draw our own)
  scale_x_discrete(labels = NULL) +
  
  # Circle labels (0.02, 0.04, …)
  annotate("text",
           x = 0.5, y = circle_breaks,
           label = circle_breaks,
           color = "red", size = 4.5, fontface = "bold") +
  
  # Custom x labels drawn inside the circle
  geom_text(
    data = your_data,
    aes(x = factor(x), y = r_max + 0.01, label = Description_wrapped),
    inherit.aes = FALSE,
    size = 5,              # text size you want
    fontface = "bold",
    lineheight = 0.9
  ) +
  
  # Theme: keep fonts large, tighten panel, add margins so nothing is cut
  theme_minimal(base_size = 16) +
  theme(
    panel.border       = element_blank(),
    axis.line          = element_blank(),
    axis.ticks         = element_blank(),
    axis.text.y        = element_blank(),
    axis.title         = element_blank(),
    panel.grid.major.y = element_line(color = "gray30", size = 0.3, linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray20", size = 0.3, linetype = "dashed"),
    legend.position    = "bottom",
    legend.title       = element_text(size = 16, face = "bold"),
    legend.text        = element_text(size = 14),
    plot.margin        = margin(10, 24, 10, 24)
  )

ggsave(file.path(plot_dir, "radar_small6.pdf"), p,
       width = 10, height = 10, dpi = 300, device = cairo_pdf)


##################################################

