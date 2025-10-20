# Write functions for them
## For GO Term Enrichment
go_enrich_func <- function(diff_exp, adjp, lfc, naming, trim_val){
  
  # diff_exp is the differential expression table 
  # Important columns: gene_id: ensemble, log2FoldChange, padj values
  # adjp = padj
  # lfc = log2FoldChange
  # naming, prefix to appear at the end of the csv file
  # minsize: minimum number of terms for goenrich
  # trim_val: trimming value to remove redunadant terms
  
  plot_list = list()
  
  sig_diff_up <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange > lfc)
  sig_diff_down <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange < -lfc)
  
  goterms_up <- enrichGO(gene  = rownames(sig_diff_up),
                         OrgDb         = org.Mm.eg.db,
                         # universe = rownames(diff_exp),
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)
  
  up_trim <- clusterProfiler::simplify(goterms_up, cutoff=trim_val, by="p.adjust", select_fun=min)
  
  p_up = dotplot(up_trim, showCategory=20, title = "Upregulated GO")
  
  plot_list[[1]] = p_up
  
  write.csv(up_trim, file = file.path(out_dir, "palmitate_upregulated_goterms.csv"))
  
  goterms_down <- enrichGO(gene  = rownames(sig_diff_down),
                           OrgDb         = org.Mm.eg.db,
                           # universe = rownames(diff_exp),
                           keyType       = 'ENSEMBL',
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
  
  down_trim <- clusterProfiler::simplify(goterms_down, cutoff=trim_val, by="p.adjust", select_fun=min)
  
  p_down = dotplot(down_trim, showCategory=20, title = "Downregulated GO")
  
  plot_list[[2]] = p_down
  
  write.csv(down_trim, file = file.path(out_dir, "palmitate_downregulated_goterms.csv"))
  
  return(plot_list)
}


## For KEGG Pathways
kegg_enrich_func <- function(diff_exp, adjp, lfc, naming){
  
  # diff_exp is the differential expression table 
  # Important columns: gene_id: ensemble, log2FoldChange, padj, entrez
  # adjp = padj
  # lfc = log2FoldChange
  # naming, prefix to appear at the end of the csv file
  
  plot_list = list()
  
  sig_diff_up <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange > lfc)
  sig_diff_down <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange < -lfc)
  
  entrez_ids_up = mapIds(org.Mm.eg.db,
                         keys= rownames(sig_diff_up),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrez_ids_up = entrez_ids_up[!is.na(entrez_ids_up)]
  
  kegg_up <- enrichKEGG(gene = entrez_ids_up,
                        organism = "mmu",
                        keyType = "kegg",
                        qvalueCutoff  = 0.05,
                        pAdjustMethod = "BH")
  
  p_up = dotplot(kegg_up, showCategory=20, title = "Upregulated KEGG")
  
  plot_list[[1]] = p_up
  
  write.csv(kegg_up, file = file.path(out_dir, "palmitate_upregulated_kegg.csv"))
  
  entrez_ids_dn = mapIds(org.Mm.eg.db,
                         keys= rownames(sig_diff_down),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrez_ids_dn = entrez_ids_dn[!is.na(entrez_ids_dn)]
  
  kegg_down <- enrichKEGG(gene = entrez_ids_dn,
                          organism = "mmu",
                          keyType = "kegg",
                          qvalueCutoff  = 0.05,
                          pAdjustMethod = "BH")
  
  p_down = dotplot(kegg_down, showCategory=20, title = "Downregulated KEGG")
  
  plot_list[[2]] = p_down
  
  write.csv(kegg_down, file = file.path(out_dir, "palmitate_downregulated_Kegg.csv"))
  
  return(plot_list)
}

## For Reactome Enrichment
reactome_enrich_func <- function(diff_exp, adjp, lfc, naming){
  
  # diff_exp is the differential expression table 
  # Important columns: gene_id: ensemble, log2FoldChange, padj values
  # adjp = padj
  # lfc = log2FoldChange
  # naming, prefix to appear at the end of the csv file
  # minsize: minimum number of terms for goenrich
  
  plot_list = list()
  
  sig_diff_up <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange > lfc)
  sig_diff_down <- diff_exp %>% dplyr::filter(padj < adjp & log2FoldChange < -lfc)
  
  entrez_ids_up = mapIds(org.Mm.eg.db,
                         keys= rownames(sig_diff_up),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrez_ids_up = entrez_ids_up[!is.na(entrez_ids_up)]
  
  reactome_up <- enrichPathway(gene= entrez_ids_up,
                               qvalueCutoff = 0.05,
                               organism = "mouse",
                               pAdjustMethod = "BH",
                               readable=TRUE)
  
  write.csv(reactome_up, file = file.path(out_dir, "palmitate_upregulated_reactome.csv"))
  
  entrez_ids_dn = mapIds(org.Mm.eg.db,
                         keys= rownames(sig_diff_down),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
  entrez_ids_dn = entrez_ids_dn[!is.na(entrez_ids_dn)]
  
  reactome_dn <- enrichPathway(gene= entrez_ids_dn,
                               qvalueCutoff = 0.05,
                               organism = "mouse",
                               pAdjustMethod = "BH",
                               readable=TRUE)
  
  write.csv(reactome_dn, file = file.path(out_dir, "palmitate_downregulated_reactome.csv"))
  
  p_up = dotplot(reactome_up, showCategory=20, title = "Upregulated Reactome")
  
  plot_list[[1]] = p_up
  
  p_down = dotplot(reactome_dn, showCategory=20, title = "Downregulated Reactome")
  
  plot_list[[2]] = p_down
  
  return(plot_list)
}

# Define a function to automatically cut long labels up to 4 separate lines
# Set the default 20, change when drawing the plot
wrap_label_4line <- function(s, width = 20) {
  if (nchar(s) <= width) return(s)
  
  find_break <- function(str, w = width) {
    idx <- gregexpr(" ", str)[[1]]
    pos <- idx[idx <= w]
    if (length(pos) == 0) return(w)
    return(max(pos))
  }
  
  # Line 1
  pos1 <- find_break(s)
  line1 <- substr(s, 1, pos1)
  rest1 <- substr(s, pos1 + 1, nchar(s))
  if (nchar(rest1) <= width) return(paste0(line1, "\n", rest1))
  
  # Line 2
  pos2 <- find_break(rest1)
  line2 <- substr(rest1, 1, pos2)
  rest2 <- substr(rest1, pos2 + 1, nchar(rest1))
  if (nchar(rest2) <= width) return(paste0(line1, "\n", line2, "\n", rest2))
  
  # Line 3
  pos3 <- find_break(rest2)
  line3 <- substr(rest2, 1, pos3)
  rest3 <- substr(rest2, pos3 + 1, nchar(rest2))
  
  # Line 4 (final)
  if (nchar(rest3) <= width) {
    result <- paste0(line1, "\n", line2, "\n", line3, "\n", rest3)
  } else {
    # Truncate last line if still too long
    rest3 <- paste0(substr(rest3, 1, width - 3), "...")
    result <- paste0(line1, "\n", line2, "\n", line3, "\n", rest3)
  }
  
  return(result)
}
