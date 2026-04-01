install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# Leer los archivos de featureCounts
files <- c(
  "BTZ_23h_rep1.txt",
  "BTZ_23h_rep2.txt",
  "BTZ_47h_rep1.txt",
  "BTZ_47h_rep2.txt",
  "CTR_23h_rep1.txt",
  "CTR_23h_rep2.txt",
  "CTR_47h_rep1.txt",
  "CTR_47h_rep2.txt",
  "DOX_23h_rep1.txt",
  "DOX_23h_rep2.txt",
  "DOX_47h_rep1.txt",
  "DOX_47h_rep2.txt"
)

#Construir matriz de conteos
count_list <- lapply(files, function(f) {
  df <- read.table(f, header = TRUE)
  rownames(df) <- df[,1]
  df[,2, drop=FALSE]
})

count_matrix <- do.call(cbind, count_list)
colnames(count_matrix) <- gsub(".txt", "", files)

#Crear metadata (condiciones)
condition <- factor(c(
  "BTZ","BTZ","BTZ","BTZ",
  "CTR","CTR","CTR","CTR",
  "DOX","DOX","DOX","DOX"
))

coldata <- data.frame(
  row.names = colnames(count_matrix),
  condition = condition
)

#Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = coldata,
  design = ~ condition
)

#Filtrado básico
dds <- dds[rowSums(counts(dds)) > 10, ]

#Ejecutar DESeq2
dds <- DESeq(dds)

#Resultados 
res_dox <- results(dds, contrast = c("condition", "DOX", "CTR")) #DOX vs CTR
res_btz <- results(dds, contrast = c("condition", "BTZ", "CTR")) #BTZ vs CTR
res_tr <- results(dds, contrast = c("condition", "BTZ", "DOX")) #BTZ vs DOX


#Guardar resultados
write.csv(as.data.frame(res_dox), "DESeq2_DOX_vs_CTR.csv")
write.csv(as.data.frame(res_btz), "DESeq2_BTZ_vs_CTR.csv")
write.csv(as.data.frame(res_tr), "DESeq2_BTZ_vs_DOX.csv")

#Comprobaciones
head(res_dox) #Primeras filas
sum(res_dox$padj < 0.05, na.rm = TRUE) #Número de genes significativos
summary(res_dox$padj) #Distribución de p-values
sum(is.na(res_dox$padj)) #NAs

###############################################################################
#Analisis de los resultados 

library(ggplot2)
library(ggrepel)

analizar_contraste <- function(res, nombre) {
  
  cat("\n=============================\n")
  cat("Análisis:", nombre, "\n")
  cat("=============================\n")
  
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Filtrar genes significativos
  res_sig <- res_df[which(res_df$padj < 0.05), ]
  up   <- res_sig[res_sig$log2FoldChange >  1, ]
  down <- res_sig[res_sig$log2FoldChange < -1, ]
  
  cat("Genes UP:", nrow(up), "\n")
  cat("Genes DOWN:", nrow(down), "\n")
  
  top_genes <- res_sig[order(res_sig$padj), ]
  print(head(top_genes, 10))
  
  res_strict <- res_sig[abs(res_sig$log2FoldChange) > 1, ]
  cat("Genes con cambio fuerte (|log2FC| > 1):", nrow(res_strict), "\n")
  
  # Clasificar puntos
  res_df$categoria <- "No significativo"
  res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange >  1] <- "UP"
  res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "DOWN"
  res_df$categoria <- factor(res_df$categoria,
                             levels = c("UP", "DOWN", "No significativo"))
  
  # Límite eje Y para valores extremos — ANTES de crear top_label
  res_df$padj_plot <- res_df$padj
  res_df$padj_plot[!is.na(res_df$padj_plot) & res_df$padj_plot < 1e-100] <- 1e-100
  
  # Top 15 genes para etiquetar — DESPUÉS de crear padj_plot
  top_label <- res_df[res_df$categoria != "No significativo", ]
  top_label <- top_label[order(top_label$padj), ]
  top_label <- head(top_label, 15)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj_plot),
                          color = categoria)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_point(data = subset(res_df, categoria != "No significativo"),
               alpha = 0.8, size = 2) +
    geom_text_repel(data = top_label,
                    aes(label = gene),
                    size = 3, max.overlaps = 20,
                    box.padding = 0.4) +
    geom_vline(xintercept = c(-1, 1),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_color_manual(values = c(
      "UP"               = "#E05C5C",
      "DOWN"             = "#4E9AF1",
      "No significativo" = "grey70"
    )) +
    annotate("text", x = max(res_df$log2FoldChange, na.rm = TRUE) * 0.8,
             y = max(-log10(res_df$padj_plot), na.rm = TRUE) * 0.95,
             label = paste("UP:", nrow(up)), color = "#E05C5C",
             fontface = "bold", size = 4) +
    annotate("text", x = min(res_df$log2FoldChange, na.rm = TRUE) * 0.8,
             y = max(-log10(res_df$padj_plot), na.rm = TRUE) * 0.95,
             label = paste("DOWN:", nrow(down)), color = "#4E9AF1",
             fontface = "bold", size = 4) +
    xlab("log2 Fold Change") +
    ylab("-log10 (p-valor ajustado)") +
    ggtitle(paste("Volcano plot:", gsub("_", " ", nombre)),
            subtitle = "Umbral: |log2FC| > 1 y p.adj < 0.05") +
    theme_classic(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "right",
      legend.title  = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)
    ) +
    labs(color = "Expresión")
  
  print(p)
  ggsave(paste0("Volcano_", nombre, ".png"), plot = p,
         width = 8, height = 6, dpi = 300)
  
  write.csv(res_df, paste0("DESeq2_", nombre, ".csv"))
}

# Ejecutar para las 3 comparaciones
analizar_contraste(res_dox, "DOX_vs_CTR")
analizar_contraste(res_btz, "BTZ_vs_CTR")
analizar_contraste(res_tr,  "BTZ_vs_DOX")

# Calcular PCA
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Añadir información de tiempo a los metadatos
pca_data$time <- factor(c("23h","23h","47h","47h",
                          "23h","23h","47h","47h",
                          "23h","23h","47h","47h"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = time)) +
  geom_point(size = 5, alpha = 0.9) +
  scale_color_manual(values = c(
    "CTR" = "#4E9AF1",
    "DOX" = "#E05C5C",
    "BTZ" = "#57B85A"
  )) +
  scale_shape_manual(values = c("23h" = 16, "47h" = 17)) +
  xlab(paste0("PC1: ", percentVar[1], "% varianza")) +
  ylab(paste0("PC2: ", percentVar[2], "% varianza")) +
  ggtitle("Análisis de componentes principales (PCA)",
          subtitle = "Muestras coloreadas por condición, forma por tiempo") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.3)
  ) +
  labs(color = "Condición", shape = "Tiempo")

ggsave("PCA_condicion_tiempo.png", width = 8, height = 6, dpi = 300)

#Crear heatmap
library(pheatmap)

# Transformación (importante)
vsd <- vst(dds, blind = FALSE)

# Seleccionar genes significativos (ejemplo DOX vs CTR)
res_sig <- res_dox[which(res_dox$padj < 0.05), ]

# Ordenar por significancia
top_genes <- rownames(res_sig[order(res_sig$padj), ])

# Quedarte con los top 50
top_genes <- top_genes[1:50]

# Extraer matriz de expresión
mat <- assay(vsd)[top_genes, ]

# Escalar por gen (muy importante para visualizar)
mat <- t(scale(t(mat)))

pheatmap(mat,
         annotation_col = as.data.frame(colData(dds)[, "condition", drop=FALSE]),
         show_rownames = FALSE,
         cluster_cols = TRUE,
         cluster_rows = TRUE)

#Enriquecimiento funcional
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

#Preparar los genes
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)  


# Quitar NA
res <- res_dox[!is.na(res_dox$padj), ]

# Filtrar significativos
res_sig <- res[res$padj < 0.05, ]

# Extraer nombres de genes
genes <- rownames(res_sig)

#Convertir a Entrez
genes_entrez <- bitr(genes,
                     fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

# GO plot
p_go <- dotplot(ego, showCategory = 20) +
  ggtitle("GO enrichment (Biological Process)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  scale_color_gradient(low = "blue", high = "red")

# KEGG plot
p_kegg <- dotplot(ekegg, showCategory = 20) +
  ggtitle("KEGG pathway enrichment") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  ) +
  scale_color_gradient(low = "blue", high = "red")

# Mostrar
p_go
p_kegg

#Enriquecimiento GSEA

preparar_gsea <- function(res, nombre) {
  
  # Convertir resultados a data.frame
  res_df <- as.data.frame(res)
  
  # Eliminar NAs en padj
  res_df <- res_df[!is.na(res_df$padj), ]
  
  # Ordenar por log2FoldChange de mayor a menor
  # GSEA necesita genes ordenados por estadístico (usamos log2FC)
  gene_list <- res_df$log2FoldChange
  names(gene_list) <- rownames(res_df)
  
  # Ordenar de mayor a menor (genes más sobreexpresados al principio)
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  return(list(gene_list = gene_list, res_df = res_df))
}


ejecutar_gsea <- function(gene_list, nombre_contraste, ontologias = c("GO", "KEGG")) {
  
  cat("\n========================================\n")
  cat("GSEA para:", nombre_contraste, "\n")
  cat("========================================\n")
  
  resultados <- list()
  
  # GSEA con GO (Gene Ontology - Biological Process)
  
  if ("GO" %in% ontologias) {
    cat("\n>>> Ejecutando GSEA con GO (Biological Process)...\n")
    
    gsea_go <- gseGO(
      geneList = gene_list,
      OrgDb = org.Hs.eg.db,
      ont = "BP",           # Biological Process
      keyType = "SYMBOL",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
      resultados$GO <- gsea_go
      cat("  -> Términos GO enriquecidos:", nrow(gsea_go), "\n")
    } else {
      cat("  -> No se encontraron términos GO significativos\n")
    }
  }
  
  
  # GSEA con KEGG
  if ("KEGG" %in% ontologias) {
    cat("\n>>> Ejecutando GSEA con KEGG...\n")
    
    # Convertir símbolos a ENTREZID para KEGG
    gene_list_entrez <- gene_list
    
    # Bitr para convertir nombres
    genes_symbol <- names(gene_list_entrez)
    genes_entrez_map <- bitr(genes_symbol, 
                             fromType = "SYMBOL", 
                             toType = "ENTREZID", 
                             OrgDb = org.Hs.eg.db)
    
    # Filtrar y mantener solo genes con conversión exitosa
    gene_list_entrez <- gene_list_entrez[names(gene_list_entrez) %in% genes_entrez_map$SYMBOL]
    names(gene_list_entrez) <- genes_entrez_map$ENTREZID[match(names(gene_list_entrez), 
                                                               genes_entrez_map$SYMBOL)]
    
    gsea_kegg <- gseKEGG(
      geneList = gene_list_entrez,
      organism = "hsa",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      verbose = FALSE
    )
    
    if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
      resultados$KEGG <- gsea_kegg
      cat("  -> Rutas KEGG enriquecidas:", nrow(gsea_kegg), "\n")
    } else {
      cat("  -> No se encontraron rutas KEGG significativas\n")
    }
  }
  
  
  return(resultados)
}


graficar_gsea <- function(resultados, nombre_contraste) {
  
  if (is.null(resultados) || length(resultados) == 0) {
    cat("No hay resultados GSEA para graficar\n")
    return()
  }
  
  # Graficar GO results
  if (!is.null(resultados$GO) && nrow(resultados$GO) > 0) {
    p_go_ridge <- ridgeplot(resultados$GO, showCategory = 20) +
      ggtitle(paste("GSEA - GO Biological Process:", nombre_contraste)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_go_ridge)
    ggsave(paste0("GSEA_GO_ridgeplot_", nombre_contraste, ".png"), 
           p_go_ridge, width = 10, height = 8, dpi = 300)
    
    # Dotplot para GO
    p_go_dot <- dotplot(resultados$GO, showCategory = 15) +
      ggtitle(paste("GSEA - GO:", nombre_contraste)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_go_dot)
    ggsave(paste0("GSEA_GO_dotplot_", nombre_contraste, ".png"), 
           p_go_dot, width = 10, height = 8, dpi = 300)
  }
  
  # Graficar KEGG results
  if (!is.null(resultados$KEGG) && nrow(resultados$KEGG) > 0) {
    p_kegg_ridge <- ridgeplot(resultados$KEGG, showCategory = 20) +
      ggtitle(paste("GSEA - KEGG Pathways:", nombre_contraste)) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    print(p_kegg_ridge)
    ggsave(paste0("GSEA_KEGG_ridgeplot_", nombre_contraste, ".png"), 
           p_kegg_ridge, width = 10, height = 8, dpi = 300)
  }
  
}


# Preparar datos
gsea_dox <- preparar_gsea(res_dox, "DOX_vs_CTR")
gsea_btz <- preparar_gsea(res_btz, "BTZ_vs_CTR")
gsea_tr <- preparar_gsea(res_tr, "BTZ_vs_DOX")

# Ejecutar GSEA para cada contraste
resultados_gsea_dox <- ejecutar_gsea(gsea_dox$gene_list, "DOX_vs_CTR", ontologias = c("GO", "KEGG"))
resultados_gsea_btz <- ejecutar_gsea(gsea_btz$gene_list, "BTZ_vs_CTR", ontologias = c("GO", "KEGG"))
resultados_gsea_tr <- ejecutar_gsea(gsea_tr$gene_list, "BTZ_vs_DOX", ontologias = c("GO", "KEGG"))

# Graficar resultados
graficar_gsea(resultados_gsea_dox, "DOX_vs_CTR")
graficar_gsea(resultados_gsea_btz, "BTZ_vs_CTR")
graficar_gsea(resultados_gsea_tr, "BTZ_vs_DOX")


guardar_resultados_gsea <- function(resultados, nombre_contraste) {
  if (!is.null(resultados$GO)) {
    write.csv(as.data.frame(resultados$GO), 
              paste0("GSEA_GO_", nombre_contraste, ".csv"))
  }
  if (!is.null(resultados$KEGG)) {
    write.csv(as.data.frame(resultados$KEGG), 
              paste0("GSEA_KEGG_", nombre_contraste, ".csv"))
  }
}

guardar_resultados_gsea(resultados_gsea_dox, "DOX_vs_CTR")
guardar_resultados_gsea(resultados_gsea_btz, "BTZ_vs_CTR")
guardar_resultados_gsea(resultados_gsea_tr, "BTZ_vs_DOX")


