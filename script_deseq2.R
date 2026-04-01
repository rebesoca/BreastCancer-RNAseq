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
