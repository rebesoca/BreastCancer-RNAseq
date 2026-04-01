#install.packages("BiocManager")
#BiocManager::install("DESeq2")
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

#Crear funcion general de analisis
analizar_contraste <- function(res, nombre) {
  
  cat("\n=============================\n")
  cat("Análisis:", nombre, "\n")
  cat("=============================\n")
  
  res_df <- as.data.frame(res)
  
  # Filtrar genes significativos
  res_sig <- res_df[which(res_df$padj < 0.05), ]
  
  # UP / DOWN
  up <- res_sig[res_sig$log2FoldChange > 0, ]
  down <- res_sig[res_sig$log2FoldChange < 0, ]
  
  cat("Genes UP:", nrow(up), "\n")
  cat("Genes DOWN:", nrow(down), "\n")
  
  # Top genes
  top_genes <- res_sig[order(res_sig$padj), ]
  print(head(top_genes, 10))
  
  # Filtro biológico
  res_strict <- res_sig[abs(res_sig$log2FoldChange) > 1, ]
  cat("Genes con cambio fuerte (|log2FC| > 1):", nrow(res_strict), "\n")
  
  # Volcano plot
  library(ggplot2)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(alpha = 0.5) +
    ggtitle(paste("Volcano plot:", nombre)) +
    theme_minimal()
  
  print(p)
  
  # Guardar resultados
  write.csv(res_df, paste0("DESeq2_", nombre, ".csv"))
}


#Ejecutar para las 3 comparaciones
analizar_contraste(res_dox, "DOX_vs_CTR")
analizar_contraste(res_btz, "BTZ_vs_CTR")
analizar_contraste(res_tr,  "BTZ_vs_DOX")


#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="condition")

