# A quick R script used to extract a list of mouse IFNs uoregulated in macrophages 
# Can be extended to other cells by downloading .RDS files from immune gene directory 

library(Seurat)


celltype <- "Macrophage"

ref_data <- readRDS("/home/s2759112/Project/ref_data_Macrophage.RDS")   # REPLACE PATH
seurat_obj <- readRDS("/home/s2759112/Project/ref_data_Macrophage.RDS") # REPLACE PATH
head(seurat_obj$sample)

Idents(seurat_obj) <- "sample"
markers_ifnb <- FindMarkers(seurat_obj, ident.1 = "IFNb", ident.2 = "PBS", min.pct = 0.25)
up_ifnb <- subset(markers_ifnb, avg_log2FC > 0 & p_val_adj < 0.05) # Can use >0 or >2 - be consistent

#ifns <- c("IFNa1", "IFNb", "IFNg", "IFNe", "IFNk", "IFNl2") #Comment out if not using full list 
ifns <- c("IFNa1", "IFNb")
upregulated_lists <- list()

for (cytokine in ifns) {
  markers <- FindMarkers(seurat_obj, ident.1 = cytokine, ident.2 = "PBS", min.pct = 0.25)
  up_genes <- subset(markers, avg_log2FC > 0 & p_val_adj < 0.05)
  upregulated_lists[[cytokine]] <- rownames(up_genes)
}

all_up_genes <- unique(unlist(upregulated_lists))
writeLines(all_up_genes, "upregulated_genes_IFNa_IFNb.txt")
