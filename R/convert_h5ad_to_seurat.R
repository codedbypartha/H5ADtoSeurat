#' Convert H5AD to Seurat
#'
#' This function reads a .h5ad file (or SCE) and converts it into a Seurat object,
#' optionally mapping Ensembl IDs to gene symbols and preserving all assays.
#'
#' @param h5ad_file Path to the .h5ad file
#' @param species Either "hsapiens" or "mmusculus" for gene mapping
#' @param counts_assay Name of counts assay in SCE (default: "raw_counts")
#' @param data_assay Name of normalized data assay in SCE (default: "X")
#' @param new_main_assay Name for the Seurat main assay (default: "Xenium")
#' @param gene_ids_are_symbols TRUE if rownames are already gene symbols
#' @return Seurat object
#' @export


# -------------------------------
# Wrapper function: fully automated AnnData → Seurat with all assays
# -------------------------------
convert_h5ad_to_seurat <- function(
  h5ad_file,                    # Path to h5ad file
  species = c("hsapiens", "mmusculus"),  # Species for gene mapping
  counts_assay = "raw_counts",  # Counts assay in SCE
  data_assay = "X",             # Normalized data assay in SCE
  new_main_assay = "Xenium",    # Name for Seurat main assay
  gene_ids_are_symbols = FALSE  # TRUE if rownames are already HGNC/MGI symbols
) {
  species <- match.arg(species)
  
  # -------------------------------
  # 1. Load AnnData object
  # -------------------------------
  sce <- zellkonverter::readH5AD(h5ad_file)
  
  # -------------------------------
  # 2. Gene ID mapping if needed
  # -------------------------------
  if (!gene_ids_are_symbols) {
    # Choose dataset
    mart_dataset <- switch(species,
                           hsapiens = "hsapiens_gene_ensembl",
                           mmusculus = "mmusculus_gene_ensembl")
    mart <- biomaRt::useEnsembl(biomart = "genes", dataset = mart_dataset)
    
    gene_ids <- rownames(sce)
    
    # Get gene symbols
    attrs <- switch(species,
                    hsapiens = c("ensembl_gene_id", "hgnc_symbol"),
                    mmusculus = c("ensembl_gene_id", "mgi_symbol"))
    
    mapping <- biomaRt::getBM(
      attributes = attrs,
      filters = "ensembl_gene_id",
      values = gene_ids,
      mart = mart
    )
    
    # Standardize column name
    colnames(mapping)[2] <- "symbol"
    
    # Map rownames
    symbol_map <- mapping$symbol
    names(symbol_map) <- mapping$ensembl_gene_id
    
    new_symbols <- symbol_map[rownames(sce)]
    new_symbols[is.na(new_symbols) | new_symbols == ""] <- rownames(sce)[is.na(new_symbols) | new_symbols == ""]
    
    rownames(sce) <- new_symbols
  }
  
  # -------------------------------
  # 3. Convert main assay to Seurat
  # -------------------------------
  seu <- Seurat::as.Seurat(sce, counts = counts_assay, data = data_assay)
  
  # Rename main assay
  seu <- Seurat::RenameAssays(seu, originalexp = new_main_assay)
  Seurat::DefaultAssay(seu) <- new_main_assay
  
  
  # -------------------------------
  # 4. Convert remaining SCE assays automatically
  # -------------------------------
  all_assays <- SingleCellExperiment::assayNames(sce)
  
  # Skip counts/data assays already converted
  remaining <- setdiff(all_assays, c(counts_assay, data_assay))
  
  for (a in remaining) {
    # Add each assay as its own Seurat assay object
    seu[[a]] <- Seurat::CreateAssayObject(counts = SingleCellExperiment::assays(sce)[[a]])
  }
  
  return(seu)
}
