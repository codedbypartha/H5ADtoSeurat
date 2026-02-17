# H5ADtoSeurat

# Easily Convert AnnData (.h5ad) Object to Seurat Object

This R script provides a flexible wrapper function to convert **AnnData (adata)** objects stored in **h5ad files** (from Python/Scanpy/Squidpy) to **Seurat** objects in R. It preserves all assays, allows gene ID mapping, and lets you rename the main Seurat assay.

---

## Features

* Converts a h5ad file to a Seurat object.
* Handles both **human (HGNC)** and **mouse (MGI)** gene IDs.
* Optional automatic gene ID mapping if rownames are Ensembl IDs.
* Preserves **all other assays** in the SCE as separate Seurat assays.
* Allows **custom naming** of the main Seurat assay.
* Flexible handling of counts and normalized data assays.

---

## Installation

### 1. Install the dependencies

```r
# Install required packages if not already installed
install.packages(c("Seurat", "SingleCellExperiment", "biomaRt", "devtools"))
BiocManager::install("zellkonverter")
BiocManager::install("SeuratDisk")
```

### 2. Install the H5ADtoSeurat package from GitHub

```r
library(devtools)
install_github("codedbypartha/H5ADtoSeurat")
```

---

## Usage

### 1. Load the wrapper function

```r
library(H5ADtoSeurat)
ls("package:H5ADtoSeurat")
# It should show "convert_h5ad_to_seurat"
```

### 2. Convert a human h5ad file

```r
seu_human <- convert_h5ad_to_seurat(
  h5ad_file = "test_data.h5ad",
  species = "hsapiens",
  counts_assay = "raw_counts",   # Name of raw counts assay in SCE
  data_assay = "X",              # Name of normalized data assay in SCE
  new_main_assay = "Xenium",     # Name for main Seurat assay
  gene_ids_are_symbols = FALSE   # FALSE if rownames are Ensembl IDs
)
```

### 3. Convert a mouse h5ad file with pre-mapped gene symbols

```r
seu_mouse <- convert_h5ad_to_seurat(
  h5ad_file = "mouse_data.h5ad",
  species = "mmusculus",
  counts_assay = "counts",
  data_assay = "logcounts",
  new_main_assay = "RNA",
  gene_ids_are_symbols = TRUE
)
```

---

## Function Arguments

| Argument               | Description                                                                      |
| ---------------------- | -------------------------------------------------------------------------------- |
| `h5ad_file`            | Path to the `.h5ad` file.                                                        |
| `species`              | `"hsapiens"` for human, `"mmusculus"` for mouse. Determines gene mapping.        |
| `counts_assay`         | Name of the raw counts assay in the SCE object.                                  |
| `data_assay`           | Name of the normalized/processed data assay in the SCE object.                   |
| `new_main_assay`       | Name for the main Seurat assay. Defaults to `"Xenium"`.                          |
| `gene_ids_are_symbols` | Set `TRUE` if the rownames are already gene symbols (HGNC/MGI). Default `FALSE`. |

---

## Output

* Returns a **Seurat object** with:

  * Main assay containing counts and normalized data.
  * All remaining SCE assays preserved as Seurat assays.
  * Default assay set to the `new_main_assay` name.

* Example inspection:

```r
seu_human
Assays(seu_human)
DefaultAssay(seu_human)
```

---

## Notes

* You must provide the correct counts and data assay names if they differ from `"raw_counts"` and `"X"`.
* All other SCE assays will be added automatically.
* If some genes cannot be mapped to symbols, the original rownames are retained.
* In case you have only raw or normalized data, use the same assay name in both data_assay and counts_assay.

---

## License

This script is open for academic and research use.

---
