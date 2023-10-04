# 10/4/2023 Single cell analysis demo for BME 4550/6550
# Part 1: Seurat analysis of Amadei et al. dataset

# Test add a new comment

library(Seurat)
library(data.table)


# Location of data files
setwd("G:\\My Drive/Working/ETX/")

fix_data <- function(txt_to_fix) {
  
  # Read in dataset
  data_in <- suppressWarnings(fread(txt_to_fix, data.table=FALSE, showProgress=FALSE))
  
  # Clean up by converting gene names (first column) to row names
  first_col_name <- colnames(data_in[1])
  data_in <- data_in %>% column_to_rownames(first_col_name)
  
  colnames(data_in) <- gsub("_(?=[^.]*_)", "-", colnames(data_in), perl=TRUE)
  
  # Downsample if desired
  if (DS_BOOL) {
    if (ncol(data_in) > DOWNSAMPLE) {
      data_in <- data_in[,sample(ncol(data_in), DOWNSAMPLE)]
    }
  }
  
  return(data_in)
}

# Function to build Seurat object and add Doublet Scores (calculated with Scrublet in Python)
# This will be used in a for loop to add all datasets to a single Seurat object. . .
seurat_in <- function(in_file) {
  # Create Seurat object with cleaned up/fixed data from a single .txt.gz file
  seurat_obj <- CreateSeuratObject(counts = fix_data(in_file), names.field = 1,
                                   names.delim = "_", project = PROJECT)
  # Read in doublet scores for the corresponding dataset
  ds <- fread(sub(".inex.txt.gz",".scr.txt",in_file), data.table=FALSE)
  rownames(ds) <- colnames(seurat_obj)
  
  # Add doublet scores to Seurat object
  seurat_obj = AddMetaData(object=seurat_obj, metadata=ds, col.name="DoubletScores")
  
  return(seurat_obj)
}

# Get all sample files to input
in_files <- list.files(pattern=".txt.gz")

# Initialize Seurat object with the first sample
cat("loading sample 1\n")
seurat_obj <- seurat_in(in_files[1])

# Add the rest of the sample files to the Seurat object one at a time
for (i in 2:length(in_files)) {
  cat("loading sample", i, "\n")
  seurat_obj <- merge(seurat_obj, y = seurat_in(in_files[i]), project = PROJECT)
}

# Export to save time later
saveRDS(seurat_obj, file = "full_dataset.rds")

