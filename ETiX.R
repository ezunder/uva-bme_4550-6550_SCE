# 10/4/2023 Single cell analysis demo for BME 4550/6550
# Part 1: Seurat analysis of Amadei et al. dataset

# Test add a new comment

library(Seurat)
library(data.table)


# Location of data files
setwd("G:\\My Drive/Working/ETX/")

# Read in dataset
data_in <- suppressWarnings(fread("./GSM5701522_ETiX_d5_1.inex.txt.gz", data.table=FALSE, showProgress=FALSE))

# Clean up by converting gene names (first column) to row names
first_col_name <- colnames(data_in[1])
data_in <- data_in %>% column_to_rownames(first_col_name)

colnames(data_in) <- gsub("_(?=[^.]*_)", "-", colnames(data_in), perl=TRUE)



remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

