library(GEOquery)
library(tidyverse)

# Set working directory to Output directory created
setwd("path/to/Output")

# Create the folder
dir.create("GSEA/Input", recursive=TRUE, showWarnings = FALSE)

# Loading the downloaded series_matrix file
my_id = "GSE64217"
gse = getGEO(my_id)
gse = gse[[1]]

# Expression data from Series_matrix
expression.data = exprs(gse)
summary(expression.data)

#expression.data = log2(exprs(gse))
#summary(expression.data)

pData(gse)$data_processing[1]

# Writing expression data into a .gct file
expression.gct <- as.data.frame(expression.data) %>%
  tibble::rownames_to_column(var = "Name") %>%
  mutate(Description = "na") %>%  # Duplicate gene name as description
  select(Name, Description, everything())

# Get dimensions
num_genes <- nrow(expression.gct)
num_samples <- ncol(expression.gct) - 2  # Exclude 'Name' and 'Description'

# Output file path
output_file <- "GSEA/Input/GSE64217expression_data.gct"

# Writing GCT file
writeLines("#1.2", output_file)

# Writing the dimensions (genes x samples)
write(paste(num_genes, num_samples, sep = "\t"), output_file, append = TRUE)


# Writing the actual data (gene expression values)
write.table(expression.gct[1, ], output_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE, append = TRUE)
write.table(expression.gct[-1, ], output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
