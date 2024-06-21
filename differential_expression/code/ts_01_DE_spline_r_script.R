library(splineTimeR)
library(Biobase)

# Function to read and preprocess expression data
preprocess_data <- function(filename) {
  data <- read.csv(filename, row.names = 1)
  exprsData <- t(as.matrix(data[, 5:ncol(data)]))
  set.seed(123)
  exprsData <- exprsData + matrix(rnorm(n = nrow(exprsData) * ncol(exprsData), mean = 0, sd = 0.05), nrow = nrow(exprsData), ncol = ncol(exprsData))
  apply(exprsData, 2, function(x) (x - min(x)) / (max(x) - min(x)))
}

# Function for differential expression analysis and enrichment analysis
perform_analysis <- function(eset, treatment, geneSets, mapping_dict, file_prefix) {
  selected_samples <- pData(eset)$Treatment %in% treatment
  eset_selected <- eset[, selected_samples]
  diffExprs <- splineDiffExprs(eSetObject = eset_selected, df = 3, cutoff.adj.pVal = 1, reference = treatment[1], intercept = TRUE)
  write.csv(diffExprs, paste0("genes/", file_prefix, ".csv"), row.names = TRUE)
  top_genes <- rownames(diffExprs)[1:15]
  
  mapped_diffExprs <- head(diffExprs, 500)
  new_row_names <- sapply(rownames(mapped_diffExprs), function(gene) {
    if (gene %in% names(mapping_dict)) {
      return(mapping_dict[[gene]])
    } else {
      return(gene)  # Keep original if no mapping found
    }
  })
  rownames(mapped_diffExprs) <- new_row_names
  
  for (geneSet in names(geneSets)) {
    enrichPath <- pathEnrich(geneList = rownames(mapped_diffExprs), geneSets = geneSets[[geneSet]], universe = 6536)
    if (!dir.exists("gsea")) {
      dir.create("gsea")
    }
    write.csv(enrichPath, paste0("gsea/", file_prefix, "_", geneSet, "_gsea.csv"), row.names = TRUE)
  }
}


calculate_regulation <- function(eset, treatment1, treatment2, file_name) {
  avgExprs_Treatment1 <- rowMeans(exprs(eset[, pData(eset)$Treatment == treatment1]))
  avgExprs_Treatment2 <- rowMeans(exprs(eset[, pData(eset)$Treatment == treatment2]))
  
  regulation <- ifelse(avgExprs_Treatment1 > avgExprs_Treatment2, "Upregulated", "Downregulated")
  regulation_df <- data.frame(row_IDs = rownames(eset), Regulation = regulation)
  
  existing_data <- read.csv(file_name)
  merged_data <- merge(existing_data, regulation_df, by = "row_IDs", all.x = TRUE)
  write.csv(merged_data, file_name, row.names = FALSE)
}


# Set working directory and read data
setwd("/Users/lenatrnovec/FAKS/BIOLAB/RNAseq_dicty/differential_expression/results/time_series") # Change to your working directory
exprsData <- preprocess_data("../combined_df.csv")

# Load gene sets
geneSets <- list(
  BP = getGmt("../genesets/GO-biological_process_mod.gmt"),
  MF = getGmt("../genesets/GO-molecular_function_mod.gmt"),
  DB = getGmt("../genesets/DictyBase-Phenotypes_mod.gmt")
)

# Create ExpressionSet
data <- read.csv("../combined_df.csv", row.names = 1)
data$sampleName <- rownames(data)
phenoData <- new("AnnotatedDataFrame", data = data.frame(SampleName = data$sampleName, replicate = data$replicate, Time = data$hour, Treatment = data$strain))
colnames(exprsData) <- rownames(phenoData)
eset <- new("ExpressionSet", exprs = exprsData, phenoData = phenoData)

# Read the gene mapping file
gene_mapping <- read.csv("../entrez_annot.csv")
mapping_dict <- setNames(gene_mapping$`Entrez.ID`, gene_mapping$Gene)

# Analysis for different treatments
treatments <- list(c("AX4", "B1"), c("AX4", "C1"), c("AX4", "rgB"), c("AX4", "B1_rgB"), c("AX4", "C1_rgB"), c("AX4", "AX4L846F"), c("AX4", "B1_L846F"))
file_prefixes <- c("DE_B1", "DE_C1", "DE_rgB", "DE_B1_rgB", "DE_C1_rgB", "DE_AX4L846F", "DE_B1_L846F")

for (i in seq_along(treatments)) {
  file_prefix <- file_prefixes[i]
  perform_analysis(eset, treatments[[i]], geneSets, mapping_dict, file_prefix)
  if (!dir.exists("DE")) {
    dir.create("DE")
  }
  calculate_regulation(eset, treatments[[i]][1], treatments[[i]][2], paste0("genes/", file_prefix, ".csv"))
}

