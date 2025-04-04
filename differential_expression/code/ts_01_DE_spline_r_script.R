library(splineTimeR)
library(Biobase)

# Function to read and preprocess expression data
data <- function(filename) {
  data <- read.csv(filename, row.names = 1)
  exprsData <- t(as.matrix(data[, 5:ncol(data)]))
  return(exprsData)
}

# Function for differential expression analysis and enrichment analysis
perform_analysis <- function(eset, treatment, geneSets, mapping_dict, file_prefix) {
  selected_samples <- pData(eset)$Treatment %in% treatment
  eset_selected <- eset[, selected_samples]
  diffExprs <- splineDiffExprs(eSetObject = eset_selected, df = 3, cutoff.adj.pVal = 1, reference = treatment[1], intercept = TRUE)
  
  diffExprs <- diffExprs[order(diffExprs$adj.P.Val), ]
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
    enrichPath <- pathEnrich(geneList = rownames(mapped_diffExprs), geneSets = geneSets[[geneSet]], universe = 12827)
    if (!dir.exists("gsea")) {
      dir.create("gsea")
    }
    write.csv(enrichPath, paste0("gsea/", file_prefix, "_", geneSet, "_gsea.csv"), row.names = TRUE)
  }
}


calculate_regulation <- function(eset, treatment1, treatment2, file_name) {
  avgExprs_Treatment1 <- rowMeans(exprs(eset[, pData(eset)$Treatment == treatment1]))
  avgExprs_Treatment2 <- rowMeans(exprs(eset[, pData(eset)$Treatment == treatment2]))
  regulation <- ifelse(avgExprs_Treatment2 > avgExprs_Treatment1, "Upregulated", "Downregulated")
  regulation_df <- data.frame(row_IDs = rownames(eset), Regulation = regulation)
  
  existing_data <- read.csv(file_name)
  merged_data <- merge(existing_data, regulation_df, by = "row_IDs", all.x = TRUE)
  merged_data <- merged_data[order(merged_data$adj.P.Val), ]
  write.csv(merged_data, file_name, row.names = FALSE)
}


# Set working directory and read data
setwd("/Users/lenatrnovec/FAKS/BIOLAB/RNAseq_dicty/differential_expression/results/time_series") # Change to your working directory
exprsData <- preprocess_data("../../data/combined_df.csv")

# Load gene sets
geneSets <- list(
  BP = getGmt("../../data/genesets/GO-biological_process_mod.gmt"),
  MF = getGmt("../../data/genesets/GO-molecular_function_mod.gmt"),
  DB = getGmt("../../data/genesets/DictyBase-Phenotypes_mod.gmt")
)

# Create ExpressionSet
data <- read.csv("../../data/combined_df.csv", row.names = 1)
data$sampleName <- rownames(data)
phenoData <- new("AnnotatedDataFrame", data = data.frame(SampleName = data$sampleName, replicate = data$replicate, Time = data$hour, Treatment = data$strain))
colnames(exprsData) <- rownames(phenoData)
eset <- new("ExpressionSet", exprs = exprsData, phenoData = phenoData)

# Read the gene mapping file
gene_mapping <- read.csv("../../data/entrez_annot.csv")
mapping_dict <- setNames(gene_mapping$`Entrez.ID`, gene_mapping$Gene)

# Analysis for different treatments
treatments <- list(c("AX4", "B1"), c("AX4", "C1"), c("AX4", "rgB"), c("AX4", "B1_rgB"), c("AX4", "AX4L846F"))
# treatments <- list(c("B1_rgB", "B1"), c("B1_rgB", "C1"), c("B1_rgB", "rgB"), c("B1_rgB", "AX4"), c("B1_rgB", "AX4L846F"))
# treatments <- list(c("AX4L846F", "B1"), c("AX4L846F", "C1"), c("AX4L846F", "rgB"), c("AX4L846F", "B1_rgB"), c("AX4L846F", "AX4"))
# treatments <- list(c("B1", "AX4"), c("B1", "C1"), c("B1", "rgB"), c("B1", "B1_rgB"), c("B1", "AX4L846F"))
# treatments <- list(c("rgB", "B1"), c("rgB", "C1"), c("rgB", "AX4"), c("rgB", "B1_rgB"), c("rgB", "AX4L846F"))
# treatments <- list(c("C1", "B1"), c("C1", "AX4"), c("C1", "rgB"), c("C1", "B1_rgB"), c("C1", "AX4L846F"))

# treatments <- list(c("AX4", "B1"), c("AX4", "C1"), c("AX4", "rgB"), c("AX4", "B1_rgB"), c("AX4", "AX4L846F"), c("B1_rgB", "B1"), c("B1_rgB", "C1"), c("B1_rgB", "rgB"), c("B1_rgB", "AX4"), c("B1_rgB", "AX4L846F"), c("AX4L846F", "B1"), c("AX4L846F", "C1"), c("AX4L846F", "rgB"), c("AX4L846F", "B1_rgB"), c("AX4L846F", "AX4"), c("B1", "AX4"), c("B1", "C1"), c("B1", "rgB"), c("B1", "B1_rgB"), c("B1", "AX4L846F"), c("rgB", "B1"), c("rgB", "C1"), c("rgB", "AX4"), c("rgB", "B1_rgB"), c("rgB", "AX4L846F"), c("C1", "B1"), c("C1", "AX4"), c("C1", "rgB"), c("C1", "B1_rgB"), c("C1", "AX4L846F"))

file_prefixes <- c("DE_B1", "DE_C1", "DE_rgB","DE_B1_rgB", "DE_AX4L846F")
# file_prefixes <- c("DE_B1", "DE_C1", "DE_rgB","DE_AX4", "DE_AX4L846F")
# file_prefixes <- c("DE_B1", "DE_C1", "DE_rgB","DE_B1_rgB", "DE_AX4")
# file_prefixes <- c("DE_AX4", "DE_C1", "DE_rgB","DE_B1_rgB", "DE_AX4L846F")
# file_prefixes <- c("DE_B1", "DE_C1", "DE_AX4","DE_B1_rgB", "DE_AX4L846F")
# file_prefixes <- c("DE_B1", "DE_AX4", "DE_rgB","DE_B1_rgB", "DE_AX4L846F")

# file_prefixes <- c("DE_AX4-B1", "DE_AX4-C1", "DE_AX4-rgB","DE_AX4-B1_rgB", "DE_AX4-AX4L846F", "DE_B1_rgB-B1", "DE_B1_rgB-C1", "DE_B1_rgB-rgB","DE_B1_rgB-AX4", "DE_B1_rgB-AX4L846F", "DE_AX4L846F-B1", "DE_AX4L846F-C1", "DE_AX4L846F-rgB","DE_AX4L846F-B1_rgB", "DE_AX4L846F-AX4", "DE_B1-AX4", "DE_B1-C1", "DE_B1-rgB","DE_B1-B1_rgB", "DE_B1-AX4L846F", "DE_rgB-B1", "DE_rgB-C1", "DE_rgB-AX4","DE_rgB-B1_rgB", "DE_rgB-AX4L846F", "DE_C1-B1", "DE_C1-AX4", "DE_C1-rgB","DE_C1-B1_rgB", "DE_C1-AX4L846F")

if (!dir.exists("genes")) {
  dir.create("genes")
}

for (i in seq_along(treatments)) {
  file_prefix <- file_prefixes[i]
  perform_analysis(eset, treatments[[i]], geneSets, mapping_dict, file_prefix)
  calculate_regulation(eset, treatments[[i]][1], treatments[[i]][2], paste0("genes/", file_prefix, ".csv"))
}

