# ===========================================
# GDSC2 gemcitabine and RNA-seq data processing
# Pancreatic target and four candidate sources
# ===========================================

# Packages
source("00_requirements.R")

# Input and output files
response_file = "data/raw/GDSC2_fitted_dose_response_27Oct23.xlsx"
model_file = "data/raw/model_list_20261013.csv"
tpm_file = "data/raw/rnaseq_merged_rsem_tpm_20260323.csv"

# Keep the gemcitabine drug responses.
response = as.data.frame(read_excel(response_file))
response = response[response$DATASET == "GDSC2" &  response$DRUG_NAME == "Gemcitabine", c("SANGER_MODEL_ID", "LN_IC50")]
response = response[!is.na(response$SANGER_MODEL_ID) &
                    !is.na(response$LN_IC50), ]
stopifnot(!anyDuplicated(response$SANGER_MODEL_ID))

# Identify the cancer groups using the model list.
model_list = read.csv(model_file, stringsAsFactors = FALSE)
model_list = model_list[, c("model_id", "model_name", "cancer_type")]

samples = merge(response, model_list, by.x = "SANGER_MODEL_ID", by.y = "model_id")

cancer_types = c("Pancreatic Carcinoma", "Colorectal Carcinoma",
                 "Gastric Carcinoma", "Non-Small Cell Lung Carcinoma",
                 "Ovarian Carcinoma")
group_names = c("Pancreatic", "Colorectal", "Gastric", "Non-small-cell lung", "Ovarian")

samples = samples[samples$cancer_type %in% cancer_types, ]
samples$group = group_names[match(samples$cancer_type, cancer_types)]
samples$role = ifelse(samples$group == "Pancreatic", "Target", "Source")

# Match each model ID to a column in the TPM file.
# The first TPM row lists model IDs. The next three rows are metadata.
tpm_header = strsplit(readLines(tpm_file, n = 1), ",", fixed = TRUE)[[1]]
samples$in_expression = samples$SANGER_MODEL_ID %in% tpm_header

# Save the matching decisions, including excluded models.
audit = samples[, c("SANGER_MODEL_ID", "model_name", "cancer_type",
                    "group", "role", "LN_IC50", "in_expression")]
names(audit)[1] = "model_id"
write.csv(audit, "data/gdsc2_gemcitabine_merge_audit.csv", row.names = FALSE)

counts = data.frame()
for (g in group_names) 
{
  group_rows = samples[samples$group == g, ]
  counts = rbind(counts, data.frame(
    group = g,
    role = ifelse(g == "Pancreatic", "Target", "Source"),
    gemcitabine_responses = nrow(group_rows),
    matched_tpm_models = sum(group_rows$in_expression),
    missing_tpm = sum(!group_rows$in_expression)))
}
write.csv(counts, "data/gdsc2_merge_counts.csv", row.names = FALSE)

samples = samples[samples$in_expression, ]
samples$group = factor(samples$group, levels = group_names)
samples = samples[order(samples$group, samples$SANGER_MODEL_ID), ]
samples$group = as.character(samples$group)

sample_info = samples[, c("SANGER_MODEL_ID", "model_name", "cancer_type",
                          "group", "role", "LN_IC50")]
names(sample_info)[1] = "model_id"
write.csv(sample_info, "data/gdsc2_gemcitabine_samples.csv", row.names = FALSE)

# Read only the matching model columns from the large TPM file.
positions = match(samples$SANGER_MODEL_ID, tpm_header)
columns_to_read = sort(c(1, 2, 3, positions))
expr = fread(tpm_file, skip = 4, header = FALSE, select = columns_to_read, na.strings = c("", "NA"))

# Name the first three gene columns and the selected cell-model columns.
selected_header = tpm_header[columns_to_read]
selected_header[1:3] = c("gene_symbol", "ensembl_gene_id", "gene_id")
setnames(expr, selected_header)

tpm = as.matrix(expr[, -c(1, 2, 3), with = FALSE])
storage.mode(tpm) = "double"
tpm = tpm[, samples$SANGER_MODEL_ID, drop = FALSE]

# Remove genes with missing, non-finite, or negative TPM values.
n_all_genes = nrow(expr)
valid_genes = rowSums(!is.finite(tpm)) == 0 &
              rowSums(tpm < 0, na.rm = TRUE) == 0 &
              !is.na(expr$ensembl_gene_id) & expr$ensembl_gene_id != "" &
              !is.na(expr$gene_id) & expr$gene_id != ""

n_valid = sum(valid_genes)
expr = expr[valid_genes, ]
tpm = tpm[valid_genes, , drop = FALSE]

# Select features using only the external source samples.
source_cols = which(samples$role == "Source")
min_samples = ceiling(0.20 * length(source_cols))
expressed = rowSums(tpm[, source_cols, drop = FALSE] > 1) >= min_samples

n_expressed = sum(expressed)
expr = expr[expressed, ]
tpm = tpm[expressed, , drop = FALSE]
log_tpm = log2(tpm + 1)

# Rank genes by variance of log2(TPM + 1) in source samples only.
source_variance = apply(log_tpm[, source_cols, drop = FALSE], 1, var)
ordered_genes = order(-source_variance, expr$ensembl_gene_id)
top_genes = ordered_genes[1:2000]

selected_genes = data.frame(
  gene_symbol = expr$gene_symbol[top_genes],
  ensembl_gene_id = expr$ensembl_gene_id[top_genes],
  gene_id = expr$gene_id[top_genes],
  external_log2_tpm_variance = source_variance[top_genes],
  rank = seq_along(top_genes))
write.csv(selected_genes, "data/gdsc2_selected_genes.csv", row.names = FALSE)

# Save a single matrix containing responses and 2,000 predictors.
rownames(log_tpm) = expr$ensembl_gene_id
X = t(log_tpm[top_genes, , drop = FALSE])
cleaned = data.frame(model_id = samples$SANGER_MODEL_ID,
                     group = samples$group,
                     role = samples$role,
                     LN_IC50 = samples$LN_IC50,
                     X, check.names = FALSE)
write.csv(cleaned, "data/gdsc2_gemcitabine_2000genes.csv", row.names = FALSE)

# Save the gene-filtering totals.
filter_counts = data.frame(
  step = c("all_genes", "valid_TPM_in_all_229_models",
           "TPM_above_1_in_at_least_40_sources", "selected_genes"),
  n = c(n_all_genes, n_valid, n_expressed, nrow(selected_genes)))
write.csv(filter_counts, "data/gdsc2_gene_filter_counts.csv", row.names = FALSE)

print(counts)
print(filter_counts)
cat("Saved", nrow(cleaned), "models with", ncol(cleaned) - 4, "gene-expression predictors.\n")
