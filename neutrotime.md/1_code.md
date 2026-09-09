# neutrotime preprocessing and QC

- run with module load R4.6.1 on frontenac CAC HPC
- https://github.com/GITC2025/neutrophils_2026/tree/main/Renv.md 

```sh
mkdir neutrotime
cd neutrotime
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165276/suppl/GSE165276_RAW.tar"
tar -xvf GSE165276_RAW.tar 
```
- we only focus on the healthy (non inflammation) datasets now 
```
[hpc6297@frnt153 neutrotime]$ ls -lht
total 44M
-rw-r-----. 1 hpc6297 hpc6297   22M Feb 12  2021 GSE165276_RAW.tar
-rw-r-----. 1 hpc6297 hpc6297  137K Jan 21  2021 GSM5029341_inflammation_dataset3_readme.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  7.4M Jan 21  2021 GSM5029341_inflammation_dataset3.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  4.2M Jan 21  2021 GSM5029340_SP_dataset2.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  4.0M Jan 21  2021 GSM5029339_BM_dataset2.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  3.0M Jan 21  2021 GSM5029338_BL_dataset2.txt.gz
-rw-r-----. 1 hpc6297 hpc6297 1008K Jan 21  2021 GSM5029337_SP_dataset1.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  1.4M Jan 21  2021 GSM5029336_BM_dataset1.txt.gz
-rw-r-----. 1 hpc6297 hpc6297  1.1M Jan 21  2021 GSM5029335_BL_dataset1.txt.gz
```

# check txtgz formats
```r
options(width = 800)
input_dir <- "/global/scratch/hpc6297/original_neutrotime"

files_to_check <- c(
"GSM5029335_BL_dataset1.txt.gz",
"GSM5029338_BL_dataset2.txt.gz"
)

for (fn in files_to_check) {
con <- gzfile(file.path(input_dir, fn), "rt")
header_line <- readLines(con, n = 1)
close(con)

raw_headers <- unlist(strsplit(header_line, "\t"))

cat("file:", fn, "\n")
cat("total columns (including gene col):", length(raw_headers), "\n")
cat("first 10 column names:\n")
print(head(raw_headers, 10))
cat("unique column names in header:", length(unique(raw_headers[-1])), "\n")
cat("summary table of header values:\n")
print(table(head(raw_headers[-1], 20)))
cat("\n")
}
```

- problem with dataset 2 headers, lack the "", has additional -1 suffix
```
total columns (including gene col): 1589 
first 10 column names:
 [1] ""                 "AAACCTGAGTCCATAC" "AAACCTGAGTCTCCTC" "AAACCTGAGTGTGAAT" "AAACGGGCAAATTGCC" "AAACGGGCAGGTCGTC" "AAACGGGTCGGAAATA" "AAACGGGTCTTAGAGC" "AAACGGGTCTTGTCAT" "AAAGATGCAACGCACC"
unique column names in header: 1588 
summary table of header values:

AAACCTGAGTCCATAC AAACCTGAGTCTCCTC AAACCTGAGTGTGAAT AAACGGGCAAATTGCC AAACGGGCAGGTCGTC AAACGGGTCGGAAATA AAACGGGTCTTAGAGC AAACGGGTCTTGTCAT AAAGATGCAACGCACC AAAGATGCAGTCAGAG AAAGCAAAGTAGGTGC AAAGCAAAGTTAAGTG AAAGCAACATTGAGCT AAAGCAAGTAACGCGA AAAGCAAGTACAGTTC AAAGCAATCATCGGAT AAAGTAGGTAAGTGTA AAAGTAGGTCTAAAGA AAATGCCCAGGCTGAA AAATGCCTCAAACGGG 
               1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1 

file: GSM5029338_BL_dataset2.txt.gz 
total columns (including gene col): 4653 
first 10 column names:
 [1] "\"AAACCTGCAATGAAAC-1\"" "\"AAACCTGCAGTCAGAG-1\"" "\"AAACCTGCATGGGAAC-1\"" "\"AAACCTGGTCGGATCC-1\"" "\"AAACCTGGTCTCTCTG-1\"" "\"AAACCTGTCGGCGGTT-1\"" "\"AAACGGGCACATTTCT-1\"" "\"AAACGGGCAGGTGCCT-1\"" "\"AAACGGGCATAAGACA-1\"" "\"AAACGGGCATGGGAAC-1\""
unique column names in header: 4652 
summary table of header values:

"AAACCTGCAGTCAGAG-1" "AAACCTGCATGGGAAC-1" "AAACCTGGTCGGATCC-1" "AAACCTGGTCTCTCTG-1" "AAACCTGTCGGCGGTT-1" "AAACGGGCACATTTCT-1" "AAACGGGCAGGTGCCT-1" "AAACGGGCATAAGACA-1" "AAACGGGCATGGGAAC-1" "AAACGGGCATTAGCCA-1" "AAACGGGGTTCTGAAC-1" "AAACGGGTCAGTTGAC-1" "AAACGGGTCTGACCTC-1" "AAAGATGAGGCTAGCA-1" "AAAGATGAGGCTCATT-1" "AAAGATGAGTATGACA-1" "AAAGATGCATATACGC-1" "AAAGATGGTCCGTGAC-1" "AAAGATGGTCTAGTGT-1" "AAAGATGGTGTTTGTG-1" 
                   1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1                    1 
```

# fix dataset 2 barcode format
```r
options(width = 800)
data_dir <- "/global/scratch/hpc6297/original_neutrotime"

# read reference header from dataset 1
ref_fn <- "GSM5029335_BL_dataset1.txt.gz"
con_ref <- gzfile(file.path(data_dir, ref_fn), "rt")
ref_header_line <- readLines(con_ref, n = 1)
close(con_ref)
ref_tokens <- unlist(strsplit(ref_header_line, "\t"))

dataset2_files <- c(
"GSM5029338_BL_dataset2.txt.gz",
"GSM5029339_BM_dataset2.txt.gz",
"GSM5029340_SP_dataset2.txt.gz"
)

for (fn in dataset2_files) {
orig_path <- file.path(data_dir, fn)
bak_path <- file.path(data_dir, sub("\\.txt\\.gz$", "_rawbackup.txt.gz", fn))

# create backup before modifying
if (!file.exists(bak_path)) {
file.copy(orig_path, bak_path)
cat("backup created:", basename(bak_path), "\n")
}

# inspect current header directly from backup file
con_in <- gzfile(bak_path, "rt")
orig_header_line <- readLines(con_in, n = 1)

orig_tokens <- unlist(strsplit(orig_header_line, "\t"))

# clean quotes, whitespace, and trailing -1
clean_barcodes <- sub("-1$", "", gsub('^["\']|["\']$', "", trimws(orig_tokens)))

# remove leading empty element if file was already partially formatted
if (clean_barcodes[1] == "") {
clean_barcodes <- clean_barcodes[-1]
}

# prepend empty tab delimiter for gene column alignment
new_header_line <- paste(c("", clean_barcodes), collapse = "\t")

# stream write fixed file via temporary file
tmp_out <- tempfile(pattern = "fix_", tmpdir = data_dir, fileext = ".txt.gz")
con_out <- gzfile(tmp_out, "wt")

# write corrected header
writeLines(new_header_line, con_out)

# copy remaining count lines directly
while (length(lines <- readLines(con_in, n = 10000)) > 0) {
writeLines(lines, con_out)
}
close(con_in)
close(con_out)

# overwrite original target file atomically
file.rename(tmp_out, orig_path)

# inspect new header for verification against dataset 1
con_new <- gzfile(orig_path, "rt")
new_header_line_read <- readLines(con_new, n = 1)
close(con_new)
new_tokens <- unlist(strsplit(new_header_line_read, "\t"))

cat("\nfile:", fn, "\n")
cat("reference file (dataset 1)   :", ref_fn, "\n")
cat("dataset 1 first 4 tokens     :", paste(head(ref_tokens, 4), collapse = " | "), "\n")
cat("dataset 2 fixed first 4 tokens:", paste(head(new_tokens, 4), collapse = " | "), "\n")
cat("leading empty token matched  :", (ref_tokens[1] == new_tokens[1]), "\n")
cat("barcode length matched (16bp):", (nchar(new_tokens[2]) == 16 && nchar(ref_tokens[2]) == 16), "\n\n")
}
```

```
backup created: GSM5029338_BL_dataset2_rawbackup.txt.gz 

file: GSM5029338_BL_dataset2.txt.gz 
reference file (dataset 1)   : GSM5029335_BL_dataset1.txt.gz 
dataset 1 first 4 tokens     :  | AAACCTGAGTCCATAC | AAACCTGAGTCTCCTC | AAACCTGAGTGTGAAT 
dataset 2 fixed first 4 tokens:  | AAACCTGCAATGAAAC | AAACCTGCAGTCAGAG | AAACCTGCATGGGAAC 
leading empty token matched  : TRUE 
barcode length matched (16bp): TRUE 

backup created: GSM5029339_BM_dataset2_rawbackup.txt.gz 

file: GSM5029339_BM_dataset2.txt.gz 
reference file (dataset 1)   : GSM5029335_BL_dataset1.txt.gz 
dataset 1 first 4 tokens     :  | AAACCTGAGTCCATAC | AAACCTGAGTCTCCTC | AAACCTGAGTGTGAAT 
dataset 2 fixed first 4 tokens:  | AAACCTGAGAGGACGG | AAACCTGAGCCACGCT | AAACCTGCATCAGTAC 
leading empty token matched  : TRUE 
barcode length matched (16bp): TRUE 

backup created: GSM5029340_SP_dataset2_rawbackup.txt.gz 

file: GSM5029340_SP_dataset2.txt.gz 
reference file (dataset 1)   : GSM5029335_BL_dataset1.txt.gz 
dataset 1 first 4 tokens     :  | AAACCTGAGTCCATAC | AAACCTGAGTCTCCTC | AAACCTGAGTGTGAAT 
dataset 2 fixed first 4 tokens:  | AAACCTGAGATCCGAG | AAACCTGAGCTGTCTA | AAACCTGCACAACGCC 
leading empty token matched  : TRUE 
barcode length matched (16bp): TRUE 
```

# compare headers for dataset 1 and 2 again
```
file: GSM5029335_BL_dataset1.txt.gz 
total columns (including gene col): 1589 
first 10 column names:
 [1] ""                 "AAACCTGAGTCCATAC" "AAACCTGAGTCTCCTC" "AAACCTGAGTGTGAAT" "AAACGGGCAAATTGCC" "AAACGGGCAGGTCGTC" "AAACGGGTCGGAAATA" "AAACGGGTCTTAGAGC" "AAACGGGTCTTGTCAT" "AAAGATGCAACGCACC"
unique column names in header: 1588 
summary table of header values:

AAACCTGAGTCCATAC AAACCTGAGTCTCCTC AAACCTGAGTGTGAAT AAACGGGCAAATTGCC AAACGGGCAGGTCGTC AAACGGGTCGGAAATA AAACGGGTCTTAGAGC AAACGGGTCTTGTCAT AAAGATGCAACGCACC AAAGATGCAGTCAGAG AAAGCAAAGTAGGTGC AAAGCAAAGTTAAGTG AAAGCAACATTGAGCT AAAGCAAGTAACGCGA AAAGCAAGTACAGTTC AAAGCAATCATCGGAT AAAGTAGGTAAGTGTA AAAGTAGGTCTAAAGA AAATGCCCAGGCTGAA AAATGCCTCAAACGGG 
               1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1 

file: GSM5029338_BL_dataset2.txt.gz 
total columns (including gene col): 4654 
first 10 column names:
 [1] ""                 "AAACCTGCAATGAAAC" "AAACCTGCAGTCAGAG" "AAACCTGCATGGGAAC" "AAACCTGGTCGGATCC" "AAACCTGGTCTCTCTG" "AAACCTGTCGGCGGTT" "AAACGGGCACATTTCT" "AAACGGGCAGGTGCCT" "AAACGGGCATAAGACA"
unique column names in header: 4653 
summary table of header values:

AAACCTGCAATGAAAC AAACCTGCAGTCAGAG AAACCTGCATGGGAAC AAACCTGGTCGGATCC AAACCTGGTCTCTCTG AAACCTGTCGGCGGTT AAACGGGCACATTTCT AAACGGGCAGGTGCCT AAACGGGCATAAGACA AAACGGGCATGGGAAC AAACGGGCATTAGCCA AAACGGGGTTCTGAAC AAACGGGTCAGTTGAC AAACGGGTCTGACCTC AAAGATGAGGCTAGCA AAAGATGAGGCTCATT AAAGATGAGTATGACA AAAGATGCATATACGC AAAGATGGTCCGTGAC AAAGATGGTCTAGTGT 
               1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1                1 
```


# merge all matrices with metadata into merged preQC rds
```r
options(width = 800)
input_dir <- "/global/scratch/hpc6297/original_neutrotime"
output_dir <- "/global/scratch/hpc6297/neutrotime_output"

if (!dir.exists(output_dir)) {
dir.create(output_dir, recursive = TRUE)
}

library(data.table)
library(Matrix)
library(BiocParallel)
library(Seurat)

n_cores <- 8
bpparam <- MulticoreParam(workers = n_cores)

files_preqc <- c(
"GSM5029335_BL_dataset1.txt.gz",
"GSM5029336_BM_dataset1.txt.gz",
"GSM5029337_SP_dataset1.txt.gz",
"GSM5029338_BL_dataset2.txt.gz",
"GSM5029339_BM_dataset2.txt.gz",
"GSM5029340_SP_dataset2.txt.gz"
)

# read individual files, capture gene rows, construct cell IDs
read_preqc_sample <- function(fn) {
file_path <- file.path(input_dir, fn)
dt <- fread(file_path, header = TRUE, sep = "\t", check.names = FALSE)
gene_names <- dt[[1]]
raw_barcodes <- colnames(dt)[-1]
mat <- as.matrix(dt[, -1, with = FALSE])
rownames(mat) <- gene_names

prefix <- sub("\\.txt\\.gz$", "", fn)
cell_ids <- paste0(prefix, "_", raw_barcodes)
colnames(mat) <- cell_ids

sparse_mat <- as(mat, "CsparseMatrix")
list(mat = sparse_mat, genes = gene_names, raw_barcodes = raw_barcodes, sample_id = fn)
}

sample_data <- bplapply(files_preqc, read_preqc_sample, BPPARAM = bpparam)
names(sample_data) <- files_preqc

# build master gene universe across all datasets
master_genes <- sort(unique(unlist(lapply(sample_data, function(x) x$genes))))

# align rows to master gene set with sparse zero-padding
align_matrix_genes <- function(item, all_genes) {
mat <- item$mat
missing_genes <- setdiff(all_genes, rownames(mat))
if (length(missing_genes) > 0) {
zero_mat <- Matrix(0, nrow = length(missing_genes), ncol = ncol(mat), sparse = TRUE)
rownames(zero_mat) <- missing_genes
colnames(zero_mat) <- colnames(mat)
mat <- rbind(mat, zero_mat)
}
mat[all_genes, ]
}

aligned_matrices <- bplapply(sample_data, align_matrix_genes, all_genes = master_genes, BPPARAM = bpparam)
merged_counts <- do.call(cbind, aligned_matrices)

# assemble metadata aligned 1:1 to matrix columns
meta_dt_list <- lapply(names(sample_data), function(fn) {
item <- sample_data[[fn]]
mat <- item$mat

site_val <- fcase(
grepl("_SP_", fn), "spleen",
grepl("_BM_", fn), "bone marrow",
grepl("_BL_", fn), "peripheral blood",
default = NA_character_
)

dataset_val <- fcase(
grepl("dataset1", fn), "1",
grepl("dataset2", fn), "2",
default = NA_character_
)

data.table(
cell_id = colnames(mat),
raw_barcode = item$raw_barcodes,
sample_id = fn,
dataset = dataset_val,
site = site_val,
strain = "C57BL/6J",
sex = "male",
health_status = "healthy"
)
})

complete_metadata <- as.data.frame(rbindlist(meta_dt_list))
rownames(complete_metadata) <- complete_metadata$cell_id

# verify 1:1 alignment
stopifnot(identical(rownames(complete_metadata), colnames(merged_counts)))

# instantiate standard seurat container
seu <- CreateSeuratObject(
counts = merged_counts,
meta.data = complete_metadata,
project = "neutrotime"
)

# compute feature-level metadata (gene-level stats)
seu[["RNA"]] <- AddMetaData(
seu[["RNA"]],
metadata = data.frame(
gene_symbol = rownames(seu),
n_cells_expressing = rowSums(GetAssayData(seu, layer = "counts") > 0),
total_counts = rowSums(GetAssayData(seu, layer = "counts")),
row.names = rownames(seu),
stringsAsFactors = FALSE
)
)

# verify canonical mouse neutrophil markers
key_markers <- c("Ly6g", "S100a8", "S100a9", "Itgam", "Mpo", "Elane", "Actb", "Gapdh")
present_markers <- key_markers[key_markers %in% rownames(seu)]

cat("marker validation:\n")
cat("expected markers found:", paste(present_markers, collapse = ", "), "\n\n")

# atomic save uncompressed rds helper
atomic_save_uncompressed_rds <- function(object, file) {
temp_file <- paste0(file, ".tmp.", Sys.getpid())
saveRDS(object, file = temp_file, compress = FALSE)
if (!file.rename(temp_file, file)) {
file.copy(temp_file, file, overwrite = TRUE)
unlink(temp_file)
}
}

out_seurat_rds <- file.path(output_dir, "neutrotime_preQC_merged.rds")
atomic_save_uncompressed_rds(seu, file = out_seurat_rds)

cat("saved uncompressed seurat rds:\n")
cat(out_seurat_rds, "\n\n")

cat("seurat object summary:\n")
print(seu)
cat("\n")

cat("metadata (first 6 rows):\n")
print(head(seu@meta.data))
cat("\n")

cat("feature metadata (first 6 rows):\n")
print(head(seu[["RNA"]][[]]))
```

```
marker validation:
expected markers found: Ly6g, S100a8, S100a9, Itgam, Mpo, Elane, Actb, Gapdh 

saved uncompressed seurat rds:
/global/scratch/hpc6297/neutrotime_output/neutrotime_preQC_merged.rds 

seurat object summary:
An object of class Seurat 
28692 features across 16907 samples within 1 assay 
Active assay: RNA (28692 features, 0 variable features)
 1 layer present: counts

metadata (first 6 rows):
                                        orig.ident nCount_RNA nFeature_RNA                                 cell_id      raw_barcode                     sample_id dataset             site   strain  sex health_status
GSM5029335_BL_dataset1_AAACCTGAGTCCATAC GSM5029335        835          349 GSM5029335_BL_dataset1_AAACCTGAGTCCATAC AAACCTGAGTCCATAC GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy
GSM5029335_BL_dataset1_AAACCTGAGTCTCCTC GSM5029335       3093          691 GSM5029335_BL_dataset1_AAACCTGAGTCTCCTC AAACCTGAGTCTCCTC GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy
GSM5029335_BL_dataset1_AAACCTGAGTGTGAAT GSM5029335       2857          785 GSM5029335_BL_dataset1_AAACCTGAGTGTGAAT AAACCTGAGTGTGAAT GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy
GSM5029335_BL_dataset1_AAACGGGCAAATTGCC GSM5029335        548          306 GSM5029335_BL_dataset1_AAACGGGCAAATTGCC AAACGGGCAAATTGCC GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy
GSM5029335_BL_dataset1_AAACGGGCAGGTCGTC GSM5029335        748          342 GSM5029335_BL_dataset1_AAACGGGCAGGTCGTC AAACGGGCAGGTCGTC GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy
GSM5029335_BL_dataset1_AAACGGGTCGGAAATA GSM5029335       1235          588 GSM5029335_BL_dataset1_AAACGGGTCGGAAATA AAACGGGTCGGAAATA GSM5029335_BL_dataset1.txt.gz       1 peripheral blood C57BL/6J male       healthy

feature metadata (canonical markers):
       gene_symbol n_cells_expressing total_counts
Ly6g          Ly6g               4874        17622
S100a8      S100a8              16440      4496585
S100a9      S100a9              16292      2538791
Itgam        Itgam               4354         6355
Mpo            Mpo                126         2614
Elane        Elane                118         4504
Actb          Actb              16286       403430
Gapdh        Gapdh               6664        13572
```

# QC view on merged preQC rds
```r
options(width = 800)
rds_path <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_preQC_merged.rds"
output_dir <- "/global/scratch/hpc6297/neutrotime_output"

library(data.table)
library(Matrix)
library(Seurat)

seu <- readRDS(rds_path)
counts_layer <- GetAssayData(seu, layer = "counts")

# identify mitochondrial and ribosomal features globally
mito_genes <- grep("^mt-", rownames(counts_layer), value = TRUE, ignore.case = TRUE)
ribo_genes <- grep("^(rpl|rps|mrpl|mrps)", rownames(counts_layer), value = TRUE, ignore.case = TRUE)

# construct cell qc data frame from seurat object
qc_df <- data.frame(
sample_id = seu$sample_id,
cell_id = colnames(seu),
ncount = seu$nCount_RNA,
nfeature = seu$nFeature_RNA,
percent_mito = as.numeric((colSums(counts_layer[mito_genes, , drop = FALSE]) / pmax(seu$nCount_RNA, 1)) * 100),
percent_ribo = as.numeric((colSums(counts_layer[ribo_genes, , drop = FALSE]) / pmax(seu$nCount_RNA, 1)) * 100),
stringsAsFactors = FALSE
)

hi_res_probs <- c(0.01, 0.02, 0.03, 0.04, 0.96, 0.97, 0.98, 0.99)
hi_res_cols <- paste0("q_", c("01", "02", "03", "04", "96", "97", "98", "99"))

# helper to build full stats row
calc_full_stats <- function(df, scope_name) {
metrics <- c("ncount", "nfeature", "percent_mito", "percent_ribo")
out_list <- list()
for (m in metrics) {
vals <- df[[m]]
sum_vals <- summary(vals)
q_vals <- quantile(vals, probs = hi_res_probs, na.rm = TRUE)
row_dt <- data.table(
scope = scope_name,
metric = m,
cells = length(vals),
min = as.numeric(sum_vals["Min."]),
q25 = as.numeric(sum_vals["1st Qu."]),
median = as.numeric(sum_vals["Median"]),
mean = as.numeric(sum_vals["Mean"]),
q75 = as.numeric(sum_vals["3rd Qu."]),
max = as.numeric(sum_vals["Max."])
)
for (i in seq_along(hi_res_cols)) {
set(row_dt, j = hi_res_cols[i], value = unname(q_vals[i]))
}
out_list[[m]] <- row_dt
}
rbindlist(out_list)
}

# helper to print formatted console block
print_block <- function(label, n_genes, n_cells, mito_len, ribo_len, df) {
cat(label, "\n\n")
if (!is.na(n_genes)) {
cat("genes:", n_genes, "\n")
}
cat("cells:", n_cells, "\n\n")
if (!is.na(mito_len)) {
cat("identified mitochondrial genes:", mito_len, "\n")
cat("identified ribosomal genes:", ribo_len, "\n\n")
}
cat("qc metric distributions across cells:\n")
print(summary(df[, c("ncount", "nfeature", "percent_mito", "percent_ribo")]))
cat("\n")
cat("quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):\n")
q_mat <- apply(df[, c("ncount", "nfeature", "percent_mito", "percent_ribo")], 2, quantile, probs = hi_res_probs, na.rm = TRUE)
print(round(q_mat, 3))
cat("\n\n")
}

samples <- sort(unique(qc_df$sample_id))
summary_list <- list()
dataset_meta <- list()

for (s in samples) {
sub_df <- qc_df[qc_df$sample_id == s, ]
sub_cells <- sub_df$cell_id
sub_counts <- counts_layer[, sub_cells, drop = FALSE]

# count expressed genes in sample
expressed_genes <- sum(rowSums(sub_counts > 0) > 0)

summary_list[[s]] <- calc_full_stats(sub_df, s)
dataset_meta[[s]] <- list(
genes = expressed_genes,
cells = nrow(sub_df),
mito = length(mito_genes),
ribo = length(ribo_genes)
)
}

# global dataset stats
global_summary <- calc_full_stats(qc_df, "global_all_datasets")

# combine and export summary tsv
complete_summary <- rbindlist(c(list(global_summary), summary_list))
num_cols <- c("min", "q25", "median", "mean", "q75", "max", hi_res_cols)
complete_summary[, (num_cols) := lapply(.SD, function(x) round(x, 3)), .SDcols = num_cols]

out_tsv <- file.path(output_dir, "neutrotime_preQC_metrics_summary.tsv")
fwrite(complete_summary, file = out_tsv, sep = "\t", quote = FALSE)
cat("saved qc metrics summary:", out_tsv, "\n\n")

# print global console block
print_block(
label = "global (all samples combined)",
n_genes = nrow(seu),
n_cells = ncol(seu),
mito_len = length(mito_genes),
ribo_len = length(ribo_genes),
df = qc_df
)

# print per-sample console blocks
for (s in samples) {
meta <- dataset_meta[[s]]
print_block(
label = paste("sample_id:", s),
n_genes = meta$genes,
n_cells = meta$cells,
mito_len = meta$mito,
ribo_len = meta$ribo,
df = qc_df[qc_df$sample_id == s, ]
)
}
```

```
saved qc metrics summary: /global/scratch/hpc6297/neutrotime_output/neutrotime_preQC_metrics_summary.tsv 

global (all samples combined) 

genes: 28692 
cells: 16907 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount           nfeature     percent_mito      percent_ribo   
 Min.   :  486.0   Min.   :  22   Min.   : 0.0000   Min.   : 0.000  
 1st Qu.:  821.5   1st Qu.: 377   1st Qu.: 0.1005   1st Qu.: 2.438  
 Median : 1288.0   Median : 482   Median : 0.2465   Median : 4.016  
 Mean   : 1875.3   Mean   : 585   Mean   : 0.7795   Mean   : 8.303  
 3rd Qu.: 2251.0   3rd Qu.: 707   3rd Qu.: 0.5843   3rd Qu.: 7.273  
 Max.   :33654.0   Max.   :4981   Max.   :30.8642   Max.   :59.914  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   543.06    63.00        0.000        0.150
2%   592.00    88.00        0.000        0.223
3%   607.00   119.18        0.000        0.312
4%   616.00   175.24        0.000        0.525
96% 5464.80  1317.52        3.850       35.254
97% 6207.00  1494.82        4.068       38.024
98% 7204.40  1677.88        4.341       40.803
99% 9166.74  1929.94        5.007       44.157


sample_id: GSM5029335_BL_dataset1.txt.gz 

genes: 7578 
cells: 1588 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount         nfeature       percent_mito     percent_ribo   
 Min.   :  486   Min.   :  22.0   Min.   :0.0000   Min.   : 0.000  
 1st Qu.:  634   1st Qu.: 320.0   1st Qu.:0.0000   1st Qu.: 1.969  
 Median :  841   Median : 391.0   Median :0.1096   Median : 2.930  
 Mean   : 1240   Mean   : 473.4   Mean   :0.2884   Mean   : 3.913  
 3rd Qu.: 1272   3rd Qu.: 502.2   3rd Qu.:0.2567   3rd Qu.: 4.095  
 Max.   :11097   Max.   :2912.0   Max.   :9.4178   Max.   :26.801  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   497.00    73.74        0.000        0.082
2%   500.74   162.62        0.000        0.265
3%   506.00   238.00        0.000        0.735
4%   512.48   248.00        0.000        0.918
96% 3964.84  1236.16        1.491       13.460
97% 4395.78  1380.97        1.839       14.697
98% 4791.64  1588.98        2.242       16.137
99% 5395.91  1794.82        3.072       19.218


sample_id: GSM5029336_BM_dataset1.txt.gz 

genes: 8657 
cells: 1271 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount         nfeature       percent_mito      percent_ribo   
 Min.   : 1191   Min.   :  44.0   Min.   :0.00000   Min.   : 0.000  
 1st Qu.: 1833   1st Qu.: 530.5   1st Qu.:0.04068   1st Qu.: 1.311  
 Median : 2743   Median : 688.0   Median :0.08529   Median : 1.607  
 Mean   : 3703   Mean   : 846.7   Mean   :0.15036   Mean   : 2.220  
 3rd Qu.: 4441   3rd Qu.: 951.0   3rd Qu.:0.15217   3rd Qu.: 1.998  
 Max.   :33654   Max.   :4981.0   Max.   :5.65410   Max.   :24.977  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   1218.1    356.7        0.000        0.584
2%   1240.0    372.8        0.000        0.773
3%   1260.1    389.0        0.000        0.866
4%   1287.8    396.6        0.000        0.909
96% 10652.0   2092.6        0.657        9.190
97% 11167.1   2205.8        0.939       11.449
98% 12221.2   2521.2        1.134       13.435
99% 14971.8   3224.1        1.399       16.136


sample_id: GSM5029337_SP_dataset1.txt.gz 

genes: 7346 
cells: 1219 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount         nfeature       percent_mito     percent_ribo     
 Min.   :  703   Min.   : 195.0   Min.   :0.0000   Min.   : 0.09166  
 1st Qu.:  903   1st Qu.: 411.0   1st Qu.:0.0000   1st Qu.: 1.58428  
 Median : 1355   Median : 511.0   Median :0.1036   Median : 2.15569  
 Mean   : 2360   Mean   : 627.8   Mean   :0.1714   Mean   : 3.13218  
 3rd Qu.: 3411   3rd Qu.: 785.0   3rd Qu.:0.1821   3rd Qu.: 3.10016  
 Max.   :16204   Max.   :3046.0   Max.   :2.8465   Max.   :34.47239  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   721.18   291.90        0.000        0.760
2%   728.00   310.36        0.000        0.890
3%   734.54   320.62        0.000        0.985
4%   743.00   327.44        0.000        1.023
96% 6455.44  1248.96        0.933       10.019
97% 6815.54  1319.92        1.136       14.695
98% 7450.40  1487.04        1.412       18.906
99% 9066.32  1685.58        1.804       25.322


sample_id: GSM5029338_BL_dataset2.txt.gz 

genes: 12311 
cells: 4653 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount        nfeature       percent_mito      percent_ribo   
 Min.   : 588   Min.   :  22.0   Min.   :0.00000   Min.   : 0.000  
 1st Qu.: 717   1st Qu.: 326.0   1st Qu.:0.04226   1st Qu.: 3.318  
 Median : 885   Median : 379.0   Median :0.15674   Median : 4.888  
 Mean   :1276   Mean   : 383.4   Mean   :0.28531   Mean   : 5.229  
 3rd Qu.:1248   3rd Qu.: 447.0   3rd Qu.:0.32733   3rd Qu.: 6.421  
 Max.   :8756   Max.   :1891.0   Max.   :9.46970   Max.   :50.436  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   593.52    43.00        0.000        0.091
2%   599.00    51.00        0.000        0.130
3%   604.00    59.00        0.000        0.151
4%   608.00    67.00        0.000        0.173
96% 3992.84   699.92        0.933       10.672
97% 4307.04   758.00        1.322       14.416
98% 4728.56   867.96        2.052       20.218
99% 5427.96  1025.40        2.885       26.946


sample_id: GSM5029339_BM_dataset2.txt.gz 

genes: 13724 
cells: 3823 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount         nfeature       percent_mito      percent_ribo    
 Min.   : 1054   Min.   :  80.0   Min.   : 0.0000   Min.   : 0.1522  
 1st Qu.: 1469   1st Qu.: 492.0   1st Qu.: 0.1949   1st Qu.: 2.4549  
 Median : 2101   Median : 625.0   Median : 0.3319   Median : 2.9341  
 Mean   : 2644   Mean   : 728.6   Mean   : 0.6466   Mean   : 5.1521  
 3rd Qu.: 3084   3rd Qu.: 812.0   3rd Qu.: 0.5388   3rd Qu.: 3.7116  
 Max.   :19736   Max.   :4092.0   Max.   :10.3448   Max.   :50.1345  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
      ncount nfeature percent_mito percent_ribo
1%   1067.00   353.00        0.000        1.515
2%   1081.44   368.00        0.000        1.665
3%   1098.00   378.00        0.048        1.774
4%   1112.00   385.88        0.059        1.848
96%  7189.12  1693.12        3.260       24.109
97%  7718.04  1781.68        3.651       26.523
98%  8879.80  1919.00        4.012       28.937
99% 10130.12  2111.84        4.434       32.940


sample_id: GSM5029340_SP_dataset2.txt.gz 

genes: 14162 
cells: 4353 

identified mitochondrial genes: 13 
identified ribosomal genes: 190 

qc metric distributions across cells:
     ncount         nfeature       percent_mito      percent_ribo   
 Min.   :  608   Min.   : 233.0   Min.   : 0.0000   Min.   : 1.737  
 1st Qu.:  766   1st Qu.: 388.0   1st Qu.: 0.3096   1st Qu.: 5.866  
 Median : 1160   Median : 587.0   Median : 1.8537   Median :17.965  
 Mean   : 1404   Mean   : 626.7   Mean   : 1.9577   Mean   :19.182  
 3rd Qu.: 1675   3rd Qu.: 779.0   3rd Qu.: 3.2588   3rd Qu.:30.743  
 Max.   :17964   Max.   :2963.0   Max.   :30.8642   Max.   :59.914  

quantiles (1%, 2%, 3%, 4%, 96%, 97%, 98%, 99%):
     ncount nfeature percent_mito percent_ribo
1%   612.00   293.00        0.000        2.856
2%   618.00   304.00        0.000        3.165
3%   623.00   311.00        0.000        3.372
4%   628.00   317.00        0.000        3.529
96% 3162.92  1182.76        4.738       43.961
97% 3535.00  1284.76        5.037       44.941
98% 4063.88  1426.00        5.710       45.954
99% 5417.12  1600.92        7.306       47.190
```



# sessioninfo()
