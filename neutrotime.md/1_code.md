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

# singleR celldex immgen LABELS MAIN 
- label.main Neutrophils
- 4 finer labels have inflamm conditions and healthy PBL
- collapsed from original 6 labels

```r
options(width = 800)
library(Seurat)
library(SingleR)
library(celldex)
library(BiocParallel)

target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

# set biocparallel to 8 workers
bpparam <- MulticoreParam(workers = 8)
register(bpparam)

# load pre-qc merged seurat object
rds_path <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_preQC_merged.rds"
seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

# extract counts matrix for singler
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

# fetch mouse immgen reference (2024-02-26 snapshot)
ref_immgen <- celldex::fetchReference("immgen", "2024-02-26")

# run singler annotation using label.main with 8 workers
pred_immgen <- SingleR(
test = counts_mat,
ref = ref_immgen,
labels = ref_immgen$label.main,
BPPARAM = bpparam
)

# append annotation metadata
seurat_obj$SingleR_ImmGen_main <- pred_immgen$labels
seurat_obj$SingleR_ImmGen_pruned <- pred_immgen$pruned.labels
seurat_obj$is_neutrophil <- ifelse(!is.na(pred_immgen$pruned.labels) & pred_immgen$pruned.labels == "Neutrophils", TRUE, FALSE)

# atomic save to preQC merged rds
tmp_rds <- paste0(rds_path, ".tmp.", Sys.getpid())
saveRDS(seurat_obj, file = tmp_rds)
file.rename(tmp_rds, rds_path)

# calculate global metrics
global_before <- ncol(seurat_obj)
global_after <- sum(!is.na(seurat_obj$SingleR_ImmGen_pruned))
global_neutro_count <- sum(seurat_obj$is_neutrophil, na.rm = TRUE)
global_neutro_pct <- (global_neutro_count / global_after) * 100

global_df <- data.frame(
sample_id = "GLOBAL",
cells_before = global_before,
cells_after = global_after,
neutrophil_count = global_neutro_count,
neutrophil_pct = round(global_neutro_pct, 2)
)

# calculate per dataset metrics
samples <- unique(seurat_obj$sample_id)
sample_df_list <- lapply(samples, function(s) {
cells_s <- colnames(seurat_obj)[seurat_obj$sample_id == s]
cnt_before <- length(cells_s)
pruned_s <- seurat_obj$SingleR_ImmGen_pruned[cells_s]
cnt_after <- sum(!is.na(pruned_s))
cnt_neutro <- sum(seurat_obj$is_neutrophil[cells_s], na.rm = TRUE)
pct_neutro <- ifelse(cnt_after > 0, (cnt_neutro / cnt_after) * 100, 0)
data.frame(
sample_id = as.character(s),
cells_before = cnt_before,
cells_after = cnt_after,
neutrophil_count = cnt_neutro,
neutrophil_pct = round(pct_neutro, 2)
)
})

metrics_df <- rbind(global_df, do.call(rbind, sample_df_list))

# print metrics to console
print(metrics_df, row.names = FALSE)

# save metrics to tsv
out_tsv <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_singler_metrics.tsv"
write.table(metrics_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
```

## Immgen neutrophils FINE LABELS
```R
options(width = 800)
library(Seurat)
library(SingleR)
library(celldex)
library(BiocParallel)

target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

# set biocparallel to 8 workers
bpparam <- MulticoreParam(workers = 8)
register(bpparam)

# load pre-qc merged seurat object
rds_path <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_preQC_merged.rds"
seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

# extract counts matrix for singler
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")

# fetch mouse immgen reference (2024-02-26 snapshot)
ref_immgen <- celldex::fetchReference("immgen", "2024-02-26")

# run singler annotation using label.fine with 8 workers
pred_immgen_fine <- SingleR(
test = counts_mat,
ref = ref_immgen,
labels = ref_immgen$label.fine,
BPPARAM = bpparam
)

# define target fine labels for neutrophils
neutro_fine_labels <- c(
"Neutrophils (GN)",
"Neutrophils (GN.ARTH)",
"Neutrophils (GN.Thio)",
"Neutrophils (GN.URAC)"
)

# append annotation metadata
seurat_obj$SingleR_ImmGen_fine <- pred_immgen_fine$labels
seurat_obj$SingleR_ImmGen_fine_pruned <- pred_immgen_fine$pruned.labels
seurat_obj$is_neutrophil_fine <- ifelse(
!is.na(pred_immgen_fine$pruned.labels) & pred_immgen_fine$pruned.labels %in% neutro_fine_labels,
TRUE,
FALSE
)

# atomic save to preQC merged rds
tmp_rds <- paste0(rds_path, ".tmp.", Sys.getpid())
saveRDS(seurat_obj, file = tmp_rds)
file.rename(tmp_rds, rds_path)

# calculate global metrics
global_before <- ncol(seurat_obj)
global_after <- sum(!is.na(seurat_obj$SingleR_ImmGen_fine_pruned))
global_neutro_count <- sum(seurat_obj$is_neutrophil_fine, na.rm = TRUE)
global_neutro_pct <- ifelse(global_after > 0, (global_neutro_count / global_after) * 100, 0)

global_df <- data.frame(
sample_id = "GLOBAL",
cells_before = global_before,
cells_after = global_after,
neutrophil_count = global_neutro_count,
neutrophil_pct = round(global_neutro_pct, 2),
stringsAsFactors = FALSE
)

# calculate per dataset metrics
samples <- unique(seurat_obj$sample_id)
sample_df_list <- lapply(samples, function(s) {
cells_s <- colnames(seurat_obj)[seurat_obj$sample_id == s]
cnt_before <- length(cells_s)
pruned_s <- seurat_obj$SingleR_ImmGen_fine_pruned[cells_s]
cnt_after <- sum(!is.na(pruned_s))
cnt_neutro <- sum(seurat_obj$is_neutrophil_fine[cells_s], na.rm = TRUE)
pct_neutro <- ifelse(cnt_after > 0, (cnt_neutro / cnt_after) * 100, 0)
data.frame(
sample_id = as.character(s),
cells_before = cnt_before,
cells_after = cnt_after,
neutrophil_count = cnt_neutro,
neutrophil_pct = round(pct_neutro, 2),
stringsAsFactors = FALSE
)
})

sample_df <- do.call(rbind, sample_df_list)

# calculate dataset 1 and dataset 2 summary block
d1_rows <- sample_df[grepl("dataset1", sample_df$sample_id), ]
d2_rows <- sample_df[grepl("dataset2", sample_df$sample_id), ]

agg_list <- list()

if (nrow(d1_rows) > 0) {
b_d1 <- sum(d1_rows$cells_before)
a_d1 <- sum(d1_rows$cells_after)
n_d1 <- sum(d1_rows$neutrophil_count)
p_d1 <- ifelse(a_d1 > 0, round((n_d1 / a_d1) * 100, 2), 0)
agg_list[[length(agg_list) + 1]] <- data.frame(
sample_id = "SUM_dataset1",
cells_before = b_d1,
cells_after = a_d1,
neutrophil_count = n_d1,
neutrophil_pct = p_d1,
stringsAsFactors = FALSE
)
}

if (nrow(d2_rows) > 0) {
b_d2 <- sum(d2_rows$cells_before)
a_d2 <- sum(d2_rows$cells_after)
n_d2 <- sum(d2_rows$neutrophil_count)
p_d2 <- ifelse(a_d2 > 0, round((n_d2 / a_d2) * 100, 2), 0)
agg_list[[length(agg_list) + 1]] <- data.frame(
sample_id = "SUM_dataset2",
cells_before = b_d2,
cells_after = a_d2,
neutrophil_count = n_d2,
neutrophil_pct = p_d2,
stringsAsFactors = FALSE
)
}

summary_block <- do.call(rbind, agg_list)
full_metrics_df <- rbind(global_df, sample_df, summary_block)

# save metrics to tsv with labelsfine appended
out_tsv <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_singler_metrics_labelsfine.tsv"
write.table(full_metrics_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# print complete tsv
print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```

# main vs fine labels concordance 
```r
out_tsv <- "/global/scratch/hpc6297/neutrotime_output/neutrotime_labels_concordance.tsv"

neutro_fine_labels <- c(
"Neutrophils (GN)",
"Neutrophils (GN.ARTH)",
"Neutrophils (GN.Thio)",
"Neutrophils (GN.URAC)"
)

cells_main <- colnames(seurat_obj)[!is.na(seurat_obj$SingleR_ImmGen_pruned) & seurat_obj$SingleR_ImmGen_pruned == "Neutrophils"]
cells_fine <- colnames(seurat_obj)[!is.na(seurat_obj$SingleR_ImmGen_fine_pruned) & seurat_obj$SingleR_ImmGen_fine_pruned %in% neutro_fine_labels]

shared_cells <- intersect(cells_main, cells_fine)
only_main <- setdiff(cells_main, cells_fine)
only_fine <- setdiff(cells_fine, cells_main)
union_cells <- union(cells_main, cells_fine)

df_overlap <- data.frame(
barcode = union_cells,
sample_id = seurat_obj$sample_id[union_cells],
in_main = union_cells %in% cells_main,
in_fine = union_cells %in% cells_fine,
stringsAsFactors = FALSE
)

overlap_by_sample <- do.call(rbind, lapply(split(df_overlap, df_overlap$sample_id), function(sub) {
n_main <- sum(sub$in_main)
n_fine <- sum(sub$in_fine)
n_shared <- sum(sub$in_main & sub$in_fine)
n_only_main <- sum(sub$in_main & !sub$in_fine)
n_only_fine <- sum(!sub$in_main & sub$in_fine)
data.frame(
sample_id = unique(sub$sample_id),
neutrophils_main = n_main,
neutrophils_fine = n_fine,
retained_in_both = n_shared,
lost_in_fine = n_only_main,
gained_in_fine = n_only_fine,
pct_fine_retained_in_main = ifelse(n_fine > 0, round((n_shared / n_fine) * 100, 2), 0),
pct_main_retained_in_fine = ifelse(n_main > 0, round((n_shared / n_main) * 100, 2), 0),
stringsAsFactors = FALSE
)
}))

global_row <- data.frame(
sample_id = "GLOBAL",
neutrophils_main = length(cells_main),
neutrophils_fine = length(cells_fine),
retained_in_both = length(shared_cells),
lost_in_fine = length(only_main),
gained_in_fine = length(only_fine),
pct_fine_retained_in_main = ifelse(length(cells_fine) > 0, round((length(shared_cells) / length(cells_fine)) * 100, 2), 0),
pct_main_retained_in_fine = ifelse(length(cells_main) > 0, round((length(shared_cells) / length(cells_main)) * 100, 2), 0),
stringsAsFactors = FALSE
)

concordance_df <- rbind(global_row, overlap_by_sample)

write.table(concordance_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("saved concordance metrics to", out_tsv, "\n\n")
print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)

# format discordance tables as two-column data frames
main_only_df <- as.data.frame(table(label_fine = seurat_obj$SingleR_ImmGen_fine_pruned[only_main], useNA = "ifany"), responseName = "cell_count")
main_only_df <- main_only_df[order(-main_only_df$cell_count), ]

fine_only_df <- as.data.frame(table(label_main = seurat_obj$SingleR_ImmGen_pruned[only_fine], useNA = "ifany"), responseName = "cell_count")
fine_only_df <- fine_only_df[order(-fine_only_df$cell_count), ]

cat("\nidentity in label.fine of cells kept only by label.main:\n")
print(main_only_df, row.names = FALSE)

cat("\nidentity in label.main of cells kept only by label.fine:\n")
print(fine_only_df, row.names = FALSE)
```

```
saved concordance metrics to /global/scratch/hpc6297/neutrotime_output/neutrotime_labels_concordance.tsv 

                     sample_id neutrophils_main neutrophils_fine retained_in_both lost_in_fine gained_in_fine pct_fine_retained_in_main pct_main_retained_in_fine
                        GLOBAL            13392            13129            13111          281             18                     99.86                     97.90
 GSM5029335_BL_dataset1.txt.gz             1394             1383             1383           11              0                    100.00                     99.21
 GSM5029336_BM_dataset1.txt.gz             1212             1198             1198           14              0                    100.00                     98.84
 GSM5029337_SP_dataset1.txt.gz             1165             1163             1163            2              0                    100.00                     99.83
 GSM5029338_BL_dataset2.txt.gz             4068             3903             3891          177             12                     99.69                     95.65
 GSM5029339_BM_dataset2.txt.gz             3427             3404             3402           25              2                     99.94                     99.27
 GSM5029340_SP_dataset2.txt.gz             2126             2078             2074           52              4                     99.81                     97.55

identity in label.fine of cells kept only by label.main:
                   label_fine cell_count
                         <NA>        150
        Monocytes (MO.6C-II+)         44
           B cells (proB.FrA)         36
              B cells (B.FrF)          8
             Macrophages (MF)          8
        Monocytes (MO.6C+II-)          8
      Monocytes (MO.6C-IIINT)          7
           B cells (proB.CLP)          6
             Stem cells (GMP)          3
           Stem cells (LTHSC)          3
         Stem cells (SC.STSL)          3
             Stem cells (MLP)          2
 Stromal cells (ST.31-38-44-)          2
         NK cells (NK.DAP10-)          1

identity in label.main of cells kept only by label.fine:
  label_main cell_count
        <NA>         15
   Basophils          1
 Macrophages          1
   Monocytes          1
```

# subset to labels fine
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

options(width = 800)
library(Seurat)
library(SingleR)
library(celldex)
library(BiocParallel)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
setwd(work_dir)

bpparam <- MulticoreParam(workers = 8)
register(bpparam)

rds_path <- file.path(work_dir, "neutrotime_preQC_merged.rds")
seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
ref_immgen <- celldex::fetchReference("immgen", "2024-02-26")

pred_immgen_fine <- SingleR(
test = counts_mat,
ref = ref_immgen,
labels = ref_immgen$label.fine,
BPPARAM = bpparam
)

neutro_fine_labels <- c(
"Neutrophils (GN)",
"Neutrophils (GN.ARTH)",
"Neutrophils (GN.Thio)",
"Neutrophils (GN.URAC)"
)

seurat_obj$SingleR_ImmGen_fine <- pred_immgen_fine$labels
seurat_obj$SingleR_ImmGen_fine_pruned <- pred_immgen_fine$pruned.labels
seurat_obj$is_neutrophil_fine <- ifelse(
!is.na(pred_immgen_fine$pruned.labels) & pred_immgen_fine$pruned.labels %in% neutro_fine_labels,
TRUE,
FALSE
)

tmp_merged_rds <- paste0(rds_path, ".tmp.", Sys.getpid())
saveRDS(seurat_obj, file = tmp_merged_rds)
file.rename(tmp_merged_rds, rds_path)

neutro_subset <- subset(seurat_obj, subset = is_neutrophil_fine == TRUE)

subset_rds_path <- file.path(work_dir, "neutrotime_preQC_neutrophils_fine.rds")
tmp_subset_rds <- paste0(subset_rds_path, ".tmp.", Sys.getpid())
saveRDS(neutro_subset, file = tmp_subset_rds)
file.rename(tmp_subset_rds, subset_rds_path)

global_before <- ncol(seurat_obj)
global_after <- sum(!is.na(seurat_obj$SingleR_ImmGen_fine_pruned))
global_neutro_count <- ncol(neutro_subset)
global_neutro_pct <- ifelse(global_after > 0, (global_neutro_count / global_after) * 100, 0)

global_df <- data.frame(
sample_id = "GLOBAL",
cells_before = global_before,
cells_after = global_after,
neutrophil_count = global_neutro_count,
neutrophil_pct = round(global_neutro_pct, 2),
stringsAsFactors = FALSE
)

samples <- unique(seurat_obj$sample_id)
sample_df_list <- lapply(samples, function(s) {
cells_s <- colnames(seurat_obj)[seurat_obj$sample_id == s]
cnt_before <- length(cells_s)
pruned_s <- seurat_obj$SingleR_ImmGen_fine_pruned[cells_s]
cnt_after <- sum(!is.na(pruned_s))
cnt_neutro <- sum(neutro_subset$sample_id == s)
pct_neutro <- ifelse(cnt_after > 0, (cnt_neutro / cnt_after) * 100, 0)
data.frame(
sample_id = as.character(s),
cells_before = cnt_before,
cells_after = cnt_after,
neutrophil_count = cnt_neutro,
neutrophil_pct = round(pct_neutro, 2),
stringsAsFactors = FALSE
)
})

sample_df <- do.call(rbind, sample_df_list)

d1_rows <- sample_df[grepl("dataset1", sample_df$sample_id), ]
d2_rows <- sample_df[grepl("dataset2", sample_df$sample_id), ]

agg_list <- list()

if (nrow(d1_rows) > 0) {
b_d1 <- sum(d1_rows$cells_before)
a_d1 <- sum(d1_rows$cells_after)
n_d1 <- sum(d1_rows$neutrophil_count)
p_d1 <- ifelse(a_d1 > 0, round((n_d1 / a_d1) * 100, 2), 0)
agg_list[[length(agg_list) + 1]] <- data.frame(
sample_id = "SUM_dataset1",
cells_before = b_d1,
cells_after = a_d1,
neutrophil_count = n_d1,
neutrophil_pct = p_d1,
stringsAsFactors = FALSE
)
}

if (nrow(d2_rows) > 0) {
b_d2 <- sum(d2_rows$cells_before)
a_d2 <- sum(d2_rows$cells_after)
n_d2 <- sum(d2_rows$neutrophil_count)
p_d2 <- ifelse(a_d2 > 0, round((n_d2 / a_d2) * 100, 2), 0)
agg_list[[length(agg_list) + 1]] <- data.frame(
sample_id = "SUM_dataset2",
cells_before = b_d2,
cells_after = a_d2,
neutrophil_count = n_d2,
neutrophil_pct = p_d2,
stringsAsFactors = FALSE
)
}

summary_block <- do.call(rbind, agg_list)
full_metrics_df <- rbind(global_df, sample_df, summary_block)

print(full_metrics_df, row.names = FALSE)
```

# preQC check on neutrophil subset
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_path <- file.path(work_dir, "neutrotime_preQC_neutrophils_fine.rds")
out_tsv <- file.path(work_dir, "neutrophil_subset_preQC_metrics.tsv")

seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

cat("metadata headers:\n")
print(colnames(seurat_obj@meta.data))
cat("\n")

if (!"percent_mito" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_mito <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
}
if (!"percent_ribo" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_ribo <- PercentageFeatureSet(seurat_obj, pattern = "^(Rps|Rpl|RPS|RPL)")
}

qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
qc_df <- seurat_obj@meta.data[, qc_cols]
colnames(qc_df) <- c("ncount", "nfeature", "percent_mito", "percent_ribo")

mito_genes <- grep("^(mt-|MT-)", rownames(seurat_obj), value = TRUE)
ribo_genes <- grep("^(Rps|Rpl|RPS|RPL)", rownames(seurat_obj), value = TRUE)

stats_list <- list(
"total_genes" = c(nrow(seurat_obj), NA, NA, NA),
"total_cells" = c(ncol(seurat_obj), NA, NA, NA),
"mito_genes" = c(length(mito_genes), NA, NA, NA),
"ribo_genes" = c(length(ribo_genes), NA, NA, NA),
"Min." = sapply(qc_df, function(x) min(x, na.rm = TRUE)),
"1%" = sapply(qc_df, function(x) quantile(x, 0.01, na.rm = TRUE)),
"2%" = sapply(qc_df, function(x) quantile(x, 0.02, na.rm = TRUE)),
"3%" = sapply(qc_df, function(x) quantile(x, 0.03, na.rm = TRUE)),
"4%" = sapply(qc_df, function(x) quantile(x, 0.04, na.rm = TRUE)),
"1st Qu." = sapply(qc_df, function(x) quantile(x, 0.25, na.rm = TRUE)),
"Median" = sapply(qc_df, function(x) median(x, na.rm = TRUE)),
"Mean" = sapply(qc_df, function(x) mean(x, na.rm = TRUE)),
"3rd Qu." = sapply(qc_df, function(x) quantile(x, 0.75, na.rm = TRUE)),
"96%" = sapply(qc_df, function(x) quantile(x, 0.96, na.rm = TRUE)),
"97%" = sapply(qc_df, function(x) quantile(x, 0.97, na.rm = TRUE)),
"98%" = sapply(qc_df, function(x) quantile(x, 0.98, na.rm = TRUE)),
"99%" = sapply(qc_df, function(x) quantile(x, 0.99, na.rm = TRUE)),
"Max." = sapply(qc_df, function(x) max(x, na.rm = TRUE))
)

mat_out <- do.call(rbind, stats_list)
res_df <- data.frame(
statistic = names(stats_list),
round(mat_out, 4),
row.names = NULL,
stringsAsFactors = FALSE
)
colnames(res_df) <- c("statistic", "ncount", "nfeature", "percent_mito", "percent_ribo")

write.table(res_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("saved metrics to", out_tsv, "\n\n")
print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```

# mito 5%
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_path <- file.path(work_dir, "neutrotime_preQC_neutrophils_fine.rds")
out_rds_path <- file.path(work_dir, "neutrotime_postmito5_neutrophils_fine.rds")

seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

if (!"percent_mito" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_mito <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
}

cells_before <- ncol(seurat_obj)
neutro_postmito <- subset(seurat_obj, subset = percent_mito < 5)
cells_after <- ncol(neutro_postmito)
cells_removed <- cells_before - cells_after

tmp_out_rds <- paste0(out_rds_path, ".tmp.", Sys.getpid())
saveRDS(neutro_postmito, file = tmp_out_rds)
file.rename(tmp_out_rds, out_rds_path)

cat(sprintf("cells before mito filter : %d\n", cells_before))
cat(sprintf("cells removed (mito >= 5%%): %d\n", cells_removed))
cat(sprintf("cells retained (< 5%% mito): %d\n", cells_after))
cat("saved atomic rds to", out_rds_path, "\n")
```
neutrophil_postmito5_neutrophils_metrics.tsv

# QC inflection graphs 
```r
options(width = 800)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_path <- file.path(work_dir, "neutrotime_postmito5_neutrophils_fine.rds")
seurat_obj <- readRDS(rds_path)
DefaultAssay(seurat_obj) <- "RNA"

if (!"percent_mito" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_mito <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
}
if (!"percent_ribo" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_ribo <- PercentageFeatureSet(seurat_obj, pattern = "^(Rps|Rpl|RPS|RPL)")
}

cfg <- list(
dpi = 300,
line_color = "#001F5B",
global_dim = c(w = 10, h = 10),
panel_dim = c(w = 10, h = 8),
metrics = list(
nCount_RNA = list(col = "ncount", label = "nCount_RNA", log = TRUE, pct = FALSE),
nFeature_RNA = list(col = "nfeature", label = "nFeature_RNA", log = TRUE, pct = FALSE),
percent_mito = list(col = "percent_mito", label = "percent_mito", log = FALSE, pct = TRUE),
percent_ribo = list(col = "percent_ribo", label = "percent_ribo", log = FALSE, pct = TRUE)
)
)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
qc_data <- seurat_obj@meta.data[, names(cfg$metrics)]

# 1. global continuous percentile inflection plot
plot_percentile_curve <- function(values, m_cfg) {
sorted_vals <- sort(values, decreasing = FALSE)
df <- data.frame(Percentile = (seq_along(sorted_vals) / length(sorted_vals)) * 100, Value = sorted_vals)
p <- ggplot(df, aes(x = Percentile, y = Value)) +
geom_line(color = cfg$line_color, linewidth = 1.1) +
labs(title = paste(m_cfg$label, "Inflection"), x = "Percentile (%)", y = m_cfg$label) +
theme_minimal(base_size = 12) +
theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13), panel.grid.minor = element_blank()) +
scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100))
if (m_cfg$log) {
p <- p + scale_y_log10(labels = comma)
} else {
p <- p + scale_y_continuous(labels = comma)
}
return(p)
}

global_plots <- lapply(names(cfg$metrics), function(m) {
plot_percentile_curve(qc_data[[m]], cfg$metrics[[m]])
})

plot_global <- (global_plots[[1]] | global_plots[[2]]) / (global_plots[[3]] | global_plots[[4]]) +
plot_annotation(
title = "Neutrotime QC Outliers (Neutrophil Subset)",
subtitle = "After 5% Mito Ceiling",
theme = theme(
plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
plot.subtitle = element_text(face = "plain", size = 12, hjust = 0.5)
)
)

out_png_plot <- file.path(work_dir, paste0("neutrotime_qc_percentile_inflections_", timestamp, ".png"))
ggsave(out_png_plot, plot_global, width = cfg$global_dim["w"], height = cfg$global_dim["h"], dpi = cfg$dpi)

# 2. sequenced percentile calculations and tsv output
probs_sequenced <- sort(unique(c(
seq(0, 0.01, by = 0.001),
c(0.02, 0.03, 0.04, 0.05),
c(0.95, 0.96, 0.97, 0.98, 0.99),
seq(0.99, 1.0, by = 0.001)
)))

q_mat <- t(sapply(probs_sequenced, function(p) sapply(qc_data, function(x) quantile(x, probs = p, na.rm = TRUE))))

percentile_labels <- ifelse(
probs_sequenced == 0, "min",
ifelse(probs_sequenced == 1, "max", paste0(round(probs_sequenced * 100, 2), "%"))
)

combined_df <- data.frame(
prob = probs_sequenced,
percentile = percentile_labels,
round(q_mat, 4),
row.names = NULL,
stringsAsFactors = FALSE
)
colnames(combined_df) <- c("prob", "percentile", "ncount", "nfeature", "percent_mito", "percent_ribo")

out_tsv <- file.path(work_dir, paste0("neutrotime_sequenced_percentiles_", timestamp, ".tsv"))
write.table(combined_df[, -1], file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# 3. panel plotting for low and high inflection windows
plot_quantile_curve <- function(df, m_cfg, range_label, x_breaks) {
p <- ggplot(df, aes(x = prob * 100, y = .data[[m_cfg$col]])) +
geom_line(color = cfg$line_color, linewidth = 1.0) +
geom_point(color = cfg$line_color, size = 2) +
labs(title = paste0(m_cfg$label, " (", range_label, ")"), x = "Percentile (%)", y = m_cfg$label) +
theme_minimal(base_size = 11) +
theme(
plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
axis.text.x = element_text(angle = 45, hjust = 1)
) +
scale_x_continuous(breaks = x_breaks, minor_breaks = NULL, labels = scales::number_format(accuracy = 0.1))
if (m_cfg$pct) {
p <- p + scale_y_continuous(labels = function(x) paste0(x, "%"))
} else {
p <- p + scale_y_continuous(labels = comma)
}
return(p)
}

render_panel <- function(df, range_label, title, out_name, x_breaks, sub_text = "Inflection Detail") {
plots <- lapply(cfg$metrics, function(m) plot_quantile_curve(df, m, range_label, x_breaks))
panel <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
plot_annotation(
title = title,
subtitle = sub_text,
theme = theme(
plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
plot.subtitle = element_text(face = "plain", size = 11, hjust = 0.5)
)
)
out_path <- file.path(work_dir, paste0(out_name, "_", timestamp, ".png"))
ggsave(out_path, panel, width = cfg$panel_dim["w"], height = cfg$panel_dim["h"], dpi = cfg$dpi)
return(out_path)
}

df_fine_low <- combined_df[combined_df$prob <= 0.01, ]
df_fine_high <- combined_df[combined_df$prob >= 0.99, ]
df_tail_low5 <- combined_df[combined_df$prob <= 0.05, ]
df_tail_high5 <- combined_df[combined_df$prob >= 0.95, ]

out_png_low <- render_panel(df_fine_low, "0% - 1.0%", "Neutrotime Fine Low Percentiles (0% to 1.0%)", "neutrotime_qc_fine_low_percentiles", seq(0, 1.0, by = 0.1), "0.1% Increments")
out_png_high <- render_panel(df_fine_high, "99.0% - 100%", "Neutrotime Fine High Percentiles (99.0% to 100%)", "neutrotime_qc_fine_high_percentiles", seq(99.0, 100.0, by = 0.1), "0.1% Increments")
out_png_low5 <- render_panel(df_tail_low5, "0% - 5.0%", "Neutrotime Low Percentiles (Min to 5%)", "neutrotime_qc_low_5pct_percentiles", seq(0, 5.0, by = 0.5), "0.5% Increments")
out_png_high5 <- render_panel(df_tail_high5, "95.0% - 100%", "Neutrotime High Percentiles (95% to Max)", "neutrotime_qc_high_5pct_percentiles", seq(95.0, 100.0, by = 0.5), "0.5% Increments")

cat(out_png_plot, "\n")
cat(out_png_low, "\n")
cat(out_png_high, "\n")
cat(out_png_low5, "\n")
cat(out_png_high5, "\n")
cat(out_tsv, "\n\n")

print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```
110516.tsv 

## QC at 1% and 99% features
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_in <- file.path(work_dir, "neutrotime_postmito5_neutrophils_fine.rds")
rds_out <- file.path(work_dir, "neutrotime_postquantile_features1_99_neutrophils_fine.rds")
out_tsv <- file.path(work_dir, "neutrophil_postquantile_features1_99_neutrophils_metrics.tsv")

seurat_obj <- readRDS(rds_in)
DefaultAssay(seurat_obj) <- "RNA"

if (!"percent_mito" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_mito <- PercentageFeatureSet(seurat_obj, pattern = "^(mt-|MT-)")
}
if (!"percent_ribo" %in% colnames(seurat_obj@meta.data)) {
seurat_obj$percent_ribo <- PercentageFeatureSet(seurat_obj, pattern = "^(Rps|Rpl|RPS|RPL)")
}

cells_before <- ncol(seurat_obj)

feature_floor <- unname(quantile(seurat_obj$nFeature_RNA, probs = 0.01, na.rm = TRUE))
feature_ceiling <- unname(quantile(seurat_obj$nFeature_RNA, probs = 0.99, na.rm = TRUE))

seurat_filtered <- subset(
seurat_obj,
subset = nFeature_RNA >= feature_floor &
nFeature_RNA <= feature_ceiling
)

tmp_rds <- paste0(rds_out, ".tmp.", Sys.getpid())
saveRDS(seurat_filtered, file = tmp_rds)
file.rename(tmp_rds, rds_out)

cells_after <- ncol(seurat_filtered)
cells_removed <- cells_before - cells_after
attrition_rate <- (cells_removed / cells_before) * 100

qc_cols <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")
qc_df <- seurat_filtered@meta.data[, qc_cols]
colnames(qc_df) <- c("ncount", "nfeature", "percent_mito", "percent_ribo")

mito_genes <- grep("^(mt-|MT-)", rownames(seurat_filtered), value = TRUE)
ribo_genes <- grep("^(Rps|Rpl|RPS|RPL)", rownames(seurat_filtered), value = TRUE)

stats_list <- list(
"total_genes" = c(nrow(seurat_filtered), NA, NA, NA),
"cells_before" = c(cells_before, NA, NA, NA),
"cells_retained" = c(cells_after, NA, NA, NA),
"cells_removed" = c(cells_removed, NA, NA, NA),
"attrition_rate_pct" = c(round(attrition_rate, 4), NA, NA, NA),
"feature_floor_1pct" = c(round(feature_floor, 4), NA, NA, NA),
"feature_ceiling_99pct" = c(round(feature_ceiling, 4), NA, NA, NA),
"mito_genes" = c(length(mito_genes), NA, NA, NA),
"ribo_genes" = c(length(ribo_genes), NA, NA, NA),
"Min." = sapply(qc_df, function(x) min(x, na.rm = TRUE)),
"1%" = sapply(qc_df, function(x) quantile(x, 0.01, na.rm = TRUE)),
"2%" = sapply(qc_df, function(x) quantile(x, 0.02, na.rm = TRUE)),
"3%" = sapply(qc_df, function(x) quantile(x, 0.03, na.rm = TRUE)),
"1st Qu." = sapply(qc_df, function(x) quantile(x, 0.25, na.rm = TRUE)),
"Median" = sapply(qc_df, function(x) median(x, na.rm = TRUE)),
"Mean" = sapply(qc_df, function(x) mean(x, na.rm = TRUE)),
"3rd Qu." = sapply(qc_df, function(x) quantile(x, 0.75, na.rm = TRUE)),
"97%" = sapply(qc_df, function(x) quantile(x, 0.97, na.rm = TRUE)),
"98%" = sapply(qc_df, function(x) quantile(x, 0.98, na.rm = TRUE)),
"99%" = sapply(qc_df, function(x) quantile(x, 0.99, na.rm = TRUE)),
"Max." = sapply(qc_df, function(x) max(x, na.rm = TRUE))
)

mat_out <- do.call(rbind, stats_list)
res_df <- data.frame(
statistic = names(stats_list),
round(mat_out, 4),
row.names = NULL,
stringsAsFactors = FALSE
)
colnames(res_df) <- c("statistic", "ncount", "nfeature", "percent_mito", "percent_ribo")

write.table(res_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(rds_out, "\n")
cat(out_tsv, "\n\n")

print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```
neutrotime_postquantile_features1_99_neutrophils_fine.rds 

# gene attrition check
```r
options(width = 800)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)
library(Matrix)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_in <- file.path(work_dir, "neutrotime_postquantile_features1_99_neutrophils_fine.rds")
out_tsv <- file.path(work_dir, "neutrophil_gene_attrition_thresholds_features1_99.tsv")

seurat_obj <- readRDS(rds_in)
DefaultAssay(seurat_obj) <- "RNA"

counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
total_genes <- nrow(counts_mat)
total_cells <- ncol(counts_mat)

# calculate nonzero cells per gene
num_cells_per_gene <- Matrix::rowSums(counts_mat > 0)

cutoffs <- c(3, 5, 10)

attrition_list <- lapply(cutoffs, function(k) {
genes_retained <- sum(num_cells_per_gene >= k)
genes_removed <- total_genes - genes_retained
attrition_pct <- (genes_removed / total_genes) * 100
data.frame(
min_cells_threshold = paste0(">=", k, " cells"),
total_genes = total_genes,
genes_retained = genes_retained,
genes_removed = genes_removed,
gene_attrition_pct = round(attrition_pct, 4),
stringsAsFactors = FALSE
)
})

res_df <- do.call(rbind, attrition_list)
write.table(res_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(out_tsv, "\n\n")
print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```

# apply min 10 cells floor
```r
options(width = 800)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)
library(Matrix)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_in <- file.path(work_dir, "neutrotime_postquantile_features1_99_neutrophils_fine.rds")
rds_out <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")
out_tsv <- file.path(work_dir, "neutrophil_postQC_final_summary.tsv")

seurat_obj <- readRDS(rds_in)
DefaultAssay(seurat_obj) <- "RNA"

counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
total_genes_before <- nrow(counts_mat)
total_cells <- ncol(counts_mat)

# identify genes expressed in at least 10 cells
num_cells_per_gene <- Matrix::rowSums(counts_mat > 0)
genes_to_keep <- names(num_cells_per_gene[num_cells_per_gene >= 10])

# subset features on seurat object
seurat_filtered <- seurat_obj[genes_to_keep, ]

tmp_rds <- paste0(rds_out, ".tmp.", Sys.getpid())
saveRDS(seurat_filtered, file = tmp_rds)
file.rename(tmp_rds, rds_out)

genes_retained <- nrow(seurat_filtered)
genes_removed <- total_genes_before - genes_retained
gene_attrition_pct <- (genes_removed / total_genes_before) * 100

summary_df <- data.frame(
metric = c(
"total_cells_retained",
"genes_before_filter",
"genes_retained_ge10",
"genes_removed",
"gene_attrition_pct"
),
value = c(
total_cells,
total_genes_before,
genes_retained,
genes_removed,
round(gene_attrition_pct, 4)
),
stringsAsFactors = FALSE
)

write.table(summary_df, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat(rds_out, "\n")
cat(out_tsv, "\n\n")

print(read.table(out_tsv, header = TRUE, sep = "\t"), row.names = FALSE)
```
neutrotime_postQC_neutrophils.rds

# check rds metadata
```r
options(width = 120)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_in <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")

seurat_obj <- readRDS(rds_in)

cat("dimensions (cells x genes):\n")
cat(ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n\n")

meta_cols <- colnames(seurat_obj@meta.data)
cat("total metadata columns:", length(meta_cols), "\n\n")

print(data.frame(
index = seq_along(meta_cols),
column_name = meta_cols,
data_type = sapply(seurat_obj@meta.data, class),
stringsAsFactors = FALSE
), row.names = FALSE)
```

12845 cells x 9735 genes

# print key metadata
```r
options(width = 800)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_in <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")

seurat_obj <- readRDS(rds_in)

target_cols <- c("sample_id", "dataset", "site", "strain", "sex", "health_status", "is_neutrophil_fine")
present_cols <- intersect(target_cols, colnames(seurat_obj@meta.data))

unique_samples <- unique(seurat_obj@meta.data[, present_cols, drop = FALSE])
unique_samples <- unique_samples[order(unique_samples$sample_id), , drop = FALSE]

print(unique_samples, row.names = FALSE)
```

# clean out main labels immgen 3 columns
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
rds_path <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")

seurat_obj <- readRDS(rds_path)

cols_to_remove <- c("SingleR_ImmGen_main", "SingleR_ImmGen_pruned", "is_neutrophil")

for (col in cols_to_remove) {
if (col %in% colnames(seurat_obj@meta.data)) {
seurat_obj@meta.data[[col]] <- NULL
cat("removed:", col, "\n")
} else {
cat("skipped (not found):", col, "\n")
}
}

tmp_rds <- paste0(rds_path, ".tmp.", Sys.getpid())
saveRDS(seurat_obj, file = tmp_rds)
file.rename(tmp_rds, rds_path)

cat("\nsaved atomically to:", rds_path, "\n\n")

print(data.frame(
index = seq_along(colnames(seurat_obj@meta.data)),
column_name = colnames(seurat_obj@meta.data),
data_type = sapply(seurat_obj@meta.data, class),
stringsAsFactors = FALSE
), row.names = FALSE)
```

```
saved atomically to: /global/scratch/hpc6297/neutrotime_output/neutrotime_postQC_neutrophils.rds 

 index                column_name data_type
     1                 orig.ident    factor
     2                 nCount_RNA   numeric
     3               nFeature_RNA   integer
     4                    cell_id character
     5                raw_barcode character
     6                  sample_id character
     7                    dataset character
     8                       site character
     9                     strain character
    10                        sex character
    11              health_status character
    12        SingleR_ImmGen_fine character
    13 SingleR_ImmGen_fine_pruned character
    14         is_neutrophil_fine   logical
    15               percent_mito   numeric
    16               percent_ribo   numeric
```

# postQC plots
```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
setwd(work_dir)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
rds_path <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")
seu <- readRDS(rds_path)

# harmonize column names if needed
if (!"percent_mito" %in% colnames(seu@meta.data) && "percent.mt" %in% colnames(seu@meta.data)) {
seu$percent_mito <- seu$percent.mt
}
if (!"percent_ribo" %in% colnames(seu@meta.data) && "percent.ribo" %in% colnames(seu@meta.data)) {
seu$percent_ribo <- seu$percent.ribo
}

qc_dt <- as.data.table(seu@meta.data)
qc_dt[, dataset := fifelse(grepl("1", as.character(dataset)), "1", "2")]
n_cells_str <- format(nrow(qc_dt), big.mark = ",")

set.seed(42)
qc_dt_shuffled <- qc_dt[sample(nrow(qc_dt))]

palettes <- list(
dataset = c("1" = "#00C5CD", "2" = "#FF69B4"),
site = c("bone marrow" = "#D55E00", "peripheral blood" = "#0072B2", "spleen" = "#009E73")
)
dark_green <- "#1b4d3e"
metrics <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo")

theme_cfg <- list(
dpi = 300,
title_size_violin = 25,
subtitle_size_violin = 17,
title_size_scatter = 27,
subtitle_size_scatter = 18,
legend_title_size = 15,
legend_text_size = 15,
legend_circle_size = 5
)

base_theme <- theme_classic(base_size = 14) +
theme(
plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
plot.subtitle = element_text(hjust = 0.5, size = 12),
axis.title = element_text(face = "bold", size = 12),
axis.text = element_text(color = "black", size = 12)
)

apply_shared_legend <- function(layout_obj, group_var, pos = "right") {
if (!is.null(group_var)) {
layout_obj <- layout_obj +
plot_layout(guides = "collect") &
theme(
legend.position = pos,
legend.title = element_text(face = "bold", size = theme_cfg$legend_title_size),
legend.text = element_text(size = theme_cfg$legend_text_size)
) &
guides(color = guide_legend(title = group_var, override.aes = list(size = theme_cfg$legend_circle_size, alpha = 1)))
} else {
layout_obj <- layout_obj & theme(legend.position = "none")
}
return(layout_obj)
}

make_violin_strip <- function(df, group_var = NULL, pal = NULL, row_title = "") {
plots <- lapply(metrics, function(m) {
p <- ggplot(df, aes(x = "All Cells", y = .data[[m]]))
if (is.null(group_var)) {
p <- p + geom_violin(fill = dark_green, color = "black", trim = FALSE)
} else {
p <- p +
geom_jitter(aes(color = .data[[group_var]]), width = 0.22, size = 0.40, alpha = 0.25, stroke = 0) +
geom_violin(color = "grey30", fill = NA, width = 0.85, trim = FALSE, linewidth = 0.55) +
scale_color_manual(values = pal)
}
p + labs(title = m, x = NULL, y = NULL) +
base_theme +
theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
})

strip <- wrap_plots(plots, nrow = 1)
strip <- apply_shared_legend(strip, group_var, pos = "right")
strip + plot_annotation(title = row_title, theme = theme(plot.title = element_text(face = "bold", size = 19, hjust = 0)))
}

row_v_global <- make_violin_strip(qc_dt, NULL, NULL, "Global")
row_v_dataset <- make_violin_strip(qc_dt_shuffled, "dataset", palettes$dataset, "By Dataset")
row_v_site <- make_violin_strip(qc_dt_shuffled, "site", palettes$site, "By Site")

combined_violin_plot <- (row_v_global / row_v_dataset / row_v_site) +
plot_annotation(
title = "Comprehensive QC Violin Distributions",
subtitle = paste0("PostQC Neutrophils (n = ", n_cells_str, ")"),
theme = theme(
plot.title = element_text(face = "bold", hjust = 0.5, size = theme_cfg$title_size_violin),
plot.subtitle = element_text(hjust = 0.5, size = theme_cfg$subtitle_size_violin)
)
)

fn_combined_violin <- file.path(work_dir, paste0("neutrotime_qc_combined_violins_", timestamp, ".png"))
ggsave(fn_combined_violin, combined_violin_plot, width = 16, height = 18, dpi = theme_cfg$dpi)

pairs_list <- list(
c("nCount_RNA", "nFeature_RNA"),
c("nCount_RNA", "percent_mito"),
c("nCount_RNA", "percent_ribo"),
c("nFeature_RNA", "percent_mito"),
c("nFeature_RNA", "percent_ribo"),
c("percent_ribo", "percent_mito")
)

make_scatter_column <- function(df, group_var = NULL, pal = NULL, col_title = "") {
sub_plots <- lapply(pairs_list, function(pr) {
x_val <- df[[pr[1]]]
y_val <- df[[pr[2]]]
ct <- cor.test(x_val, y_val, method = "pearson")
stat_lbl <- paste0("r = ", round(ct$estimate, 2), ", ", fifelse(ct$p.value < 0.01, "p < 0.01", paste("p =", round(ct$p.value, 3))))

p <- ggplot(df, aes(x = .data[[pr[1]]], y = .data[[pr[2]]]))
if (is.null(group_var)) {
p <- p + geom_point(color = dark_green, alpha = 0.20, size = 0.70, stroke = 0)
} else {
p <- p + geom_point(aes(color = .data[[group_var]]), alpha = 0.30, size = 0.70, stroke = 0) +
scale_color_manual(values = pal)
}
p + labs(title = stat_lbl, x = pr[1], y = pr[2]) +
base_theme +
theme(
axis.title = element_text(face = "bold", size = 14),
axis.text = element_text(color = "black", size = 12),
axis.text.x = element_text(angle = 45, hjust = 1),
plot.title = element_text(face = "plain", size = 16, hjust = 0.5)
)
})

col_layout <- wrap_plots(sub_plots, ncol = 1)
col_layout <- apply_shared_legend(col_layout, group_var, pos = "bottom")
col_layout + plot_annotation(title = col_title, theme = theme(plot.title = element_text(face = "bold", size = 19, hjust = 0.5)))
}

col_sc_global <- make_scatter_column(qc_dt, NULL, NULL, "Global")
col_sc_dataset <- make_scatter_column(qc_dt_shuffled, "dataset", palettes$dataset, "By Dataset")
col_sc_site <- make_scatter_column(qc_dt_shuffled, "site", palettes$site, "By Site")

combined_scatter_plot <- (col_sc_global | col_sc_dataset | col_sc_site) +
plot_annotation(
title = "Comprehensive Pairwise QC Metrics",
subtitle = paste0("PostQC Neutrophils (n = ", n_cells_str, ")"),
theme = theme(
plot.title = element_text(face = "bold", hjust = 0.5, size = theme_cfg$title_size_scatter),
plot.subtitle = element_text(hjust = 0.5, size = theme_cfg$subtitle_size_scatter)
)
)

fn_combined_scatter <- file.path(work_dir, paste0("neutrotime_qc_combined_pairwise_scatters_", timestamp, ".png"))
ggsave(fn_combined_scatter, combined_scatter_plot, width = 18, height = 22, dpi = theme_cfg$dpi)

all_saved_plots <- c(fn_combined_violin, fn_combined_scatter)
cat("saved combined plots:\n", paste(all_saved_plots, collapse = "\n"), "\n")
```

20260908_105717.png    

# SCTransform
- don't set multithreading to avoid freezing

```r
options(width = 800)
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
.libPaths(c(target_lib, .libPaths()))

library(Seurat)
library(future)

# disable forked futures to prevent BLAS/OpenMP mutex deadlocks
plan(sequential)

work_dir <- "/global/scratch/hpc6297/neutrotime_output"
setwd(work_dir)

rds_path <- file.path(work_dir, "neutrotime_postQC_neutrophils.rds")
neu <- readRDS(rds_path)
DefaultAssay(neu) <- "RNA"

# harmonize sample and library columns
if ("sample_id" %in% colnames(neu@meta.data)) {
neu$lib_ID <- factor(neu$sample_id)
} else if (!"lib_ID" %in% colnames(neu@meta.data)) {
neu$lib_ID <- factor(neu$orig.ident)
}

if ("site" %in% colnames(neu@meta.data)) {
neu$lib <- factor(neu$site)
} else if (!"lib" %in% colnames(neu@meta.data)) {
neu$lib <- neu$lib_ID
}

neu$lib_ID <- droplevels(neu$lib_ID)
neu$lib <- droplevels(neu$lib)

cat("retained libraries:", levels(neu$lib_ID), "\n")
cat("retained cells:", ncol(neu), "\n\n")

cat("running SCTransform v2...\n")
neu <- SCTransform(
neu,
assay = "RNA",
new.assay.name = "SCT",
vst.flavor = "v2",
min_cells = 5,
vars.to.regress = NULL,
verbose = TRUE
)

rna_feats <- nrow(GetAssayData(neu, assay = "RNA", layer = "counts"))
sct_counts_feats <- nrow(GetAssayData(neu, assay = "SCT", layer = "counts"))
sct_data_feats <- nrow(GetAssayData(neu, assay = "SCT", layer = "data"))
sct_scale_feats <- nrow(GetAssayData(neu, assay = "SCT", layer = "scale.data"))
sct_hvgs <- length(VariableFeatures(neu, assay = "SCT"))

cat("\nRNA assay features:", rna_feats, "\n")
cat("SCT assay counts features:", sct_counts_feats, "\n")
cat("SCT assay data features:", sct_data_feats, "\n")
cat("SCT assay scale.data features:", sct_scale_feats, "\n")
cat("SCT top HVGs selected:", sct_hvgs, "\n\n")

out_rds <- file.path(work_dir, "neutrotime_postQC_neutrophils_SCT.rds")
tmp_rds <- paste0(out_rds, ".tmp.", Sys.getpid())
saveRDS(neu, tmp_rds)
file.rename(tmp_rds, out_rds)
cat("sctransformed neutrophil object saved to:", out_rds, "\n\n")

scale_mat <- GetAssayData(neu, assay = "SCT", layer = "scale.data")
hvg_n <- nrow(scale_mat)

global_vals <- as.numeric(scale_mat)
global_res <- data.frame(
Min = round(min(global_vals), 5),
Q1 = round(as.numeric(quantile(global_vals, 0.25)), 5),
Median = round(median(global_vals), 5),
Mean = round(mean(global_vals), 5),
Q3 = round(as.numeric(quantile(global_vals, 0.75)), 5),
Max = round(max(global_vals), 5)
)

cat("overall sct scale.data residual summary:\n")
print(global_res, row.names = FALSE)

lib_ids <- as.character(neu$lib_ID)
lib_names <- as.character(neu$lib)

if ("nCount_SCT" %in% colnames(neu@meta.data)) {
sct_umi <- as.numeric(neu$nCount_SCT)
sct_feat <- as.numeric(neu$nFeature_SCT)
} else {
sct_counts <- GetAssayData(neu, assay = "SCT", layer = "counts")
sct_umi <- as.numeric(colSums(sct_counts))
sct_feat <- as.numeric(colSums(sct_counts > 0))
}

unique_libs <- unique(lib_ids)

summary_rows <- lapply(unique_libs, function(l) {
cell_idx <- which(lib_ids == l)
vals <- as.numeric(scale_mat[, cell_idx, drop = FALSE])

data.frame(
lib_ID = l,
sample = lib_names[cell_idx[1]],
total_cells = length(cell_idx),
hvg_in_scaledata = hvg_n,
median_corrected_umi = round(median(sct_umi[cell_idx]), 1),
median_corrected_features = round(median(sct_feat[cell_idx]), 1),
res_Min = round(min(vals), 5),
res_Q1 = round(as.numeric(quantile(vals, 0.25)), 5),
res_Median = round(median(vals), 5),
res_Mean = round(mean(vals), 5),
res_Q3 = round(as.numeric(quantile(vals, 0.75)), 5),
res_Max = round(max(vals), 5),
stringsAsFactors = FALSE
)
})

combined_summary <- do.call(rbind, summary_rows)

cat("\nconsolidated sct diagnostic summary:\n")
options(max.print = 500)
print(combined_summary)

out_tsv <- file.path(work_dir, "neutrotime_neutrophils_SCT_combined_summary.tsv")
write.table(combined_summary, file = out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\ncombined summary saved to:", out_tsv, "\n")
```

```
retained libraries: GSM5029335_BL_dataset1.txt.gz GSM5029336_BM_dataset1.txt.gz GSM5029337_SP_dataset1.txt.gz GSM5029338_BL_dataset2.txt.gz GSM5029339_BM_dataset2.txt.gz GSM5029340_SP_dataset2.txt.gz 
retained cells: 12845 

running SCTransform v2...
Running SCTransform on assay: RNA
vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
Calculating cell attributes from input UMI matrix: log_umi
Variance stabilizing transformation of count matrix of size 9735 by 12845
Model formula is y ~ log_umi
Get Negative Binomial regression parameters per gene
Using 2000 genes, 5000 cells
Found 358 outliers - those will be ignored in fitting/regularization step

Second step: Get residuals using fitted parameters for 9735 genes
Computing corrected count matrix for 9735 genes
Calculating gene attributes
Wall clock passed: Time difference of 35.13214 secs
Determine variable features
Centering data matrix

Set default assay to SCT

RNA assay features: 9735 
SCT assay counts features: 9735 
SCT assay data features: 9735 
SCT assay scale.data features: 3000 
SCT top HVGs selected: 3000 

[1] TRUE
sctransformed neutrophil object saved to: /global/scratch/hpc6297/neutrotime_output/neutrotime_postQC_neutrophils_SCT.rds 

overall sct scale.data residual summary:
      Min       Q1   Median Mean       Q3     Max
 -4.48685 -0.36132 -0.23887    0 -0.12964 20.7146

consolidated sct diagnostic summary:
                         lib_ID           sample total_cells hvg_in_scaledata median_corrected_umi median_corrected_features  res_Min   res_Q1 res_Median res_Mean   res_Q3  res_Max
1 GSM5029335_BL_dataset1.txt.gz peripheral blood        1343             3000                949.0                       364 -4.47493 -0.31385   -0.20456  0.04217 -0.11130 20.68878
2 GSM5029336_BM_dataset1.txt.gz      bone marrow        1151             3000               1198.0                       442 -4.32677 -0.45460   -0.31741 -0.09446 -0.21734 20.70423
3 GSM5029337_SP_dataset1.txt.gz           spleen        1154             3000               1072.5                       393 -4.32156 -0.40413   -0.26670 -0.02639 -0.14691 20.70423
4 GSM5029338_BL_dataset2.txt.gz peripheral blood        3819             3000                956.0                       376 -4.48685 -0.30363   -0.20024  0.05051 -0.11550 20.70142
5 GSM5029339_BM_dataset2.txt.gz      bone marrow        3321             3000               1202.0                       464 -4.38417 -0.41928   -0.29233 -0.07166 -0.20419 20.71460
6 GSM5029340_SP_dataset2.txt.gz           spleen        2057             3000                934.0                       370 -4.36107 -0.29643   -0.19179  0.06205 -0.10989 20.69828

combined summary saved to: /global/scratch/hpc6297/neutrotime_output/neutrotime_neutrophils_SCT_combined_summary.tsv 
```



