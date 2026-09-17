# R4.6.1 and dependencies
- for Neutrotime dataset onwards
- also can be used to rerun Liang 2022
- module R4.6.1 on frontenac CAC has internet access, allowing downloading on the fly for reference sets
- versions are locked as updates are not automatically installed
- freeze into a container at the end of analysis 

# check existing R4.5.3 dependencies and ext library
[R4.5.3_v2container_deps_20260907_150357.txt](https://github.com/user-attachments/files/31919634/R4.5.3_v2container_deps_20260907_150357.txt)

```
> sessionInfo()
R version 4.5.3 (2026-03-11)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.4 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
[1] C

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.5.3    cli_3.6.6         tools_4.5.3       sessioninfo_1.2.4
```

```
# extLib

Package              Version
BPCells              0.3.1  
RcppProgress         0.4.2  
SpatialArtifacts     1.1.0  
biglm                0.9-3  
harmony              2.0.5  
lisi                 1.0    
monocle3             1.4.29 
presto               1.0.0  
scIntegrationMetrics 1.2.0  
speedglm             0.3-5
```

# install all deps/ext lib from R4.5.3 into R4.6.1

```sh
cat << 'EOF' > install_R4.6.1_deps.slurm
#!/bin/bash
#SBATCH --job-name=install_r_deps
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/global/scratch/hpc6297/install_R4.6.1_deps_%j_progress.log
#SBATCH --error=/global/scratch/hpc6297/install_R4.6.1_deps_%j_error.log
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --mail-type=END,FAIL

module purge
module load StdEnv/2023 gcc udunits gdal geos r/4.6.1

Rscript -e '
options(Ncpus = 8)

target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
dir.create(target_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(target_lib, .libPaths()))

old_deps_file <- "R4.5.3_v2container_deps_20260907_150357.txt"
old_deps <- read.delim(old_deps_file, stringsAsFactors = FALSE)
old_pkgs <- unique(old_deps$Name)

base_recommended <- c(
  "base", "compiler", "datasets", "graphics", "grDevices", "grid",
  "methods", "parallel", "splines", "stats", "stats4", "tcltk",
  "tools", "utils"
)
old_pkgs <- setdiff(old_pkgs, base_recommended)

current_installed <- installed.packages(lib.loc = .libPaths())[, "Package"]
missing_pkgs <- setdiff(old_pkgs, current_installed)

message(sprintf("Found %d packages missing from current environment.", length(missing_pkgs)))

if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = target_lib, repos = "https://cloud.r-project.org")
  }
  BiocManager::install(
    missing_pkgs,
    lib = target_lib,
    update = FALSE,
    ask = FALSE,
    checkBuilt = FALSE
  )
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = target_lib, repos = "https://cloud.r-project.org")
}

github_targets <- c(
  "monocle3" = "cole-trapnell-lab/monocle3",
  "scIntegrationMetrics" = "carmonalab/scIntegrationMetrics",
  "SpatialArtifacts" = "MarioniLab/SpatialArtifacts"
)

for (pkg in names(github_targets)) {
  if (!pkg %in% installed.packages(lib.loc = .libPaths())[, "Package"]) {
    tryCatch(
      remotes::install_github(github_targets[[pkg]], lib = target_lib, upgrade = "never"),
      error = function(e) message(sprintf("%s installation failed: %s", pkg, conditionMessage(e)))
    )
  }
}

post_installed <- installed.packages(lib.loc = .libPaths())[, "Package"]
failed_pkgs <- setdiff(c(missing_pkgs, names(github_targets)), post_installed)

failed_log <- "/global/scratch/hpc6297/failed_deps_R4.6.1.txt"
if (length(failed_pkgs) > 0) {
  writeLines(failed_pkgs, con = failed_log)
  message(sprintf("%d packages failed to install. Logged to %s", length(failed_pkgs), failed_log))
} else {
  writeLines(character(0), con = failed_log)
  message("All packages installed successfully. Empty log created.")
}
'
EOF

sbatch install_R4.6.1_deps.slurm
```

```
16 packages failed to install. Logged to /global/scratch/hpc6297/failed_deps_R4.6.1.txt

arrow
CellChat
copykat
DoubletFinder
Giotto
GiottoClass
GiottoUtils
GiottoVisuals
infercnv
kBET
liana
RcisTarget
rjags
SCENIC
SeuratWrappers
SpatialArtifacts
```

# part 2
```r
cat << 'EOF' > install_R4.6.1_deps_part2.slurm
#!/bin/bash
#SBATCH --job-name=install_r_deps_part2
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=/global/scratch/hpc6297/install_R4.6.1_deps_part2_%j_progress.log
#SBATCH --error=/global/scratch/hpc6297/install_R4.6.1_deps_part2_%j_error.log
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --mail-type=END,FAIL

module purge
module load StdEnv/2023 gcc udunits gdal geos jags r/4.6.1

Rscript -e '
options(Ncpus = 8)

target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
dir.create(target_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(target_lib, .libPaths()))

is_installed <- function(pkg) {
  pkg %in% installed.packages(lib.loc = .libPaths())[, "Package"]
}

# Read the remaining 16 failed packages from part 1 log
failed_file <- "/global/scratch/hpc6297/failed_deps_R4.6.1.txt"
if (file.exists(failed_file)) {
  target_pkgs <- readLines(failed_file)
  target_pkgs <- target_pkgs[nzchar(target_pkgs)]
} else {
  message("Failed log from part 1 not found. Aborting part 2.")
  quit(status = 1)
}

# 1. arrow with static thrift/c++ engine
if ("arrow" %in% target_pkgs && !is_installed("arrow")) {
  message(">>> Installing arrow...")
  Sys.setenv(ARROW_WITH_THRIFT = "ON", LIBARROW_BUILD = "true")
  install.packages("arrow", lib = target_lib, repos = "https://cloud.r-project.org")
} else if ("arrow" %in% target_pkgs) {
  message(">>> Already installed: arrow")
}

# 2. rjags (system jags loaded via module)
if ("rjags" %in% target_pkgs && !is_installed("rjags")) {
  message(">>> Installing rjags...")
  install.packages("rjags", lib = target_lib, repos = "https://cloud.r-project.org")
} else if ("rjags" %in% target_pkgs) {
  message(">>> Already installed: rjags")
}

# 3. Bioconductor targets
bioc_pkgs <- c("RcisTarget", "infercnv")
bioc_targets <- intersect(bioc_pkgs, target_pkgs)

if (length(bioc_targets) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE, lib.loc = .libPaths())) {
    install.packages("BiocManager", lib = target_lib, repos = "https://cloud.r-project.org")
  }
  for (bpkg in bioc_targets) {
    if (!is_installed(bpkg)) {
      message(sprintf(">>> Installing Bioc package: %s", bpkg))
      BiocManager::install(bpkg, lib = target_lib, update = FALSE, ask = FALSE)
    } else {
      message(sprintf(">>> Already installed: %s", bpkg))
    }
  }
}

# 4. GitHub targets
github_pkgs <- c(
  "CellChat"         = "sqjin/CellChat",
  "copykat"          = "navinlabcode/copykat",
  "DoubletFinder"    = "chris-mcginnis-ucsf/DoubletFinder",
  "kBET"             = "theislab/kBET",
  "liana"            = "saezlab/liana",
  "SCENIC"           = "aertslab/SCENIC",
  "SeuratWrappers"   = "satijalab/seurat-wrappers",
  "SpatialArtifacts" = "CambridgeCat13/SpatialArtifacts",
  "GiottoClass"      = "giotto-suite/GiottoClass",
  "GiottoUtils"      = "giotto-suite/GiottoUtils",
  "GiottoVisuals"    = "giotto-suite/GiottoVisuals",
  "Giotto"           = "giotto-suite/Giotto"
)

github_targets <- intersect(names(github_pkgs), target_pkgs)

if (length(github_targets) > 0) {
  if (!requireNamespace("remotes", quietly = TRUE, lib.loc = .libPaths())) {
    install.packages("remotes", lib = target_lib, repos = "https://cloud.r-project.org")
  }
  for (pkg in github_targets) {
    if (!is_installed(pkg)) {
      message(sprintf(">>> Installing GitHub package: %s (%s)", pkg, github_pkgs[[pkg]]))
      tryCatch(
        remotes::install_github(github_pkgs[[pkg]], lib = target_lib, upgrade = "never"),
        error = function(e) message(sprintf("FAILED: %s -> %s", pkg, conditionMessage(e)))
      )
    } else {
      message(sprintf(">>> Already installed: %s", pkg))
    }
  }
}

# 5. verify residual failures
missing_final <- target_pkgs[!sapply(target_pkgs, is_installed)]
failed_log_part2 <- "/global/scratch/hpc6297/failed_deps_R4.6.1_part2.txt"

if (length(missing_final) > 0) {
  writeLines(missing_final, con = failed_log_part2)
  message(sprintf("Finished. %d packages failed again, logged to %s", length(missing_final), failed_log_part2))
} else {
  writeLines(character(0), con = failed_log_part2)
  message("All 16 missing packages installed successfully. Empty log created.")
}
'
EOF

sbatch install_R4.6.1_deps_part2.slurm
```

```
# failed_deps_R4.6.1_part2.txt
arrow
RcisTarget
SCENIC
```

# install last libraries on the fly within R

```r
target_lib <- "/global/home/hpc6297/R/x86_64-pc-linux-gnu-library/4.6"
options(Ncpus = 8)

# 1. rcistarget
BiocManager::install("RcisTarget", lib = target_lib, update = FALSE, ask = FALSE)

# 2. scenic
remotes::install_github("aertslab/SCENIC", lib = target_lib, upgrade = "never")

# 3. arrow (easiest way for HPC)
Sys.setenv(NOT_CRAN = "true")
install.packages("arrow", lib = target_lib, repos = "https://cloud.r-project.org")

# verification
sapply(c("RcisTarget", "SCENIC", "arrow"), library, character.only = TRUE, lib.loc = target_lib)
```

- [arrow installation troubleshoot](https://arrow.apache.org/docs/r/articles/install.html)
