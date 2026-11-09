# scRNAseq_basics

Basics of single-cell RNA-seq (scRNA-seq) data analysis with **Seurat 5 (R)**, run from Jupyter notebooks on **WSL2 or Linux**.

| Notebook | What it covers | Input |
|---|---|---|
| `scripts/seurat_analysis/analysis_stage1.ipynb` | Tutorial: count matrix → ambient-RNA check (SoupX) → doublets (scDblFinder) → QC metrics → MAD-based filtering → normalization → highly variable genes → scaling, with extensive explanations | 10x PBMC dataset |
| `scripts/seurat_analysis/preprocess_preInt.ipynb` | Two NSCLC datasets: merging, the same preprocessing, and PCA-based tests of whether the structure is driven by batch (dataset / sample) or by biology, before integration | 10x CellPlex NSCLC dataset + a public GEO dataset |
| `scripts/scanpy_analysis/*.ipynb` | Python / Scanpy versions of earlier sessions (not covered by `install.sh`) | |

Both Seurat notebooks end with a **runtime-per-stage** table and plot. They only **estimate** ambient RNA contamination with SoupX and do not correct the counts. Each notebook has a "Full analysis: include the correction" box with the commands for a complete run.

---

## 1. Repository layout

```
scRNAseq_basics/                 <- project root (the notebooks find it by looking for install.sh + scripts/)
├── install.sh                   <- creates the conda environment "seurat_env" and the Jupyter kernel
├── scripts/
│   ├── seurat_analysis/         <- R / Seurat notebooks
│   └── scanpy_analysis/         <- Python / Scanpy notebooks
├── NSCLC_7pt/                   <- input: NSCLC 10x data, one folder per sample (git-ignored, download yourself)
│   ├── P1/  barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
│   └── P2/ … P4/
└── data/                        <- input PBMC / GEO data + ALL results (git-ignored)
    ├── pbmc10_data/             <- input: 10x PBMC filtered matrix
    ├── GSE<ID>/                 <- input: GEO count matrices
    ├── figures/  tables/        <- results written by the notebooks
    └── *.rds                    <- saved Seurat objects
```

All paths in the notebooks are **relative to the project root**, so the repository works wherever it is cloned. `data/`, `NSCLC_7pt/`, `*.rds` and slide decks are git-ignored, so **no data or results are pushed to GitHub**.

---

## 2. Installation (WSL2 or Linux)

### 2.1 Requirements

| | Minimum | Recommended |
|---|---|---|
| OS | Linux x86-64, or Windows 10/11 with WSL2 (Ubuntu 22.04 / 24.04) | |
| Disk | ~6 GB for the conda environment + ~3 GB for the data | 15 GB free |
| RAM | 8 GB (PBMC tutorial) | 16 GB or more for `preprocess_preInt.ipynb` |
| Internet | to download packages and data | |

Use a **WSL2 terminal**, not PowerShell or Git Bash, for everything below.

### 2.2 Set up WSL2 (Windows only)

In **PowerShell as Administrator**:

```powershell
wsl --install -d Ubuntu-22.04
```

Restart when asked, open "Ubuntu" from the Start menu and create a Linux user. Then update the system and install the basic tools:

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y git curl wget ca-certificates
```

> **Performance tip:** clone the repository into the Linux file system (e.g. `~/projects/`) rather than under `/mnt/c/...`. File access across the Windows / Linux boundary is much slower, and reading the large text matrices is noticeably faster on the Linux side.

### 2.3 Install conda (Miniforge)

`install.sh` needs `conda`. If you do not have it yet, install **Miniforge**, which uses the conda-forge channel and includes `mamba`:

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
exec bash                       # reload the shell so that "conda" is on PATH
conda --version
```

### 2.4 Get the code

```bash
mkdir -p ~/projects && cd ~/projects
git clone <REPOSITORY_URL> scRNAseq_basics
cd scRNAseq_basics
```

### 2.5 Create the environment

From the project root:

```bash
bash install.sh                 # create "seurat_env" (5–15 minutes)
FORCE=1 bash install.sh         # delete an existing seurat_env and rebuild it from scratch
```

`install.sh`:

1. checks that it runs on Linux / WSL2 and finds `conda` (and `mamba`, if present);
2. creates the conda environment **`seurat_env`** with **R 4.4.3** from `conda-forge` + `bioconda` only (`--override-channels`, strict channel priority). It installs Seurat 5, SeuratObject, SoupX, scDblFinder, SingleCellExperiment, BiocParallel, harmony, leidenbase, Matrix, dplyr, tidyr, ggplot2, patchwork, IRkernel, JupyterLab, Notebook and nbconvert;
3. checks that every R package loads and that R is 4.4.3;
4. registers the Jupyter kernel **"R 4.4.3 (seurat_env)"** for your user;
5. runs a one-cell test notebook on that kernel.

If bash reports `$'\r': command not found`, the script has Windows line endings. Fix them with `sed -i 's/\r$//' install.sh`.

---

## 3. Download the datasets

The data are **not** in the repository. Download them into the folders below; the names must match, because the notebooks look for them there. Replace the `<...>` placeholders with the dataset you use.

### 3.1 PBMC 10k (10x Genomics), for `analysis_stage1.ipynb`

The Cell Ranger **filtered feature-barcode matrix** (MEX format) of a ~10,000-cell PBMC dataset from the 10x Genomics dataset portal (<https://www.10xgenomics.com/datasets>).

### 3.2 NSCLC, 10x CellPlex (multiplexed), for `preprocess_preInt.ipynb`

A 10x **Cell Ranger `multi`** dataset in which several NSCLC donors were pooled with CellPlex (CMO tags). It has gene expression, TotalSeq antibodies and CMO tags. The notebook needs **one folder per donor** with that donor's `sample_filtered_feature_bc_matrix` files:

```bash
mkdir -p NSCLC_7pt
# Download the per-sample outputs of the multiplexed run from the 10x dataset portal:
#   "<10X_NSCLC_CELLPLEX_DATASET_URL>"
# For each donor, extract its sample_filtered_feature_bc_matrix into NSCLC_7pt/P<n>/, e.g.:
for n in 1 2 3 4; do
  mkdir -p NSCLC_7pt/P$n
  wget -O /tmp/P$n.tar.gz "<10X_NSCLC_SAMPLE_${n}_FILTERED_MATRIX_TAR_GZ_URL>"
  tar -xzf /tmp/P$n.tar.gz -C NSCLC_7pt/P$n --strip-components=1
done
ls NSCLC_7pt/P1                 # barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

| Folder | Donor / CellPlex tag |
|---|---|
| `NSCLC_7pt/P1` | `<CMO_TAG_1>` (e.g. CMO301) |
| `NSCLC_7pt/P2` | `<CMO_TAG_2>` |
| `NSCLC_7pt/P3` | `<CMO_TAG_3>` |
| `NSCLC_7pt/P4` | `<CMO_TAG_4>` |

To use more donors, add their folders and extend `sample_names` in the notebook's paths cell. All folders must share the same `features.tsv.gz`.

### 3.3 Public NSCLC dataset from GEO, for `preprocess_preInt.ipynb`

The second dataset is a GEO series that provides **one dense, tab-separated genes × barcodes matrix of integer UMI counts per patient**, named `GSM<id>_P<n>_exp.txt`. The header row lists only the barcodes, and the first column holds the gene symbols. GEO stores supplementary files under `series/<first digits>nnn/<GSE ID>/suppl/`:

```bash
GSE_ID="<GSE_ACCESSION>"                        # e.g. GSE123456
GSE_DIR="${GSE_ID:0:${#GSE_ID}-3}nnn"           # GSE123456 -> GSE123nnn
mkdir -p "data/$GSE_ID"
wget -O "/tmp/${GSE_ID}_RAW.tar" "https://ftp.ncbi.nlm.nih.gov/geo/series/${GSE_DIR}/${GSE_ID}/suppl/${GSE_ID}_RAW.tar"

tar -tf "/tmp/${GSE_ID}_RAW.tar"                # list the per-sample files
# Extract only the patients you want, then decompress (the notebook reads plain .txt files)
tar -xf "/tmp/${GSE_ID}_RAW.tar" -C "data/$GSE_ID" "<GSM_ID_1>_P1_exp.txt.gz" "<GSM_ID_2>_P3_exp.txt.gz"
gunzip data/"$GSE_ID"/*.txt.gz
ls -lh "data/$GSE_ID"
```

In the notebook's paths cell, set `gse_dir` and `gse_id` to the same accession. Every `GSM*_exp.txt` file in that folder is loaded, and the patients are labelled `<GSE_ID>_P<n>`. The four files used in the notebook need ~1.2 GB of disk; reading the largest one needs ~1 GB of RAM.

---

## 4. Run the analysis

### 4.1 JupyterLab

```bash
cd ~/projects/scRNAseq_basics
conda activate seurat_env
jupyter lab --no-browser        # open the printed http://localhost:8888/... link in a Windows browser
```

Open the notebooks in `scripts/seurat_analysis/` and select the kernel **"R 4.4.3 (seurat_env)"**. Do not use an existing `ir` kernel; it may point to a different R installation.

### 4.2 VS Code

Open the project folder in a **WSL window** (extension "WSL"/"Remote - WSL"), install the Jupyter extension if asked, open a notebook and pick the kernel **"R 4.4.3 (seurat_env)"**.

> If you edit a notebook outside VS Code while it is open, VS Code may save its old copy over your changes. Close the notebook (or use **Revert File**) first.

### 4.3 Headless, top to bottom (best for runtimes)

```bash
conda run -n seurat_env jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=seurat_env --ExecutePreprocessor.timeout=-1 \
  scripts/seurat_analysis/analysis_stage1.ipynb

conda run -n seurat_env jupyter nbconvert --to notebook --execute --inplace \
  --ExecutePreprocessor.kernel_name=seurat_env --ExecutePreprocessor.timeout=-1 \
  scripts/seurat_analysis/preprocess_preInt.ipynb
```

The runtime-per-stage table at the end of each notebook is only meaningful for an uninterrupted run like this.

### 4.4 Outputs

| Notebook | Figures | Tables | Objects |
|---|---|---|---|
| `analysis_stage1.ipynb` | `data/figures/pbmc10k_stage1/` | `data/tables/pbmc10k_stage1/` | `data/pbmc10k_stage1_preprocessed.rds` |
| `preprocess_preInt.ipynb` | `data/figures/NSCLC_combined/` | `data/tables/NSCLC_combined/` | `data/nsclc7pt_merged.rds`, `data/nsclc_combined_merged.rds`, `data/nsclc_combined_preInt.rds` |

---

## 5. A reproducible environment in a clean WSL / Linux system

`install.sh` pins **R 4.4.3** but takes the newest compatible versions of all other packages, so two installs made months apart can differ. To reproduce an analysis exactly, **freeze** the environment once it works and rebuild from that file.

### 5.1 Freeze a working environment

```bash
conda activate seurat_env
# Exact package URLs + checksums (linux-64); the most reproducible option
conda list -n seurat_env --explicit --md5 > envs/seurat_env.linux-64.lock.txt
# Human-readable version list (portable, but re-solved on install)
conda env export -n seurat_env --no-builds | grep -v '^prefix:' > envs/seurat_env.yml
# R-side record: R, Seurat and every loaded package version
Rscript -e 'sessionInfo()' > envs/seurat_env.sessionInfo.txt
```

Commit the files in `envs/` (create the folder first with `mkdir -p envs`). They are small text files.

### 5.2 Rebuild it on a clean machine

On a fresh WSL2 Ubuntu or Linux install: install Miniforge (Section 2.3), clone the repository (Section 2.4), then:

```bash
cd ~/projects/scRNAseq_basics
# Exact rebuild from the lock file (no dependency solving; same platform: linux-64)
conda create -n seurat_env --file envs/seurat_env.linux-64.lock.txt
# ... or, if no lock file exists, the slower re-solved version:
# conda env create -n seurat_env -f envs/seurat_env.yml

# Register the Jupyter kernel (install.sh does this automatically on a normal install)
conda run -n seurat_env Rscript -e 'IRkernel::installspec(name = "seurat_env", displayname = "R 4.4.3 (seurat_env)", user = TRUE)'
jupyter kernelspec list         # seurat_env should be listed
```

Then download the data (Section 3) and run the notebooks headless (Section 4.3).

### 5.3 What else keeps results reproducible

* **Random seeds:** every random step (scDblFinder, module scores, silhouettes, integer rounding) uses the seed set in the setup cell (`seed <- 198716`).
* **Single-threaded doublet detection:** scDblFinder runs with `BiocParallel::SerialParam(RNGseed = seed)`, so results do not depend on the number of cores.
* **Versions in the output:** each notebook ends with `sessionInfo()`. Keep the executed notebook (or its `sessionInfo()` output) together with the results.
* **Same input files:** record the dataset URLs / GEO accessions you used (replace the placeholders above). If possible, also record checksums: `md5sum data/pbmc10_data/* NSCLC_7pt/*/* data/GSE*/* > envs/input_md5.txt`.
* **BLAS threads:** for bit-identical PCA results across machines, fix the number of BLAS threads before starting Jupyter: `export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.

---

## 6. Troubleshooting

| Problem | Fix |
|---|---|
| `conda: command not found` | Run `exec bash` after `conda init`, or call `~/miniforge3/bin/conda` directly. `install.sh` also searches the usual install locations. |
| `$'\r': command not found` | `sed -i 's/\r$//' install.sh` (Windows line endings) |
| `A conda environment named 'seurat_env' already exists` | `FORCE=1 bash install.sh` |
| Kernel "R 4.4.3 (seurat_env)" not listed | `conda run -n seurat_env Rscript -e 'IRkernel::installspec(name = "seurat_env", displayname = "R 4.4.3 (seurat_env)", user = TRUE)'`, then restart Jupyter / VS Code |
| `Project root (folder with install.sh and scripts/) not found` | Start Jupyter inside the repository, and keep `install.sh` and `scripts/` in the project root |
| `missing barcodes.tsv.gz ...` / `No GSM*_exp.txt files` | Check the folder names and file names against Section 3 |
| R is killed / runs out of memory | Close other programs, or give WSL2 more memory in `%UserProfile%\.wslconfig` (`[wsl2]` then `memory=16GB`), then run `wsl --shutdown` in PowerShell |

---

## License

See [LICENSE](LICENSE).
