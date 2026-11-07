#!/usr/bin/env bash
# =============================================================================================
# install.sh: create the "seurat_env" conda environment for the DeCodeLife scRNA-seq notebooks
#
# Target: WSL2 (Ubuntu or another Linux distribution). Plain Linux works too.
# Installs: R 4.4.3, Seurat 5, SeuratObject, SoupX, scDblFinder, SingleCellExperiment, harmony,
#           Matrix, dplyr, tidyr, ggplot2, patchwork, plus JupyterLab, Notebook, nbconvert and IRkernel.
# Registers the Jupyter kernel "seurat_env" (shown as "R 4.4.3 (seurat_env)"); the notebooks use this kernel.
#
# Usage (in a WSL2 terminal, not PowerShell or Git Bash):
#   bash install.sh              # create the environment
#   FORCE=1 bash install.sh      # delete an existing seurat_env first, then rebuild it
#
# Afterwards, from the project root (the folder that holds this script); the notebooks are in scripts/seurat_analysis/:
#   conda activate seurat_env && jupyter lab --no-browser
#       -> open the printed http://localhost:8888/... link in a Windows browser
#   VS Code (WSL window): open a notebook and select the kernel "R 4.4.3 (seurat_env)"
#   conda run -n seurat_env jupyter nbconvert --to notebook --execute --inplace scripts/seurat_analysis/preprocess_preInt.ipynb
#       -> runs a notebook without opening it
#
# If bash reports "$'\r': command not found", the file got Windows line endings: sed -i 's/\r$//' install.sh
# =============================================================================================
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # this script lives in the project root
ENV_NAME="seurat_env"
R_VERSION="4.4.3"
KERNEL_NAME="seurat_env"                        # must match the kernelspec "name" in the notebooks
KERNEL_DISPLAY="R ${R_VERSION} (seurat_env)"

# conda-forge first, then bioconda (Bioconductor packages), with strict priority. --override-channels ignores the
# channels in ~/.condarc, so channels such as "r", "defaults" or "pytorch" cannot mix incompatible R builds in.
CHANNELS=(--override-channels -c conda-forge -c bioconda --strict-channel-priority)

PACKAGES=(
  "r-base=${R_VERSION}"
  # single-cell analysis
  "r-seurat>=5.0"
  "r-seuratobject>=5.0"
  r-soupx                              # ambient RNA correction (needs Cell Ranger's raw_feature_bc_matrix)
  bioconductor-scdblfinder             # doublet detection
  bioconductor-singlecellexperiment
  bioconductor-summarizedexperiment
  bioconductor-biocparallel
  r-harmony                            # batch integration (next analysis step)
  # data handling and plotting
  r-matrix
  r-dplyr
  r-tidyr
  r-ggplot2
  r-patchwork
  # Jupyter: notebooks, headless execution, and the R kernel
  r-irkernel
  jupyterlab
  notebook
  nbconvert
)

info() { printf '\n\033[1;34m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[1;33mWARNING: %s\033[0m\n' "$*" >&2; }
die()  { printf '\033[1;31mERROR: %s\033[0m\n' "$*" >&2; exit 1; }

# ---------------------------------------------------------------------------------------------
# 1. Platform: Linux / WSL2
# ---------------------------------------------------------------------------------------------
info "Checking the platform"
[[ "$(uname -s)" == "Linux" ]] || die "Run this script inside WSL2 (or Linux), not from Windows PowerShell / Git Bash."
if grep -qi microsoft /proc/version 2>/dev/null; then
  echo "WSL detected: $(uname -r)"
else
  warn "This does not look like WSL2; continuing as a plain Linux install."
fi

# ---------------------------------------------------------------------------------------------
# 2. Find conda. A non-interactive shell does not read ~/.bashrc, so "conda init" may not have put conda on PATH.
# ---------------------------------------------------------------------------------------------
info "Looking for conda"
CONDA=""
for c in "${CONDA_EXE:-}" "$(command -v conda 2>/dev/null || true)" \
         "$HOME/miniforge3/bin/conda" "$HOME/mambaforge/bin/conda" "$HOME/miniconda3/bin/conda" \
         "$HOME/anaconda3/bin/conda" "/opt/conda/bin/conda"; do
  if [[ -n "$c" && -x "$c" ]]; then CONDA="$c"; break; fi
done
if [[ -z "$CONDA" ]]; then
  die "conda not found. Install Miniforge in WSL2 first, then re-run this script:
    curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
    bash Miniforge3-Linux-x86_64.sh -b -p \$HOME/miniforge3 && \$HOME/miniforge3/bin/conda init bash"
fi
CONDA_BASE="$("$CONDA" info --base)"
echo "conda: $CONDA ($("$CONDA" --version))"

# mamba solves and downloads faster; plain conda (libmamba solver) works too
SOLVER="$CONDA"
if [[ -x "$CONDA_BASE/bin/mamba" ]]; then SOLVER="$CONDA_BASE/bin/mamba"; fi
echo "solver: $SOLVER"

# Rough disk-space check (the environment needs about 3-4 GB)
free_gb="$(df -Pk "$CONDA_BASE" | awk 'NR == 2 {printf "%d", $4 / 1024 / 1024}')"
if (( free_gb < 6 )); then warn "Only ${free_gb} GB free under $CONDA_BASE; the environment needs about 3-4 GB."; fi

# ---------------------------------------------------------------------------------------------
# 3. Existing environment
# ---------------------------------------------------------------------------------------------
env_names="$("$CONDA" env list | awk '!/^#/ && NF {print $1}')"
if grep -qx "$ENV_NAME" <<< "$env_names"; then
  if [[ "${FORCE:-0}" == "1" ]]; then
    info "Removing the existing '$ENV_NAME' environment (FORCE=1)"
    "$CONDA" remove -n "$ENV_NAME" --all --yes
  else
    die "A conda environment named '$ENV_NAME' already exists. Re-run with FORCE=1 to rebuild it, or remove it with:
    conda remove -n $ENV_NAME --all"
  fi
fi

# ---------------------------------------------------------------------------------------------
# 4. Create the environment
# ---------------------------------------------------------------------------------------------
info "Creating '$ENV_NAME' with R $R_VERSION (this can take 5-15 minutes)"
"$SOLVER" create -n "$ENV_NAME" --yes "${CHANNELS[@]}" "${PACKAGES[@]}"

run_in_env() { "$CONDA" run -n "$ENV_NAME" --no-capture-output "$@"; }
ENV_PREFIX="$(run_in_env printenv CONDA_PREFIX)"
echo "environment: $ENV_PREFIX"

# While seurat_env is activated, later "conda install" calls also prefer conda-forge / bioconda with strict priority
# (from another environment, pass the channels explicitly; see the hint printed at the end)
"$CONDA" config --file "$ENV_PREFIX/.condarc" --add channels bioconda
"$CONDA" config --file "$ENV_PREFIX/.condarc" --add channels conda-forge
"$CONDA" config --file "$ENV_PREFIX/.condarc" --set channel_priority strict

# ---------------------------------------------------------------------------------------------
# 5. Check the R version and that every package loads
# ---------------------------------------------------------------------------------------------
info "Checking R and the R packages"
run_in_env Rscript --vanilla -e '
pkgs <- c("Seurat", "SeuratObject", "SoupX", "scDblFinder", "SingleCellExperiment", "SummarizedExperiment",
          "BiocParallel", "harmony", "Matrix", "dplyr", "tidyr", "ggplot2", "patchwork", "IRkernel")
failed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed) > 0) stop("These packages failed to load: ", paste(failed, collapse = ", "))
cat(R.version.string, "\n")
for (p in pkgs) cat(sprintf("  %-22s %s\n", p, as.character(packageVersion(p))))
if (getRversion() != "'"$R_VERSION"'") stop("Expected R '"$R_VERSION"', got ", as.character(getRversion()))
'

# ---------------------------------------------------------------------------------------------
# 6. Register the Jupyter kernel for the current user (~/.local/share/jupyter/kernels/seurat_env).
#    A user-level kernel is visible to every Jupyter install and to VS Code in a WSL window. It has its own name,
#    so an existing "ir" kernel (e.g. the system R in /usr/lib/R) is left untouched.
# ---------------------------------------------------------------------------------------------
info "Registering the Jupyter kernel '$KERNEL_DISPLAY'"
run_in_env Rscript --vanilla -e "IRkernel::installspec(name = '$KERNEL_NAME', displayname = '$KERNEL_DISPLAY', user = TRUE)"
run_in_env jupyter kernelspec list

# ---------------------------------------------------------------------------------------------
# 7. Run a one-cell test notebook on the new kernel (the same way nbconvert would run the analysis notebooks)
# ---------------------------------------------------------------------------------------------
info "Testing notebook execution with the '$KERNEL_NAME' kernel"
TEST_DIR="$(mktemp -d)"
trap 'rm -rf "$TEST_DIR"' EXIT
cat > "$TEST_DIR/kernel_test.ipynb" <<EOF
{
 "cells": [
  {"cell_type": "code", "execution_count": null, "metadata": {}, "outputs": [],
   "source": ["suppressPackageStartupMessages(library(Seurat))\n",
              "cat(R.version.string, '| Seurat', as.character(packageVersion('Seurat')), '\\\\n')"]}
 ],
 "metadata": {"kernelspec": {"display_name": "$KERNEL_DISPLAY", "language": "R", "name": "$KERNEL_NAME"}},
 "nbformat": 4,
 "nbformat_minor": 4
}
EOF
run_in_env jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=600 \
  --output-dir "$TEST_DIR" --output kernel_test_out.ipynb "$TEST_DIR/kernel_test.ipynb"
if grep -q "R version ${R_VERSION}" "$TEST_DIR/kernel_test_out.ipynb"; then
  echo "Notebook test passed: $(grep -o "R version ${R_VERSION}[^\"]*" "$TEST_DIR/kernel_test_out.ipynb" | head -1 | sed 's/\\n$//')"
else
  die "The test notebook ran, but did not report R ${R_VERSION}. Check: jupyter kernelspec list"
fi

# ---------------------------------------------------------------------------------------------
info "Done"
cat <<EOF
Environment : $ENV_NAME  ($ENV_PREFIX)
Jupyter     : kernel "$KERNEL_DISPLAY" (name: $KERNEL_NAME)

Next steps (the notebooks are in scripts/seurat_analysis/ and find the project root themselves):
  cd "$PROJECT_ROOT"
  conda activate $ENV_NAME
  jupyter lab --no-browser        # then open the printed http://localhost:8888/... link in a Windows browser

  VS Code: open the project folder in a WSL window (Remote - WSL), open a notebook in scripts/seurat_analysis/
           and pick "$KERNEL_DISPLAY" as the kernel. Install the Jupyter extension if VS Code asks for it.

  Run a notebook without opening it (outputs are written back into the .ipynb):
  conda run -n $ENV_NAME jupyter nbconvert --to notebook --execute --inplace scripts/seurat_analysis/preprocess_preInt.ipynb

To add packages later, keep the same channels:
  conda install -n $ENV_NAME --override-channels -c conda-forge -c bioconda <package>
EOF
