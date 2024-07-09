# SmartQC

## Installation

1. First, install the required dependencies:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("preprocessCore", "BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "HDF5Array",
                       "terra", "ggrastr", "densvis", "biomaRt"))

# Install GitHub dependencies
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github("arc85/singleseqgset")
devtools::install_github('cole-trapnell-lab/monocle3')

# Install SmartQC
devtools::install_github('your_github_username/smartqc')
```

2. After installation, add the SmartQC executable to your PATH:

```R
# Run this in R to get the installation path
smartqc_path <- system.file("bin", package = "smartqc")
cat(paste0("\nAdd the following line to your ~/.bashrc or ~/.bash_profile file:\n\n",
           "export PATH=\"$PATH:", smartqc_path, "\"\n\n",
           "Then, run 'source ~/.bashrc' or 'source ~/.bash_profile' in your terminal.\n"))
```

3. Copy the output line and add it to your `~/.bashrc` or `~/.bash_profile` file.

4. Reload your shell configuration:

```bash
source ~/.bashrc  # or source ~/.bash_profile
```

## Usage

After installation and adding SmartQC to your PATH, you can use it from the command line:

```bash
smartqc -r /path/to/root -s StudyName -g MOUSE
```

or

```bash
smartqc -i /path/to/input -o /path/to/output.rds -g HUMAN
```

For more options, run:

```bash
smartqc --help
```

## Troubleshooting

If you encounter any issues with the installation or usage of SmartQC, please check the following:

1. Ensure all dependencies are correctly installed.
2. Verify that the SmartQC bin directory is in your PATH.
3. If you're using RStudio or a specific R environment, make sure it's aware of the updated PATH.

For further assistance, please open an issue on the GitHub repository.
