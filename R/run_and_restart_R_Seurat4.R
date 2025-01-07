if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
options(repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("SeuratObject", "4.1.4")
remotes::install_version("Seurat", "4.4.0", upgrade = FALSE) 