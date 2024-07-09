#' Run SmartQC pipeline
#'
#' @param input Input directory with raw data (for single sample processing)
#' @param output Output file name (for single sample processing)
#' @param root Root directory containing Studies folder (for multi-sample processing)
#' @param study Study name (for multi-sample processing)
#' @param organism Organism (MOUSE, HUMAN, or RAT)
#'
#' @import Seurat
#' @import DoubletFinder
#' @import ggplot2
#' @import dplyr
#' @import cluster
#' @import ggstatsplot
#' @import SCINA
#' @import plyr
#' @import monocle3
#' @import singleseqgset
#' @import SeuratWrappers
#' @import scCustomize
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggExtra ggMarginal
#' @importFrom cowplot plot_grid
#' @importFrom biomaRt useMart getBM
#' @importFrom densvis densmap
#' @importFrom pheatmap pheatmap
#' @importFrom GSA GSA
#' @importFrom preprocessCore normalize.quantiles
#'
#' @export
run_smartqc <- function(input = NULL, output = NULL, root = NULL, study = NULL, organism = "MOUSE") {
  #silence package start up and warnings noise but ##TODO: find fix
  old_options <- options(
    warn = -1, 
    readr.show_progress = FALSE,  
    dplyr.summarise.inform = FALSE 
  )
  on.exit(options(old_options))  #rstore options when function exits

  #get the markers and cell cycle genesw
  data <- load_necessary_data(organism)
  cellcyclegenes <- data$cell_cycle_genes
  
if (!is.null(input)) {
  #single sample processing
  message(sprintf("Processing single sample: %s", input))
  output_dir <- input
  sample_name <- basename(input)
  output_file <- file.path(output_dir, paste0("SeuratObject_", sample_name, ".rds"))
  
  if (file.exists(output_file)) {
    message(sprintf("Processed file already exists: %s. Skipping.", output_file))
  } else {
    seur <- process_sample(input, organism, cellcyclegenes, output_dir)
    if (!is.null(seur)) {
      saveRDS(seur, file = output_file)
      message(sprintf("Processed data saved to: %s", output_file))
    }
  }
} else if (!is.null(root) && !is.null(study)) {
  #multi-sample processing
  message(sprintf("Processing study: %s in %s", study, root))
  
  study_dir <- file.path(root, study)
  results_dir <- file.path(study_dir, "results")
  message(sprintf("Study directory: %s", study_dir))
  message(sprintf("Results directory: %s", results_dir))
  
  #check if study directory exists
  if (!dir.exists(study_dir)) {
    stop(sprintf("Study directory does not exist: %s", study_dir))
  }
  
  #list all files and directories in the study directory
  all_contents <- list.files(study_dir, full.names = TRUE)
  message(sprintf("All contents of study directory: %s", paste(all_contents, collapse = ", ")))
  
  #list only directories
  dir.list <- list.dirs(study_dir, full.names = TRUE, recursive = FALSE)
  message(sprintf("Directories in study directory: %s", paste(dir.list, collapse = ", ")))
  
  #filter out the results directory if it exists
  dir.list <- dir.list[!grepl("results$", dir.list)]
 #message(sprintf("Filtered directories: %s", paste(dir.list, collapse = ", ")))
  
  #exract organs
  organs <- unique(sapply(strsplit(basename(dir.list), '-'), `[`, 1))
  message(sprintf("Identified organs: %s", paste(organs, collapse = ", ")))
  
  if (length(organs) == 0) {
    stop("No organs identified. Check the directory structure.")
  }
  
  for (organ in organs) {
    message(sprintf("Processing organ: %s", organ))
    
    organ_dirs <- dir.list[grepl(paste0("^", organ), basename(dir.list))]
    message(sprintf("Organ directories: %s", paste(organ_dirs, collapse = ", ")))
    
    if (length(organ_dirs) == 0) {
      warning(sprintf("No directories found for organ: %s", organ))
      next
    }
    
    for (inp.dir in organ_dirs) {
      sample <- basename(inp.dir)
      message(sprintf("Processing sample: %s", sample))
      message(sprintf("Input directory: %s", inp.dir))
      
      #ceck if input directory exists and has expected files
      if (!dir.exists(inp.dir)) {
        warning(sprintf("Input directory does not exist: %s", inp.dir))
        next
      }
      ##TO DO: RESTORE THE NOT ONLY CELLRANGER FUNCIONALITY
      ##DONE, SEE IF REMOVE THIS
     # expected_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
      #existing_files <- list.files(inp.dir)
      #message(sprintf("Files in input directory: %s", paste(existing_files, collapse = ", ")))
      
      #if (!all(expected_files %in% existing_files)) {
       # warning(sprintf("Input directory %s does not contain all expected files", inp.dir))
       # next
      #}
      
      organ_results_dir <- file.path(results_dir, organ)
      dir.create(organ_results_dir, showWarnings = FALSE, recursive = TRUE)
      sample_dir <- file.path(organ_results_dir, sample)
      dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
      output_file <- file.path(sample_dir, paste0("SeuratObject_", sample, ".rds"))
      
      message(sprintf("Output file: %s", output_file))
      
      if (file.exists(output_file)) {
        message(sprintf("Sample %s already processed. Skipping.", sample))
        next
      }
      
      seur <- process_sample(inp.dir, organism, cellcyclegenes, sample_dir)
      if (!is.null(seur)) {
        saveRDS(seur, file = output_file)
        message(sprintf("Processed data saved to: %s", output_file))
      }
    }
  }
}
}
#' Process a single sample
#'
#' This function processes a single sample of single-cell RNA-seq data.
#'
#' @param input_dir Directory containing the input data
#' @param organism Organism (MOUSE, HUMAN, or RAT)
#' @param cellcyclegenes Cell cycle genes data
#' @param output_dir Directory to save output files (optional)
#' @param run_umap Whether to run UMAP (default: TRUE)
#' @param run_densmap Whether to run DensMAP (default: FALSE)
#' @param run_scina Whether to run SCINA (default: FALSE)
#' @param umap_neighbors Number of neighbors for UMAP (default: 20)
#'
#' @return A Seurat object containing the processed data
#'
#' @export
process_sample <- function(input_dir, 
                           organism, 
                           cellcyclegenes, 
                           output_dir = NULL,
                           run_umap = TRUE,
                           run_densmap = FALSE,
                           run_scina = FALSE,
                           umap_neighbors = 20) {
  tryCatch({
  
  if (length(list.files(input_dir)) == 1) {
    file_name <- list.files(input_dir)[1]
    if (grepl("\\.h5$", file_name)) {
      seur <- Read10X_h5(file.path(input_dir, file_name))
      seur <- CreateSeuratObject(counts = seur, project = "Seurat", min.cells = 3, min.features = 200)
    } else {
      seur <- readRDS(file.path(input_dir, file_name))
    }
  } else {
    seur <- Read10X(data.dir = input_dir)
    seur <- CreateSeuratObject(counts = seur, project = "Seurat", min.cells = 3, min.features = 200)
  }
    #Quality control plots
    p1 <- ggplot(seur@meta.data, aes(x=nCount_RNA, y=nFeature_RNA)) + geom_point() + geom_smooth(method="lm")
    p1 <- ggMarginal(p1, type = "histogram", fill="lightgrey")
    
    p2 <- ggplot(seur@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) + geom_point() + geom_smooth(method="lm")
    p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
    
    p <- plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
    
    if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "Count_Feature_Marginal_1stPass.png"), plot = p)
    }
      
    seur <- filterCells(seur, mad.coeff = 3, pass = 1, org = organism, output_path = output_dir)
  
  # Normalization and feature selection
  seur <- NormalizeData(seur)
  seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 3000)
  
  # Scaling
  seur <- ScaleData(seur, features = rownames(seur))
  
  # Run PCA
  seur <- RunPCA(seur, features = VariableFeatures(object = seur), npcs = 100)
  
  # Cell cycle scoring
  g2m_genes <- cellcyclegenes$symbol[cellcyclegenes$phase == "G2/M"]
  g2m_genes <- intersect(g2m_genes, rownames(seur))
  s_genes <- cellcyclegenes$symbol[cellcyclegenes$phase == "S"]
  s_genes <- intersect(s_genes, rownames(seur))
  
  seur <- CellCycleScoring(seur, 
						   s.features = s_genes, 
						   g2m.features = g2m_genes, 
						   set.ident = F)
  
  if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "PCA_CellCycleGenes.png"), DimPlot(seur, group.by = "Phase"))
    }
  
  seur$CC.Difference <- seur$S.Score - seur$G2M.Score
  
  # SCTransform and PCA
  seur <- SCTransform(seur, vars.to.regress = "CC.Difference", vst.flavor = "v2", verbose = T)
  seur <- RunPCA(seur, npcs = 100, verbose = T)
	
  if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "ElbowPlot.png"), ElbowPlot(seur, ndims = 100))
    }
	
  # Determine number of PCs
  numPCs <- 50  # Hard-coded
  
  #Doublet detection
  sweep.list <- paramSweep(seur, PCs = 1:numPCs,sct = T)
    
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    
  bcmvn <- find.pK(sweep.stats)
    
  ##Estimate expected percentage of doublets from 10X Genomics estimates from 3' Gene Expression v3.1 assay##
  estDoublets <- c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
  numCellsRec <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
  lm_doublets <- lm(estDoublets ~ numCellsRec)
  summary(lm_doublets) #Perfect linear relationship r2 = 1
    
  nExp <- round(ncol(seur) * (unname(predict(lm_doublets,data.frame(numCellsRec = ncol(seur))))/100))
  pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
  seur <- doubletFinder(seu = seur,
                                    PCs = 1:numPCs,
                                    pN = 0.25, #default
                                    pK = pK,
                                    nExp = nExp,
                                    reuse.pANN = FALSE,
                                    sct = TRUE,
                                    annotations = NULL)
    
  seur$Doublets <- "Singlet"
  seur$Doublets[seur[[paste0("pANN_0.25_",pK,"_",nExp)]] >= 0.5] <- "Doublet"
  seur$Doublets <- factor(seur$Doublets,levels = c("Doublet","Singlet"))
	
  p <- ggplot(seur@meta.data, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
      geom_point(aes(colour = Doublets), fill = "black",pch=21) +
      scale_color_manual(breaks = c("Singlet", "Doublet"),
                         values = c("black","firebrick1")) +
      geom_smooth(method="lm") +
      theme(legend.position="none") +
      annotate(geom = "text", label = paste0(as.numeric(table(seur@meta.data$Doublets)[2]), " Singlets\n",
                                             as.numeric(table(seur@meta.data$Doublets)[1]), " Doublets"), x = 4, y = 3.8)
    
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "Doublets.png"), plot = p)
  } 
  # Filter doublets
  seur <- subset(seur, subset = Doublets == "Singlet")
  
  # Run UMAP and t-SNE
  seur <- RunTSNE(seur, dims = 1:numPCs)
  # Run UMAP (optional)
  if (run_umap) {
   seur <- RunUMAP(seur, dims = 1:numPCs, n.neighbors = umap_neighbors)
   if (!is.null(output_dir)) {
     ggsave(file.path(output_dir, "UMAP.png"), DimPlot(seur, reduction = "umap"))
   }
  }
	  
  # Clustering
  seur <- findOptimalResolution(seur, pcs = 1:numPCs, output_path = output_dir)
  if (!is.null(output_dir)) {
      ggsave(file.path(output_dir, "tSNE.png"), DimPlot(seur, reduction = "tsne"))
    }
        # Run DensMAP (optional)
    if (run_densmap) {
      seur <- runDensUMAP(seur, pcs = 1:numPCs)
      if (!is.null(output_dir)) {
        ggsave(file.path(output_dir, "densMap.png"), DimPlot(seur, reduction = "densumap"))
      }
    }

  # Run SCINA (optional)
  if (run_scina) {
    tryCatch({
      seur <- annotateCells(seur, org = organism, id.type = "Symbol", tissue = "ALL")
      if (!is.null(output_dir)) {
        ggsave(file.path(output_dir, "UMAP_SCINA.png"), 
               DimPlot(seur, reduction = "umap", group.by = "SCINA_annot"), 
               width = 15)
      }
    }, error = function(e) {
      warning(paste("SCINA annotation failed:", e$message))
      if (!is.null(output_dir)) {
        write("SCINA failed", file = file.path(output_dir, "failed_SCINA.txt"))
      }
    })
  }

  return(seur)
}, error = function(e) {
    message(sprintf("Error processing sample: %s", e$message))
    return(NULL)
  })
}
