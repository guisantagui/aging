################################################################################
# Aging: given the folder where the preprocessed scRNAseq files are, normalizes#
# and integrates them in a single Seurat object, and saves the result in a RDS #
# file.                                                                        #
################################################################################

if(!require("argparser")){
        install.packages("argparser",
                         repos = 'https://pbil.univ-lyon1.fr/CRAN/')
}

if(!require("Seurat", quietly = T)) remotes::install_github("satijalab/seurat",
                                                            "seurat5",
                                                            quiet = TRUE,
                                                            upgrade = "never")

if(!require("SeuratData", quietly = T)){
        remotes::install_github("satijalab/seurat-data",
                                "seurat5",
                                quiet = TRUE)
}

if(!require("patchwork", quietly = T)){
        install.packages("patchwork")
}

if(!require("dplyr", quietly = T)){
        install.packages("dplyr")
}

if(!require("ggplot2", quietly = T)){
        install.packages("ggplot2")
}


library(Seurat)
library(patchwork)
library(SeuratData)
library(dplyr)
library(ggplot2)

# Terminal argument parser
################################################################################

parser <- arg_parser("scRNAseq preprocessing")

parser <- add_argument(parser = parser,
                       arg = "--preProcDir",
                       help = "Directory where folders containing the RDS file for each sample are stored.",
                       flag = F)

parsed <- parse_args(parser)

# Directory stuff
################################################################################

resuPrepDir <- parsed$preProcDir


plotIntgDir <- resuPrepDir %>%
        dirname() %>%
        dirname() %>%
        dirname() %>%
        paste0("/plots/integration/")

resuIntgDir <- resuPrepDir %>%
        dirname() %>%
        dirname() %>%
        paste0("/integration/")

outFile <- sprintf("%s%s_int.rds", resuIntgDir, basename(resuPrepDir))

mkdir_ifNot <- function(path){
        if(!dir.exists(path)){
                dir.create(path, recursive = T)
        }
}

mkdir_ifNot(resuIntgDir)
mkdir_ifNot(plotIntgDir)

# Functions
################################################################################

# Obtains paths to all the preprocessed files given the directory where they're 
# stored
getPreprocFiles <- function(prepDir){
        prepFiles <- list.files(prepDir, full.names = T)# %>%
                #list.files(full.names = T)
        RDataFiles <- prepFiles[grep("RData", prepFiles)]
        rdsFiles <- prepFiles[grep("rds", prepFiles)]
        out <- list(RData = RDataFiles,
                    rds = rdsFiles)
        return(out)
}

# Integration
################################################################################

# Obtain paths to preprocessed files
prepFiles <- getPreprocFiles(resuPrepDir)

# Create a list with each sample's scTransformed seurat object
seurList <- list()
for(i in seq_along(prepFiles$rds)){
        f <- prepFiles$rds[i]
        samp <- gsub(".rds", "", basename(f))
        seur <- readRDS(f)
        seur <- SCTransform(seur, vst.flavor = "v2", verbose = F) %>%
                RunPCA(npcs = 100, verbose = F)
        seurList[[samp]] <- seur
}

# Merge all the objectds in the list in a single seurat object
seur <- merge(seurList[[1]], seurList[2:length(seurList)])

seur <- IntegrateLayers(object = seur,
                        method = CCAIntegration,
                        orig.reduction = "PCA",
                        new.reduction = "integrated.cca")



# Perform an integrated analysis
################################################################################

integrated <- RunPCA(integrated, verbose = F, npcs = 100)
integrated <- RunUMAP(integrated,
                      reduction = "pca",
                      dims = 1:100,
                      verbose = F)

# Find k.param nearest neighbours using as input the PCA dimensionality
# reduction
integrated <- FindNeighbors(integrated,
                            reduction = "pca",
                            dims = 1:100)

# Identify clusters of cells based on the shared nearest neighbour obtained in
# the previous step. The resolution parameter gives the size of communities 
# obtained. The bigger the resolution the smaller the communities, so more 
# groups will be obtained
integrated <- FindClusters(integrated,
                           resolution = .3)

# The split.by argument allows us to obtain separate plots for each stim 
# factor (either control or stim). DimPlot graphs the output of a dimensional
# reduction technique on a 2D scatterplot. Each point is a cell and it's 
# positioned based on the cell embeddings determined by the reduction technique.
# By default, cells are colored by their identity class. This can be changed
# with group.by argument.

sampVec <- colnames(integrated@assays$SCT$data)
names(sampVec) <- names(integrated$SCINA_annot)
integrated$sample <- sampVec

intUMAP_plt <- DimPlot(integrated,
                       reduction = "umap",
                       split.by = "sample")

#intUMAP_scina_plt <- DimPlot(integrated,
#                             reduction = "umap",
#                             split.by = "sample",
#                             group.by = "SCINA_annot")

ggsave(paste0(plotIntgDir, "intUMAP.pdf"),
       intUMAP_plt, width = 50, height = 10, limitsize = F)

ggsave(paste0(plotIntgDir, "intUMAP_scina.pdf"),
       intUMAP_scina_plt, width = 50, height = 50, limitsize = F)

# Save integrated RDS file
saveRDS(integrated, file = outFile)

print(sprintf("%s saved at %s.",
              basename(outFile),
              dirname(outFile)))