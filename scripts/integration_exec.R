################################################################################
# Aging: given the folder where the preprocessed scRNAseq files are, normalizes#
# and integrates them in a single Seurat object, and saves the result in a RDS #
# file. This integration pipe is based on the RPCA vignette from Seurat, where #
# they not use oneliner(https://satijalab.org/seurat/articles/integration_rpca)#
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
                       arg = c("--preProcDir", "--method"),
                       help = c("Directory where folders containing the RDS file for each sample are stored.", "Integration method, either CCA or RPCA."),
                       flag = c(F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################

resuPrepDir <- parsed$preProcDir
intMethod <- parsed$method
tissue <- basename(resuPrepDir)

plotIntgDir <- resuPrepDir %>%
        dirname() %>%
        dirname() %>%
        dirname() %>%
        paste0("/plots/integration/")

resuIntgDir <- resuPrepDir %>%
        dirname() %>%
        dirname() %>%
        paste0("/integration/")

outFile <- sprintf("%s%s_int_%s.rds",
                   resuIntgDir,
                   basename(resuPrepDir),
                   intMethod)

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

# Selects features to use for integrating multiple datasets. In this case ctrl
# and stim. We give the number of features. It ranks the features by the number
# of datasets where these features are deemed variable in.
features <- SelectIntegrationFeatures(object.list = seurList,
                                      nfeatures = 3000)

# This function does some steps for preparing the two objects in the list to be
# integrated

seurList <- PrepSCTIntegration(object.list = seurList,
                               anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = seurList,
                                         normalization.method = "SCT",
                                         anchor.features = features,
                                         reduction = intMethod)

#if(tissue == "brain"){
#        kWeight <- 32
#}else{
        kWeight <- 100
#}

integrated <- IntegrateData(anchorset = immune.anchors,
                            normalization.method = "SCT",
                            k.weight = kWeight)

# Save integrated RDS file
saveRDS(integrated, file = outFile)

print(sprintf("%s saved at %s.",
              basename(outFile),
              dirname(outFile)))