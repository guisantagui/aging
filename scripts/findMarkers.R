################################################################################
# Aging: given a Seurat object computes unsupervised clustering (with a given  #
# resolution) and obtains the markers for the clusters (all vs all).           #
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

# Terminal argument parser
################################################################################

parser <- arg_parser("scRNAseq markers finding")

parser <- add_argument(parser = parser,
                       arg = c("--seuratFile", "--resolution", "--nPCs"),
                       help = c("Seurat file to be clustered and markers find in.",
                                "Resolution for the clustering.",
                                "Number of PCs to be used in the clustering."),
                       flag = c(F, F, F))

parsed <- parse_args(parser)

# Directory stuff
################################################################################
inFile <- parsed$seuratFile
reso <- as.numeric(parsed$resolution)
nPCs <- as.numeric(parsed$nPCs)

resuAnnoDir <- sprintf("%s/annotation/", dirname(inFile))
tissue <- gsub("\\_.*", "", basename(inFile))
intMethod <- gsub(".rds",
                  "",
                  gsub(".*_",
                       "", 
                       basename(inFile)),
                  fixed = T)

outFile <- sprintf("%s%s_%s_reso_%s_markers.rds",
                   resuAnnoDir,
                   tissue,
                   intMethod,
                   as.character(reso))

mkdir_ifNot <- function(path){
        if(!dir.exists(path)){
                dir.create(path, recursive = T)
        }
}

mkdir_ifNot(resuAnnoDir)

# Do the clustering
################################################################################
integrated <- readRDS(inFile)

integrated <- RunPCA(integrated, verbose = F, npcs = 100)

integrated <- FindNeighbors(integrated,
                            reduction = "pca",
                            dims = 1:nPCs,
                            assay = "integrated")

integrated <- FindClusters(integrated,
                           resolution = reso)

# Find the markers
################################################################################

# Find markers of each cluster to the other ones (univariate test) and save 
# result. This will be used to refine SCINA annotation.
if(!file.exists(outFile)){
        clustMarkers <- FindAllMarkers(integrated)
        saveRDS(clustMarkers, file = outFile)
}#else{
 #       clustMarkers <- readRDS(markFile)
#}