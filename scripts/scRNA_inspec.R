if(!require("remotes", quietly = T)) install.packages("remotes")
if(!require("BiocManager", quietly = T)) install.packages("BiocManager")
if(!require("Seurat", quietly = T)) remotes::install_github("satijalab/seurat",
                                                            "seurat5",
                                                            quiet = TRUE,
                                                            upgrade = "never")
library(Seurat)

if(!require("SCINA", quietly = T)){
        install.packages("SCINA",
                         repos='http://cran.us.r-project.org')
}
library(SCINA)

if(!require("plyr", quietly = T)){
        install.packages("plyr",
                         repos='http://cran.us.r-project.org')
}
library(SCINA)
if(!require("ggplot2", quietly = T)){
        install.packages("ggplot2",
                         repos='http://cran.us.r-project.org')
}
library(ggplot2)
if(!require("ggExtra", quietly = T)){
        install.packages("ggExtra",
                         repos='http://cran.us.r-project.org')
}
library(ggExtra)
if(!require("cowplot", quietly = T)){
        install.packages("cowplot",
                         repos='http://cran.us.r-project.org')
}
library(cowplot)
if(!require("biomaRt", quietly = T)){
        install.packages("biomaRt",
                         repos='http://cran.us.r-project.org')
}
library(biomaRt)
if(!require(dplyr)) install.packages("dplyr",
                                     repos = 'http://cran.us.r-project.org')
if(!require("DoubletFinder", quietly = T)){
        remotes::install_github("chris-mcginnis-ucsf/DoubletFinder",
                                quiet = T,
                                upgrade = "never")
}
library(DoubletFinder)
if(!require("monocle3", quietly = T)){
        remotes::install_github("cole-trapnell-lab/monocle3",
                                quiet = T,
                                upgrade = "never")
}
library(monocle3)
if(!require("SeuratWrappers", quietly = T)){
        devtools::install_github('satijalab/seurat-wrappers')
}
library(SeuratWrappers)

#library(SeuratObject)

# Directory stuff
################################################################################

rootDir <- "/Users/guillem.santamaria/Documents/postdoc/comput/aging/"
plotPrepDir <- paste0(rootDir, "plots/preprocessing/")
resuPrepDir <- paste0(rootDir, "results/preprocessing/")
dataDir <- paste0(rootDir, "data/")
agAnDir <- paste0(dataDir, "ageAnno/")
markerFile <- paste0(agAnDir, "AgeAnno-main/scRNA/scRNAmarker.txt")
metDatDir <- paste0(dataDir, "gkac847_supplemental_files/")
utilDir <- paste0(dataDir, "utility_files/")

metDatFile <- paste0(metDatDir, "Supplementary Tables.xlsx")
cellCycGenesFile <- paste0(utilDir, "CellCycleGenes_Human.csv")

if(!dir.exists(plotPrepDir)){
        dir.create(plotPrepDir, recursive = T)
}
if(!dir.exists(resuPrepDir)){
        dir.create(resuPrepDir, recursive = T)
}

# User parameters
################################################################################
id.type = "Symbol"

# Functions
################################################################################

# For parsing metadata. There are some ages that are ranges and others that are
# higher than something, expressed with a plus. Let's convert it to numeric. For
# ranges calculate average. For the plus just remove it.
ageToNum <- function(age_years){
        convAge <- function(x){
                if(grepl("+", x, fixed = T)){
                        x <- as.numeric(gsub("+", "", x, fixed = T))
                }else if(grepl("-", x, fixed = T)){
                        splt <- as.numeric(strsplit(x, split = "-")[[1]])
                        x <- mean(splt)
                }else{
                        x <- as.numeric(x)
                }
                return(x)
        }
        age_num <- c()
        for(a in age_years){
                numA <- convAge(a)
                age_num <- c(age_num, numA)
        }
        return(age_num)
}

# For parsing metadata. Fills the tissue column, which is partially empty
fillOrgVec <- function(organs){
        notNAIdxs <- which(!is.na(organs))
        notNAIdxs <- c(notNAIdxs, (length(organs) + 1))
        idxVec <- 1:length(organs)
        for(i in 1:(length(notNAIdxs) - 1)){
                idx1 <- notNAIdxs[i]
                idx2 <- notNAIdxs[i + 1]
                organ <- organs[idx1]
                logVec <- is.na(organs) & idxVec < idx2 & idxVec > idx1
                organs[logVec] <- organ
        }
        return(organs)
}


# seur: A Seurat object
# mad.coeff: median absolute deviation (MAD) coefficient. MAD is the absolute
# deviation of each value from the median, so the coefficient determines 
# the maximums and minimums we are allowing.
# pass: An identifier of the number of times this function is running. Only for plotting.
# org: The organism the data is from. Only relevant if not all metadata info is there.
filterCells <- function(seur,mad.coeff = 3,pass = -1,org = "HUMAN",plot.path = "./"){
        seur@meta.data$Barcodes <- rownames(seur@meta.data)
        #Calculate percent.mito and percent.ribo metadata columns if they are not there
        if(!any(colnames(seur@meta.data) == "percent.mito")){
                if(org == "HUMAN"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^MT-")
                } else if(org == "MOUSE"){
                        seur[["percent.mito"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "^mt-")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        if(!any(colnames(seur@meta.data) == "percent.ribo")){
                if(org == "HUMAN"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^RPL|^RPS|^MRP)")
                } else if(org == "MOUSE"){
                        seur[["percent.ribo"]] <- PercentageFeatureSet(seur,
                                                                       pattern = "(^Rpl|^Rps|^Mrp)")
                } else {
                        stop("The specified organism is not supported")
                }
        }
        # Filtering cells based on percentage of mitochondrial transcripts
        cell.QC.stat <- seur@meta.data
        
        max.mito.thr <- median(cell.QC.stat$percent.mito) + mad.coeff*mad(cell.QC.stat$percent.mito)
        min.mito.thr <- median(cell.QC.stat$percent.mito) - mad.coeff*mad(cell.QC.stat$percent.mito)
        
        p1 <- ggplot(cell.QC.stat, aes(x=nFeature_RNA, y=percent.mito)) +
                geom_point() +
                geom_hline(aes(yintercept = max.mito.thr), colour = "red", linetype = 2) +
                geom_hline(aes(yintercept = min.mito.thr), colour = "red", linetype = 2) +
                annotate(geom = "text",
                         label = paste0(as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[2]),
                                        " cells removed\n",
                                        as.numeric(table(cell.QC.stat$percent.mito > max.mito.thr | cell.QC.stat$percent.mito < min.mito.thr)[1]),
                                        " cells remain"),
                         x = 6000,
                         y = -10)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey", bins=100) 
        ggsave(paste0(plot.path,"Mitofilter_Marginal_Pass",pass,".png"),plot = p)
        
        cell.QC.stat <- cell.QC.stat %>%
                dplyr::filter(percent.mito <= max.mito.thr) %>% 
                dplyr::filter(percent.mito >= min.mito.thr)
        
        # Filtering cells based on number of genes and transcripts detected
        # Set low and hight thresholds on the number of detected genes
        min.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) - mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        max.features.thr <- median(log10(cell.QC.stat$nFeature_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nFeature_RNA))
        
        # Set hight threshold on the number of transcripts
        max.count.thr <- median(log10(cell.QC.stat$nCount_RNA)) + mad.coeff*mad(log10(cell.QC.stat$nCount_RNA))
        
        p1 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2)
        
        p <- ggMarginal(p1, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountFilter_1_Pass",pass,".png"),plot = p)
        
        # Filter cells base on both metrics
        cell.QC.stat <- cell.QC.stat %>% 
                dplyr::filter(log10(nFeature_RNA) > min.features.thr) %>%
                dplyr::filter(log10(nFeature_RNA) < max.features.thr) %>%
                dplyr::filter(log10(nCount_RNA) < max.count.thr)
        
        lm.model <- lm(data = cell.QC.stat, formula = log10(nFeature_RNA) ~ log10(nCount_RNA))
        
        p2 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point() +
                geom_smooth(method="lm") +
                geom_hline(aes(yintercept = min.features.thr), colour = "green", linetype = 2) +
                geom_hline(aes(yintercept = max.features.thr), colour = "green", linetype = 2) +
                geom_vline(aes(xintercept = max.count.thr), colour = "red", linetype = 2) +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") +
                annotate(geom = "text", label = paste0(dim(cell.QC.stat)[1], " QC passed cells"), x = 4, y = 3.8)
        
        p <- ggMarginal(p2, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_2_Pass",pass,".png"),plot = p)
        
        # Cells to exclude lie below an intercept offset of -0.09
        cell.QC.stat$validCells <- log10(cell.QC.stat$nFeature_RNA) > (log10(cell.QC.stat$nCount_RNA) * lm.model$coefficients[2] + (lm.model$coefficients[1] - 0.09))
        
        p3 <- ggplot(cell.QC.stat, aes(x=log10(nCount_RNA), y=log10(nFeature_RNA))) +
                geom_point(aes(colour = validCells), fill = "black",pch=21) +
                scale_color_manual(breaks = c("TRUE", "FALSE"), 
                                   values = c("black","firebrick1")) +
                geom_smooth(method="lm") +
                geom_abline(intercept = lm.model$coefficients[1] - 0.09 , slope = lm.model$coefficients[2], color="orange") + 
                theme(legend.position="none") +
                annotate(geom = "text", label = paste0(as.numeric(table(cell.QC.stat$validCells)[2]), " QC passed cells\n",
                                                       as.numeric(table(cell.QC.stat$validCells)[1]), " QC filtered"), x = 4, y = 3.8)
        
        p <- ggMarginal(p3, type = "histogram", fill="lightgrey")
        ggsave(paste0(plot.path,"FeatureAndCountOutlier_3_Pass",pass,".png"),plot = p)
        
        # Remove invalid cells
        cell.QC.stat <- cell.QC.stat %>% dplyr::filter(validCells)
        
        seur <- subset(seur, subset = Barcodes %in% cell.QC.stat$Barcodes)
        return(seur)
}

# Find the optimal number of PCs

findNumPCs <- function(seur){
        pct <- seur[["pca"]]@stdev / sum(seur[["pca"]]@stdev) * 100
        cumu <- cumsum(pct)
        
        co1 <- which(cumu > 90 & pct < 5)[1]
        
        co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05),
                    decreasing = T)[1] + 1
        
        return(min(co1,co2))
}

# Annotate cells using SCINA
annotateCells <- function(seur,
                          plot.path = "./",
                          id.type = "Symbol",
                          tissue = "Brain",
                          useVarFeats = F,
                          markDF = NULL){
        require(preprocessCore)
        #inp_df <- as.matrix(seur@assays$RNA[,])
        if(!useVarFeats){
                inp_df <- as.matrix(seur@assays$RNA@data)
                inp_df <- log(inp_df+1)
                inp_df[] = normalize.quantiles(as.matrix(inp_df))
        }else{
                inp_df <- as.matrix(seur@assays$RNA@scale.data)
        }
        
        if(!is.null(markDF)){
                cm_human <- markDF
                cm_human <- cm_human[grepl(tissue, tolower(cm_human$tissue)),]
                cm_human <- cm_human[,c("cell_name","Symbol")]
                cm_human <- cm_human[cm_human$Symbol != "",]
        }else{
                cm_human <- read.table(sprintf("%sCell_marker_Human.txt",
                                               utilDir),
                                       header = T,
                                       sep = "\t",
                                       stringsAsFactors = F,
                                       quote = "",
                                       comment.char = "",
                                       encoding = "UTF-8")
                
                cm_human$tissue_class <- gsub("Ê",
                                              " ",
                                              iconv(cm_human$tissue_class,
                                                    from = "ISO-8859-1",
                                                    to = "UTF-8"))
                
                cm_human$tissue_class <- gsub("Bone marrow",
                                              "bone-marrow",
                                              cm_human$tissue_class)
                cm_human$tissue_class <- gsub("Skeletal muscle",
                                              "skeletal-muscle",
                                              cm_human$tissue_class)
                cm_human$tissue_class[cm_human$tissue_type == "Retina"] <- "Retina"
                cm_human$tissue_class <- tolower(cm_human$tissue_class)
                
                cm_human <- cm_human[which(grepl(tissue,cm_human$tissue_class)),]
                cm_human <- cm_human[cm_human$cancer_type == "Normal",]
                cm_human <- cm_human[cm_human$cell_type != "Cancer cell",]
                cm_human <- cm_human[,c("cell_name","Symbol")]
                #cm_human <- ddply(cm_human,.(cell_name,Symbol),nrow)
                cm_human <- cm_human[cm_human$Symbol != "",]
                #cm_human <- cm_human[cm_human$V1 > 1,]
                #cm_human <- cm_human[,c("cell_name","Symbol")]
        }
        
        if(id.type == "Ensembl"){
                require(biomaRt)
                mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
                symbToEns <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=cm_human$Symbol,mart= mart)
                cm_human <- merge(cm_human,symbToEns,by.x = 2, by.y = 2, all.x = F)
                cm_human <- cm_human[,-1]
                colnames(cm_human) <- c("cell_name","Symbol")
        }
        print(cm_human)
        sigs <- lapply(unique(cm_human$cell_name),function(x){cm_human$Symbol[cm_human$cell_name == x]})
        names(sigs) <- unique(cm_human$cell_name)
        sigs <- lapply(sigs,function(x){intersect(x,rownames(inp_df))})
        #print(sigs)
        if(!is.null(markDF)){
                sigs <- sigs[sapply(sigs,length) >= 1]
        }else{
                sigs <- sigs[sapply(sigs,length) > 3]
        }
        print(sigs)
        
        scina.results <- SCINA(inp_df,
                               sigs, 
                               max_iter = 100, 
                               convergence_n = 10, 
                               convergence_rate = 0.99, 
                               sensitivity_cutoff = 0.9, 
                               rm_overlap=F, 
                               allow_unknown=T
        )
        
        seur$SCINA_annot <- scina.results$cell_labels
        ggsave(paste0(plot.path,"UMAP_SCINA.png"),
               DimPlot(seur, reduction = "umap", group.by = "SCINA_annot"),
               width = 20, height = 10)
        return(seur)
}



# Load and parse data
################################################################################

# Load and parse metadata file
metDat <- data.frame(readxl::read_xlsx(metDatFile, sheet = 1))

colnames(metDat) <- unlist(metDat[2, ])
colnames(metDat) <- gsub(" ", "_", colnames(metDat))
colnames(metDat) <- gsub("(", "", colnames(metDat), fixed = T)
colnames(metDat) <- gsub(")", "s", colnames(metDat), fixed = T)

metDat <- metDat[3:(which(metDat$organs == "scATAC") - 2), ]

metDat$age_years_num <- ageToNum(metDat$age_years)
metDat$organs <- fillOrgVec(metDat$organs)

# marker file from ageAnno paper
markers <- read.table(markerFile, sep = "\t", header = T)
colnames(markers) <- markers[1, ]
colnames(markers) <- tolower(colnames(markers))
markers <- markers[2:nrow(markers), ]
colnames(markers)[2:3] <- c("cell_name","Symbol")

#markers$cell_name <- gsub("Epithelial cells",
#                          "Epithelial cell",
#                          markers$cell_name)
#markers$cell_name <- gsub("  ", " ", markers$cell_name)
#markers$cell_name <- gsub("Endothelial cells",
#                          "Endothelial cell",
#                          markers$cell_name)

# cellMarker (10.1093/nar/gky900) marker file.
cellMrkr <- read.table(paste0(utilDir, "Human_cell_markers.txt"), sep = "\t", header = T)
cellMrkr <- cellMrkr[cellMrkr$cellType == "Normal cell", ]
colnames(cellMrkr) <- gsub("tissueType", "tissue", colnames(cellMrkr))
cellMrkr_parsed <- data.frame(matrix(nrow = 0,
                                     ncol = 3,
                                     dimnames = list(NULL,
                                                     c("tissue",
                                                       "cell_name",
                                                       "Symbol"))))
for(i in 1:nrow(cellMrkr)){
        mrkrs <- strsplit(cellMrkr$geneSymbol[i], split = ", ")[[1]]
        if(sum(is.na(mrkrs)) > 0){
                print(i)
        }
        gsub("[|]", "", mrkrs)
        tissVec <- rep(cellMrkr$tissue[i], length(mrkrs))
        cell <- rep(cellMrkr$cellName[i], length(mrkrs))
        toBind <- data.frame(tissue = tissVec,
                             cell_name = cell,
                             Symbol = mrkrs)
        cellMrkr_parsed <- rbind.data.frame(cellMrkr_parsed, toBind)
}

cellMrkr_parsed$Symbol <- gsub("[", "", cellMrkr_parsed$Symbol, fixed = T)
cellMrkr_parsed$Symbol <- gsub("]", "", cellMrkr_parsed$Symbol, fixed = T)
cellMrkr_parsed <- cellMrkr_parsed[!is.na(cellMrkr_parsed$Symbol), ]

cellMrkr_parsed$tissue <- gsub("Bone marrow", "bone-marrow", cellMrkr_parsed$tissue)
cellMrkr_parsed$tissue <- gsub("Skeletal muscle", "skeletal-muscle", cellMrkr_parsed$tissue)


# Load tissue scRNAseq data in a list
tissFiles <- list.files(agAnDir, full.names = T)
tissFiles <- tissFiles[grep(".rds", tissFiles)]
organs <- gsub(".rds", "", basename(tissFiles))


# Get H. sapiens ensemble dataset
mart <- useDataset("hsapiens_gene_ensembl",useMart(biomart = "ensembl"))

alrDoneOrgs <- list.files(resuPrepDir)[grep(".rds",
                                            list.files(resuPrepDir))]

alrDoneOrgs <- unique(gsub("_preproc.rds|_preproc_ageAnnot.rds",
                           "",
                           alrDoneOrgs))

startedOrgs <- list.files(resuPrepDir)[!grepl(".rds", list.files(resuPrepDir))]
# Iterate over tissues and samples to preprocess each one individually, write
# RDS files of each sample (we will use these files for integration later on)
# and merge samples of same tissue in a single file (we will use them later on
# for statistical analysis, after we perform annotation with the integrated
# files).
orgList <- list()
for(i in seq_along(tissFiles)){
        #i <- 2
        f <- tissFiles[i]
        org <- organs[i]
        
        plotTissDir <- sprintf("%s%s/", plotPrepDir, org)
        resuTissDir <- sprintf("%s%s/", resuPrepDir, org)
        if(!dir.exists(resuTissDir)){
                dir.create(resuTissDir)
        }
        
        tissList <- list()
        if(!org %in% alrDoneOrgs){
                tiss <- readRDS(f)
                count <- tiss@assays$RNA$counts
                
                samps <- unique(gsub(".*_", "", colnames(count)))
                for(s in samps){
                        #s <- samps[1]
                        sampTissFiles <- list.files(resuTissDir, full.names = T)
                        alrDoneSamps <- gsub("_preproc.rds",
                                             "",
                                             basename(sampTissFiles))
                        if(!s %in% alrDoneSamps){
                                #s <- "mid4"
                                print(sprintf("Preprocessing %s sample %s...",
                                              org, s))
                                
                                contSamp <- count[, grep(sprintf("_%s", s),
                                                         colnames(count))]
                                
                                minCells <- 3
                                minFeats <- 700
                                nfeatures <- Matrix::colSums(x = contSamp > 0)
                                nCellsPassed <- length(which(x = nfeatures >= minFeats))
                                #Matrix::rowSums(contSamp > 0)
                                #sum(Matrix::rowSums(contSamp) >= minFeats)
                                #sum(Matrix::rowSums(contSamp) >= minCells)
                                #nFeatsFilt <- sum(rowSums(contSamp) >= minFeats)
                                #nCellsFilt <- sum(colSums(contSamp) >= minCells)
                                if(nCellsPassed > 1){
                                        tissSamp <- CreateSeuratObject(contSamp,
                                                                       min.cells = minCells,
                                                                       min.features = minFeats)
                                        plot.path <- sprintf("%s%s/", plotTissDir, s)
                                        if(!dir.exists(plot.path)){
                                                dir.create(plot.path, recursive = T)
                                        }
                                        numCellsRecovered <- ncol(tissSamp)
                                        if(id.type == "Ensembl"){
                                                require(biomaRt)
                                                mart <- useDataset("hsapiens_gene_ensembl",
                                                                   useMart("ensembl"))
                                                ensToSymb <- getBM(filters= "ensembl_gene_id",
                                                                   attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                                   values=rownames(tissSamp),
                                                                   mart= mart)
                                                mt.features <- ensToSymb$ensembl_gene_id[which(grepl("^MT-",
                                                                                                     ensToSymb$hgnc_symbol))]
                                                ribo.features <- ensToSymb$ensembl_gene_id[which(grepl("(^RPL|^RPS|^MRP)",
                                                                                                       ensToSymb$hgnc_symbol))]
                                                
                                                # Percent of mitochondrial counts
                                                tissSamp[["percent.mito"]] <- PercentageFeatureSet(tissSamp,
                                                                                                   features = mt.features)
                                                
                                                # Percent of mitochondrial ribosomal
                                                tissSamp[["percent.ribo"]] <- PercentageFeatureSet(tissSamp,
                                                                                                   features = ribo.features)
                                                
                                        }else if(id.type == "Symbol"){
                                                # Percent of mitochondrial counts
                                                tissSamp[["percent.mito"]] <- PercentageFeatureSet(tissSamp,
                                                                                                   pattern = "^MT-")
                                                
                                                # Percent of mitochondrial ribosomal
                                                tissSamp[["percent.ribo"]] <- PercentageFeatureSet(tissSamp,
                                                                                                   pattern = "(^RPL|^RPS|^MRP)")
                                        }
                                        
                                        p1 <- ggplot(tissSamp@meta.data,
                                                     aes(x=nCount_RNA,
                                                         y=nFeature_RNA)) +
                                                geom_point() +
                                                geom_smooth(method="lm")
                                        p1 <- ggMarginal(p1,
                                                         type = "histogram",
                                                         fill="lightgrey")
                                        
                                        p2 <- ggplot(tissSamp@meta.data,
                                                     aes(x=log10(nCount_RNA),
                                                         y=log10(nFeature_RNA))) +
                                                geom_point() +
                                                geom_smooth(method="lm")
                                        p2 <- ggMarginal(p2, type = "histogram", fill="lightgrey")
                                        
                                        p <- plot_grid(plotlist = list(p1,p2), ncol=2, align='h', rel_widths = c(1, 1))
                                        
                                        ggsave(paste0(plot.path,
                                                      "Count_Feature_Marginal_1stPass.png"),
                                               plot = p)
                                        
                                        # Filter cells based on percent.mito, percent.ribo and count/feature relationships.
                                        # Filter cells based on percent.mito, percent.ribo and count/feature relationships.
                                        tissSampFilt <- filterCells(tissSamp,
                                                                    mad.coeff = 5,
                                                                    pass = 1,
                                                                    org = "HUMAN",
                                                                    plot.path = plot.path)
                                        
                                        # Log-normalize data
                                        tissSampFilt <- NormalizeData(tissSampFilt)
                                        
                                        # Find variable features --> Identify features that are outliers on a mean
                                        # variability plot
                                        tissSampFilt <- FindVariableFeatures(tissSampFilt,
                                                                             selection.method = "vst",
                                                                             nfeatures = 3000)
                                        
                                        # Scale data --> Scales and centers features in the dataset
                                        tissSampFilt <- ScaleData(tissSampFilt,
                                                                  features = VariableFeatures(object = tissSampFilt))
                                        
                                        # Run PCA
                                        if(ncol(tissSampFilt) < 100){
                                                nPCs <- ncol(tissSampFilt) - 1
                                        }else{
                                                nPCs <- 100
                                        }
                                        
                                        tissSampFilt <- RunPCA(tissSampFilt,
                                                               features = VariableFeatures(object = tissSampFilt),
                                                               npcs = nPCs)
                                        
                                        # Assign cell cycle phase to each cell by looking at the expression of
                                        # genes associated with S and G2M phases.
                                        
                                        cellcyclegenes <- read.csv(cellCycGenesFile)
                                        
                                        ensToSymb <- getBM(filters= "ensembl_gene_id",
                                                           attributes= c("ensembl_gene_id","hgnc_symbol"),
                                                           values=cellcyclegenes$geneID,mart= mart)
                                        cellcyclegenes <- merge(cellcyclegenes,ensToSymb,by.x = 2, by.y = 1, all.x = T)
                                        g2m_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "G2/M"]
                                        g2m_genes <- intersect(g2m_genes,rownames(tissSampFilt))
                                        s_genes <- cellcyclegenes$hgnc_symbol[cellcyclegenes$phase == "S"]
                                        s_genes <- intersect(s_genes,rownames(tissSampFilt))
                                        
                                        if(length(s_genes) > 1 & length(g2m_genes) > 1){
                                                tissSampFilt <- CellCycleScoring(tissSampFilt,
                                                                                 s.features = s_genes,
                                                                                 g2m.features = g2m_genes,
                                                                                 set.ident = F)
                                                
                                                ggsave(paste0(plot.path,
                                                              "PCA_CellCycleGenes.png"),
                                                       DimPlot(tissSampFilt,
                                                               group.by = "Phase"))
                                                
                                                # Regress out the expression of the genes in the cell cycle
                                                tissSampFilt$CC.Difference <- tissSampFilt$S.Score - tissSampFilt$G2M.Score
                                                
                                                tissSampFilt <- SCTransform(tissSampFilt,
                                                                            vars.to.regress = "CC.Difference",
                                                                            vst.flavor = "v2",
                                                                            verbose = T)
                                                
                                                tissSampFilt <- RunPCA(tissSampFilt,
                                                                       assay = "SCT",
                                                                       npcs = nPCs,
                                                                       verbose = T)
                                                
                                                numPCs <- findNumPCs(tissSampFilt) #+
                                                #round(findNumPCs(tissSampFilt)/2)
                                                
                                                elbPlot <- ElbowPlot(tissSampFilt,
                                                                     ndims = nPCs) +
                                                        geom_vline(aes(xintercept = numPCs),
                                                                   colour = "red",
                                                                   linetype = "dashed")
                                                
                                                ggsave(paste0(plot.path,
                                                              "ElbowPlot.png"),
                                                       elbPlot)
                                                
                                                # paramSweep_v3 finds pN and pK values for removing the doublets
                                                # later on
                                                sweep.list <- paramSweep_v3(tissSampFilt,
                                                                            PCs = 1:numPCs,
                                                                            sct = T)
                                                sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
                                                bcmvn <- find.pK(sweep.stats)
                                                
                                                ## Estimate expected percentage of doublets from 10X Genomics 
                                                # estimates from 3' Gene Expression v3.1 assay##
                                                estDoublets <- c(0.4,0.8,1.6,2.4,3.2,4,4.8,5.6,6.4,7.2,8)
                                                numCellsRec <- c(500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
                                                
                                                scatter.smooth(numCellsRec,estDoublets) #Looks linear
                                                lm_doublets <- lm(estDoublets ~ numCellsRec)
                                                summary(lm_doublets) #Perfect linear relationship r2 = 1
                                                
                                                nExp <- round(ncol(tissSampFilt) * (unname(predict(lm_doublets,
                                                                                                   data.frame(numCellsRec = ncol(tissSampFilt))))/100))
                                                
                                                pK <- as.numeric(levels(bcmvn$pK)[bcmvn$BCmetric == max(bcmvn$BCmetric)])
                                                tissSampFilt <- doubletFinder_v3(seu = tissSampFilt,
                                                                                 PCs = 1:numPCs,
                                                                                 pN = 0.25, #default
                                                                                 pK = pK,
                                                                                 nExp = nExp,
                                                                                 reuse.pANN = FALSE,
                                                                                 sct = TRUE,
                                                                                 annotations = NULL)
                                                
                                                tissSampFilt$Doublets <- "Singlet"
                                                
                                                tissSampFilt$Doublets[tissSampFilt[[paste0("pANN_0.25_",
                                                                                           pK,
                                                                                           "_",
                                                                                           nExp)]] >= 0.5] <- "Doublet"
                                                
                                                tissSampFilt$Doublets <- factor(tissSampFilt$Doublets,
                                                                                levels = c("Doublet",
                                                                                           "Singlet"))
                                                
                                                p <- ggplot(tissSampFilt@meta.data,
                                                            aes(x=log10(nCount_RNA),
                                                                y=log10(nFeature_RNA))) +
                                                        geom_point(aes(colour = Doublets), fill = "black",pch=21) +
                                                        scale_color_manual(breaks = c("Singlet", "Doublet"), 
                                                                           values = c("black","firebrick1")) +
                                                        geom_smooth(method="lm") +
                                                        theme(legend.position="none") +
                                                        annotate(geom = "text",
                                                                 label = paste0(as.numeric(table(tissSampFilt@meta.data$Doublets)[2]),
                                                                                " Singlets\n",
                                                                                as.numeric(table(tissSampFilt@meta.data$Doublets)[1]),
                                                                                " Doublets"),
                                                                 x = 4,
                                                                 y = 3.8)
                                                
                                                ggsave(paste0(plot.path,"Doublets.png"),plot = p)
                                                
                                                tissSampFilt_nodups <- subset(tissSampFilt,
                                                                              subset = Doublets == "Singlet")
                                                
                                                # Run UMAP and tSNE
                                                tissSampFilt_nodups <- RunUMAP(tissSampFilt_nodups,
                                                                               dims = 1:numPCs,
                                                                               n.neighbors = 20)
                                                
                                                if(ncol(tissSampFilt_nodups) < 90){
                                                        perp <- floor(ncol(tissSampFilt_nodups)/3)
                                                }else{
                                                        perp <- 30
                                                }
                                                
                                                
                                                tissSampFilt_nodups <- RunTSNE(tissSampFilt_nodups,
                                                                               dims = 1:numPCs,
                                                                               perplexity = perp)
                                                
                                                cds <- as.cell_data_set(tissSampFilt_nodups)
                                                cds <- cluster_cells(cds,
                                                                     k = 10,
                                                                     reduction_method = "UMAP",
                                                                     cluster_method = "leiden",
                                                                     #resolution = 0.0001,
                                                                     num_iter = 5
                                                )
                                                
                                                tissSampFilt_nodups$seurat_clusters <- clusters(cds)
                                                Idents(tissSampFilt_nodups) <- clusters(cds)
                                                
                                                #Let's see how it looks
                                                ggsave(paste0(plot.path,"UMAP.png"),
                                                       DimPlot(tissSampFilt_nodups,
                                                               reduction = "umap"))
                                                ggsave(paste0(plot.path,"tSNE.png"),
                                                       DimPlot(tissSampFilt_nodups,
                                                               reduction = "tsne"))
                                                
                                                #Annotate
                                                #tissSampFilt_nodups <- annotateCells(tissSampFilt_nodups,
                                                #                                     plot.path = plot.path,
                                                #                                     id.type = "Symbol",
                                                #                                     tissue = org,
                                                #                                     useVarFeats = T,
                                                #                                     markDF = markers
                                                #                                     #, markDF = cellMrkr_parsed
                                                #)
                                                
                                                #tissSampFilt_nodups_ageAnnot <- annotateCells(tissSampFilt_nodups,
                                                #                                              plot.path = plot.path,
                                                #                                              id.type = "Symbol",
                                                #                                              tissue = org,
                                                #                                              useVarFeats = T,
                                                #                                              markDF = markers
                                                #                                              #, markDF = cellMrkr_parsed
                                                #                                              )
                                                
                                                #tissSampFilt_nodups_cmAnnot <- annotateCells(tissSampFilt_nodups,
                                                #                                              plot.path = plot.path,
                                                #                                             id.type = "Symbol",
                                                #                                             tissue = org,
                                                #                                              useVarFeats = T#,
                                                #                                              #markDF = markers
                                                #                                              #, markDF = cellMrkr_parsed
                                                #                                             )
                                                
                                                #FindMarkers(tissSampFilt_nodups, ident.1 = "2")
                                                #ggsave(paste0(plot.path,"UMAP_SCINA.png"),
                                                #       DimPlot(tissSampFilt_nodups,
                                                #               reduction = "umap",
                                                #               group.by = "SCINA_annot"),
                                                #       width = 15)
                                                tissList[[s]] <- tissSampFilt_nodups
                                                # Save RDS file of the individual sample
                                                saveRDS(tissSampFilt_nodups,
                                                        file = sprintf("%s%s_preproc.rds",
                                                                       resuTissDir,
                                                                       s))
                                        }else{
                                                unlink(plot.path, recursive = T)
                                        }
                                }
                        }else{
                                sampTissF <- sampTissFiles[grep(sprintf("%s_",
                                                                        s),
                                                                sampTissFiles)]
                                print(sprintf("%s sample %s has already been preprocessed. Loading RDS file %s",
                                              org,
                                              s,
                                              sampTissF))
                                tissSampFilt_nodups <- readRDS(sampTissF)
                                tissList[[s]] <- tissSampFilt_nodups
                        }
                }
                # Merge Seurat objects of the samples of the same tissue and
                # save it in a RDS file.
                tissPreproc <- merge(tissList[[1]], tissList[2:length(tissList)])
                resuFile <- sprintf("%s%s_preproc.rds", resuPrepDir, org)
                saveRDS(tissPreproc, file = resuFile)
        }else{
                print(sprintf("%s samples have already been preprocessed. Loading RDS files from %s...",
                              org, resuTissDir))
                sampTissFiles <- list.files(resuTissDir, full.names = T)
                for(j in seq_along(sampTissFiles)){
                        sampTissF <- sampTissFiles[j]
                        s <- gsub("_preproc.rds", "", basename(sampTissF))
                        tissList[[s]] <- readRDS(sampTissF)
                }
        }
        #orgList[[org]] <- tissList
}




tiss <- readRDS(tissFiles[1])

tiss@assays$RNA

unique(gsub(".*_", "", colnames(tiss@assays$RNA$counts)))

tissue <- "bladder"
cm_human <- read.table(sprintf("%sCell_marker_Human.txt",
                               utilDir),
                       header = T,
                       sep = "\t",
                       stringsAsFactors = F,
                       quote = "",
                       comment.char = "",
                       encoding = "UTF-8")

cm_human$tissue_class <- gsub("Ê",
                              " ",
                              iconv(cm_human$tissue_class,
                                    from = "ISO-8859-1",
                                    to = "UTF-8"))

cm_human$tissue_class <- gsub("Bone marrow",
                              "bone-marrow",
                              cm_human$tissue_class)
cm_human$tissue_class <- gsub("Skeletal muscle",
                              "skeletal-muscle",
                              cm_human$tissue_class)
cm_human$tissue_class[cm_human$tissue_type == "Retina"] <- "Retina"
cm_human$tissue_class <- tolower(cm_human$tissue_class)



cm_human <- cm_human[which(grepl(tissue,cm_human$tissue_class)),]
cm_human <- cm_human[cm_human$cancer_type == "Normal",]
cm_human <- cm_human[cm_human$cell_type != "Cancer cell",]
cm_human <- cm_human[,c("cell_name","Symbol")]
#cm_human <- ddply(cm_human,.(cell_name,Symbol),nrow)
cm_human <- cm_human[cm_human$Symbol != "",]
#cm_human <- cm_human[cm_human$V1 > 1,]
#cm_human <- cm_human[,c("cell_name","Symbol")]

sigs <- lapply(unique(cm_human$cell_name),function(x){cm_human$Symbol[cm_human$cell_name == x]})
names(sigs) <- unique(cm_human$cell_name)
#sigs <- lapply(sigs,function(x){intersect(x,rownames(inp_df))})






CreateAssayObject <- function(
                counts,
                data,
                min.cells = 0,
                min.features = 0,
                key = NULL,
                check.matrix = FALSE,
                ...
) {
        if (missing(x = counts) && missing(x = data)) {
                abort(message = "Must provide either 'counts' or 'data'")
        } else if (!missing(x = counts) && !missing(x = data)) {
                abort(message = "Either 'counts' or 'data' must be missing; both cannot be provided")
        } else if (!missing(x = counts)) {
                # check that dimnames of input counts are unique
                if (anyDuplicated(x = rownames(x = counts))) {
                        warn(
                                message = "Non-unique features (rownames) present in the input matrix, making unique"
                        )
                        rownames(x = counts) <- make.unique(names = rownames(x = counts))
                }
                if (anyDuplicated(x = colnames(x = counts))) {
                        warn(
                                message = "Non-unique cell names (colnames) present in the input matrix, making unique"
                        )
                        colnames(x = counts) <- make.unique(names = colnames(x = counts))
                }
                if (is.null(x = colnames(x = counts))) {
                        abort(message = "No cell names (colnames) names present in the input matrix")
                }
                if (any(rownames(x = counts) == '')) {
                        abort(message = "Feature names of counts matrix cannot be empty")
                }
                if (nrow(x = counts) > 0 && is.null(x = rownames(x = counts))) {
                        abort(message = "No feature names (rownames) names present in the input matrix")
                }
                if (!inherits(x = counts, what = 'dgCMatrix')) {
                        if (inherits(x = counts, what = "data.frame")) {
                                counts <- as.sparse(x = counts, ...)
                        } else {
                                counts <- as.sparse(x = counts)
                        }
                }
                if (isTRUE(x = check.matrix)) {
                        CheckMatrix(object = counts)
                }
                print(class(counts))
                # Filter based on min.features
                if (min.features > 0) {
                        nfeatures <- Matrix::colSums(x = counts > 0)
                        counts <- counts[, which(x = nfeatures >= min.features)]
                }
                # filter genes on the number of cells expressing
                if (min.cells > 0) {
                        print("jeje")
                        print(counts)
                        print(dim(counts))
                        print(class(counts))
                        num.cells <- Matrix::rowSums(x = counts > 0)
                        print("jojo")
                        counts <- counts[which(x = num.cells >= min.cells), ]
                }
                data <- counts
        } else if (!missing(x = data)) {
                # check that dimnames of input data are unique
                if (anyDuplicated(x = rownames(x = data))) {
                        warn(
                                message = "Non-unique features (rownames) present in the input matrix, making unique"
                        )
                        rownames(x = data) <- make.unique(names = rownames(x = data))
                }
                if (anyDuplicated(x = colnames(x = data))) {
                        warn(
                                message = "Non-unique cell names (colnames) present in the input matrix, making unique"
                        )
                        colnames(x = data) <- make.unique(names = colnames(x = data))
                }
                if (is.null(x = colnames(x = data))) {
                        abort(message = "No cell names (colnames) names present in the input matrix")
                }
                if (any(rownames(x = data) == '')) {
                        abort(message = "Feature names of data matrix cannot be empty", call. = FALSE)
                }
                if (nrow(x = data) > 0 && is.null(x = rownames(x = data))) {
                        abort(message = "No feature names (rownames) names present in the input matrix")
                }
                if (min.cells != 0 | min.features != 0) {
                        warn(
                                message = "No filtering performed if passing to data rather than counts"
                        )
                }
                counts <- new(Class = 'matrix')
        }
        # Ensure row- and column-names are vectors, not arrays
        if (!is.vector(x = rownames(x = counts))) {
                rownames(x = counts) <- as.vector(x = rownames(x = counts))
        }
        if (!is.vector(x = colnames(x = counts))) {
                colnames(x = counts) <- as.vector(x = colnames(x = counts))
        }
        if (!is.vector(x = rownames(x = data))) {
                rownames(x = data) <- as.vector(x = rownames(x = data))
        }
        if (!is.vector(x = colnames(x = data))) {
                colnames(x = data) <- as.vector(x = colnames(x = data))
        }
        counts <- CheckFeaturesNames(data = counts)
        data <- CheckFeaturesNames(data = data)
        # Initialize meta.features
        init.meta.features <- data.frame(row.names = rownames(x = data))
        misc <- if (.GetSeuratCompat() < '5.0.0') {
                list()
        } else {
                calcN_option <- getOption(
                        x = 'Seurat.object.assay.calcn',
                        default =  Seurat.options$Seurat.object.assay.calcn
                )
                list(calcN = calcN_option %||% TRUE)
        }
        assay <- new(
                Class = 'Assay',
                counts = counts,
                data = data,
                scale.data = new(Class = 'matrix'),
                key = Key(object = key)[1L] %||% '',
                meta.features = init.meta.features,
                misc = misc
        )
        return(assay)
}

CreateAssayObject(contSamp, min.cells = minCells, min.features = minFeats)

min.features <- minFeats
counts <- contSamp
nfeatures <- Matrix::colSums(x = counts > 0)
counts <- counts[, which(x = nfeatures >= min.features)]
