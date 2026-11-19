# Load necessary libraries
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(monocle3)
library(SeuratWrappers)
library(EnhancedVolcano)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(harmony)
library(igraph)
library(magrittr)
library(JASPAR2020)
library(JASPAR2024)
library(motifmatchr)
library(TFBSTools)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(xgboost)
library(presto)



# Define file paths for matrices right now I am excluding the first patient in this paper because the basic read quality metrics were way off

#devtools::install_github('smorabit/hdWGCNA', ref='dev')


# ---------------------------
#  Paste GSM table omg I love chatgpt, it just rewrote everything formatted correctly from the geo site...
#nvm, it made a mistake, be very careful
# ---------------------------
gsm_data <- tribble(
  ~GSM, ~full_label,


  "GSM6186022", "ALS4 GM MA001",
  "GSM6186023", "ALS4 GM MA002",
  "GSM6186024", "ALS4 FN MA003",
  "GSM6186025", "ALS4 SC MA004",
  "GSM6186046", "ALS1 GM PM036",
  "GSM6186047", "ALS1 GM PM037",
  "GSM6186048", "ALS1 FN PM038",
  "GSM6186049", "ALS1 SC PM039",
  "GSM6186050", "ALS2 GM PM040",
  "GSM6186051", "ALS2 GM PM041",
  "GSM6186052", "ALS2 FN PM042",
  "GSM6186053", "ALS2 SC PM043",
  "GSM6186054", "ALS3 GM PM045",
  "GSM6186055", "ALS3 GM PM046",
  "GSM6186056", "ALS3 FN PM047",
  "GSM6186057", "ALS3 SC PM048",
  "GSM6186058", "ALS5 GM PM051",
  "GSM6186059", "ALS5 SC PM052",
  "GSM6186065", "ALS6 BA9+BA4+SC PM058",
  "GSM6186069", "ALS7 BA9+BA4+SC PM062",
  "GSM6186070", "ALS/FTD1 BA9+BA4+SC PM063",
  "GSM6186076", "ALS8 BA9+BA4+SC PM069",
  "GSM6186082", "ALS9 BA9+BA4+SC PM080",
)


#"GSM6186081", "LOAD28 BA9+BA20+H PM078",
#"GSM6186083", "LOAD31 BA9+BA20+H PM082",
#"GSM6186084", "LOAD32 BA9+BA20+H PM083",
#"GSM6186085", "LOAD33 BA9+BA20+H PM084",
#"GSM6186086", "LOAD34 BA9+AWS PM085",
#"GSM6186087", "LOAD35 BA9+BA20+H;HD1 BA9+BA20+H PM086",
#"GSM6186078", "LOAD26 GM PM075",
#"GSM6186079", "LOAD27 BA9+H PM076",
#"GSM6186072", "LOAD20 BA9+BA20+H PM065",
#"GSM6186073", "EOAD2 BA9+BA20+H PM066",
#"GSM6186074", "LOAD21 BA9+BA20+H PM067",
#"GSM6186066", "LOAD1 BA9+BA20+H PM059",
#"GSM6186067", "LOAD2 BA9+BA20+H PM060",
#"GSM6186068", "EOAD1 BA9+BA20+H PM061",
#"GSM6186030", "LOAD12 GM PM013",
#"GSM6186031", "LOAD12 AH PM014",
#"GSM6186032", "LOAD13 GM PM016",
#"GSM6186033", "LOAD13 H PM017",
#"GSM6186034", "LOAD13 SN PM018",
#"GSM6186035", "LOAD14 GM PM020",
#"GSM6186036", "LOAD14 H PM021",
#"GSM6186037", "LOAD14 SN PM022",
#"GSM6186038", "LOAD15 GM PM023",
#"GSM6186039", "LOAD15 H PM024",
#"GSM6186040", "LOAD15 SN PM025",
#"GSM6186021", "LOAD11 GM DeJager_MO_LOAD11_GM",




#  "GSM6186077", "CNTRL1 BA9+BA4+H+SN PM073",
#"GSM6186081", "LOAD28 BA9+BA20+H PM078",
#"GSM6186082", "ALS9 BA9+BA4+SC PM080",
#"GSM6186083", "LOAD31 BA9+BA20+H PM082",
#"GSM6186084", "LOAD32 BA9+BA20+H PM083",
#"GSM6186085", "LOAD33 BA9+BA20+H PM084",
#"GSM6186086", "LOAD34 BA9+AWS PM085",
#"GSM6186087", "LOAD35 BA9+BA20+H;HD1 BA9+BA20+H PM086"
#"GSM6186077", "CNTRL1 BA9+BA4+H+SN PM073",
#"GSM6186078", "LOAD26 GM PM075",
#"GSM6186079", "LOAD27 BA9+H PM076",
#"GSM6186071", "MS2 BA9+TH+Lesion+AWS PM064",
#"GSM6186072", "LOAD20 BA9+BA20+H PM065",
#"GSM6186073", "EOAD2 BA9+BA20+H PM066",
#"GSM6186074", "LOAD21 BA9+BA20+H PM067",
#"GSM6186075", "DLBD-PD4 BA9+SN PM068",
#"GSM6186066", "LOAD1 BA9+BA20+H PM059",
#"GSM6186067", "LOAD2 BA9+BA20+H PM060",
#"GSM6186068", "EOAD1 BA9+BA20+H PM061",
#"GSM6186060", "DLBD-PD3 GM PM053",
#"GSM6186061", "DLBD-PD3 SN PM054",
#"GSM6186062", "MS1 AWS PM055",
#"GSM6186063", "MS1 GM PM056",
#"GSM6186064", "MS1 TH PM057",
#"GSM6186026", "DLBD-PD2 GM PM004",
#"GSM6186027", "DLBD-PD2 SN PM005",
#"GSM6186028", "DLBD-PD1 GM PM006",
#"GSM6186029", "DLBD-PD1 SN PM007",
#"GSM6186030", "LOAD12 GM PM013",
#"GSM6186031", "LOAD12 AH PM014",
#"GSM6186032", "LOAD13 GM PM016",
#"GSM6186033", "LOAD13 H PM017",
#"GSM6186034", "LOAD13 SN PM018",
#"GSM6186035", "LOAD14 GM PM020",
#"GSM6186036", "LOAD14 H PM021",
#"GSM6186037", "LOAD14 SN PM022",
#"GSM6186038", "LOAD15 GM PM023",
#"GSM6186039", "LOAD15 H PM024",
#"GSM6186040", "LOAD15 SN PM025",
#"GSM6186041", "PSP2 H PM026",
#"GSM6186042", "PSP2 SN PM027",
#"GSM6186043", "PSP3 GM PM029",
#"GSM6186044", "FTD1 GM PM030",
#"GSM6186045", "FTD1 H PM031",
#"GSM6186021", "LOAD11 GM DeJager_MO_LOAD11_GM",






# ---------------------------
#  Extract metadata (all diagnoses, exact brain regions)
# ---------------------------
sample_metadata <- gsm_data %>%
  mutate(
    diagnosis    = stringr::word(full_label, 1),                     # just take first word as diagnosis
    brain_region = stringr::word(full_label, 2) |>                   # second word is always brain area
      stringr::str_replace_all("[;,:]+$", ""),          # drop trailing punctuation
    sample_id    = stringr::word(full_label, -1)                     # last token (PM#### / MA### / DeJager_...)
  ) %>%

  dplyr::select(GSM, sample_id, diagnosis, brain_region)




#load all data
data_dir  <- "/Volumes/PC60/scRNASeqAnalyses/ALS_scRNASEQ/GSE204702_RAW"
h5_files  <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)
seurat_list <- list() # seurat object to load

#gsm_id <- str_extract(basename(file), "GSM\\d+")
#meta_row <- sample_metadata %>% filter(GSM == gsm_id)


for (file in h5_files) {
  base <- basename(file)

  # Try to identify by PM/MA first, else GSM
  key_pm_ma <- str_extract(base, "(PM\\d+|MA\\d+)")
  key_gsm   <- str_extract(base, "GSM\\d+")

  meta_row <- NULL
  if (!is.na(key_pm_ma)) {
    meta_row <- sample_metadata %>% dplyr::filter(sample_id == key_pm_ma)

  } else if (!is.na(key_gsm)) {
    meta_row <- sample_metadata %>% dplyr::filter(GSM == key_gsm)

  } else {
    message("Skipping (no PM/MA/GSM match in filename): ", base)
    next
  }

  if (nrow(meta_row) != 1) {
    message("Skipping (no unique metadata match): ", base)
    next
  }

  message("Processing: ", base, "  ->  ",
          meta_row$diagnosis, " | ", meta_row$brain_region, " | ", meta_row$sample_id)

  counts <- Read10X_h5(file)
  if (is.list(counts)) counts <- counts[["Gene Expression"]]

  seu <- CreateSeuratObject(
    counts = counts,
    project = meta_row$sample_id,
    min.features = 200,
    min.cells = 3
  )

  seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
  ribo_genes            <- grep("^RP[SL]", rownames(seu), value = TRUE) # I found during an earlier analysis that the ribo genes were a significant source of variation, I will consider regressing them out
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, features = ribo_genes)

  # Attach metadata (keep both GSM and sample_id)
  seu$GSM         <- meta_row$GSM
  seu$sample_id   <- meta_row$sample_id
  seu$diagnosis   <- meta_row$diagnosis
  seu$brain_region<- meta_row$brain_region

  # Apply filters Here I am removing cells with less than 200 or greater than 2500 individual transcripts (should get rid of many doublets and other non-descriptive low count cells)
  seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)

  seurat_list[[meta_row$sample_id]] <- seu
}

#  Merge everything into one Seurat object
if (length(seurat_list) > 1) {
  combined <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list))
} else if (length(seurat_list) == 1) {
  combined <- seurat_list[[1]]
} else {
  stop("No Seurat objects were created—check filename/ID matching.")
}



#checking to make sure it looks kosher

p <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, raster=FALSE)
print(p)







#normalizing data (this is all the default currently): for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p
allCells <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)

# let's start by just nabbing highly differentially expressed genes
allCells <- FindVariableFeatures(allCells, selection.method = "vst", nfeatures = 2000)

# Identify the 20 most highly variable genes
top100 <- head(VariableFeatures(allCells), 100)






# Scale the data
allCells <- ScaleData(allCells, vars.to.regress = c("percent.mt",   "nCount_RNA")) #I'm not sure if this is typical but I find that ribosomal genes dominate differential signature

#In the tuddenham paper they just remove all ribo and mitochondrial genes.... I think this is preferable
allCells <- RunPCA(allCells, features = VariableFeatures(object = allCells))

allCells <- RunHarmony(allCells, group.by.vars = "sample_id") #running harmony to standardize

#run PCA


#checking out loadings
VizDimLoadings(allCells, dims = 1:4, reduction = "harmony")

#basic PCA plot w/ 2000 most diff expressed genes
DimPlot(allCells, reduction = "harmony") + NoLegend()

#now looking at cells from across the spectrum to really check out how things are divided

DimHeatmap(allCells, dims = 1:5, cells = 500,  balanced = TRUE)

#now an elbow plot to see which PCAs get fed into UMAP
ElbowPlot(allCells)



# based on elbow plot I will use PCA 1:20

dims_use <- 1:20
set.seed(42)

allCells <- FindNeighbors(allCells, reduction = "harmony", dims = dims_use, k.param = 30)
allCells <- FindClusters(allCells, resolution = 0.5)  # sweep later if needed

allCells <- RunUMAP(
  allCells,
  reduction   = "harmony",
  dims        = dims_use,
  n.neighbors = 30,
  min.dist    = 0.3,
  umap.method = "uwot",   # default
  metric      = "cosine"  # default for PCA/harmony I just like being explicit
)



DimPlot(allCells, reduction = "umap",label=TRUE) #plot it all

FeaturePlot(allCells, features = c("CHIT1","JAML","GPNMB","PLA2G7","S100A4","CD74"))


FeaturePlot(allCells, features = c("CHIT1","JAML","GPNMB","PLA2G7","S100A4","CD74",
                                 "AKR1B1","CAPG","CD48","CHI3L1", "CLEC5A","CPB2",
                                "F11R","HAVCR2","SERPINA4","TGFBR1","TNFSF8" ))



microglia_only <- subset(allCells, idents = c(0,1,2,3,5,7))

DimPlot(microglia_only, reduction = "umap")+
  ggtitle("Microglial Subclusters in ALS") #plot it all





genes <- c("CHIT1","JAML","GPNMB","PLA2G7","S100A4")

genes <- genes[genes %in% rownames(microglia_only)]

# Pull expression as FeaturePlot does
expr_df <- FetchData(microglia_only, vars = genes)  # cells x genes

# Average log1p expression across genes
module_score <- rowMeans(log1p(as.matrix(expr_df)))  # log(1 + x)

# =z-score for nicer spread if ya want
module_score_z <- scale(module_score)[, 1]

# Store in metadata
microglia_only$DAM_score <- module_score

# Plot
FeaturePlot(microglia_only, features = "DAM_score")+
ggtitle("Mean Expression of CHIT1 JAML GPNMB PLA2G7 and S100A4 in Microglia")


FeaturePlot(microglia_only, features = c("GPNMB"))+
  ggtitle("GPNMB Expression in Microglia")




FeaturePlot(microglia_only, features= c("GPNMB","CD74","P2RY12","LGALS3"))

FeaturePlot(allCells, features = c("CHI3L2","FCGR2A","FCGR3A","SERPINA5","SERPINA6","VSIG4"))


FeaturePlot(allCells, features = c("CD8","CD2","CD74","GZMK","TREM2","CX3CR1","P2RY12"))


allCells$MajorCellType <- "Other"

allCells$MajorCellType[ Idents(allCells) %in% c( 5) ] <- "DAM Microglia" #unfortunately in this dataset the GPNMB high cells are too sparse and get knocked out thus just focusing on this one cluster

Idents(allCells) <- "MajorCellType"

allCells <- JoinLayers(allCells)






#just microglia excludign t cells and whatnot
microglia_only <- JoinLayers(microglia_only)
GPNMBCluster.markers <- FindMarkers(microglia_only, ident.1 = 5, ident.2=c(0,1,2,3,7) )



#saveRDS(microglia_only, file = "/Volumes/PC60/scRNASeqAnalyses/ALS_scRNASEQ/MicrogliaTuddenhamEtal.rds")



EnhancedVolcano(GPNMBCluster.markers ,
                rownames(GPNMBCluster.markers ),
                x ="avg_log2FC",
                y ="p_val",
                title = 'GPNMB High Cluster DEGs in ALS',
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)




write.csv(GPNMBCluster.markers, "/Volumes/PC60/scRNASeqAnalyses/adSCRNASeq/TuddenhamEtAl_GPNMBHighCluster_DEGs_updated.csv")


EnhancedVolcano(clust1.markers  ,
                rownames(clust1.markers ),
                x ="avg_log2FC",
                y ="p_val",
                title = 'Putative DAM cluster vs Putative Transition Cluster',
                pCutoff = 10e-6,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 3.0)






DefaultAssay(microglia_only)

expr <- GetAssayData(microglia_only, slot = "data")["GPNMB", ]

sum(expr > 0)/ sum(expr>-1)  #13.4% are positive for GPNMB in ALS




microglia_only<-readRDS("/Volumes/PC60/scRNASeqAnalyses/ALS_scRNASEQ/MicrogliaTuddenhamEtal.rds")



toTest <- subset(microglia_only, idents=5)

exprinCluster <- GetAssayData(toTest, slot = "data")["GPNMB", ]
exprinTot <- GetAssayData(microglia_only, slot = "data")["GPNMB", ]


sum(exprinCluster > 0)/ sum(exprinTot>-1)  # 6.4% are found in cluster for ALS



P2RY12expr <- GetAssayData(microglia_only, slot="data")["P2RY12",]
GPNMBexpr <- GetAssayData(microglia_only, slot="data")["GPNMB",]

sum(GPNMBexpr > 0 & P2RY12expr==0  )/ sum(P2RY12expr > 0 & GPNMBexpr==0      )  # .24 is the ratio for ALS




microglia_only<-readRDS("/Volumes/PC60/scRNASeqAnalyses/ALS_scRNASEQ/MicrogliaTuddenhamEtal.rds")





selected_cells <- CellSelector(plot = DimPlot(microglia_only, reduction = "umap"))

microglia_only <- subset(microglia_only, cells = selected_cells)



FeaturePlot3(
  seu       = microglia_only,
  feature.1 = "GPNMB",
  feature.2 = "P2RY12",
  color     = "rgb",
  reduction = "umap",
  pt.size   = 0.2,
  order     = TRUE,
  color.range = c(0, 1)
) +
  plot_annotation(title = "ALS Microglia")



FeaturePlot(microglia_only, features = c("GPNMB", "P2RY12"), blend = TRUE)


FeaturePlot(microglia_only, features='GPNMB')+
  ggtitle( "GPNMB Expression in ALS Microglia")





allCellsAtPlay <-subset(microglia_only, idents = 5)


allCellsAtPlay <- JoinLayers(allCellsAtPlay)


#ok now we are going to start with wcgna


# start by setting up

# this was all taken from the tutorial  https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html


#by default this filters for the variable genes identified in findvariable genes above, thus is 2000 righ tnow
allCellsAtPlay <- SetupForWGCNA(
  allCellsAtPlay,
  wgcna_name = "alsMicroglia" # the name of the hdWGCNA experiment
)

allCellsAtPlay$WGCNA_group <- "DAM_GPNMB_high"



allCellsAtPlay <- MetacellsByGroups(
  seurat_obj  = allCellsAtPlay,
  group.by    = "WGCNA_group",  # all DAM pooled
  reduction   = "harmony",
  k           = 25,             #
  max_shared  = 10,
  min_cells   = 30,
  ident.group = "WGCNA_group"
)


# normalize metacell expression matrix:
allCellsAtPlay <- NormalizeMetacells(allCellsAtPlay)


#ok nice now on to frank analysis



allCellsAtPlay <- SetDatExpr(
  allCellsAtPlay,
  group_name = "DAM_GPNMB_high", # the name of the group of interest in the group.by column
  group.by='WGCNA_group', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)


#selecting soft power threshold


# Test different soft powers:
allCellsAtPlay <- TestSoftPowers(
  allCellsAtPlay,
  networkType = 'signed' #  can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(allCellsAtPlay)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(allCellsAtPlay)
head(power_table)


#ok based on this we should use 4 as our soft power threshold, not sure how automatic thresholding works, should test this


# construct co-expression network:
allCellsAtPlay <- ConstructNetwork(
  allCellsAtPlay,
  soft_power = 4,
  tom_name = 'Microglia', # name of the topoligical overlap matrix written to disk,
  overwrite_tom = TRUE,
  minModuleSize = 25,
  mergeCutHeight = 0.2,
  detectCutHeight = 0.95,
  deepSplit = 4

)



# look at dendrogram
PlotDendrogram(allCellsAtPlay, main='GPNMB High Microglial Cluster  ALS ')



# compute all MEs in the full single-cell dataset
allCellsAtPlay <- ModuleEigengenes(
  allCellsAtPlay,
  group.by.vars="sample_id"
)


# compute eigengene-based connectivity (kME):
allCellsAtPlay <- ModuleConnectivity(
  allCellsAtPlay,
  group.by.vars="sample_id"
)


allCellsAtPlay$metacell_grouping


p <- PlotKMEs(allCellsAtPlay)

p




# get the module assignment table:
modules <- GetModules(allCellsAtPlay) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])


# get hub genes
hub_df <- GetHubGenes(allCellsAtPlay, n_hubs = 15)


#write.csv(hub_df, "/Volumes/PC60/ALS_scRNASEQ/HubGenes.csv")




blue_hubs <- subset(hub_df, module == "blue")

library(ggplot2)

ggplot(blue_hubs,
       aes(x = reorder(gene_name, kME), y = kME)) +
  geom_col(fill = "blue") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "kME values for Brown Module Hub Genes",
    x = "Gene",
    y = "kME (connectivity)"
  )




green_hubs <- subset(hub_df, module == "green")

library(ggplot2)

ggplot(green_hubs,
       aes(x = reorder(gene_name, kME), y = kME)) +
  geom_col(fill = "forestgreen") +
  coord_flip() +
  theme_bw() +
  labs(
    title = "kME values for Green Module Hub Genes",
    x = "Gene",
    y = "kME (connectivity)"
  )





library(UCell)
allCellsAtPlay <- ModuleExprScore(
  allCellsAtPlay,
  n_genes = 25,
  method='UCell'
)





# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  allCellsAtPlay,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)


# stitch together with patchwork
wrap_plots(plot_list, ncol=6)





# make a featureplot of hub scores for each module
plot_list <- ModuleFeaturePlot(
  allCells,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE
)



## okay now onto plotting these network relationships




ModuleNetworkPlot(
  allCellsAtPlay,
  outdir='ModuleNetworks2', # new folder name
  n_inner = 15, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)


# get the list of modules:
modules <- GetModules(allCellsAtPlay)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
HubGeneNetworkPlot(
  allCellsAtPlay,
  n_hubs = 5, n_other=10,
  edge_prop = 0.75,
  edge.alpha=0.25,
  mods = mods[1:7] # only select 5 modules
)




allCellsAtPlay <- RunModuleUMAP(
  allCellsAtPlay,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)



# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(allCellsAtPlay)

# plot with ggplot
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()




ModuleUMAPPlot(
  allCellsAtPlay,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (10% here)
  label_hubs=3 ,# how many hub genes to plot per module I think?
  keep_grey_edges=FALSE

)

## and transcription factor network analysis




#ok first we will output position weight matrices for each individual TF thesw will then be used to scan all of our module inhabitants


JASPAR2024 <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

pfm_core <- TFBSTools::getMatrixSet(
  x = sq24,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)



# now we run the motif scan
allCellsAtPlay <- MotifScan(
  allCellsAtPlay,
  species_genome = 'hg38',
  pfm = pfm_core,
  EnsDb = EnsDb.Hsapiens.v86
)


# get the motif df:
motif_df <- GetMotifs(allCellsAtPlay)



# Overlap TF motif targets with modules
allCellsAtPlay <- OverlapModulesMotifs(allCellsAtPlay)

#  Get the overlap table
overlap_df <- GetMotifOverlap(allCellsAtPlay)
head(overlap_df)





# keep all TFs, and then remove all genes from the grey module
tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(allCellsAtPlay)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

# update the gene list and re-run SetDatExpr
allCellsAtPlay <- SetWGCNAGenes(allCellsAtPlay, genes_use)
allCellsAtPlay <- SetDatExpr(allCellsAtPlay, group.by = 'WGCNA_group', group_name='DAM_GPNMB_high')



# define model params:
model_params <- list(
  objective = 'reg:squarederror',
  max_depth = 1,
  eta = 0.1,
  nthread=16,
  alpha=0.5
)

# construct the TF network
allCellsAtPlay <- ConstructTFNetwork(allCellsAtPlay, model_params)



results <- GetTFNetwork(allCellsAtPlay)
head(results)





allCellsAtPlay <- AssignTFRegulons(
  allCellsAtPlay,
  strategy = "A",   # Strategy “A” selects the top TFs for each gene.kind of funny they did it this way lol
  reg_thresh = 0.01,
  n_tfs = 10
)




p2 <- RegulonBarPlot(allCellsAtPlay, selected_tf='MITF')

 p2


 brown_hub_genes <- brown_hubs$gene_name





 results <- GetTFNetwork(allCellsAtPlay)

 # Keep only TF→target edges where the target is a brown hub gene
 brown_tf_net <- results %>%
   dplyr::filter(gene %in% brown_hub_genes)


 top_n_tfs <- 10  #

 top_tf_per_hub <- brown_tf_net %>%
   group_by(gene) %>%
   slice_max(order_by = Gain, n = top_n_tfs, with_ties = FALSE) %>%
   arrange(gene, desc(Gain)) %>%
   ungroup()

 top_tf_per_hub <- top_tf_per_hub %>%
   left_join(
     brown_hubs %>% select(gene_name, kME),
     by = c("gene" = "gene_name")
   ) %>%
   rename(
     target_gene = gene,
     tf_gene     = tf
   )

 write.csv(top_tf_per_hub, "/Volumes/PC60/ALS_scRNASEQ/TFAssociatedWithGPNMBHubGenes.csv")





green_hubs <- subset(hub_df, module == "green")


green_hub_genes <- green_hubs$gene_name


green_tf_net <- results %>%
  filter(gene %in% green_hub_genes)


top_n_tfs <- 10  #

top_tf_per_hub <- green_tf_net %>%
  group_by(gene) %>%
  slice_max(order_by = Gain, n = top_n_tfs, with_ties = FALSE) %>%
  arrange(gene, desc(Gain)) %>%
  ungroup()

top_tf_per_hub <- top_tf_per_hub %>%
  left_join(
    green_hubs %>% select(gene_name, kME),
    by = c("gene" = "gene_name")
  ) %>%
  rename(
    target_gene = gene,
    tf_gene     = tf
  )


write.csv(top_tf_per_hub, "/Volumes/PC60/ALS_scRNASEQ/TFAssociatedWithLPLCluster.csv")



brown_df <- read.csv("/Volumes/PC60/ALS_scRNASEQ/TFAssociatedWithGPNMBHubGenes.csv")
green_df <- read.csv("/Volumes/PC60/ALS_scRNASEQ/TFAssociatedWithLPLCluster.csv")


brown_genes <- unique(brown_df$target_gene)
green_genes <- unique(green_df$target_gene)



# How many targets per TF?
brown_tf_table <- table(brown_df$tf_gene)

# Keep TFs with >= 5 targets
brown_tf_ge5 <- names(brown_tf_table[brown_tf_table >= 5])
brown_tf_ge1 <- names(brown_tf_table[brown_tf_table > 1])

brown_member_tfs <- intersect(unique(brown_df$tf_gene), brown_genes)
brown_member_tfs
brown_consistent_tfs <- sort(unique(c(brown_tf_ge5, brown_member_tfs)))
brown_consistent_tfs




green_tf_table <- table(green_df$tf_gene)

green_tf_ge5 <- names(green_tf_table[green_tf_table >= 5])
green_tf_ge1 <- names(green_tf_table[green_tf_table > 1])

green_member_tfs <- intersect(unique(green_df$tf_gene), green_genes)
green_member_tfs
green_consistent_tfs <- sort(unique(c(green_tf_ge5, green_member_tfs)))
green_consistent_tfs



shared_tfs <- intersect(brown_consistent_tfs, green_consistent_tfs)
shared_tfs

brown_only_tfs <- setdiff(brown_consistent_tfs, green_tf_ge1)
brown_only_tfs


green_only_tfs <- setdiff(green_consistent_tfs, brown_tf_ge1)
green_only_tfs








