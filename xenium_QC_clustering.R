
#-------------------------------------------------------------------------
# Set up CLI options with optparse
#-------------------------------------------------------------------------

library(optparse)

option_list <- list(
    make_option(
        c('--indir'), type='character', default=NULL,
        help='Directory with input files, should be the output of the 10X software.', metavar='character'
    ),
    make_option(
        c('--panel'), type='character', default=NULL,
        help='Path to .csv file containing the Xenium panel information', metavar='character'
    ),
    make_option(
        c('--sampleid'), type='character', default=NULL,
        help='Identifier for the biological sample', metavar='character'
    ),
    make_option(
        c('--slideid'), type='character', default=NULL,
        help='Identifier for the Xenium slide', metavar='character'
    ),
    make_option(
        c('--outdir'), type='character', default='./',
        help='Directory to place output files', metavar='character'
    ),
    make_option(
        c('--nCount-thresh'), type='numeric', default=10,
        help='Minimum number of counts required to pass QC', metavar='character'
    ),
    make_option(
        c('--nGene-thresh'), type='numeric', default=10,
        help='Minimum number of genes required to pass QC', metavar='character'
    ),
    make_option(
        c('--negcontrol-thresh'), type='numeric', default=5,
        help='Max number of negative control counts to allow', metavar='character'
    ),
    make_option(
        c('--cluster_res'), type='numeric', default=0.7,
        help='Resolution parameter', metavar='character'
    )
)

#------------------------------------------------------------------------
# Load R Packages
#------------------------------------------------------------------------

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if(!dir.exists(opt$outdir)){
  tt <- tryCatch(
    dir.create(opt$outdir),
    error=function(e) e,
    warning=function(w) w
  )
  if(is(tt,"warning")){
    stop(paste0("Permission denied, cannot create output directory ", opt$outdir))
  }
}
setwd(opt$outdir)

fig_dir <- 'figures/'
data_dir <- 'data/'

dir.create(fig_dir)
dir.create(data_dir)

#------------------------------------------------------------------------
# Load R Packages
#------------------------------------------------------------------------

library(Voyager)
library(SpatialFeatureExperiment)
library(rjson)
library(Matrix)
library(SFEData)
library(RBioFormats)

library(SingleCellExperiment)
library(SpatialExperiment)
library(ggplot2)
library(stringr)
library(scater) 
library(scuttle)
library(BiocParallel)
library(BiocSingular)
library(bluster)
library(scran)
library(patchwork)

library(tidyverse)
library(cowplot)
library(Seurat)
theme_set(theme_cowplot())

#------------------------------------------------------------------------
# Load the SFE object, Xenium panel
#------------------------------------------------------------------------

# load the SFE object
try(sfe <- readXenium(opt$indir, add_molecules=TRUE))
sfe <- readXenium(opt$indir, add_molecules=TRUE)

# load the Xenium panel metadata
panel_meta <- read.csv(opt$panel)
rownames(panel_meta) <- panel_meta$Gene
panel_meta <- subset(panel_meta, Gene %in% rownames(sfe))
rowData(sfe)[panel_meta$Gene, 'Annotation'] <- panel_meta$Annotation

# add sample name and slide name 
cur_sample <- opt$sampleid
colData(sfe)$SampleID <- cur_sample
colData(sfe)$slideID <- opt$slideid

# add centroid coordinates to the colData
coords <- spatialCoords(sfe)
colData(sfe)$x_centroid <- coords[,1]
colData(sfe)$y_centroid <- coords[,2]


#------------------------------------------------------------------------
# Cell Segmentation plots
#------------------------------------------------------------------------

# plot the geometry 
p <- plotGeometry(sfe, "cellSeg")
pdf(paste0(fig_dir, 'cell_segmentation_geometry.pdf'), width=5, height=5)
print(p)
dev.off()

# plot cell density 
p <- plotCellBin2D(sfe) + coord_equal()
pdf(paste0(fig_dir, 'cell_segmentation_density.pdf'), width=5, height=5)
print(p)
dev.off()

#------------------------------------------------------------------------
# Quality Control 
#------------------------------------------------------------------------

# calculate the number of counts and genes expressed in each cell 
tmp <- counts(sfe)
sfe$nCounts <- colSums(tmp)
tmp[tmp != 0] <- 1
sfe$nGenes <- colSums(tmp)

# calculate the proportion of nuc area / cell area 
colData(sfe)$prop_nuc <- sfe$nucleus_area / sfe$cell_area

# QC stat distributions
p <- plotColDataHistogram(sfe, c("nCounts", "nGenes", "cell_area", "nucleus_area"), ncol=2)
pdf(paste0(fig_dir, 'qc_stats_histogram.pdf'), width=5, height=5)
print(p)
dev.off()

p <- scater::plotColData(sfe, x="nCounts", y="nGenes", bins = 100)
pdf(paste0(fig_dir, 'qc_genes_vs_counts.pdf'), width=5, height=5)
print(p)
dev.off()

p <- plotColData(sfe, x='cell_area', y='nucleus_area', bins=100)
pdf(paste0(fig_dir, 'qc_cell_vs_nuc.pdf'), width=5, height=5)
print(p) 
dev.off()

p <- plotSpatialFeature(sfe, "prop_nuc", colGeometryName = "cellSeg")
pdf(paste0(fig_dir, 'qc_proportion_nuclei_per_cell.pdf'), width=5, height=5)
print(p) 
dev.off()


#------------------------------------------------------------------------
# Quality Control with Neg Control Features
#------------------------------------------------------------------------

sfe <- addPerCellQCMetrics(
    sfe, subsets = list(
        "negCodeword" = rowData(sfe)$Type == 'Negative Control Codeword',
        "negProbe" = rowData(sfe)$Type == 'Negative Control Probe',
        "unassigned" = rowData(sfe)$Type == 'Unassigned Codeword'
    )
)

cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]

p <- plotColDataHistogram(sfe, cols_use, bins = 100, ncol = 3) + 
    scale_x_log10() +
    annotation_logticks(sides = "b")

pdf(paste0(fig_dir, 'qc_negControl_histogram.pdf'), width=6, height=3)
print(p)
dev.off()


nc_thresh <- opt[["negcontrol-thresh"]]
colData(sfe)$control_percent <- colData(sfe)$subsets_negCodeword_percent + colData(sfe)$subsets_negProbe_percent + colData(sfe)$subsets_unassigned_percent
colData(sfe)$is_nc_outlier <- colData(sfe)$subsets_negCodeword_percent > nc_thresh | 
    colData(sfe)$subsets_negProbe_percent > nc_thresh | 
    colData(sfe)$subsets_unassigned_percent > nc_thresh 


#------------------------------------------------------------------------
# Quality Control of Genes
#------------------------------------------------------------------------

rowData(sfe)$means <- rowMeans(counts(sfe))
rowData(sfe)$vars <- rowVars(counts(sfe))

p <- plotRowData(sfe, x = "means", y = "Type") +
    scale_y_log10() +
    annotation_logticks(sides = "b")

pdf(paste0(fig_dir, 'qc_Genes_mean_expression.pdf'), width=6, height=3)
print(p)
dev.off()


p <- plotRowData(sfe, x="means", y="vars", color_by = "Type") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    scale_x_log10() + scale_y_log10() +
    annotation_logticks() +
    coord_equal() +
    labs(color = "Negative control")

pdf(paste0(fig_dir, 'qc_Genes_mean_var_expression.pdf'), width=6, height=6)
print(p)
dev.off()


#------------------------------------------------------------------------
# QC Filtering of low-quality cells
# * No nucleus found 
# * Expression of negative controls 
# * Very low transcript counts
#------------------------------------------------------------------------

count_thresh <- opt[["nCount-thresh"]]
gene_thresh <- opt[["nGene-thresh"]]

colData(sfe)$low_quality <- is.na(colData(sfe)$nucleus_area) | colData(sfe)$nCounts < count_thresh | colData(sfe)$nGenes < gene_thresh | colData(sfe)$is_nc_outlier

print('number of LQ cells:')
print(table(colData(sfe)$low_quality))

# where are the low quality cells ?
p <- plotSpatialFeature(
    sfe, "low_quality", 
    colGeometryName = "cellSeg",
    dark=TRUE, image_id='morphology_focus',
    channel = 2,
    maxcell=10e+06
)
pdf(paste0(fig_dir, 'qc_low_quality_spatial.pdf'), width=6, height=6)
print(p) 
dev.off()

# filter out the lowquality cells:
sfe <- sfe[,!colData(sfe)$low_quality]
dim(sfe)


#------------------------------------------------------------------------
# Moran's I of QC features
# not super informative here
#------------------------------------------------------------------------

colGraph(sfe, "knn5") <- findSpatialNeighbors(
    sfe, 
    method = "knearneigh", 
    dist_type = "idw", k = 5, 
    style = "W"
)

sfe <- colDataMoransI(sfe, c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                      colGraphName = "knn5")

sfe <- colDataUnivariate(sfe, type = "localmoran", 
                         features = c("nCounts", "nGenes", "cell_area", 
                                      "nucleus_area"),
                         colGraphName = "knn5", BPPARAM = MulticoreParam(2))

p <- plotLocalResult(sfe, "localmoran",
                features = c("nCounts", "nGenes", "cell_area", "nucleus_area"),
                colGeometryName = "centroids", scattermore = TRUE,
                divergent = TRUE, diverge_center = 0, pointsize = 1)


pdf(paste0(fig_dir, 'qc_autocorrelation.pdf'), width=10, height=5)
print(p) 
dev.off()

sfe <- colDataUnivariate(sfe, "moran.plot", "nCounts", colGraphName = "knn5")

p1 <- moranPlot(sfe, "nCounts", binned = TRUE, plot_influential = FALSE) 
p2 <- moranPlot(sfe, "nCounts", binned = TRUE)

patch <- p1 / p2 + plot_layout(guides = "collect")

pdf(paste0(fig_dir, 'qc_autocorrelation_lineplot.pdf'), width=12, height=6)
print(p)
dev.off()


#-------------------------------------------------------------------------
# Moran's I Gene Expression
#-------------------------------------------------------------------------

# remove the negative controls
sfe <- sfe[rowData(sfe)$Type == 'Gene Expression',]

# log normalize the counts
sfe <- logNormCounts(sfe)

sfe <- runMoransI(sfe, colGraphName = "knn5", BPPARAM = MulticoreParam(8))

p <- plotRowData(sfe, x = "means", y = "moran_sample01") + 
    scale_x_log10() 


pdf(paste0(fig_dir, 'qc_mean_vs_moranI.pdf'), width=6, height=6)
print(p)
dev.off()


#-------------------------------------------------------------------------
# Dim Reduction and clustering
#-------------------------------------------------------------------------

set.seed(2024)
sfe <- runPCA(sfe, ncomponents = 30, scale = TRUE, BSPARAM = IrlbaParam())

# leiden clustering
colData(sfe)$cluster <- clusterRows(
    reducedDim(sfe, "PCA")[,1:15],
    BLUSPARAM = SNNGraphParam(
        cluster.fun = "leiden",
        cluster.args = list(
            resolution_parameter = opt$cluster_res,
            objective_function = "modularity"
        )
    )
)

# UMAP 
sfe <- runUMAP(sfe, dimred='PCA')


# plot the clusters in space
p <- plotSpatialFeature(sfe, "cluster", colGeometryName = "cellSeg") + facet_wrap(~cluster)
pdf(paste0(fig_dir, 'leiden_spatial_split.pdf'), width=8, height=6)
print(p)
dev.off()

p <- plotSpatialFeature(sfe, "cluster", colGeometryName = "cellSeg", dark=TRUE,
    image_id='morphology_focus',
    channel = 2,
    maxcell=10e+06) 
pdf(paste0(fig_dir, 'leiden_spatial_dark.pdf'), width=8, height=6)
print(p)
dev.off()

# plot the clusters in space
p <- plotSpatialFeature(sfe, "cluster", colGeometryName = "cellSeg", dark=TRUE) + facet_wrap(~cluster)
pdf(paste0(fig_dir, 'leiden_spatial_dark_split.pdf'), width=8, height=6)
print(p)
dev.off()

p <- plotSpatialFeature(sfe, "cluster", colGeometryName = "centroids") 
pdf(paste0(fig_dir, 'leiden_spatial_centroids.pdf'), width=8, height=6)
print(p)
dev.off()

p <- plotSpatialFeature(sfe, "cluster", colGeometryName = "cellSeg") 
pdf(paste0(fig_dir, 'leiden_spatial.pdf'), width=8, height=6)
print(p)
dev.off()

# get the color scheme
g <- ggplot_build(p)
color_df <- g$data[[1]] %>% dplyr::select(c('group', 'fill')) %>% dplyr::distinct()
cp <- color_df$fill 
names(cp) <- color_df$group

# plot the UMAP colored by cluster
p <- plotUMAP(sfe, colour_by = 'cluster') + scale_color_manual(values=cp)
pdf(paste0(fig_dir, 'leiden_umap.pdf'), width=6, height=6)
print(p)
dev.off()

p <- plotColData(sfe, x = "cluster", y = "nCounts")
pdf(paste0(fig_dir, 'leiden_nCounts.pdf'), width=6, height=6)
print(p)
dev.off()

#-------------------------------------------------------------------------
# Cluster marker genes (Seurat FindMarkers)
#-------------------------------------------------------------------------

X <- as(logcounts(sfe), "dgCMatrix")
rownames(X) <- rownames(sfe); colnames(X) <- colnames(sfe)
X_assay <- SeuratObject::CreateAssayObject(X)

markers <- data.frame()
groups <- as.character(unique(colData(sfe)$cluster))
for(cur_group in groups){

    print(cur_group)

    b1 <- colnames(sfe)[colData(sfe)$cluster == cur_group]
    b2 <- colnames(sfe)[colData(sfe)$cluster != cur_group] 

    # run differential test on the ME assay
    cur_markers <- FindMarkers(
      X_assay,
      cells.1 = b1,
      cells.2 = b2,
      slot='counts',
      test.use='wilcox',
      only.pos=FALSE,
      logfc.threshold=0,
      min.pct=0.0
    )
    cur_markers$gene <- rownames(cur_markers)
    cur_markers$group <- cur_group

    # add annotation info
    cur_markers$Annotation <- rowData(sfe)[cur_markers$gene, 'Annotation']
    
    markers <- rbind(markers, cur_markers)

}

# write the output table
write.table(
    markers, 
    file=paste0(data_dir, cur_sample, '_leiden_',opt$cluster_res ,'_markers.tsv'), 
    sep='\t', quote=FALSE, row.names=FALSE
)

# get list of genes to plot
genes_plot <- markers %>% subset(avg_log2FC > 2) %>% .$gene %>% unique()
genes_plot_groups <- lapply(groups, function(x){subset(markers, group == x & avg_log2FC > 2) %>% .$gene})
names(genes_plot_groups) <- groups
gene_list_lengths <- as.numeric(unlist(lapply(genes_plot_groups, function(x){length(x)})))

# set cluster factor level
markers$group <- factor(as.character(markers$group), levels=groups)
plot_df <- subset(markers, gene %in% genes_plot) %>% 
    arrange(group, -avg_log2FC)

# limit fill plotting range
fc_max <- 3
plot_df$avg_log2FC <- ifelse(abs(plot_df$avg_log2FC) > fc_max, sign(plot_df$avg_log2FC)*fc_max, plot_df$avg_log2FC )

# set gene factor level
plot_df$gene <- factor(as.character(plot_df$gene), levels=rev(genes_plot))

# create gene + annotation label
plot_df$label <- paste0(
    as.character(plot_df$gene), ', ', as.character(plot_df$Annotation)
)

# order the labels
labels_plot <- paste0(genes_plot, ', ', panel_meta[genes_plot, 'Annotation'])
plot_df$label <- factor(as.character(plot_df$label), levels=rev(labels_plot))

# make a separate heatmap for each cluster
plot_list <- list()
for(cur_group in groups){
    cur_genes <- genes_plot_groups[[cur_group]]
    p <- plot_df %>% 
        subset(gene %in% cur_genes) %>% 
         ggplot(aes(x=group, y=gene, fill=avg_log2FC)) + 
        geom_tile() + 
        scale_fill_gradient2(high='yellow', low='purple', mid='black') + 
        labs(fill='log2FC') +
        xlab('Cluster') + ylab('') + theme(
            axis.text.y = element_text(size=3),
            plot.margin = margin(0,0,0,0)
        ) 

    if(cur_group != groups[[length(groups)]]){
        p <- p + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank()
        )
    }
    plot_list[[cur_group]] <- p

}

# plot with patchwork
p <- wrap_plots(plot_list, ncol=1) + plot_layout(guides='collect', heights=gene_list_lengths)
pdf(paste0(fig_dir, 'leiden_marker_heatmap.pdf'), width=4, height=10)
p
dev.off()



# make a separate heatmap for each cluster
plot_list <- list()
for(cur_group in groups){
    cur_genes <- genes_plot_groups[[cur_group]]
    p <- plot_df %>% 
        subset(gene %in% cur_genes) %>% 
         ggplot(aes(x=group, y=label, fill=avg_log2FC)) + 
        geom_tile() + 
        scale_fill_gradient2(high='yellow', low='purple', mid='black') + 
        labs(fill='log2FC') +
        xlab('Cluster') + ylab('') + theme(
            axis.text.y = element_text(size=3),
            plot.margin = margin(0,0,0,0)
        ) 

    if(cur_group != groups[[length(groups)]]){
        p <- p + theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.line.x = element_blank()
        )
    }
    plot_list[[cur_group]] <- p

}

# plot with patchwork
p <- wrap_plots(plot_list, ncol=1) + plot_layout(guides='collect', heights=gene_list_lengths)
pdf(paste0(fig_dir, 'leiden_marker_heatmap_annotation.pdf'), width=10, height=10)
p
dev.off()

#-------------------------------------------------------------------------
# save the SFE object
#-------------------------------------------------------------------------
saveRDS(sfe, file= paste0(data_dir, 'voyager_sfe_', cur_sample, '.rds'))

#-------------------------------------------------------------------------
# Save the meta-data and the CxG matrix 
# will help with the downstream integrated analysis
#-------------------------------------------------------------------------

meta_df <- colData(sfe)


# write the output table
write.table(
    meta_df, 
    file=paste0(data_dir, cur_sample, '_metadata.tsv'), 
    sep='\t', quote=FALSE, row.names=FALSE
)


X <- as(counts(sfe), "dgCMatrix")
writeMM(X, file=paste0(data_dir, cur_sample, '_counts.mtx'))
