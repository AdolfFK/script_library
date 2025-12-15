#!/usr/bin/env Rscript
#==============================================================================
# 脚本名: run_MultiSamples.V5_gpt_v3.R
# 功能:  单细胞的第一步脚本，处理cellranger的结果文件，可以做分组。结果最好还是人工审核一下。把图上的点变的更小
# 作者:  adolf1
# 邮箱:  wangjyafk@126.com
# 创建:  2025-06-01
# 更新:  2025-12-15 adolf1
# 依赖:  R≥4.0 data.table ggplot2
# 用法:  Rscript run_MultiSamples.V5_gpt_v3.R --参数
# 仓库:  https://github.com/AdolfFK/script_library
#==============================================================================


######################################################################################
#
# Fixed the issue where points turned white after processing more than 100,000 cells.
#
######################################################################################

# Suppress package loading messages
suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(argparser)
  library(Seurat)
  library(corrplot)
  library(grid)
  library(cowplot)
  library(viridis)
  library(Matrix)
  library(rlang) # Added for check_installed
})

# Argument parsing
argv <- arg_parser('debug:hulongfei@singleronbio.com!')
argv <- add_argument(argv, "--expM", help="Cell gene expression file list, split by ,")
argv <- add_argument(argv, "--datatype", help="Cell gene expression file type: 1=EXPM, 0=10XV2, 2=10XV3, split by ,", default="1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1")
argv <- add_argument(argv, "--ugname", help="Sample group names, split by ,")
argv <- add_argument(argv, "--spname", help="Sample names, split by ,")
argv <- add_argument(argv, "--sptype", help="Sample types, split by ,", default="F")
argv <- add_argument(argv, "--prefix", help="Group name prefix")
argv <- add_argument(argv, "--group", help="Compare group names: G1vsG2,G3vsG4, split by ,", default="F")
argv <- add_argument(argv, "--tsne_method", help="tSNE method: Rtsne or FIt-SNE", default="Rtsne")
argv <- add_argument(argv, "--dim", help="Number of clusters (1:20)", default=20)
argv <- add_argument(argv, "--mtfilter", help="Max mitochondrial gene percentage (%)", default=20)
argv <- add_argument(argv, "--ngenefiltermin", help="Minimum number of genes", default=200)
argv <- add_argument(argv, "--ngenefiltermax", help="Max number of genes; if >1, use cutoff; if <=1, use quantile", default=5000)
argv <- add_argument(argv, "--umifiltermin", help="Minimum UMI count", default=0)
argv <- add_argument(argv, "--umifiltermax", help="Max UMI count; if >1, use cutoff; if <=1, use quantile", default=30000)
argv <- add_argument(argv, "--resolution", help="Clustering resolution", default=1.2)
argv <- add_argument(argv, "--logfcthreshold", help="Min log-fold change for gene testing", default=0.25)
argv <- add_argument(argv, "--minpct", help="Min fraction of cells expressing genes", default=0.1)
argv <- add_argument(argv, "--outdir", help="Output directory")
argv <- add_argument(argv, "--object", help="Analysis objects: 1=matrix/10X and basic analysis, 2=diffgene, 3=markersgene plot, 4=average expression, 5=trajectories, 6=compare group, split by ,", default="1,1,0,0,0,0")
argv <- parse_args(argv)

# Parse and validate inputs
print(argv$expM)
expfile <- unlist(strsplit(argv$expM, split=","))
name <- unlist(strsplit(argv$spname, split=","))
groupname <- unlist(strsplit(argv$ugname, split=","))
sampletype <- unlist(strsplit(argv$sptype, split=","))
gnumber <- unique(groupname)
group <- unlist(strsplit(argv$group, split=","))
object <- as.numeric(unlist(strsplit(argv$object, split=",")))
datatype <- as.numeric(unlist(strsplit(argv$datatype, split=",")))
cnum <- as.numeric(argv$dim)
tsne_method <- argv$tsne_method
resol <- as.numeric(argv$resolution)
outdir <- argv$outdir
dir.create(outdir, showWarnings=FALSE)
compare <- argv$prefix
mtfilter <- as.numeric(argv$mtfilter)
ngenefiltermin <- as.numeric(argv$ngenefiltermin)
ngenefiltermax <- as.numeric(argv$ngenefiltermax)
umifiltermin <- as.numeric(argv$umifiltermin)
umifiltermax <- as.numeric(argv$umifiltermax)
logfc <- as.numeric(argv$logfcthreshold)
minpct_u <- as.numeric(argv$minpct)
print(compare)
print(object)
print(groupname)
print(gnumber)

# Input validation
if (length(expfile) != length(name) || length(expfile) != length(groupname) || (argv$sptype != "F" && length(sampletype) != length(expfile))) {
  quit()
}

# 增强的颜色定义 - 使用更鲜艳的颜色
clustcol <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
              "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
              "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
              "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA",
              "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
              "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1B9E77",
              "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
              "#A6761D", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A",
              "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
              "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#8DD3C7")

col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "cyan", "#007FFF", "blue"))
corrcol <- colorRampPalette(c("red", "orange", "blue", "white", "white"))
col2 <- colorRampPalette(clustcol)

# 自定义绘图函数 - 禁用点阵化并增强颜色
enhanced_dim_plot <- function(seurat_obj, reduction_type = "umap", pt_size = 0.1, 
                             group_by = "ident", label = FALSE, cols = clustcol, 
                             label_size = 4) {
  DimPlot(seurat_obj, 
          reduction = reduction_type,
          pt.size = pt_size,
          group.by = group_by,
          label = label,
          label.size = label_size,
          cols = cols,
          raster = FALSE)  # 关键：禁用点阵化
}

enhanced_feature_plot <- function(seurat_obj, features, reduction_type = "umap", 
                                 pt_size = 0.1, cols = c("lightgrey", "#FF0000")) {
  FeaturePlot(seurat_obj,
              features = features,
              reduction = reduction_type,
              pt.size = pt_size,
              cols = cols,
              min.cutoff = "q5",   # 调整分位数范围
              max.cutoff = "q95",
              order = TRUE,        # 高表达细胞在上层
              raster = FALSE)      # 关键：禁用点阵化
}

enhanced_tsne_plot <- function(seurat_obj, pt_size = 0.1, group_by = "ident", 
                              label = FALSE, cols = clustcol, label_size = 4) {
  DimPlot(seurat_obj,
          reduction = "tsne",
          pt.size = pt_size,
          group.by = group_by,
          label = label,
          label.size = label_size,
          cols = cols,
          raster = FALSE)  # 关键：禁用点阵化
}

# Read and process data
read_data <- function(file, type) {
  if (type == 1) {
    data <- read.table(file, sep="\t")
    return(as(data, "dgCMatrix")) # Explicitly convert to dgCMatrix
  } else if (type == 0) {
    return(Read10X(data.dir=file))
  } else if (type == 2) {
    matrix_dir <- file
    barcode.path <- paste0(matrix_dir, "/barcodes.tsv.gz")
    features.path <- paste0(matrix_dir, "/features.tsv.gz")
    matrix.path <- paste0(matrix_dir, "/matrix.mtx.gz")
    RAW <- readMM(file=matrix.path)
    feature.names <- read.delim(features.path, header=FALSE, stringsAsFactors=FALSE)
    barcode.names <- read.delim(barcode.path, header=FALSE, stringsAsFactors=FALSE)
    colnames(RAW) <- barcode.names$V1
    rownames(RAW) <- make.unique(feature.names$V2, sep='-') # Use dash instead of underscore
    rm(feature.names, barcode.names)
    return(as(RAW, "dgCMatrix")) # Explicitly convert to dgCMatrix
  }
}

# Process first sample
RAW1 <- read_data(expfile[1], datatype[1])
sampproj <- CreateSeuratObject(counts=RAW1, project=name[1], min.cells=5)
sampproj <- RenameCells(sampproj, add.cell.id=name[1])
sampproj$sample <- groupname[1]
if (argv$sptype != "F") {
  sampproj$sampletype <- sampletype[1]
}
mt.genes <- rownames(sampproj)[grep("^MT-", rownames(sampproj), ignore.case=TRUE)]
sampproj[["percent.mt"]] <- PercentageFeatureSet(sampproj, features=mt.genes)
ngenemax <- if (ngenefiltermax <= 1) unname(quantile(sampproj$nFeature_RNA, ngenefiltermax)) else ngenefiltermax
umimax <- if (umifiltermax <= 1) unname(quantile(sampproj$nCount_RNA, umifiltermax)) else umifiltermax
print(ngenemax)
print(umimax)
sampproj <- subset(sampproj, subset=nFeature_RNA > ngenefiltermin & nFeature_RNA < ngenemax & percent.mt < mtfilter & nCount_RNA > umifiltermin & nCount_RNA < umimax)
sampproj <- NormalizeData(object=sampproj)

# Process additional samples
if (length(expfile) > 1) {
  for (i in 2:length(expfile)) {
    RAW1 <- read_data(expfile[i], datatype[i])
    sampproj_cc <- CreateSeuratObject(counts=RAW1, project=name[i], min.cells=5)
    rm(RAW1)
    sampproj_cc <- RenameCells(sampproj_cc, add.cell.id=name[i])
    sampproj_cc$sample <- groupname[i]
    if (argv$sptype != "F") {
      sampproj_cc$sampletype <- sampletype[i]
    }
    mt.genes <- rownames(sampproj_cc)[grep("^MT-", rownames(sampproj_cc), ignore.case=TRUE)]
    sampproj_cc[["percent.mt"]] <- PercentageFeatureSet(sampproj_cc, features=mt.genes)
    ngenemax <- if (ngenefiltermax <= 1) unname(quantile(sampproj_cc$nFeature_RNA, ngenefiltermax)) else ngenefiltermax
    umimax <- if (umifiltermax <= 1) unname(quantile(sampproj_cc$nCount_RNA, umifiltermax)) else umifiltermax
    print(ngenemax)
    print(umimax)
    sampproj_cc <- subset(sampproj_cc, subset=nFeature_RNA > ngenefiltermin & nFeature_RNA < ngenemax & percent.mt < mtfilter & nCount_RNA > umifiltermin & nCount_RNA < umimax)
    sampproj_cc <- NormalizeData(object=sampproj_cc)
    sampproj <- merge(sampproj, y=sampproj_cc, project="PRO")
    rm(sampproj_cc)
  }
}

# Main analysis
PRO <- FindVariableFeatures(sampproj, selection.method="vst", nfeatures=2000)
rm(sampproj)
top10 <- head(VariableFeatures(PRO), 10)
all.genes <- rownames(PRO)
PRO <- ScaleData(PRO, features=all.genes)
PRO <- RunPCA(PRO, features=VariableFeatures(object=PRO))

# PCA diagnostics
PRO <- JackStraw(PRO, num.replicate=100)
PRO <- ScoreJackStraw(PRO, dims=1:20)
pdf(paste0(outdir, '/', compare, '.PCJackStrawPlot.pdf'))
JackStrawPlot(PRO, dims=1:20)
dev.off()
pdf(paste0(outdir, '/', compare, '.PCElbowPlot.pdf'))
ElbowPlot(PRO)
dev.off()

# Dimensionality reduction and clustering
PRO <- RunUMAP(PRO, reduction="pca", dims=1:cnum)
PRO <- FindNeighbors(PRO, reduction="pca", dims=1:cnum)
PRO <- FindClusters(PRO, resolution=resol)
if (tsne_method == "Rtsne") {
  PRO <- RunTSNE(PRO, reduction="pca", dims=1:cnum, check_duplicates=FALSE)
} else {
  print(tsne_method)
  PRO <- RunTSNE(PRO, reduction="pca", dims.use=1:cnum, tsne.method=tsne_method, nthreads=2, max_iter=2000, check_duplicates=FALSE)
}

# Adjust point size based on cell number - 优化点大小设置
clust <- summary(PRO@active.ident)
cluster_cell <- as.data.frame(clust)
cell_number <- sum(cluster_cell$clust)
print(cell_number)
pt_use <- 0.6
if (cell_number > 1000) pt_use <- 0.4
if (cell_number > 2500) pt_use <- 0.3
if (cell_number > 4000) pt_use <- 0.2
if (cell_number > 5500) pt_use <- 0.15
if (cell_number > 6500) pt_use <- 0.1
if (cell_number > 100000) pt_use <- 0.08
if (cell_number > 500000) pt_use <- 0.04
if (cell_number > 800000) pt_use <- 0.02

# Visualization with NA check - 使用增强的绘图函数
p1 <- enhanced_tsne_plot(PRO, pt_size=pt_use, group_by="sample", cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="tsne")))) {
  warning("NA values detected in tSNE embeddings, skipping tSNE plot")
} else {
  pdf(paste0(outdir, '/', compare, '.tSNE_samples.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p1)
  dev.off()
  png(paste0(outdir, '/', compare, '.tSNE_samples.png'), width=2000, height=1600, res=300)
  print(p1)
  dev.off()
}

p2 <- enhanced_dim_plot(PRO, reduction_type="umap", pt_size=pt_use, group_by="sample", cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="umap")))) {
  warning("NA values detected in UMAP embeddings, skipping UMAP plot")
} else {
  pdf(paste0(outdir, '/', compare, '.UMAP_samples.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p2)
  dev.off()
  png(paste0(outdir, '/', compare, '.UMAP_samples.png'), width=2000, height=1600, res=300)
  print(p2)
  dev.off()
}

# Save dimensional reduction data
pca_data <- Embeddings(object=PRO, reduction="pca")[,1:5]
tsne_data <- Embeddings(object=PRO, reduction="tsne")[,1:2]
umap_data <- Embeddings(object=PRO, reduction="umap")[,1:2]
write.table(data.frame(cellID=rownames(pca_data), pca_data), file=paste0(outdir, '/', compare, '.pca_file.xls'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(data.frame(cellID=rownames(tsne_data), tsne_data), file=paste0(outdir, '/', compare, '.tsne_file.xls'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(data.frame(cellID=rownames(umap_data), umap_data), file=paste0(outdir, '/', compare, '.umap_file.xls'), sep='\t', quote=FALSE, row.names=FALSE)

# Cluster analysis
ident <- as.numeric(levels(x=PRO))
new <- ident + 1
names(new) <- levels(PRO)
PRO <- RenameIdents(PRO, new)
clust <- summary(Idents(object=PRO))
cluster_num <- table(PRO@active.ident, PRO@meta.data$sample)
ident <- rownames(cluster_num)
freq_table <- prop.table(x=table(PRO@active.ident, PRO@meta.data$sample), margin=2)
labels <- paste('cluster', rownames(freq_table), sep='')
rownames(freq_table) <- labels
rownames(cluster_num) <- labels
write.table(cluster_num, file=paste0(outdir, '/', compare, '.CellsPerCluster.xls'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
write.table(freq_table, file=paste0(outdir, '/', compare, '.PercentPerCluster.xls'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)

# Barplot visualization
mix <- 30 / length(gnumber)
if (mix < 1.5) mix <- 1.5
print(mix)
pdf(paste0(outdir, '/', compare, '.PercentPerCell.pdf'), width=10, height=7)
barplot(height=freq_table, width=mix, xlim=c(1,60), col=clustcol, legend=rownames(freq_table), args.legend=list(x="right", ncol=3), las=2, xlab="")
mtext(text="Samples", side=1, line=6.5)
dev.off()
png(paste0(outdir, '/', compare, '.PercentPerCell.png'), width=1024, height=640)
barplot(height=freq_table, width=mix, xlim=c(1,60), col=clustcol, legend=rownames(freq_table), args.legend=list(x="right", ncol=3), las=2, xlab="")
mtext(text="Samples", side=1, line=6.5)
dev.off()

# Save cluster information
write.table(data.frame(CellBarcode=rownames(PRO@meta.data), PRO@meta.data), file=paste0(outdir, '/', compare, '.cell_cluster_info.xls'), sep='\t', quote=FALSE, row.names=FALSE)

# tSNE and UMAP visualization with labels
subc <- levels(PRO)
print(subc)
cell.label.size <- 4.5
if (length(subc) > 10) cell.label.size <- 4.2
if (length(subc) > 25) cell.label.size <- 4
if (length(subc) > 30) cell.label.size <- 3.5

p3 <- enhanced_tsne_plot(PRO, pt_size=pt_use, cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="tsne")))) {
  warning("NA values detected in tSNE embeddings, skipping tSNE plot")
} else {
  pdf(paste0(outdir, '/', compare, '.tsne.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p3)
  dev.off()
  png(paste0(outdir, '/', compare, '.tsne.png'), width=2000, height=1600, res=300)
  print(p3)
  dev.off()
}

p5 <- enhanced_tsne_plot(PRO, pt_size=pt_use, cols=clustcol, label=TRUE, label_size=cell.label.size) + NoLegend()
if (any(is.na(Embeddings(PRO, reduction="tsne")))) {
  warning("NA values detected in tSNE embeddings, skipping labeled tSNE plot")
} else {
  pdf(paste0(outdir, '/', compare, '.labtsne.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p5)
  dev.off()
  png(paste0(outdir, '/', compare, '.labtsne.png'), width=2000, height=1600, res=300)
  print(p5)
  dev.off()
}

p4 <- enhanced_dim_plot(PRO, reduction_type="umap", pt_size=pt_use, cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="umap")))) {
  warning("NA values detected in UMAP embeddings, skipping UMAP plot")
} else {
  pdf(paste0(outdir, '/', compare, '.umap.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p4)
  dev.off()
  png(paste0(outdir, '/', compare, '.umap.png'), width=2000, height=1600, res=300)
  print(p4)
  dev.off()
}

p5 <- enhanced_dim_plot(PRO, reduction_type="umap", pt_size=pt_use, cols=clustcol, label=TRUE, label_size=cell.label.size) + NoLegend()
if (any(is.na(Embeddings(PRO, reduction="umap")))) {
  warning("NA values detected in UMAP embeddings, skipping labeled UMAP plot")
} else {
  pdf(paste0(outdir, '/', compare, '.labumap.pdf'), width=10, height=8, useDingbats=FALSE)
  print(p5)
  dev.off()
  png(paste0(outdir, '/', compare, '.labumap.png'), width=2000, height=1600, res=300)
  print(p5)
  dev.off()
}

# Save Seurat object
PRO <- DietSeurat(PRO, layers=c("counts", "data"), dimreducs=c('pca', 'tsne', 'umap'))
saveRDS(PRO, file=paste0(outdir, '/', compare, '.diff_PRO.rds'))

# Sample type visualization
if (argv$sptype != "F") {
  p6 <- enhanced_dim_plot(PRO, reduction_type="umap", pt_size=pt_use, group_by="sampletype", cols=clustcol)
  if (any(is.na(Embeddings(PRO, reduction="umap")))) {
    warning("NA values detected in UMAP embeddings, skipping sample type UMAP plot")
  } else {
    pdf(paste0(outdir, '/', compare, '.umap.sampletype.pdf'), width=10, height=8, useDingbats=FALSE)
    print(p6)
    dev.off()
    png(paste0(outdir, '/', compare, '.umap.sampletype.png'), width=2000, height=1600, res=300)
    print(p6)
    dev.off()
  }
}

# Combined reduction plots
if (all(!is.na(Embeddings(PRO, reduction="tsne"))) && all(!is.na(Embeddings(PRO, reduction="umap")))) {
  pdf(paste0(outdir, '/', compare, '.reduction.pdf'), width=12, height=10, useDingbats=FALSE)
  plot_grid(p1, p2, p3, p4, ncol=2)
  dev.off()
  png(paste0(outdir, '/', compare, '.reduction.png'), width=2400, height=2000, res=300)
  plot_grid(p1, p2, p3, p4, ncol=2)
  dev.off()
}

# Split plots
p_s1 <- enhanced_dim_plot(PRO, reduction_type="umap", pt_size=pt_use, group_by="sample", cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="umap")))) {
  warning("NA values detected in UMAP embeddings, skipping sample split UMAP plot")
} else {
  pdf(paste0(outdir, '/', compare, '.sample_split_umap.pdf'), width=12, height=8, useDingbats=FALSE)
  print(p_s1)
  dev.off()
  png(paste0(outdir, '/', compare, '.sample_split_umap.png'), width=2400, height=1600, res=300)
  print(p_s1)
  dev.off()
}

p_s1 <- enhanced_tsne_plot(PRO, pt_size=pt_use, group_by="sample", cols=clustcol)
if (any(is.na(Embeddings(PRO, reduction="tsne")))) {
  warning("NA values detected in tSNE embeddings, skipping sample split tSNE plot")
} else {
  pdf(paste0(outdir, '/', compare, '.sample_split_tsne.pdf'), width=12, height=8, useDingbats=FALSE)
  print(p_s1)
  dev.off()
  png(paste0(outdir, '/', compare, '.sample_split_tsne.png'), width=2400, height=1600, res=300)
  print(p_s1)
  dev.off()
}

# Differential expression analysis
DefaultAssay(PRO) <- "RNA"
print(ident)
ident <- as.numeric(levels(x=PRO))
PRO <- JoinLayers(PRO)
for (l in ident) {
  print(l)
  if (object[2]) {
    cluster.markers <- FindMarkers(object=PRO, ident.1=l, min.pct=minpct_u, logfc.threshold=logfc)
    write.table(data.frame(gene_id=rownames(cluster.markers), cluster.markers), file=paste0(outdir, '/', compare, '.cluster', l, '_diffgenes.xls'), sep='\t', quote=FALSE, row.names=FALSE)
  }
}

# Average expression analysis
if (object[4]) {
  cluster.averages_rep <- AverageExpression(object=PRO, layer="data", add.ident="sample")
  write.table(cluster.averages_rep[["RNA"]], file=paste0(outdir, '/', compare, '.cluster_averages_repeat.xls'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  cluster.averages1 <- AverageExpression(object=PRO, layer="data")
  write.table(cluster.averages1[["RNA"]], file=paste0(outdir, '/', compare, '.cluster_averages.xls'), sep='\t', quote=FALSE, row.names=TRUE, col.names=NA)
  M <- cor(cluster.averages1[["RNA"]])
  order.hc <- corrMatOrder(M, order="hclust")
  M.hc <- M[order.hc, order.hc]
  pdf(paste0(outdir, '/', compare, '.corrplot.pdf'))
  corrplot(M.hc, method="number", cl.lim=c(0,1), tl.col="black", col=rev(corrcol(50)))
  dev.off()
  png(paste0(outdir, '/', compare, '.corrplot.png'))
  corrplot(M.hc, method="number", cl.lim=c(0,1), tl.col="black", col=rev(corrcol(50)))
  dev.off()
}

# Marker gene analysis
if (object[3]) {
  mgp <- paste0(outdir, '/markergeneplot')
  dir.create(mgp, showWarnings=FALSE)
  PRO.marker <- FindAllMarkers(PRO, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
  markergene <- PRO.marker %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)
  markergenetop10 <- PRO.marker %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
  print(markergene$gene)
  allmarker <- as.data.frame(markergene)
  allmarkertop10 <- as.data.frame(markergenetop10)
  marker <- allmarker$gene
  heatname <- paste(compare, 'top 10 marker genes heatmap', sep=' ')
  cex_row <- 8
  cex_col <- 9
  if (length(markergene$gene) > 30) { cex_row <- 6; cex_col <- 8 }
  if (length(markergene$gene) > 60) { cex_row <- 5; cex_col <- 7 }
  if (length(markergene$gene) > 80) { cex_row <- 4; cex_col <- 6 }
  if (length(markergene$gene) > 100) { cex_row <- 3.5; cex_col <- 5 }
  print(marker)
  all.genes <- rownames(PRO)
  PRO_sub <- subset(PRO, downsample=300)
  all.genes <- rownames(PRO_sub)
  PRO_sub <- ScaleData(PRO_sub, features=all.genes)
  pdf(paste0(mgp, '/', compare, '.top10markergene_DOHeatmapplot.pdf'))
  print(DoHeatmap(PRO_sub, features=markergene$gene, size=2, hjust=0, angle=0, layer="data", draw.lines=TRUE, raster=FALSE, group.by="ident") + scale_fill_viridis() + theme(text=element_text(size=3)) + guides(color=FALSE))
  dev.off()
  png(paste0(mgp, '/', compare, '.top10markergene_DOHeatmapplot.png'))
  print(DoHeatmap(PRO_sub, features=markergene$gene, size=2, hjust=0, angle=0, layer="data", draw.lines=TRUE, raster=FALSE, group.by="ident") + scale_fill_viridis() + theme(text=element_text(size=3)) + guides(color=FALSE))
  dev.off()
  write.table(allmarker, file=paste0(mgp, '/', compare, '.ALLmarkergenelist.xls'), sep='\t', quote=FALSE, row.names=FALSE)
  write.table(allmarkertop10, file=paste0(outdir, '/', compare, '.ALLmarkerTop10genelist.xls'), sep='\t', quote=FALSE, row.names=FALSE)
  marker <- unique(marker)
  g <- 0
  pnum <- 1
  plot_g <- c()
  x_size <- 9
  title_size <- 11
  y_size <- 5
  for (i in marker) {
    if (g == 4) {
      pdf(paste0(mgp, '/', compare, '.markergene_tsneFeatureplot', pnum, '.pdf'), width=10, height=8, useDingbats=FALSE)
      print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="tsne", pt_size=pt_use))
      dev.off()
      png(paste0(mgp, '/', compare, '.markergene_tsneFeatureplot', pnum, '.png'), width=2000, height=1600, res=300)
      print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="tsne", pt_size=pt_use))
      dev.off()
      pdf(paste0(mgp, '/', compare, '.markergene_umapFeatureplot', pnum, '.pdf'), width=10, height=8, useDingbats=FALSE)
      print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="umap", pt_size=pt_use))
      dev.off()
      png(paste0(mgp, '/', compare, '.markergene_umapFeatureplot', pnum, '.png'), width=2000, height=1600, res=300)
      print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="umap", pt_size=pt_use))
      dev.off()
      pdf(paste0(mgp, '/', compare, '.markergene_vlnplot', pnum, '.pdf'), width=10, height=8)
      print(VlnPlot(PRO, features=plot_g, pt.size=0.01, cols=clustcol, ncol=2))
      dev.off()
      png(paste0(mgp, '/', compare, '.markergene_vlnplot', pnum, '.png'), width=2000, height=1600, res=300)
      print(VlnPlot(PRO, features=plot_g, pt.size=0.01, cols=clustcol, ncol=2))
      dev.off()
      pnum <- pnum + 1
      g <- 1
      plot_g <- c(i)
    } else {
      plot_g <- c(plot_g, i)
      g <- g + 1
    }
  }
  if (g > 0) {
    pdf(paste0(mgp, '/', compare, '.markergene_tsneFeatureplot', pnum, '.pdf'), width=10, height=8, useDingbats=FALSE)
    print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="tsne", pt_size=pt_use))
    dev.off()
    png(paste0(mgp, '/', compare, '.markergene_tsneFeatureplot', pnum, '.png'), width=2000, height=1600, res=300)
    print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="tsne", pt_size=pt_use))
    dev.off()
    pdf(paste0(mgp, '/', compare, '.markergene_umapFeatureplot', pnum, '.pdf'), width=10, height=8, useDingbats=FALSE)
    print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="umap", pt_size=pt_use))
    dev.off()
    png(paste0(mgp, '/', compare, '.markergene_umapFeatureplot', pnum, '.png'), width=2000, height=1600, res=300)
    print(enhanced_feature_plot(PRO, features=plot_g, reduction_type="umap", pt_size=pt_use))
    dev.off()
    pdf(paste0(mgp, '/', compare, '.markergene_vlnplot', pnum, '.pdf'), width=10, height=8)
    print(VlnPlot(PRO, features=plot_g, pt.size=0.01, cols=clustcol, ncol=2))
    dev.off()
    png(paste0(mgp, '/', compare, '.markergene_vlnplot', pnum, '.png'), width=2000, height=1600, res=300)
    print(VlnPlot(PRO, features=plot_g, pt.size=0.01, cols=clustcol, ncol=2))
    dev.off()
  }
  dot_size <- 5
  if (length(marker) > 10) dot_size <- 4.5
  if (length(marker) > 15) dot_size <- 4
  if (length(marker) > 20) dot_size <- 3.5
  if (length(marker) > 25) dot_size <- 3
  if (length(marker) > 30) dot_size <- 2.5
  pdf(paste0(mgp, '/', compare, '.markergene_Dotplot.pdf'), width=10, height=8)
  print(DotPlot(PRO, features=marker, cols=c("blue", "red"), dot.scale=dot_size) + RotatedAxis())
  dev.off()
  png(paste0(mgp, '/', compare, '.markergene_Dotplot.png'), width=2000, height=1600, res=300)
  print(DotPlot(PRO, features=marker, cols=c("blue", "red"), dot.scale=dot_size) + RotatedAxis())
  dev.off()
  sdp <- DotPlot(PRO, features=marker, cols=clustcol, dot.scale=dot_size, split.by="sample") + RotatedAxis()
  pdf(paste0(mgp, '/', compare, '.SplitDotPlotGG.pdf'), width=12, height=10)
  print(sdp)
  dev.off()
  png(paste0(mgp, '/', compare, '.SplitDotPlotGG.png'), width=2400, height=2000, res=300)
  print(sdp)
  dev.off()
}

# Differential expression between groups
if (object[6]) {
  PRO$celltype.sample <- paste(Idents(PRO), PRO$sample, sep="_")
  PRO$celltype <- Idents(PRO)
  Idents(PRO) <- "celltype.sample"
  PRO <- JoinLayers(PRO)
  for (cm in group) {
    comparegroup <- unlist(strsplit(cm, split="vs"))
    for (l in ident) {
      s1 <- paste(l, '_', comparegroup[1], sep='')
      s2 <- paste(l, '_', comparegroup[2], sep='')
      b.interferon.response <- FindMarkers(PRO, ident.1=s1, ident.2=s2, verbose=FALSE)
      write.table(b.interferon.response, file=paste0(outdir, '/', compare, '_', cm, '.cluster', l, '_diffexpressed.xls'), sep='\t', quote=FALSE, row.names=TRUE)
      g_u <- head(rownames(b.interferon.response), 4)
      pdf(paste0(outdir, '/', cm, '_', l, '.splitFeaturePlot.pdf'), width=12, height=8, useDingbats=FALSE)
      print(enhanced_feature_plot(PRO, features=g_u, pt_size=pt_use))
      dev.off()
      png(paste0(outdir, '/', cm, '_', l, '.splitFeaturePlot.png'), width=2400, height=1600, res=300)
      print(enhanced_feature_plot(PRO, features=g_u, pt_size=pt_use))
      dev.off()
      pdf(paste0(outdir, '/', cm, '_', l, '.splitFeaturePlot_tsne.pdf'), width=12, height=8, useDingbats=FALSE)
      print(enhanced_feature_plot(PRO, features=g_u, reduction_type="tsne", pt_size=pt_use))
      dev.off()
      png(paste0(outdir, '/', cm, '_', l, '.splitFeaturePlot_tsne.png'), width=2400, height=1600, res=300)
      print(enhanced_feature_plot(PRO, features=g_u, reduction_type="tsne", pt_size=pt_use))
      dev.off()
      plots <- VlnPlot(PRO, features=g_u, split.by="sample", group.by="celltype", pt.size=0, combine=FALSE)
      pdf(paste0(outdir, '/', cm, '_', l, '.splitVlnPlot.pdf'), width=10, height=12)
      print(CombinePlots(plots=plots, ncol=1))
      dev.off()
      png(paste0(outdir, '/', cm, '_', l, '.splitVlnPlot.png'), width=2000, height=2400, res=300)
      print(CombinePlots(plots=plots, ncol=1))
      dev.off()
    }
  }
}