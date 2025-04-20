############################################ figure1 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
rm(list = ls()); gc(); graphics.off()

fwd <- "~/09_paper_figure/figure1"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}


# volcano plot -------------------------------------------------------------------------
library(DESeq2)
load("~/01_bulk/data/DESeq_dds.rdata")
res <- results(dds) %>% as.data.frame() %>% na.omit() %>% mutate(gene = rownames(.))

res$threshold <- "Not(n = 17584)"  
res$threshold[res$log2FoldChange > 0 & res$padj < 0.05] <- "Up(n = 431)"  
res$threshold[res$log2FoldChange < 0 & res$padj < 0.05] <- "Down(n = 2658)"

genes <- c(
  'ACSL1', 'ANPEP', 'APLNR', 'ARL1', 'BCL2', 'BCL6', 'CCT5', 'CREB5',
  'CYTH1', 'EZR', 'FLOT2', 'GCA', 'GGA2', 'HSPA8', 'IKZF1', 'LAMTOR4',
  'LBH', 'MNDA', 'MRPS23', 'MTATP6P1', 'MTRNR2L12', 'NEAT1', 'PTPN6',
  'RAB11FIP4', 'RBM47', 'RPL3P4', 'RPL6P27', 'RPS17', 'RUNX2', 'SNX30',
  'STIP1', 'TSPYL1', 'YBEY', 'ZBTB18'); genes %>% length()
data_ggrepel <- res %>% filter(gene %in% genes) %>% arrange(-log2FoldChange) %>% head(10)
(plot <- ggplot(
  data = res,
  aes(x = log2FoldChange, y = -log10(padj), colour = threshold, fill = threshold)) +
    scale_color_manual(values = c("#3471ab", "grey", "#b23f31")) +
    geom_point(alpha = 1, size = 1.2) +
    geom_vline(xintercept = 0, lty = 4, col = "grey", lwd = 0.6) +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.6) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggtitle("cfRNA DEGs") + 
    ggrepel::geom_label_repel(
      data = data_ggrepel, 
      aes(label = gene), 
      color = "black", 
      fill = "white",
      alpha = 0.8, 
      fontface = "italic",
      min.segment.length = 0,
      box.padding = 0.5,
      seed = 42,
      size = 3) +
    labs(x = "log2FC", y = "-log10(pvalue)"))
ggsave(paste0(fwd, "/cfRNA_volcano.plot.pdf"), plot, width = 4, height = 5)


# enrich plot -------------------------------------------------------------------------
enrichr <- read.table("~/09_paper_figure/cfRNA_enrich_filter.txt", sep = "\t", header = T, check.names = F)
enrichr$p <- enrichr$`-log10(pvalue)`
enrichr <- enrichr %>% arrange(Cluster, -p)
rich.res <- list()
rich.res[["GOBP up"]] <- enrichr %>% filter(Cluster == "GOBP up")
rich.res[["KEGG up"]] <- enrichr %>% filter(Cluster == "KEGG up")
rich.res[["GOBP down"]] <- enrichr %>% filter(Cluster == "GOBP down")
rich.res[["KEGG down"]] <- enrichr %>% filter(Cluster == "KEGG down")

rich_all <- do.call(rbind, lapply(names(rich.res), function(i) {
  tmp.go <- rich.res[[i]] %>% mutate(Cluster = i)
  return(tmp.go)
}))
rich_all <- rich_all %>% group_by(Cluster) %>% top_n(n = 100, wt = p) %>% as.data.frame()

rich_all$Cluster <- factor(rich_all$Cluster, levels = c("GOBP up", "KEGG up", "GOBP down", "KEGG down"))
rich_all$Description <- make.unique(rich_all$Description)
rich_all$Description <- gsub(".1", " ", rich_all$Description)
rich_all$Description <- gsub(".2", "  ", rich_all$Description)
rich_all$Description <- factor(rich_all$Description, levels = rich_all$Description)

(p1 <- ggplot(rich_all, aes(x = `-log10(pvalue)`, y = rev(Description), fill = Cluster)) +
    geom_bar(stat = "identity", width = 0.9, alpha = 0.6) +
    geom_text(aes(x = 0.1, y = rev(Description), label = Description), size = 4, hjust = 0) +
    theme_bw() +
    scale_fill_manual(values = c("#ff6666", "#ff9999", "#3471ab", "#56b1f7")) +
    scale_x_continuous(expand = c(0.01, 0.01)) +
    scale_y_discrete(expand = c(0.02, 0.02)))
ggsave(paste0(fwd, "/cfRNA_Enrichment_barplot.pdf"), plot, width = 5, height = 6)


# cfRNA 47 genes -------------------------------------------------------------------------
feature_import <- read.csv("~/09_paper_figure/cfRNA_47genes_feature_importances.csv") %>%
  arrange(Importance) %>% {.$Feature = factor(.$Feature, levels = .$Feature); .}

tmp_color <- c("#1e466e", "#376795", "#528fad", "#72bcd5", "#aadce0", "#ffe6b7", "#ffd06f", "#F05454", "#CE1212", "#810000")
(plot <- ggplot(feature_import, aes(y = Importance, x = Feature)) +
    geom_segment(aes(y = 0, x = Feature, xend = Feature, yend = Importance), size = 1, color = "grey40") +
    scale_color_gradient(low = "#ffe6b7", high = "#1e466e") +
    geom_point(aes(color = Importance), size = 4) +
    labs(x = "", y = "Importance") +
    coord_flip() +
    scale_y_continuous(expand = expansion(add = c(0.01, 0.05))))
ggsave(paste0(fwd, "/cfRNA_47Genes_importance.pdf"), plot, width = 5, height = 8, limitsize = FALSE)



############################################ figure2 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
rm(list = ls()); gc(); graphics.off()

fwd <- "~/09_paper_figure/figure2"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}

load("~/02_scRNA/result/cols.for.Major.Cluster.rdata")
load("~/02_scRNA/data/brain.Major.Cluster.rdata") # brain.Celldefine
tmp.seuObj <- brain.Celldefine 


# DimPlot AD pie -------------------------------------------------------------------------
library(scatterpie)
library(tidydr)

tmp.seuObj <- brain.Celldefine 
tmp.seuObj$Group[which(tmp.seuObj$Group == "Aging")] <- "Control"

tmp.seuObj$Major.Cluster[which(tmp.seuObj$Major.Cluster == "Perivascular Fibroblast")] <- "PVFs"
tmp.seuObj$Major.Cluster[which(tmp.seuObj$Major.Cluster == "Endothelial Cell")] <- "Endothelial"

df <- tmp.seuObj@reductions$umap@cell.embeddings %>% as.data.frame() %>% 
  cbind(CellType = tmp.seuObj@meta.data$Major.Cluster)

(freq <- prop.table(table(tmp.seuObj$Major.Cluster, tmp.seuObj$Group), margin = 1) %>% as.data.frame())
(freq <- spread(freq, Var2, Freq))
colnames(freq)[1] <- "CellType"

label <- df %>% group_by(CellType) %>% summarise(umap_1 = median(umap_1), umap_2 = median(umap_2)) %>% as.data.frame()
cell_number <- as.data.frame(table(tmp.seuObj$Major.Cluster))
colnames(cell_number)[2] <- "cellnumber"

data <- merge(freq, label, by.x = "CellType", by.y = "CellType")
data <- merge(data, cell_number, by.x = "CellType", by.y = "Var1")

col1 <- c("#06b95c", "#B0C399", "#B8C6E9", "#00738C", "#036C00", "#FFE4B5", "#DD9FD7", "#F39800")
col1 <- setNames(col1, c("Neuron", "mOli", "Astrocyte", "OPC", "Microglia", "PVFs", "Endothelial", "Pericyte")); col1
(plot <- ggplot() +
    geom_point(data = df, aes(x = umap_1, y = umap_2, color = CellType), size = 1, shape = 16) +
    guides(color = guide_legend(override.aes = list(size = 6), keywidth = 0)) +
    scale_color_manual(values = col1) +
    geom_label(data = label, aes(x = umap_1, y = umap_2 + 1, label = CellType), fill = "white", color = "black", size = 3.5) +
    geom_scatterpie(
      data = data,
      aes(x = umap_1, y = umap_2 - 0.5, group = CellType, r = cellnumber / 10),
      cols = names(freq)[2:3]) +
    scale_fill_manual(values = c("#f8bf98", "#d3d3d2"), name = ""))   # 修改扇形图填充颜色
ggsave(paste0(fwd, "/DimPlot_AD.pie.pdf"), plot, width = 6.6, height = 5)


# AD Amount -------------------------------------------------------------------------
tmp.seuObj <- brain.Celldefine
tmp.seuObj$Group[which(tmp.seuObj$Group == "Aging")] <- "Control"

col1 <- c("#06b95c", "#B0C399", "#B8C6E9", "#00738C", "#036C00", "#FFE4B5", "#DD9FD7", "#F39800")
col1 <- setNames(col1, c("Neuron", "mOli", "Astrocyte", "OPC", "Microglia", "PVFs", "Endothelial", "Pericyte")); col1

(Cellratio <- prop.table(table(tmp.seuObj$Major.Cluster, tmp.seuObj$Group), margin = 2) %>% data.frame())
levels1 <- c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs")
colnames(Cellratio) <- c("Celltype", "sample", "amount")
Cellratio$sample <- as.factor(Cellratio$sample)
Cellratio$amount <- Cellratio$amount * 100
Cellratio$Celltype <- factor(Cellratio$Celltype, levels = rev(levels1))
Cellratio$sample <- Cellratio$sample %>% as.character()
Cellratio$sample[which(Cellratio$sample == "Aging")] <- "Control"
Cellratio$sample <- factor(Cellratio$sample, levels = c("AD", "Control"))

(plot <- ggplot(Cellratio, aes(x = sample, y = amount, fill = Celltype)) +
    xlab("Group") + 
    guides(fill = guide_legend(reverse = T)) + 
    geom_bar(stat = "identity", position = "fill", width = .95) +
    scale_fill_manual(values = col1) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))))
ggsave(paste0(fwd, "/AD_Amount.pdf"), plot, width = 2.5, height = 3.2, limitsize = FALSE)


# scRNA.signature.genes_in_cfRNA -------------------------------------------------------------------------
# - violin -------------------------------------------------------------------------
# cfRNA data
library(DESeq2)
load("~/01_bulk/data/DESeq_dds.rdata")
vsd <- vst(dds, blind = FALSE)
normalizeExp <- assay(vsd)
normalizeExp[1:4, 1:4]

load("~/02_scRNA/data/allmarkers.Major.Cluster.rdata") # allmarkers
markers <- allmarkers %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% 
  arrange(-avg_log2FC) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% data.frame()
markers$cluster %>% table()

markers$cluster <- markers$cluster %>% as.character()
markers$cluster[which(markers$cluster == "Endothelial Cell")] <- "Endothelial"
markers$cluster[which(markers$cluster == "Perivascular Fibroblast")] <- "PVFs"
markers$cluster <- factor(markers$cluster, levels = c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs"))

cfRNA_seurat <- CreateSeuratObject(counts = normalizeExp, project = "cfRNA")
# cfRNA_seurat@assays$RNA@data[1:4, 1:4]
# cfRNA_seurat@assays$RNA@counts[1:4, 1:4]

# unique(markers$cluster)
for (tmp_cluster in unique(markers$cluster)) {
  # tmp_cluster <- "Neuron"
  cfRNA_seurat <- AddModuleScore(cfRNA_seurat, features = list(markers %>% filter(cluster == tmp_cluster) %>% .$gene), name = tmp_cluster)
}

load("~/01_bulk/data/colDatas.rdata")
identical(Cells(cfRNA_seurat), rownames(colDatas))
# [1] TRUE

cfRNA_seurat@meta.data %>% head()
cfRNA_seurat$group <- colDatas$condition
cfRNA_seurat$group <- cfRNA_seurat$group %>% as.character()
cfRNA_seurat$group[which(cfRNA_seurat$group == "control")] <- "Control"
# cfRNA_seurat %>% VlnPlot(features = "mOli1", group.by = "group")
colnames(cfRNA_seurat@meta.data) <- gsub("1", "", colnames(cfRNA_seurat@meta.data))

col1 <- c("#B8C6E9", "#DD9FD7", "#036C00", "#B0C399", "#06b95c", "#00738C", "#F39800", "#FFE4B5")
col1 <- setNames(col1, c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs")); col1

tmp_Plist <- list()
for (features in c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs")) {
  tmp.seuObj <- cfRNA_seurat
  # features <- "mOli"
  values <- tmp.seuObj[[features]][, 1]
  group <- tmp.seuObj$group
  
  df <- data.frame(cluster = features, value = values, group = group)
  df$group <- as.character(df$group)
  
  (sorted_clusters <- data.frame(cluster = group, value = values) %>% 
      group_by(cluster) %>% summarize(mean_value = mean(value)) %>%
      arrange(mean_value) %>% pull(cluster) %>% rev() %>% as.character())
  
  cols <- col1 %>% as.data.frame() %>% 
    {.$cellty = rownames(.); colnames(.) = c("col", "cellty"); .} %>%
    filter(cellty == features) %>% .$col; cols
  
  plot <- ggplot(df, aes(x = group, y = value, group = group, fill = group)) +
    geom_jitter(position = position_jitter(w = 0.3), size = 0, alpha = 0.1) +
    geom_violin() +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    xlab("Age") +
    ylab("Expression") +
    facet_grid(. ~ cluster) +
    scale_fill_manual(values = c("#f9d2b6", "#e0e0e0")) +
    RotatedAxis()
  tmp_Plist[[as.character(features)]] <- plot
}
(plot <- Reduce("+", tmp_Plist) + plot_layout(ncol = 8, guides = "collect"))
ggsave(paste0(fwd, "/scRNA.signature.genes_in_cfRNA_VlnPlot.pdf"), plot, width = 14, height = 6)

# - Count -------------------------------------------------------------------------
library(DESeq2)
load("~/01_bulk/data/DESeq_dds.rdata")
res <- results(dds) %>% as.data.frame() %>% na.omit()
res$gene <- rownames(res)

res$threshold <- "Not"  
res$threshold[res$log2FoldChange > 0] <- "Up"  
res$threshold[res$log2FoldChange < 0] <- "Down"
res <- res %>% dplyr::select("gene", "threshold")

merge_data <- merge(res, markers, by.x = "gene", by.y = "gene") %>% filter(threshold %in% c("Up", "Down"))
res_df <- table(merge_data$threshold, merge_data$cluster) %>% data.frame()

(plot <- ggplot(res_df, aes(x = Var2, y = Freq, fill = Var2, group = Var2)) +
    labs(x = "", y = "") +
    geom_bar(stat = "identity", width = .95) +
    facet_grid(. ~ Var1, scales = "free") +
    scale_fill_manual(values = col1) +
    ylab("Count") +
    guides(fill = guide_legend(reverse = T)) +
    RotatedAxis() +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_x_discrete(expand = expansion(mult = c(0.09, 0.09))))
ggsave(paste0(fwd, "/scRNA.signature.genes_in_cfRNA_Count.pdf"), plot, width = 6, height = 4)




############################################ figure3 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
# rm(list = ls()); gc(); graphics.off()
source("/home/wuli/R/Rfile/load_packages.R")

fwd <- "~/09_paper_figure/figure3"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}


# multiple cluster genes -------------------------------------------------------------------------
rm(list = setdiff(ls(), "fwd"))

pseudo_all <- readRDS("scRNA_DEGs.rds"); dim(pseudo_all)
df.pb <- pseudo_all %>% select(genes, log2FoldChange, cellty, padj) %>% 
  {colnames(.)[2] = "log2FC"; colnames(.)[2] = "log2FC_pseudo"; colnames(.)[4] = "padj_pseudo"; .}

load("~/01_bulk/data/DESeq_dds.rdata") 
cfRNA_regulated <- DESeq2::results(dds) %>% as.data.frame() %>% na.omit(); dim(cfRNA_regulated)
df.cf <- cfRNA_regulated %>% 
  mutate(genes = rownames(.)) %>% select(genes, log2FoldChange, padj) %>% 
  {colnames(.)[2] = "log2FC"; colnames(.)[2] = "log2FC_cfRNA"; colnames(.)[3] = "padj_cfRNA"; .}
df.all <- merge(df.pb, df.cf, by = "genes")
saveRDS(df.all, file = "merge_cfRNA.and.scRNA_DEGs.rds")

# - up -------------------------------------------------------------------------
merge_DEGs <- readRDS("merge_cfRNA.and.scRNA_DEGs.rds")
tmp_df <- merge_DEGs %>% filter(log2FC_pseudo > 0.25 & log2FC_cfRNA > 0 & padj_pseudo < 0.05 & padj_cfRNA < 0.05)

(res <- tmp_df %>%
    group_by(genes) %>%
    summarize(count = n_distinct(cellty)) %>%
    filter(count >= 3) %>%
    arrange(-count) %>%
    as.data.frame())
order <- c(res$genes)

n <- 1.98
res.df <- tmp_df %>% filter(genes %in% res$genes) %>% mutate(per = n)

labels <- res
labels$idx <- seq(1:length(labels$gene))
angle <- 90 - 360 * (labels$idx - 0.5) / length(labels$gene)
labels$hjust <- ifelse(angle < -90, 1, 0)
labels$angle <- ifelse(angle < -90, angle + 180, angle)

(plot <- ggplot(res.df, aes(x = genes, y = per, fill = cellty)) +
    geom_bar(stat = "identity", color = "black", size = 0.4, position = "stack", width = 0.8) +
    ylim(-3, n * 9 + 0.5) +
    scale_fill_manual(values = col1) +
    coord_polar(start = 0) +
    scale_x_discrete(limits = order) +
    geom_rect(
      xmin = length(res$gene) + 0.5, 
      xmax = length(res$gene) + length(res$gene) + 0.4,
      ymin = -0.8, 
      ymax = -0.2, 
      fill = "#c86468") +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()) +
    geom_text(
      data = labels,
      aes(x = genes, y = count * n + 0.1, label = genes, angle = angle, hjust = hjust),
      color = "black",
      size = 4,
      inherit.aes = FALSE) +
    guides(fill = FALSE))
ggsave(paste0(fwd, "/multiple.cluster.genes_up.pdf"), plot, width = 8, height = 8, limitsize = FALSE)

# - down -------------------------------------------------------------------------
merge_DEGs <- readRDS("merge_cfRNA.and.scRNA_DEGs.rds")
tmp_df <- merge_DEGs %>% filter(log2FC_pseudo < -0.25 & log2FC_cfRNA < 0 & padj_pseudo < 0.05 & padj_cfRNA < 0.05)

(res <- tmp_df %>%
    group_by(genes) %>%
    summarize(count = n_distinct(cellty)) %>%
    filter(count >= 4) %>%
    arrange(-count) %>%
    as.data.frame())
order <- c(res$genes)

n <- 1.98
res.df <- tmp_df %>% filter(genes %in% res$genes) %>% mutate(per = n)

labels <- res
labels$idx <- seq(1:length(labels$gene))
angle <- 90 - 360 * (labels$idx - 0.5) / length(labels$gene)
labels$hjust <- ifelse(angle < -90, 1, 0)
labels$angle <- ifelse(angle < -90, angle + 180, angle)

(plot <- ggplot(res.df, aes(x = genes, y = per, fill = cellty)) +
    geom_bar(stat = "identity", color = "black", size = 0.4, position = "stack", width = 0.8) +
    ylim(-3, n * 9 + 0.5) +
    scale_fill_manual(values = col1) +
    coord_polar(start = 0) +
    scale_x_discrete(limits = order) +
    geom_rect(
      xmin = length(res$gene) + 0.5, 
      xmax = length(res$gene) + length(res$gene) + 0.4,
      ymin = -0.8, 
      ymax = -0.2, 
      fill = "#4d83aa") +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank()) +
    geom_text(
      data = labels,
      aes(x = genes, y = count * n + 0.1, label = genes, angle = angle, hjust = hjust),
      color = "black",
      size = 4,
      inherit.aes = FALSE) +
    guides(fill = FALSE))
ggsave(paste0(fwd, "/multiple.cluster.genes_down.pdf"), plot, width = 8, height = 8, limitsize = FALSE)



############################################ figure4 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
rm(list = ls()); gc(); graphics.off()

fwd <- "~/09_paper_figure/figure4"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}

load("~/02_scRNA/result/cols.for.Major.Cluster.rdata")
load("~/02_scRNA/data/brain.Major.Cluster.rdata") # brain.Celldefine
tmp.seuObj <- brain.Celldefine 



# 34Biomarkers Importance -------------------------------------------------------------------------
genes <- c("MNDA", "ANPEP", "CREB5", "RUNX2", "ACSL1", "MTATP6P1", "BCL6", "LBH", "LAMTOR4", "HSPA8", "APLNR", "NEAT1", "GGA2", "LRRK2", "CCT5", "IKZF1", "SELPLG", "GCA", "RPS17", "MTRNR2L12", "RPL6P27", "ECI2", "ZBTB18", "ARL1", "ZFP36L1", "RASSF2", "PTPN6", "SESN3", "FLOT2", "RBM47")
expression_values <- c(0.05198736, 0.03127601, 0.03058308, 0.02764245, 0.02441252, 0.02313158, 0.02211968, 0.02192715, 0.01954864, 0.01949003, 0.01941373, 0.01816507, 0.0179387, 0.01743533, 0.01583516, 0.01479922, 0.01472388, 0.01330699, 0.01244648, 0.01210679, 0.01196544, 0.01194029, 0.01112793, 0.0111106, 0.01103344, 0.01075214, 0.01074222, 0.01046621, 0.00984259, 0.00968171)

feature_import <- data.frame(Feature = genes, Importance = expression_values) %>% 
  arrange(-Importance) %>% {.$Feature = factor(.$Feature, levels = .$Feature); .}

tmp_color <- c("#1e466e", "#376795", "#528fad", "#72bcd5", "#aadce0", "#ffe6b7", "#ffd06f", "#F05454", "#CE1212", "#810000")
(plot <- ggplot(feature_import, aes(y = Importance, x = Feature)) +
    geom_segment(aes(y = 0, x = Feature, xend = Feature, yend = Importance), size = 1, color = "grey40") +
    scale_color_gradient(low = "#ffe6b7", high = "#1e466e") +
    geom_point(aes(color = Importance), size = 4) +
    labs(x = "", y = "Importance") +
    scale_y_continuous(expand = expansion(add = c(0.001, 0.001))))
ggsave(paste0(fwd, "/34Biomarkers_Importance.png"), plot, width = 8, height = 5, limitsize = FALSE)
ggsave(paste0(fwd, "/34Biomarkers_Importance.pdf"), plot, width = 8, height = 5, limitsize = FALSE)


# 34Biomarkers.in.scRNA_heatmap -------------------------------------------------------------------------
load("~/02_scRNA/data/pseudo_exprMat_noOverlap_cfRNA.rdata") # pseudo_exprMat

clusters <- c("Astrocyte", "Endothelial Cell", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "Perivascular Fibroblast")
genes <- c(
  'ACSL1', 'ANPEP', 'APLNR', 'ARL1', 'BCL2', 'BCL6', 'CCT5', 'CREB5',
  'CYTH1', 'EZR', 'FLOT2', 'GCA', 'GGA2', 'HSPA8', 'IKZF1', 'LAMTOR4',
  'LBH', 'MNDA', 'MRPS23', 'MTATP6P1', 'MTRNR2L12', 'NEAT1', 'PTPN6',
  'RAB11FIP4', 'RBM47', 'RPL3P4', 'RPL6P27', 'RPS17', 'RUNX2', 'SNX30',
  'STIP1', 'TSPYL1', 'YBEY', 'ZBTB18')

i <- 1
for (tmp_cluster in clusters) {
  # tmp_cluster <- "Perivascular Fibroblast"
  print(tmp_cluster)
  
  exprMat <- pseudo_exprMat[, grep(tmp_cluster, colnames(pseudo_exprMat))]
  # exprMat %>% colnames()
  (ad_len <- length(grep("AD", colnames(exprMat))))
  (case_len <- length(grep("Aging", colnames(exprMat))))
  (group <- rep(c("AD", "Control"), c(ad_len, case_len)))
  (celltype <- rep(tmp_cluster, ad_len + case_len))
  
  condition <- factor(rep(c("AD", "Control"), c(ad_len, case_len)), levels = c("AD", "Control"))
  condition <- relevel(condition, ref = "Control")
  colDatas <- data.frame(row.names = colnames(exprMat), condition)
  
  exprMat1 <- round(exprMat, 0)
  dds1 <- DESeqDataSetFromMatrix(
    countData = exprMat1,
    colData = colDatas,
    design = ~ condition)
  dds2 <- dds1[rowSums(counts(dds1)) >= 0, ]
  # dds <- DESeq(dds2)
  
  vsd <- vst(dds2, blind = FALSE)
  tmp_normalizeExp <- assay(vsd)
  diff_expr <- tmp_normalizeExp[genes, ] %>% data.frame()
  
  zscore <- (diff_expr - apply(diff_expr, 1, mean)) / apply(diff_expr, 1, sd)
  zscore[is.na(zscore)] <- 0
  
  annotation_col <- data.frame(group = group, celltype = celltype)
  row.names(annotation_col) <- colnames(zscore)
  
  # celltype
  celltypecolor <- c("#B8C6E9", "#DD9FD7", "#036C00", "#B0C399", "#06b95c", "#00738C", "#F39800", "#FFE4B5") 
  names(celltypecolor) <- clusters
  celltypecolor <- celltypecolor %>% data.frame() %>% {colnames(.) = "col";. } %>% mutate(cellty = rownames(.)) %>% 
    filter(cellty == tmp_cluster) %>% .$col
  names(celltypecolor) <- tmp_cluster
  
  # group
  groupcolor <- c("#eeb996", "#d1d2d0")
  names(groupcolor) <- c('AD', 'Control')
  
  # pdf(paste0(fwd, "/pheatmap_34ML_Genes_scRNA.pdf"), width = 14, height = 6)
  ann_colors <- list(celltype = celltypecolor, group = groupcolor)
  (plot <- pheatmap::pheatmap(
    zscore,
    # scale = "column",
    annotation_colors = ann_colors,
    annotation_col = annotation_col,
    column_split = annotation_col$celltype, 
    annotation_names_row = F, 
    annotation_names_col = F, 
    column_title = NULL,
    row_title = NULL, 
    border_color = NA))
  # dev.off()
  
  fwd1 <- "~/09_paper_figure/figure4/ML34genes"
  if (!file.exists(fwd1)) {dir.create(fwd1, recursive = T)}
  
  # i <- 1
  pdf(paste0(fwd, "/pheatmap_34ML_Genes_scRNA_", i, ".pdf"), width = (case_len * 1.2 + ad_len * 1.2), height = 4)
  print(plot)
  i <- i +1
  dev.off()
}


# 34Biomarkers AddModuleScore -------------------------------------------------------------------------
genes <- c(
  'ACSL1', 'ANPEP', 'APLNR', 'ARL1', 'BCL2', 'BCL6', 'CCT5', 'CREB5',
  'CYTH1', 'EZR', 'FLOT2', 'GCA', 'GGA2', 'HSPA8', 'IKZF1', 'LAMTOR4',
  'LBH', 'MNDA', 'MRPS23', 'MTATP6P1', 'MTRNR2L12', 'NEAT1', 'PTPN6',
  'RAB11FIP4', 'RBM47', 'RPL3P4', 'RPL6P27', 'RPS17', 'RUNX2', 'SNX30',
  'STIP1', 'TSPYL1', 'YBEY', 'ZBTB18')

tmp.seuObj <- brain.Celldefine
tmp.seuObj <- AddModuleScore(tmp.seuObj, features = list(intersect(genes, rownames(tmp.seuObj))), name = "a34genes")
tmp.seuObj$Group[which(tmp.seuObj$Group == "Aging")] <- "Control"

gg_data <- data.frame(cluster = tmp.seuObj$Major.Cluster, value = tmp.seuObj$a34genes1, group = tmp.seuObj$Group) %>% 
  filter(cluster %in% c("Microglia", "Astrocyte", "mOli", "OPC", "Neuron"))
gg_data %>% head()

### violin
(levs <- data.frame(cluster = factor(gg_data$cluster), value = gg_data$value) %>%
    group_by(cluster) %>% dplyr::summarize(mean_value = mean(value)) %>%
    arrange(mean_value) %>% pull(cluster) %>% rev())
gg_data$cluster <- factor(gg_data$cluster, levels = levs)

gg_data %>% head()
pval_df <- data.frame()
n <- 1
for (tmp_cluster in unique(gg_data$cluster)) {
  # tmp_cluster <- "Astrocyte"
  print(tmp_cluster)
  
  tmp1 <- gg_data %>% filter(cluster == tmp_cluster) %>% filter(group == "AD")
  tmp2 <- gg_data %>% filter(cluster == tmp_cluster) %>% filter(group == "Control")
  
  (tmp_p <- wilcox.test(tmp1$value, tmp2$value)$p.value)
  tmp_df <- data.frame(cluster = tmp_cluster, pvalue = tmp_p, n = n)
  pval_df <- rbind(pval_df, tmp_df)
  n <- n + 1
}

(pval_df <- pval_df %>%
    mutate(significant = case_when(
      pvalue < 0.0001 ~ "****",
      pvalue >= 0.0001 & pvalue < 0.001 ~ "***",
      pvalue >= 0.001 & pvalue < 0.01 ~ "**",
      pvalue >= 0.01 & pvalue < 0.05 ~ "*",
      pvalue >= 0.05 ~ "ns")))

pval_idx <- data.frame(cluster = factor(c("Microglia", "Neuron", "Astrocyte", "mOli", "OPC")))
pval_idx$x1 <- 1:5
pval_idx$y1 <- 0.7
pval_idx$significant <- "****"

(plot <- ggplot() +
    scale_fill_manual(values = c("#f3cfb6", "#dfdfdf")) +
    ylab("Module Score") + xlab("") +
    geom_boxplot(data = gg_data, aes(x = cluster, y = value, fill = group), width = 0.3, alpha = 0.7, outlier.shape = NA) +
    geom_text(data = pval_idx, mapping = aes(x = x1, y = y1, label = significant)))
ggsave(paste0(fwd, "/34Biomarkers_ModuleScore_in.scRNA.pdf"), plot, width = 4, height = 4, limitsize = FALSE)




############################################ figure5 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
source("/home/wuli/R/Rfile/load_packages.R")

fwd <- "~/09_paper_figure/figure5"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}

# load("~/02_scRNA/result/cols.for.Major.Cluster.rdata")
# load("~/02_scRNA/data/brain.Major.Cluster.rdata") # brain.Celldefine
# tmp.seuObj <- brain.Celldefine 


# MSBB NMF -------------------------------------------------------------------------
# https://cloud.tencent.com/developer/article/1936854
# https://mp.weixin.qq.com/s/l5GQGK3hiTjxsZOFInlJiA
library(NMF)

# - run NMF -------------------------------------------------------------------------
dds <- readRDS("~/08_RNA_vaild/data/MSBB_dds_res_24.09.17.rds")
filter_meta <- readRDS("~/08_RNA_vaild/data/MSBB_meta.rds")

genes <- c(
  'ACSL1', 'ANPEP', 'APLNR', 'ARL1', 'BCL2', 'BCL6', 'CCT5', 'CREB5',
  'CYTH1', 'EZR', 'FLOT2', 'GCA', 'GGA2', 'HSPA8', 'IKZF1', 'LAMTOR4',
  'LBH', 'MNDA', 'MRPS23', 'MTATP6P1', 'MTRNR2L12', 'NEAT1', 'PTPN6',
  'RAB11FIP4', 'RBM47', 'RPL3P4', 'RPL6P27', 'RPS17', 'RUNX2', 'SNX30',
  'STIP1', 'TSPYL1', 'YBEY', 'ZBTB18')
vsd <- vst(dds, blind = FALSE) # 标准化Counts
normalizeExp <- assay(vsd)

meta1 <- filter_meta %>% filter(AD_Status %in% c("AD")) %>% arrange(AD_Status)
diff_expr <- normalizeExp[genes, meta1$specimenID] %>% data.frame()
identical(colnames(diff_expr), meta1$specimenID)

nmf_data <- diff_expr
nmf_data[nmf_data < 0] <- 0
table(nmf_data < 0)

result2 <- nmf(
  nmf_data,
  rank = 2,
  seed = 1234)
index <- extractFeatures(result2, "max")
group <- predict(result2) 
table(group)
# group
# 1   2 
# 139 207 

filter_meta %>% data.frame() %>% head()
sample_df <- data.frame(group = group) %>% mutate(sample = rownames(.))
sample_df %>% head()
sample_df1 <- left_join(sample_df, filter_meta, by = c("sample" = "specimenID")) %>% 
  dplyr::select(sample, ageDeath, group)
sample_df1 %>% head()

cairo_pdf(paste0(fwd, "/MSBB_consensusmap.pdf"), width = 5, height = 4, pointsize = 18, bg = "white")
consensusmap(
  result2,
  labRow  = NA,
  labCol = NA,
  annCol = data.frame("cluster" = group[colnames(nmf_data)]))
dev.off()

cairo_pdf(paste0(fwd, "/MSBB_basismap_genes.pdf"), width = 4, height = 5, pointsize = 18, bg = "white")
basismap(
  result2,
  cexCol = 1,
  cexRow = 1)
dev.off()

# - plaqueMean -------------------------------------------------------------------------
NMF_result <- as.data.frame(group) %>% mutate(sample = rownames(.)) %>% left_join(., meta1, by = c("sample" = "specimenID"))
(mean1 <- NMF_result %>% filter(group == 1) %>% .$plaqueMean %>% mean())
(mean2 <- NMF_result %>% filter(group == 2) %>% .$plaqueMean %>% mean())

df_anno <- data.frame(xid = (mean1 + mean2) / 2, yid = 0.065, label = "**")
(plot <- ggplot() +
    geom_density(data = NMF_result, aes(x = plaqueMean, fill = group), linewidth = 0, alpha = 0.5, size = 0.8)+
    geom_vline(xintercept = mean1, lty = 4, col = "#5080a2", lwd = 0.6) +
    geom_vline(xintercept = mean2, lty = 4, col = "#bc6166", lwd = 0.6) +
    scale_fill_manual(values = c("#5080a2", "#bc6166")) + 
    guides(fill = guide_legend(title = "group")) +
    geom_text(data = df_anno, aes(x = xid, y = yid, label = label), size = 4, fontface = "italic") +
    labs(title = "MSBB", x = "plaque Mean", y = "density") +
    scale_y_continuous(expand = c(0.001, 0.001)) +
    scale_x_continuous(expand = c(0.01, 0.01)))
ggsave(paste0(fwd, "/MSBB_NMF.Group_plaqueMean.pdf"), plot, width = 4, height = 3, limitsize = FALSE)



############################################ figure s1 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
rm(list = ls()); gc(); graphics.off()

fwd <- "~/09_paper_figure/figure.s1"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}


# cfRNA sample -------------------------------------------------------------------------
meta <- read.csv("~/09_paper_figure/ad_sample_info.csv", sep = ",") %>% filter(Center != "Unknown")
meta %>% head()

(Center <- unique(meta$Center))
sample_cols <- c("#C66594", "#D1641F", "#3C6AA1", "#5E4286", "#25823C")
sample_cols <- setNames(sample_cols, Center); sample_cols

(res <- table(meta$Disease, meta$Center) %>% data.frame() %>% 
    filter(Var1 != "None") %>% {colnames(.) = c("group", "Center", "count"); .})
(tmp_label <- table(meta$Disease) %>% data.frame())
(plot <- ggplot() +
    geom_bar(data = res, aes(x = count, y = group, fill = Center), width = 0.9, stat = "identity") +
    ylab("# Samples") +
    scale_fill_manual(values = sample_cols) +
    geom_text(data = tmp_label, aes(x = Freq, y = Var1, label = Freq), hjust = 0.5, color = "black", size = 4) +
    xlim(c(0, 190)) +
    scale_x_continuous(limits = c(0, 180), expand = c(0.01, 0.01)) +
    scale_y_discrete(expand = c(0.29, 0.29)))
ggsave(paste0(fwd, "/cfRNA_sample.pdf"), plot, width = 8.5, height = 1.8, limitsize = FALSE, device = cairo_pdf)


# cfRNA batch -------------------------------------------------------------------------
library(DESeq2)
load("~/01_bulk/data/DESeq_dds.rdata")
load("~/01_bulk/data/colDatas.rdata")
vsd <- vst(dds, blind = FALSE)

plot1 <- tinyarray::draw_pca(exp = assay(vsd), group_list = factor(colDatas$condition))
(plot2 <- plot1 + ggtitle("Batch Correction"))
ggsave(paste0(fwd, "/batch_corr.pdf"), plot2, width = 5, height = 5, limitsize = FALSE)

load("~/01_bulk/data/exprMat_final.rdata")
plot1 <- tinyarray::draw_pca(exp = exprMat, group_list = factor(colDatas$condition))
(plot2 <- plot1 + ggtitle("Batch Effect") + theme_bw())
ggsave(paste0(fwd, "/batch_effect.pdf"), plot2, width = 5, height = 5, limitsize = FALSE)



############################################ figure s2 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
source("/home/wuli/R/Rfile/load_packages.R")

fwd <- "~/09_paper_figure/figure.s2"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}


# - paper sample -------------------------------------------------------------------------
tmp.seuObj <- brain.Celldefine 
table(tmp.seuObj$paper, tmp.seuObj$Group) %>% data.frame()
df1 <- table(tmp.seuObj$paper, tmp.seuObj$sample, tmp.seuObj$Group) %>% data.frame() %>% filter(Freq > 0)
(df2 <- table(df1$Var1, df1$Var3) %>% data.frame())

(plot <- ggplot(df2, aes(x = Var1, y = Freq, fill = Var2, group = Var2)) +
    labs(x = "", y = "") +
    geom_bar(stat = "identity", width = .95) +
    scale_fill_manual(values = c("#f8bf98", "#d3d3d2")) +
    ylab("Count") +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))))
ggsave(paste0(fwd, "/paper_sample.pdf"), plot, width = 4, height = 6.5)
ggsave(paste0(fwd, "/paper_sample.png"), plot, width = 4, height = 6.5)

# - AD celltype count -------------------------------------------------------------------------
col1 <- c("#06b95c", "#B0C399", "#B8C6E9", "#00738C", "#036C00", "#FFE4B5", "#DD9FD7", "#F39800")
col1 <- setNames(col1, c("Neuron", "mOli", "Astrocyte", "OPC", "Microglia", "PVFs", "Endothelial", "Pericyte")); col1

tmp.seuObj <- brain.Celldefine 
(df <- table(tmp.seuObj$Major.Cluster, tmp.seuObj$Group) %>% data.frame())
df$Var2 <- df$Var2 %>% as.character()
df$Var2[which(df$Var2 == "Aging")] <- "Control"
df$Var1 <- df$Var1 %>% as.character()
df$Var1[which(df$Var1 == "Endothelial Cell")] <- "Endothelial"
df$Var1[which(df$Var1 == "Perivascular Fibroblast")] <- "PVFs"
(plot <- ggplot(df, aes(x = Var1, y = Freq, fill = Var1, group = Var1)) +
    labs(x = "", y = "") +
    geom_bar(stat = "identity", width = .95) +
    facet_grid(. ~ Var2, scales = "free") +
    scale_fill_manual(values = col1) +
    ylab("Count"))
ggsave(paste0(fwd, "/AD.celltype.count.pdf"), plot, width = 8, height = 4)

# - paper cells -------------------------------------------------------------------------
tmp_color <- c("#2b79b2", "#fa7c1a", "#389e6b", "#d01f27", "#a645f8")
papers <- c("HongjunSong_Nature_2022", "Andrew_Nature_2022", "Shun_ProcNatl_2020", "Samuel_NatGenet_2021", "Alexandra_NatNeurosc_2019")
tmp_color <- setNames(tmp_color, papers); tmp_color

(res <- table(tmp.seuObj$paper) %>% data.frame())
res <- res %>% arrange(-Freq) %>% {.$Var1 = factor(.$Var1, levels = rev(.$Var1)); .}
(plot <- ggplot(res, aes(x = Freq, y = Var1)) +
    geom_bar(aes(fill = Var1), width = 0.9, stat = "identity", position = "dodge") +
    # facet_grid(. ~ Var2, scales = "free") +
    scale_fill_manual(values = tmp_color) +
    xlab("# Cells"))
ggsave(paste0(fwd, "/paper.cells.pdf"), plot, width = 7, height = 4, limitsize = FALSE)

(res <- table(tmp.seuObj$Major.Cluster, tmp.seuObj$paper) %>% data.frame())
(plot <- ggplot(res, aes(x = Var1, y = Freq)) +
    geom_bar(aes(fill = Var1), width = 0.9, stat = "identity", position = "dodge") +
    facet_grid(. ~ Var2, scales = "free") +
    scale_fill_manual(values = col1) +
    ylab("Cell number"))
ggsave(paste0(fwd, "/Cell_Count_splitpaper.pdf"), plot, width = 16, height = 5, limitsize = FALSE)



# ML genes -------------------------------------------------------------------------
# - cfRNA -------------------------------------------------------------------------
# load("~/01_bulk/data/DESeq_dds.rdata")
# res <- results(dds) %>% as.data.frame() %>% na.omit()
up_regulated <- res %>% filter(as.numeric(padj) < 0.05 & as.numeric(log2FoldChange) >    1) %>% rownames()
dw_regulated <- res %>% filter(as.numeric(padj) < 0.05 & as.numeric(log2FoldChange) < (-3)) %>% rownames()
genes1 <- c(up_regulated, dw_regulated); genes1 %>% length()

train <- read.csv("~/01_bulk/data/cfRNA_train.csv", sep = ",")
vaild <- read.csv("~/01_bulk/data/cfRNA_vaild.csv", sep = ",")
identical(colnames(train), colnames(vaild))
colnames(train)[c(1, 109)]
(genes2 <- colnames(train)[-c(1, 109)])
genes2 <- gsub("\\.", "-", genes2); genes2 %>% length()

setdiff(genes2, genes1)
# character(0)

# - cfRNA&scRNA -------------------------------------------------------------------------
train <- read.csv("~/05_ML/data/data_train.csv", sep = ",")
vaild <- read.csv("~/05_ML/data/data_vaild.csv", sep = ",")
identical(colnames(train), colnames(vaild))
colnames(train)[c(1, 114)]
(genes1 <- colnames(train)[-c(1, 114)])
genes1 <- gsub("\\.", "-", genes1); genes1 %>% length()

cfRNA_regulated <- readRDS("cfRNA_DEGs.rds")
cf_up <- cfRNA_regulated %>% filter(group == "up")   %>% rownames(.); length(cf_up)
cf_dw <- cfRNA_regulated %>% filter(group == "down") %>% rownames(.); length(cf_dw)

scRNA <- readRDS("scRNA_DEGs.rds")
scRNA_up <- scRNA %>% filter(log2FoldChange > 0.25 & padj < 0.05);  dim(scRNA_up)
scRNA_dw <- scRNA %>% filter(log2FoldChange < -0.25 & padj < 0.05); dim(scRNA_dw)

(inter_up_genes <- intersect(cf_up, unique(scRNA_up$genes))); length(inter_up_genes)
(inter_dw_genes <- intersect(cf_dw, unique(scRNA_dw$genes))); length(inter_dw_genes)
setdiff(inter_up_genes, genes1)
setdiff(inter_dw_genes, genes1)
intersect(inter_dw_genes, genes1)



############################################ figure s3 ############################################
# library -------------------------------------------------------------------------
setwd("~/09_paper_figure")
rm(list = ls()); gc(); graphics.off()

fwd <- "~/09_paper_figure/figure.s3"
if (!file.exists(fwd)) {dir.create(fwd, recursive = T)}

# load("~/02_scRNA/result/cols.for.Major.Cluster.rdata")
# load("~/02_scRNA/data/brain.Major.Cluster.rdata") # brain.Celldefine
# tmp.seuObj <- brain.Celldefine 


# Detection rate of scRNA signature gene in cfRNA -------------------------------------------------------------------------
library(DESeq2)

load("~/02_scRNA/data/allmarkers.Major.Cluster.rdata") # allmarkers
markers <- allmarkers %>% filter(avg_log2FC > 0, p_val_adj < 0.05) %>% 
  arrange(-avg_log2FC) %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC) %>% data.frame()
markers$cluster %>% table()

markers$cluster <- markers$cluster %>% as.character()
markers$cluster[which(markers$cluster == "Endothelial Cell")] <- "Endothelial"
markers$cluster[which(markers$cluster == "Perivascular Fibroblast")] <- "PVFs"
markers$cluster <- factor(markers$cluster, levels = c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs"))

load("~/01_bulk/data/DESeq_dds.rdata")
res <- results(dds) %>% as.data.frame() %>% na.omit()
res$gene <- rownames(res)

res$threshold <- "Not"
res$threshold[res$log2FoldChange > 0] <- "Up"  
res$threshold[res$log2FoldChange < 0] <- "Down"
res <- res %>% dplyr::select("gene", "threshold")

merge_data <- merge(res, markers, by.x = "gene", by.y = "gene") %>% filter(threshold %in% c("Up", "Down"))
(res_df <- table(merge_data$cluster) %>% data.frame() %>% arrange(-Freq) %>% 
    {.$Var1 = factor(.$Var1, levels = .$Var1); .})

col1 <- c("#B8C6E9", "#DD9FD7", "#036C00", "#B0C399", "#06b95c", "#00738C", "#F39800", "#FFE4B5")
col1 <- setNames(col1, c("Astrocyte", "Endothelial", "Microglia", "mOli", "Neuron", "OPC", "Pericyte", "PVFs")); col1
(plot <- ggplot(res_df, aes(x = Var1, y = Freq / 100, fill = Var1, group = Var1)) +
    labs(x = "", y = "") +
    geom_bar(stat = "identity", width = .95) +
    scale_fill_manual(values = col1) +
    geom_text(aes(label = Freq / 100), vjust = 0.5, color = "black", size = 4) +
    ylab("Detection rate") +
    guides(fill = guide_legend(reverse = T)) +
    RotatedAxis() +
    scale_x_discrete(expand = expansion(mult = c(0.09, 0.09))))
ggsave(paste0(fwd, "/scRNA.signature.genes_Detection.rate_in.cfRNA.Count.pdf"), plot, width = 5, height = 4)