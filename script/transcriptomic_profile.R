library(ggrepel)
library(ggplot2)
library(pheatmap)
library(DESeq2)
library(apeglm)
library(readxl)
library(reshape)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

## ------- Inputs and parameters -------
meta_file="../data/meta_data/metadata_for_eytan-BIDMC_PCBN_UM.xlsx"
gene_anno_file='../data/transcript_gene_id_maps.xlsx'
diff_gene_file='../data/PCBN_tumor_DESeq2_res.txt'
PCBN_tumor_tpm_file='../data/PCBN_tumor_TPM.txt'
PCBN_tumor_pathways_file='../data/PCBN_tumor_pathways.txt'
expr_cor_PCBN_tumor_benign_file='../data/Expression_cor_PCBN_tumor_benign.txt'
expr_cor_PCBN_tumor_benign_pathway_file='../data/Expression_cor_PCBN_tumor_benign_pathway.txt'

outdir='../results/transcritomic_profile'
if(!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ------- Load and process count and meta data -------
## Load meta and count data 
meta_df=read_excel(meta_file)
meta_df=meta_df[,1:18]
meta_df$tissue_type=meta_df$`PRESUMED SAMPLE PHENOTYPE`
meta_df$tissue_type[which(meta_df$`SAMPLE TYPE`=='Blood')]='BLOOD'
meta_df$`SAMPLE ID`=gsub('\\-','.',meta_df$`SAMPLE ID`)
meta_df=meta_df[!is.na(meta_df$tissue_type),]
meta_df$Cohort=meta_df$`SAMPLE SOURCE`
meta_df$Cohort[which(meta_df$Cohort=='BIDMC')]='BIDMC_UM'
meta_df$Cohort[which(meta_df$Cohort=='UM')]='BIDMC_UM'

## ------- Functions -------
# Identify up and down genes
diff_rnaseq=function(dd,cutoffpvalue,cutofflogfc){
  dd$Significant <- ifelse(dd$padj < cutoffpvalue, "asig", "nosig")
  dd$group_rnaseq='other'
  dd$group_rnaseq[which(dd$padj < cutoffpvalue & dd$log2FoldChange>=cutofflogfc)] = 'up'
  dd$group_rnaseq[which(dd$padj < cutoffpvalue & dd$log2FoldChange<=(-cutofflogfc))] = 'down'
  return(dd)
}

# Volcano plot
volcano_plot=function(dd,genes=NULL,cutoffpvalue,cutofflogfc,filename,dir='.',width = 5,height = 4.5){
  #dd$Significant <- ifelse(dd$pvalue < 0.05, "p < 0.05", "p >= 0.05")
  dd$symbol=rownames(dd)
  dd=na.omit(dd)
  p=ggplot(dd, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = group_rnaseq)) +
    scale_color_manual(values = c("blue", "grey","red")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_text_repel(
      data = dd[genes,],
      aes(label = symbol),
      size = 5,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines"),
      max.overlaps = Inf
    )+
    geom_hline(yintercept=-log10(cutoffpvalue), linetype="dashed", color = "black")+
    geom_vline(xintercept = cutofflogfc, linetype="dotted", color = "black")+
    geom_vline(xintercept = -cutofflogfc, linetype="dotted", color = "black")
  
  print(p)
  ggsave(filename = paste0(dir,'/',filename,'.pdf'),p,width = width,height = height)
}


## ------- Volcano plot -------
## PCBN tumor, rec vs non-rec
b_cell_lable=c('IGKC','CD79A','TNFRSF13B','IGLC2','BLK','IGHM')
cell_cycle=c('TTK','MCM10','POLQ','UBE2C','CDK1','MKI67')
genes=c(b_cell_lable,cell_cycle)

diff_genes_df=read.delim(diff_gene_file)

# Identify up and down genes
diff_gene_group=diff_rnaseq(diff_genes_df,cutoffpvalue = 0.1,cutofflogfc = 0.585)
# Volcano plot
volcano_plot(diff_gene_group,cutoffpvalue = 0.1,cutofflogfc = 0.585, filename = 'PCBN_tumor_Volcano',dir=outdir,genes,width = 8,height = 8)

## ------- Heatmap, label key genes -------
## PCBN tumor 
# Load data
tpm=read.delim(PCBN_tumor_tpm_file)
# group info
groups <- PCBN_tumor_meta[,c('SAMPLE ID','RECUR')]
colnames(groups)=c('sample','group')
groups$group=ifelse(groups$group=='1','Rec','Non-Rec')
groups$sample <- as.character(groups$sample)
groups$group <- factor(groups$group, levels = unique(groups$group))
ordered_samples <- groups$sample[order(groups$group)]
tpm <- tpm[, ordered_samples]

# Row label logic (only label selected genes)
b_cell_label <- c("IGKC","CD79A","TNFRSF13B","IGLC2","BLK","IGHM")
cell_cycle   <- c("TTK","MCM10","POLQ","UBE2C","CDK1","MKI67")
selected_genes <- c(b_cell_label, cell_cycle)

# Differential genes (keep ALL DEGs)
diff_genes_df=diff_genes_df[which(diff_genes_df$padj<=0.1),]
# keep padj <= 0.1
deg_genes <- rownames(diff_genes_df)

# Extract TPM of ALL differential genes
heat_data <- tpm[rownames(tpm) %in% deg_genes, , drop = FALSE]

# Z-score normalization
heat_scaled <- t(scale(t(heat_data)))

# Gene annotation (rows)
row_anno_df <- data.frame(
  B_cell    = ifelse(rownames(heat_scaled) %in% b_cell_label, "B-cell", "Other"),
  CellCycle = ifelse(rownames(heat_scaled) %in% cell_cycle, "Cell-cycle", "Other")
)

row_ha <- rowAnnotation(
  Bcell = row_anno_df$B_cell,
  CellCycle = row_anno_df$CellCycle,
  col = list(
    Bcell = c("B-cell" = "steelblue", "Other" = "grey90"),
    CellCycle = c("Cell-cycle" = "firebrick", "Other" = "grey90")
  )
)

# Create the mark annotation
label_pos <- match(selected_genes, rownames(heat_scaled))
mark_ha <- rowAnnotation(
  mark = anno_mark(
    at = label_pos,
    labels = selected_genes,
    labels_gp = gpar(fontsize = 7),
    link_width = unit(5, "mm")   # length of leader lines
  )
)


# Sample group annotation (columns)
group_levels <- levels(groups$group)

if (length(group_levels) == 1) {
  group_colors <- c("gray60")
} else if (length(group_levels) == 2) {
  group_colors <- c("#66C2A5", "#FC8D62")  # Set2 first two colors
} else {
  group_colors <- brewer.pal(n = length(group_levels), "Set2")
}
names(group_colors) <- group_levels

col_ha <- HeatmapAnnotation(
  Group = groups$group[match(ordered_samples, groups$sample)],
  col = list(Group = group_colors),
  annotation_height = unit(4, "mm")
)

# Row label logic (only label selected genes)
selected_genes <- c(b_cell_label, cell_cycle)

row_labels <- ifelse(
  rownames(heat_scaled) %in% selected_genes,
  rownames(heat_scaled),
  ""   # suppress other labels
)

# Draw the heatmap
pdf(file.path(outdir,'heatmap_PCBN_tumor.pdf'),width = 6,height = 6)
Heatmap(
  heat_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick")),
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,    
  
  right_annotation = mark_ha,
  row_names_side = "left",
  
  top_annotation = col_ha,
  
  show_row_names = F,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 8),
  
  show_column_names = F,
  column_labels = as.character(groups$group[match(ordered_samples, groups$sample)]),
  column_names_gp = gpar(fontsize = 10)
)
dev.off()

## ------- Barplot if top B cell and Cell cycle pathways -------
## PCBN tumor
PCBN_tumor_pathways <- read.delim(PCBN_tumor_pathways_file)
PCBN_tumor_pathways_up <- PCBN_tumor_pathways[which(PCBN_tumor_pathways$Direction=='Up'),]
PCBN_tumor_pathways_Down <- PCBN_tumor_pathways[which(PCBN_tumor_pathways$Direction=='Down'),]
PCBN_tumor_pathways_up <- PCBN_tumor_pathways_up[order(PCBN_tumor_pathways_up$logq,decreasing = F),]
PCBN_tumor_pathways_Down <- PCBN_tumor_pathways_Down[order(PCBN_tumor_pathways_Down$logq,decreasing = F),]
PCBN_tumor_pathways_order <- rbind(PCBN_tumor_pathways_Down,PCBN_tumor_pathways_up)
PCBN_tumor_pathways_order$Description <- factor(PCBN_tumor_pathways_order$Description,levels = PCBN_tumor_pathways_order$Description)

ggplot(PCBN_tumor_pathways_order, aes(y = logq * ifelse(Direction == "Up", 1, -1), x = Description, fill = Direction)) +
  geom_col() +
  scale_fill_manual(values = c("Up" = "#f46d43", "Down" = "#4682b4")) +
  labs(
    y = expression(-log[10](qvalue)),
    x = "",
    title = ""
  ) +
  geom_text(aes(
    x = (Description),
    y = 0,
    label = Description,
    hjust = label_hjust
  ),
  size = 7 * 0.35,
  angle = 0)+
  geom_vline(xintercept = 0, color = "black") +
  theme_classic()+
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )+coord_flip()

ggsave(file.path(outdir,'PCBN_tumor_pathways.pdf'), width = 8,height = 6) 

## ------- Correlation between tumor and benign -------
## Barplot of enriched pathways, PCBN
enrich_cohort1_b_t=read.delim(expr_cor_PCBN_tumor_benign_pathway_file)
sel_path=c('immunoglobulin complex','antigen binding','adaptive immune response','leukocyte mediated immunity','B cell mediated immunity','antigen processing and presentation','immunoglobulin mediated immune response','regulation of leukocyte activation','leukocyte mediated cytotoxicity','positive regulation of lymphocyte activation','cell killing')

sel_df=enrich_cohort1_b_t[which(enrich_cohort1_b_t$Description%in%sel_path),]
sel_df$neglog10qvalue=-log10(sel_df$qvalue)
sel_df$Direction='Down'
sel_df$Description=str_to_title(sel_df$Description)
sel_df=sel_df[order(sel_df$qvalue,decreasing = T),]
sel_df$Description=factor(sel_df$Description,levels = sel_df$Description)

p<-ggplot(sel_df, aes(y = neglog10qvalue, x = Description, fill = Direction))  +
  geom_bar(stat="identity")+scale_fill_manual(values = c("Up" = "#f46d43", "Down" = "#4682b4")) +
  coord_flip()+theme_classic()+xlab('')+theme(legend.position = 'none')+ylab(expression(-log[10]('Adjusted Pvalue')))
p
ggsave(file.path(outdir,'expr_cor_PCBN_tumor_benign_path_barplot.pdf'), p, width = 5, height = 3)

## ------- Volcano with B-cell labels -------
## PCBN tumor vs benign
cohort1_b_t <- read.delim(expr_cor_PCBN_tumor_benign_file)

user_labels <- c("IGKC","CD79A","TNFRSF13B","IGLC2","BLK","IGHM")
# Canonical B-cell markers (extend as needed)
canonical <- c(
  "CD19","MS4A1","CD79A","CD79B","CD22","CD72","BLK","BANK1","PAX5","SPIB",
  "TNFRSF13B","TNFRSF13C","CR2","FCRL1","FCRL2","FCRL3","FCRL4","FCRL5",
  "IGHM","IGHD","IGKC","IGLC1","IGLC2","IGLC3","IGLL5","JCHAIN","CXCR5",
  "BACH2","BCL6","MEF2C"
)

bcell_markers <- unique(toupper(c(user_labels, canonical)))

cohort1_b_t <- cohort1_b_t %>%
  mutate(
    gene_upper = toupper(Genes),
    is_bcell   = gene_upper %in% bcell_markers
  )

# Choose labels: top-10 B-cell by p-value + specified present
top10_bcell <- cohort1_b_t %>%
  filter(is_bcell) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  pull(gene_upper)

user_present <- intersect(toupper(user_labels), unique(cohort1_b_t$gene_upper))
label_set <- unique(c(top10_bcell, user_present))

cohort1_b_t <- cohort1_b_t %>%
  mutate(to_label = gene_upper %in% label_set)

# Volcano plot with reduced overlaps ----------------------------
theme_set(theme_classic(base_size = 12))
base_col_other <- "grey75"
base_col_bcell <- "#1f77b4"     # blue for B-cell markers
label_col      <- "#0d3b66"     # dark navy for labels

df=cohort1_b_t
df$adj=p.adjust(df$pvalue,method = 'BH')
df$neglog10p_adjust=-log10(df$adj)
volcano <- ggplot(df, aes(x = cor, y = neglog10p_adjust)) +
  # Non-B-cell genes
  geom_point(
    data = ~ filter(.x, !is_bcell),
    color = base_col_other, size = 0.8, alpha = 0.55
  ) +
  # B-cell genes
  geom_point(
    data = ~ filter(.x, is_bcell),
    color = base_col_bcell, size = 1.4, alpha = 0.95
  ) +
  # Nominal p=0.05 reference (optional)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  labs(
    x = "Correlation coefficient (r)",
    y = expression(-log[10]~p),
    title = ""
  ) +
  # Tweak panel spacing and text
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

# Label points (reduced overlaps)
volcano_labeled <- volcano +
  geom_text_repel(
    data = df %>% filter(to_label),
    aes(label = Genes),
    size = 3.1,
    color = label_col,
    max.overlaps = Inf,          # practical since we limited the label set
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = label_col,
    segment.size = 0.3,
    segment.alpha = 0.7,
    force = 3,                   # increase repulsion
    force_pull = 2,
    max.time = 2,                # keep layout time bounded (seconds)
    nudge_y = 0, nudge_x = 0,    # pure repel
    seed = 123                   # reproducible layout
  )

ggsave(file.path(outdir,'volcano_bcell_deoverlap.pdf'), volcano_labeled, width = 3.5, height = 3.5)

## ------- Plot B cell gene correlation compared with background -------
## PCBN
sel_df_b_t_plot=cohort1_b_t[which(cohort1_b_t$Genes%in%label_set),c("Genes","cor","Diff_group")]

df_wide <- sel_df_b_t_plot %>%
  mutate(Diff_group = recode(Diff_group,
                             "Benign_Tumor"     = "matched",
                             "Benign_Tumor_BG"  = "background")) %>%
  pivot_wider(names_from = Diff_group, values_from = cor)

# Order genes by matched correlation (descending)
df_wide <- df_wide %>%
  arrange(desc(matched)) %>%
  mutate(Genes = factor(Genes, levels = rev(Genes)))  # reversed so top goes first

# Plot (segments + open/filled points) 
p <- ggplot(df_wide, aes(y = Genes)) +
  # connecting segment from background to matched for each gene
  geom_segment(aes(x = background, xend = matched, yend = Genes),
               linewidth = 0.6, color = "grey35") +
  # background: open circles
  geom_point(aes(x = background),
             shape = 21, size = 2.8, stroke = 0.7,
             fill = "white", color = "grey20") +
  # matched: filled circles
  geom_point(aes(x = matched),
             shape = 16, size = 3.2, color = "#b22222") +
  scale_x_continuous(limits = c(min(df_wide$background, df_wide$matched, na.rm = TRUE) - 0.05,
                                max(df_wide$background, df_wide$matched, na.rm = TRUE) + 0.05)) +
  labs(
    title = "Matched vs background correlation (Benign Tumor vs BG)",
    subtitle = "Filled = matched correlation; Open = background correlation",
    x = "Pearson correlation (R)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(face = "bold", size = 13),
    plot.subtitle= element_text(size = 11),
    axis.text.y  = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

p

ggsave(file.path(outdir,"B_cell_genes_matched_vs_background_corr.pdf"), p, width = 5, height = 4.5)


