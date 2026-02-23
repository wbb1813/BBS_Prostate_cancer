library(ggplot2)
library(tidyr)
## ========== Inputs and outputs ==========
# Read data, staining results
df <- read.delim("../data/Slides_Bcells_intensity.tsv", header = TRUE, sep = "\t")

# Normalize Bcell.intensity by Luminal.intensity
df$normalized_Bcell <- df$Bcell.intensity / df$Luminal.intensity
df$study <- ifelse(df$study=='15c-008', 'Non-recurrence', 'Recurrence')

# 15c-008 Non-recurrence
# Moonshot Recurrence

# B cell score 
benign_sig_score_file='../data/Signature_score_benign.txt'
benign_sig_score_df=read.delim(benign_sig_score_file)

outdir <- '../results/Bcell_slides'
if (!dir.exists(outdir)){
  dir.create(outdir,T)
}

# ============================================================
# 1. Correlation between Tumor and Benign in Non-recurrence cohort
# ============================================================
cohort_15c <- df[df$study == "Non-recurrence", ]

# Pivot to wide format: each patient (pt) has Benign and Tumor values
benign_vals <- cohort_15c[cohort_15c$benign.tumor == "Benign", c("pt", "normalized_Bcell")]
tumor_vals  <- cohort_15c[cohort_15c$benign.tumor == "Tumor",  c("pt", "normalized_Bcell")]
colnames(benign_vals) <- c("pt", "Benign")
colnames(tumor_vals)  <- c("pt", "Tumor")

paired_df <- merge(benign_vals, tumor_vals, by = "pt")

# Pearson correlation
cor_test <- cor.test(paired_df$Benign, paired_df$Tumor, method = "pearson")
cat("=== Correlation (Non-recurrence: Benign vs Tumor) ===\n")
print(cor_test)

# Correlation plot
p1 <- ggplot(paired_df, aes(x = Benign, y = Tumor)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "",
    x = "Benign (normalized_Bcell)",
    y = "Tumor (normalized_Bcell)",
    subtitle = paste0("r = ", round(cor_test$estimate, 3),
                      ", p = ", signif(cor_test$p.value, 3))
  ) +
  theme_classic()
p1
ggsave("correlation_15c008.pdf", p1, path = outdir, width = 4, height = 4)

# ============================================================
# 2. Compare Non-recurrence Benign with Recurrence (boxplot + Wilcox)
# ============================================================
benign_15c <- cohort_15c[cohort_15c$benign.tumor == "Benign", ]
benign_15c$group <- "Non-recurrence adjacent"

moonshot <- df[df$study == "Recurrence", ]
moonshot$group <- "Recurrence"

comp2_df <- rbind(
  benign_15c[, c("normalized_Bcell", "group")],
  moonshot[, c("normalized_Bcell", "group")]
)

wilcox2 <- wilcox.test(normalized_Bcell ~ group, data = comp2_df)
cat("\n=== Wilcoxon test: Non-recurrence Benign vs Moonshot ===\n")
print(wilcox2)

p2 <- ggplot(comp2_df, aes(x = group, y = normalized_Bcell, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 2) +
  labs(
    title = "",
    y = "normalized_Bcell",
    x ='',
    subtitle = paste0("Wilcoxon p = ", signif(wilcox2$p.value, 3))
  ) +
  scale_fill_brewer(palette="Dark2")+
  theme_classic() +
  theme(legend.position = "none")
p2
ggsave("Non_rec_benign_vs_rec.pdf", p2,path = outdir, width = 2.5, height = 4)

# ============================================================
# 3. Compare Non-recurrence Tumor with Recurrence (boxplot + Wilcox)
# ============================================================
tumor_15c <- cohort_15c[cohort_15c$benign.tumor == "Tumor", ]
tumor_15c$group <- "Non-recurrence tumor"

moonshot2 <- moonshot
moonshot2$group <- "Recurrence"

comp3_df <- rbind(
  tumor_15c[, c("normalized_Bcell", "group")],
  moonshot2[, c("normalized_Bcell", "group")]
)

wilcox3 <- wilcox.test(normalized_Bcell ~ group, data = comp3_df)
cat("\n=== Wilcoxon test: 15c-008 Tumor vs Moonshot ===\n")
print(wilcox3)

p3 <- ggplot(comp3_df, aes(x = group, y = normalized_Bcell, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  labs(
    title = "",
    y = "normalized_Bcell",
    x = 'Tumor',
    subtitle = paste0("Wilcoxon p = ", signif(wilcox3$p.value, 3))
  ) +
  scale_fill_brewer(palette="Dark2")+
  theme_classic() +
  theme(legend.position = "none")
p3
ggsave("Non_rec_tumor_vs_rec.pdf", p3,path = outdir, width = 3, height = 4)

# ============================================================
# 4. Matched paired comparison: Benign vs Tumor in 15c-008
# ============================================================
# Reshape paired_df (already created above) to long format for boxplot

paired_long <- pivot_longer(paired_df, cols = c("Benign", "Tumor"),
                            names_to = "type", values_to = "normalized_Bcell")
paired_long$type <- factor(paired_long$type, levels = c("Benign", "Tumor"))

# Paired Wilcoxon signed-rank test
wilcox4 <- wilcox.test(paired_df$Benign, paired_df$Tumor, paired = TRUE)
cat("\n=== Paired Wilcoxon test: 15c-008 Benign vs Tumor (matched) ===\n")
print(wilcox4)

p4 <- ggplot(paired_long, aes(x = type, y = normalized_Bcell, fill = type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
  geom_point(size = 2.5) +
  geom_line(aes(group = pt), linetype = "dashed", color = "grey40") +
  labs(
    title = "",
    y = "normalized_Bcell",
    x = "",
    subtitle = paste0("Paired Wilcoxon p = ", signif(wilcox4$p.value, 3))
  ) +
  scale_fill_brewer(palette="Set1")+
  theme_classic() +
  theme(legend.position = "none")
p4
ggsave("matched_benign_vs_tumor_15c008.pdf", p4, path = outdir, width = 3, height = 4)

cat("\nDone! PDFs saved.\n")



# ============================================================
# 5. youden index 
# ============================================================
df_benign=df[which(df$benign.tumor=='Benign'),]
scores <- as.numeric(df_benign$normalized_Bcell)
labels <- df_benign$study

# youden idex
youden_coords <- coords(roc_obj, "best", best.method = "youden", ret = "threshold")
cutoff <- as.numeric(youden_coords$threshold[1])

# Generate density plot
p <- ggplot(df_benign, aes_string(x = 'normalized_Bcell', fill = 'study')) +
  geom_density(alpha = 0.7) +
  scale_fill_brewer(palette="Dark2")+
  geom_vline(xintercept = as.numeric(cutoff), color = "red", linetype = "dashed") +
  labs(title = paste0(" (Youden cutoff = ", round(cutoff, 3), ")"),
       x = "Score", y = "Density") +
  theme_classic()
p
ggsave('B_cell_density_youden.pdf', p, path = outdir,width = 6,height = 4)

## Confusion matrix
df_benign$obs=ifelse(df_benign$study=='Non-recurrence','yes','no')
colnames(df_benign)[which(colnames(df_benign)=='normalized_Bcell')]='score'
df2 <- df_benign %>%
  mutate(
    truth = ifelse(obs == "yes", 1, 0),
    pred  = ifelse(score >= cutoff, 1, 0),
    truth_lab = ifelse(truth == 1, "Truth: yes", "Truth: no"),
    pred_lab  = ifelse(pred  == 1, "Pred: yes",  "Pred: no"),
    type = case_when(
      truth == 1 & pred == 1 ~ "TP",
      truth == 0 & pred == 0 ~ "TN",
      truth == 1 & pred == 0 ~ "FN",
      truth == 0 & pred == 1 ~ "FP"
    )
  )

conf_df <- df2 %>%
  group_by(truth_lab, pred_lab, type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(label = paste0(type, "\n", n))
conf_df$truth_lab <- factor(conf_df$truth_lab, levels = c("Truth: yes", "Truth: no"))
conf_df$pred_lab  <- factor(conf_df$pred_lab,  levels = c("Pred: yes", "Pred: no"))

# Ensure these are characters (in case they are factors)
conf_df <- conf_df %>%
  mutate(
    truth_lab = as.character(truth_lab),
    pred_lab  = as.character(pred_lab)
  )

# Define the complete 2×2 grid and their types
all_cells <- tribble(
  ~truth_lab,   ~pred_lab,    ~type,
  "Truth: no",  "Pred: no",   "TN",
  "Truth: no",  "Pred: yes",  "FP",
  "Truth: yes", "Pred: no",   "FN",
  "Truth: yes", "Pred: yes",  "TP"
)

conf_df2 <- all_cells %>%
  left_join(conf_df, by = c("truth_lab","pred_lab","type")) %>%
  mutate(
    n = coalesce(n, 0L),
    label = if_else(is.na(label), paste0(type, "\n", n), label)
  ) %>%
  mutate(type = factor(type, levels = c("TN","FP","FN","TP"))) %>%
  arrange(type)

p = ggplot(conf_df2, aes(x = pred_lab, y = truth_lab, fill = n)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label), size = 6) +
  scale_fill_gradient(low = "#deebf7", high = "#08519c") +  # scientific blues
  coord_equal() +
  labs(
    x = "Predicted",
    y = "True",
    fill = "Count"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.line.x = element_blank(),   # remove x-axis line
    axis.line.y = element_blank(),   # remove y-axis line
    axis.ticks = element_blank(),    # optional: remove ticks
    panel.border = element_blank(),  # remove the outer box
    panel.grid = element_blank(),    # optional: remove grid lines
    plot.title = element_text(hjust = 0.5)
  )
p
ggsave('B_cell_density_confusion_matrix.pdf',path = outdir,width = 5,height = 4.5)

df_benign$normalized_Bcell=df_benign$Bcell.intensity/df_benign$Luminal.intensity
median(df_benign$normalized_Bcell[which(df_benign$study=='Non-recurrence')])
median(df_benign$normalized_Bcell[which(df_benign$study=='Recurrence')])
# ============================================================
# 6. Correlation between B cell density and B cell score 
# ============================================================
benign_sig_score_df$pt=as.integer(sub(".*\\.", "", benign_sig_score_df$Patient))

df_score_density=merge(benign_sig_score_df[,c('pt','B_cell_33520406')],paired_df,by='pt')

cor.test(df_score_density$Benign,df_score_density$B_cell_33520406)