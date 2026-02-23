# ===============================
# Prostate clinical scores + AUC/PRAUC evaluation
# Scores evaluated: CAPRA-S, D'Amico (ordinal), PSA, Gleason total
# ===============================
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(precrec)
library(tibble)
library(reshape)
library(ggpubr)
library(VennDiagram)
library(grid)     
library(gridExtra)
library(UpSetR)
library(patchwork)
# ---------- OUTPUT DIR ----------
# outdir <- "/data/Binbin_Kun/binbin/adam/results/recurency_prediction/clincal_score"
outdir <- "../results/clincal_score_prediction/"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---------- Meta file ----------
# infile <- "/data/Binbin_Kun/Flex_2025/data/meta_data/metadata_for_eytan-BIDMC_PCBN_UM.xlsx"
infile <- "../data/metadata_for_eytan-BIDMC_PCBN_UM.xlsx"
sheet  <- "metadata"

# Column mapping in your file
COL_PATIENT     <- "SAMPLE ID"
COL_AGE         <- "AGE AT RP (YRS)"
COL_PSA_PREOP   <- "PRE-OP PSA"
COL_STAGE_CLIN  <- "PATIENT STAGE"
COL_ISUP_GRADE   <- "PATIENT ISUP GRADE"
COL_GS_BIOPSY   <- "PATIENT GLEASON SCORE"  # e.g., "3+4" or "3+4 T5"

# Optional biopsy cores (for CAPRA strict % positive cores)
COL_CORES_POS   <- NA_character_            # e.g., "POS_CORES"
COL_CORES_TOTAL <- NA_character_            # e.g., "TOTAL_CORES"

# Optional post-op fields for CAPRA-S (compute only if present)
COL_GS_PATH     <- "SAMPLE GLEASON SCORE"
COL_MARGIN_POS  <- NA_character_
COL_EPE         <- NA_character_
COL_SVI         <- NA_character_
COL_LNI         <- NA_character_

# ---------- Functions ----------
to_num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]+", "", as.character(x))))

parse_gleason_pair <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x)) return(c(NA_integer_, NA_integer_, NA_integer_))
  s <- trimws(as.character(x))
  if (s == "" || toupper(s) %in% c("NA","N/A","UNKNOWN","UNK")) return(c(NA,NA,NA))
  m <- regexpr("(\\d)\\s*\\+\\s*(\\d)", s, perl = TRUE)
  if (is.na(m[1]) || m[1] == -1) return(c(NA,NA,NA))
  gs <- regmatches(s, m)
  parts <- strsplit(gs, "\\+")[[1]]
  parts <- trimws(parts)
  p <- suppressWarnings(as.integer(parts[1]))
  q <- suppressWarnings(as.integer(parts[2]))
  if (is.na(p) || is.na(q)) return(c(NA,NA,NA))
  c(p, q, p + q)
}

gleason_to_gradegroup <- function(primary, secondary, total = NA_integer_) {
  if (is.na(primary) || is.na(secondary)) {
    if (!is.na(total)) {
      if (total <= 6) return(1L)
      if (total == 7) return(3L)  # conservative default if pattern unknown
      if (total == 8) return(4L)
      if (total >= 9) return(5L)
    }
    return(NA_integer_)
  }
  t <- primary + secondary
  if (t <= 6) return(1L)
  if (primary == 3 && secondary == 4) return(2L)
  if (primary == 4 && secondary == 3) return(3L)
  if (t == 8) return(4L)
  if (t >= 9) return(5L)
  NA_integer_
}

map_stage <- function(stage_chr) {
  s <- toupper(gsub("\\s+", "", as.character(stage_chr)))
  s <- gsub("^T0$", "T1", s)
  s
}

yn_to01 <- function(x) {
  if (is.na(x)) return(NA_real_)
  xv <- toupper(as.character(x))
  if (xv %in% c("1","TRUE","YES","POS","POSITIVE")) return(1)
  if (xv %in% c("0","FALSE","NO","NEG","NEGATIVE")) return(0)
  suppressWarnings(as.numeric(x))  # fallback
}

# ---- D'Amico ----
damico_one <- function(psa, stage, gg) {
  stage <- map_stage(stage)
  stage_low  <- stage %in% c("T1","T1A","T1B","T1C","T2A")
  stage_int  <- stage %in% c("T2B")
  stage_high <- stage %in% c("T2C","T3","T3A","T3B","T4")
  psa_low  <- !is.na(psa) && psa <= 10
  psa_int  <- !is.na(psa) && psa > 10 && psa <= 20
  psa_high <- !is.na(psa) && psa > 20
  g_low  <- !is.na(gg) && gg == 1L
  g_int  <- !is.na(gg) && gg %in% c(2L,3L)
  g_high <- !is.na(gg) && gg %in% c(4L,5L)
  if (stage_low && psa_low && g_low) return("Low")
  if (stage_high || psa_high || g_high) return("High")
  if (stage_int || psa_int || g_int)   return("Intermediate")
  NA_character_
}

# ---- CAPRA (UCSF) ----
capra_points_components <- function(psa, stage, gg, pct_pos_cores, age_years) {
  psa_pts <- if (is.na(psa)) NA_real_ else if (psa <= 6) 0 else if (psa <= 10) 1 else if (psa <= 20) 2 else if (psa <= 30) 3 else 4
  g_pts   <- if (is.na(gg)) NA_real_ else if (gg == 1) 0 else if (gg == 2) 1 else if (gg == 3) 2 else 3
  s <- map_stage(stage)
  stg_pts <- ifelse(s %in% c("T1","T1A","T1B","T1C","T2A"), 0,
                    ifelse(s %in% c("T2B","T2C","T3","T3A","T3B","T4"), 1, NA_real_))
  cores_pts <- if (is.na(pct_pos_cores)) NA_real_ else if (pct_pos_cores <= 33) 0 else 1
  age_pts <- if (is.na(age_years)) NA_real_ else if (age_years <= 50) 0 else 1
  c(psa = psa_pts, gleason = g_pts, stage = stg_pts, cores = cores_pts, age = age_pts)
}
capra_total_strict <- function(psa, stage, gg, pct_pos_cores, age_years) {
  pts <- capra_points_components(psa, stage, gg, pct_pos_cores, age_years)
  if (any(is.na(pts))) return(NA_integer_)
  as.integer(sum(pts))
}
capra_total_no_cores <- function(psa, stage, gg, age_years) {
  pts <- capra_points_components(psa, stage, gg, NA_real_, age_years)
  if (any(is.na(pts[c("psa","gleason","stage","age")]))) return(NA_integer_)
  as.integer(sum(pts[c("psa","gleason","stage","age")]))
}

# ===================== LOAD & SCORE =====================
dat <- read_excel(infile, sheet = sheet)
dat$`PATIENT ISUP GRADE`[which(dat$`PATIENT ISUP GRADE`=='3+5')]='5'

# Patient-level DF for scoring outputs (not strictly needed for AUC/PR)
keep_cols <- unique(na.omit(c(
  COL_PATIENT, COL_AGE, COL_PSA_PREOP, COL_STAGE_CLIN, COL_ISUP_GRADE, COL_GS_BIOPSY,
  COL_CORES_POS, COL_CORES_TOTAL, COL_GS_PATH, COL_MARGIN_POS, COL_EPE, COL_SVI, COL_LNI
)))
df <- unique(dat[keep_cols])

# Parse biopsy GS
if (!is.na(COL_GS_BIOPSY) && COL_GS_BIOPSY %in% names(df)) {
  gs_list <- lapply(df[[COL_GS_BIOPSY]], parse_gleason_pair)
  gsb <- do.call(rbind, gs_list)
  colnames(gsb) <- c("bx_g_primary","bx_g_secondary","bx_g_total")
  df <- cbind(df, gsb)
} else {
  df$bx_g_primary <- df$bx_g_secondary <- df$bx_g_total <- NA_integer_
}

# Parse path GS (optional)
if (!is.na(COL_GS_PATH) && COL_GS_PATH %in% names(df)) {
  gs_list_p <- lapply(df[[COL_GS_PATH]], parse_gleason_pair)
  gsp <- do.call(rbind, gs_list_p)
  colnames(gsp) <- c("path_g_primary","path_g_secondary","path_g_total")
  df <- cbind(df, gsp)
} else {
  df$path_g_primary <- df$path_g_secondary <- df$path_g_total <- NA_integer_
}

# Numerics
df$psa_preop_num <- to_num(df[[COL_PSA_PREOP]])
df$age_num       <- to_num(df[[COL_AGE]])

# % positive cores (optional)
if (!is.na(COL_CORES_POS) && !is.na(COL_CORES_TOTAL) &&
    COL_CORES_POS %in% names(df) && COL_CORES_TOTAL %in% names(df)) {
  pos <- to_num(df[[COL_CORES_POS]])
  tot <- to_num(df[[COL_CORES_TOTAL]])
  df$pct_pos_cores <- ifelse(!is.na(pos) & !is.na(tot) & tot > 0, 100*pos/tot, NA_real_)
} else {
  df$pct_pos_cores <- NA_real_
}

# Grade groups
df$gg_bx   <- mapply(gleason_to_gradegroup, df$bx_g_primary, df$bx_g_secondary, df$bx_g_total)
df$gg_path <- mapply(gleason_to_gradegroup, df$path_g_primary, df$path_g_secondary, df$path_g_total)

# D'Amico (categorical)
df$Damico_Risk <- mapply(damico_one, psa = df$psa_preop_num,
                         stage = df[[COL_STAGE_CLIN]], gg = df$gg_bx)

# CAPRA strict + no-cores + preferred CAPRA
df$CAPRA_strict <- mapply(capra_total_strict,
                          psa = df$psa_preop_num,
                          stage = df[[COL_STAGE_CLIN]],
                          gg = df$gg_bx,
                          pct_pos_cores = df$pct_pos_cores,
                          age_years = df$age_num)

df$CAPRA_no_cores <- mapply(capra_total_no_cores,
                            psa = df$psa_preop_num,
                            stage = df[[COL_STAGE_CLIN]],
                            gg = df$gg_bx,
                            age_years = df$age_num)

df$CAPRA <- ifelse(!is.na(df$CAPRA_strict), df$CAPRA_strict, df$CAPRA_no_cores)

# -------- write scored CSV (optional)
cols_to_keep <- c(
  COL_PATIENT, COL_AGE, COL_PSA_PREOP, COL_STAGE_CLIN, COL_GS_BIOPSY,
  "bx_g_primary","bx_g_secondary","bx_g_total",
  "Damico_Risk","CAPRA_strict","CAPRA_no_cores","CAPRA",
  COL_GS_PATH,"path_g_primary","path_g_secondary","path_g_total",
  COL_CORES_POS, COL_CORES_TOTAL, "pct_pos_cores",
  "CAPRAS_Total"
)
cols_to_keep <- unique(cols_to_keep[!is.na(cols_to_keep) & cols_to_keep %in% colnames(df)])
scored <- df[, cols_to_keep, drop = FALSE]
write.csv(scored, file.path(outdir,"scored_prostate_risk.csv"), row.names = FALSE)

# ===================== AUC / PRAUC (cohort × tissue) =====================
# Build sample-level eval DF from original 'dat' (aligns with RECUR/cohort/tissue)
eval_df <- dat %>%
  mutate(
    RECUR  = to_num(RECUR),
    PSA    = to_num(!!sym(COL_PSA_PREOP)),
    Age    = to_num(!!sym(COL_AGE)),
    ISUP_GRADE = to_num(`PATIENT ISUP GRADE`),
    Cohort = `SAMPLE SOURCE`,
    Tissue = `PRESUMED SAMPLE PHENOTYPE`
  )

# Replace missing tissue with "blood"
eval_df$Tissue[eval_df$`SAMPLE TYPE`=='Blood'] <- "blood"

# Parse biopsy Gleason & Grade Group; map stage
gs_list_all <- lapply(eval_df[[COL_GS_BIOPSY]], parse_gleason_pair)
gs_mat_all  <- do.call(rbind, gs_list_all)
colnames(gs_mat_all) <- c("bx_p","bx_s","Gleason_Total")
eval_df <- bind_cols(eval_df, as.data.frame(gs_mat_all)) %>%
  mutate(
    GradeGroup = mapply(gleason_to_gradegroup, bx_p, bx_s, Gleason_Total),
    StageMap   = map_stage(!!sym(COL_STAGE_CLIN))
  )

# CAPRA (strict prefers, else no-cores) at sample level
# Try to compute %pos cores if columns exist in original dat
pct_pos_sample <- if (!is.na(COL_CORES_POS) && !is.na(COL_CORES_TOTAL) &&
                      COL_CORES_POS %in% names(dat) && COL_CORES_TOTAL %in% names(dat)) {
  100 * to_num(dat[[COL_CORES_POS]]) / pmax(1, to_num(dat[[COL_CORES_TOTAL]]))
} else { rep(NA_real_, nrow(eval_df)) }

eval_df$CAPRA_strict   <- mapply(capra_total_strict, eval_df$PSA, eval_df$StageMap, eval_df$GradeGroup, pct_pos_sample, eval_df$Age)
eval_df$CAPRA_no_cores <- mapply(capra_total_no_cores, eval_df$PSA, eval_df$StageMap, eval_df$GradeGroup, eval_df$Age)
eval_df$CAPRA          <- ifelse(!is.na(eval_df$CAPRA_strict), eval_df$CAPRA_strict, eval_df$CAPRA_no_cores)

# D'Amico numeric (Low=0, Intermediate=1, High=2)
eval_df$DAmico_cat <- mapply(damico_one, eval_df$PSA, eval_df$StageMap, eval_df$GradeGroup)
eval_df$DAmico_num <- case_when(
  eval_df$DAmico_cat == "Low"          ~ 0,
  eval_df$DAmico_cat == "Intermediate" ~ 1,
  eval_df$DAmico_cat == "High"         ~ 2,
  TRUE ~ as.numeric(NA)
)

## save evaluation score
write.csv(eval_df, file.path(outdir,"evaluation_clinical_score.csv"), row.names = FALSE)

## ------- Plot AUC and PRAUC by merging cohort 2 and 3 -------
eval_df_merged=eval_df
eval_df_merged=eval_df_merged[!is.na(eval_df_merged$Cohort),]

eval_df_merged$Cohort[which(eval_df_merged$Cohort%in%c('BIDMC','UM'))]='BIDMC_UM'

# Iterate groups
metrics <- list()
keys <- eval_df_merged %>% distinct(Cohort, Tissue)
for (i in seq_len(nrow(keys))) {
  cohort_i <- keys$Cohort[i]
  tissue_i <- keys$Tissue[i]
  sub <- eval_df_merged %>% filter(Cohort == cohort_i, Tissue == tissue_i)
  
  # need both classes and enough rows
  if (sum(!is.na(sub$RECUR)) < 3 || length(unique(na.omit(sub$RECUR))) < 2) next
  
  # Models to compare: CAPRA, D'Amico (ordinal), PSA, Gleason_Total
  sc_list <- list(
    CAPRA          = sub$CAPRA,
    `D'Amico`      = sub$DAmico_num,
    PSA            = sub$PSA,
    Gleason_Total  = sub$Gleason_Total
  )
  
  mm <- mmdata(scores = sc_list, labels = sub$RECUR, modnames = names(sc_list), na_worst = FALSE)
  ev <- evalmod(mm)
  
  # AUC table (NOTE: precrec columns are modnames + aucs)
  auc_tbl <- precrec::auc(ev) %>%
    tibble::as_tibble() %>%                                   # cols: curvetypes, modnames, dsids, aucs
    dplyr::select(Score = modnames, curvetypes, auc = aucs) %>%  # <-- rename here
    dplyr::group_by(Score, curvetypes) %>%
    dplyr::summarise(auc = base::mean(auc), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = curvetypes, values_from = auc) %>%
    dplyr::rename(AUC = ROC, PRAUC = PRC) %>%
    dplyr::mutate(Cohort = cohort_i, Tissue = tissue_i) %>%
    dplyr::relocate(Cohort, Tissue, Score, AUC, PRAUC)
  
  metrics[[length(metrics)+1]] <- auc_tbl
}

# Save metrics
metrics_df <- bind_rows(metrics) %>% arrange(Cohort, Tissue, Score)
write.csv(metrics_df, file.path(outdir, "metrics_by_group_merged_cohort.csv"), row.names = FALSE)

## ------- Plot AUC and PRAUC merged -------
## filter UM cohort due to insuficient non-recurrence samples 
metrics_df_BP_tumor=metrics_df[which(metrics_df$Tissue=='TUMOR'),]
metrics_df_BP_tumor$Cohort=factor(metrics_df_BP_tumor$Cohort,levels = c('PCBN','BIDMC_UM'))

# AUC
ggplot(metrics_df_BP_tumor, aes(x = Score, y = AUC)) +
  geom_bar(stat = "identity",fill='#89CFF0',width = 0.75) +
  facet_wrap(~Cohort, scales = "free_y") +
  #geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", AUC)),
            position = position_dodge(width = 0.7), vjust = -0.25, size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(title = "",y = "AUC",x='') 
ggsave(file.path(outdir, "AUC_merged_cohort.pdf"), width = 8,height = 4.5)

# PRAUC
ggplot(metrics_df_BP_tumor, aes(x = Score, y = PRAUC)) +
  geom_bar(stat = "identity",fill='#FF7F50',width = 0.75) +
  facet_wrap(~Cohort, scales = "free_y") +
  #geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", PRAUC)),
            position = position_dodge(width = 0.7), vjust = -0.25, size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(title = "",y = "PRAUC",x='') 
ggsave(file.path(outdir, "PRAUC_merged_cohort.pdf"), width = 8,height = 4.5)

## ------- Visualize score distribution between recurrence and non-recurrence -------
eval_df_merged_tumor=eval_df_merged[which(eval_df_merged$`PRESUMED SAMPLE PHENOTYPE`=='TUMOR'),]
eval_df_merged_tumor=eval_df_merged_tumor[!is.na(eval_df_merged_tumor$`SAMPLE ID`),]

df_merged_sel=eval_df_merged_tumor[,c("SAMPLE ID",'Cohort','RECUR',"PSA","Gleason_Total","ISUP_GRADE")]
df_merged_sel$RECUR=ifelse(df_merged_sel$RECUR==1,'Rec','Non-rec')
df_merged_sel=as.data.frame(df_merged_sel)
df_merged_sel=melt(df_merged_sel,id.vars = c("SAMPLE ID",'Cohort','RECUR'))

for (i in unique(df_merged_sel$Cohort)){
  tmp_cohort=df_merged_sel[which(df_merged_sel$Cohort==i),]
  for (j in unique(tmp_cohort$variable)){
    tmp_cohort_var=tmp_cohort[which(tmp_cohort$variable==j),]
    if (j =='PSA'){
      tmp_cohort_var$value=log2(tmp_cohort_var$value+1)
      y_lab=paste0('Log2(PSA+1)')
    }else{
      y_lab=j
    }
    
    p=ggplot(tmp_cohort_var, aes(x = RECUR, y = value, fill = RECUR)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, color = 'black') +
      scale_fill_brewer(palette="Dark2")+
      labs(title = i, y = y_lab,x='') +
      theme_classic()+stat_compare_means(label = "p.format",method = 'wilcox',label.x = 1.5)+
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")  # centered
      )
    ggsave(filename = paste0("feature_plot_", i,'_',j, ".pdf"), plot = p,path=outdir, width = 3.5, height = 4)
  }
}

## ------- Plot of matched samples of benign, and tumor -------
# PCBN67 patient has two benign samples 
raw_meta=dat
raw_meta=raw_meta[!is.na(raw_meta$`PRESUMED SAMPLE PHENOTYPE`),]
raw_meta$`PRESUMED SAMPLE PHENOTYPE`[which(raw_meta$`PRESUMED SAMPLE PHENOTYPE`=="N/A")]='BLOOD'

meta_res=raw_meta[!duplicated(raw_meta$PATIENT),]
meta_res$Response_Status=ifelse(meta_res$RECUR=='1','Recurrence','Non-recurrence')

## Function
# Pie chat
plot_single_response_pie <- function(group_name, group_ids, all_data) {
  
  current_set_df <- data.frame(PATIENT = group_ids)
  plot_data <- current_set_df %>%
    dplyr::left_join(all_data, by = "PATIENT")
  
  response_summary <- plot_data %>%
    dplyr::group_by(Response_Status) %>%
    dplyr::summarise(Count = n()) %>%
    dplyr::mutate(
      Percentage = Count / sum(Count),
      y_mid = cumsum(Count) - Count / 2, # Angular position
    )
  
  p <- ggplot2::ggplot(response_summary, ggplot2::aes(x = "", y = Count, fill = Response_Status)) +
    ggplot2::geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.5) +
    
    ggplot2::geom_text(
      ggplot2::aes( y = y_mid, label = Count, color = 'black'), 
      size = 4, 
      show.legend = FALSE, # Hide the aesthetic legend for text color
      position = ggplot2::position_nudge(x = -0.1)
    ) +
    
    ggplot2::scale_fill_brewer(palette = "Dark2") +
    # Manually assign the text colors so the geom_text color scale doesn't interfere with the fill scale
    ggplot2::scale_color_identity() + 
    
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme_void(base_size = 14) + 
    ggplot2::labs(title = paste("Sample Size (N=", sum(response_summary$Count), "):", group_name), fill = "Response") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      legend.position = "bottom"
    )
  return(p)
}

## BIDMC
tmp_meta=raw_meta[which(raw_meta$`SAMPLE SOURCE`=='BIDMC'),]

tmp_tumor=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='TUMOR'),]
tmp_benign=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='BENIGN'),]
#tmp_blood=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='BLOOD'),]

data_list <- list(
  "TUMOR" = tmp_tumor$PATIENT,
  "BENIGN" = tmp_benign$PATIENT
)

# Individual PIE plots
p_TUMOR_pie <- plot_single_response_pie("TUMOR", data_list$TUMOR, meta_res)
p_BENIGN_pie <- plot_single_response_pie("BENIGN", data_list$BENIGN, meta_res)

ggsave('BIDMC_TUMOR_pie.pdf',p_TUMOR_pie,width = 4,height = 3,path = outdir)
ggsave('BIDMC_BENIGN_pie.pdf',p_BENIGN_pie,width = 4,height = 3,path = outdir)

# upset plot
main_colors <- c(
  bars = "#0072B2",           # set size bars
  intersections = "#0072B2"   # main intersection bars
)

main_colors <- c(bars = "#B05C2F", intersections = "#1F9E89")
main_colors <- c(bars = "#4E79A7", intersections = "#E15759")
main_colors <- c(bars = "#2E7D32", intersections = "#C2185B")
ain_colors <- c(bars = "#1F9E89", intersections = "#B05C2F")

tiff(file.path(outdir, 'BIDMC_upset.tiff'), width = 3, height = 4, units = "in", res = 600)
upset(
  fromList(data_list),
  sets.bar.color = main_colors["bars"],
  order.by = "freq",
  main.bar.color = main_colors["intersections"],
  sets = names(data_list),
  # REMOVE attribute.plots to keep the figure clean
  
  # Ensure the intersection matrix is square/readable
  point.size = 3.5,
  line.size = 1.5,
  mb.ratio = c(0.65, 0.35) # Adjusts main bar vs set size bar ratio
)
dev.off()

## PCBN
tmp_meta=raw_meta[which(raw_meta$`SAMPLE SOURCE`=='PCBN'),]

tmp_tumor=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='TUMOR'),]
tmp_benign=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='BENIGN'),]

data_list <- list(
  "TUMOR" = tmp_tumor$PATIENT,
  "BENIGN" = tmp_benign$PATIENT
)

# Individual PIE plots
p_TUMOR_pie <- plot_single_response_pie("TUMOR", data_list$TUMOR, meta_res)
p_BENIGN_pie <- plot_single_response_pie("BENIGN", data_list$BENIGN, meta_res)

ggsave('PCBN_TUMOR_pie.pdf',p_TUMOR_pie,width = 4,height = 3,path = outdir)
ggsave('PCBN_BENIGN_pie.pdf',p_BENIGN_pie,width = 4,height = 3,path = outdir)

# upset plot
tiff(file.path(outdir, 'PCBN_upset.tiff'), width = 3, height = 4, units = "in", res = 600)
upset(
  fromList(data_list),
  sets.bar.color = main_colors["bars"],
  order.by = "freq",
  main.bar.color = main_colors["intersections"],
  sets = names(data_list),
  # REMOVE attribute.plots to keep the figure clean
  
  # Ensure the intersection matrix is square/readable
  point.size = 3.5,
  line.size = 1.5,
  mb.ratio = c(0.65, 0.35) # Adjusts main bar vs set size bar ratio
)
dev.off()

## UM
tmp_meta=raw_meta[which(raw_meta$`SAMPLE SOURCE`=='UM'),]

tmp_tumor=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='TUMOR'),]
tmp_benign=tmp_meta[which(tmp_meta$`PRESUMED SAMPLE PHENOTYPE`=='BENIGN'),]

data_list <- list(
  "TUMOR" = tmp_tumor$PATIENT,
  "BENIGN" = tmp_benign$PATIENT
)

# Individual PIE plots
p_TUMOR_pie <- plot_single_response_pie("TUMOR", data_list$TUMOR, meta_res)
p_BENIGN_pie <- plot_single_response_pie("BENIGN", data_list$BENIGN, meta_res)

ggsave('UM_TUMOR_pie.pdf',p_TUMOR_pie,width = 4,height = 3,path = outdir)
ggsave('UM_BENIGN_pie.pdf',p_BENIGN_pie,width = 4,height = 3,path = outdir)

# upset plot
tiff(file.path(outdir, 'UM_upset.tiff'), width = 3, height = 4, units = "in", res = 600)
upset(
  fromList(data_list),
  sets.bar.color = main_colors["bars"],
  order.by = "freq",
  main.bar.color = main_colors["intersections"],
  sets = names(data_list),
  # REMOVE attribute.plots to keep the figure clean
  
  # Ensure the intersection matrix is square/readable
  point.size = 3.5,
  line.size = 1.5,
  mb.ratio = c(0.65, 0.35) # Adjusts main bar vs set size bar ratio
)
dev.off()







## ------- Cohort Overview -------
# Set theme 
theme_set(theme_minimal(base_size = 10, base_family = "Arial"))

# Define color palette
colors <- list(
  non_recurrence = "#3498db",
  recurrence = "#e74c3c",
  discovery = "#27ae60",
  validation = "#8e44ad",
  matched = "#16a085",
  tumor = "#f1c40f",
  benign = "#95a5a6",
  rna_seq = "#0984e3",
  arrow = "#7f8c8d"
)

# ============================================================================
# COHORT DATA
# ============================================================================
cohorts <- data.frame(
  cohort = c("PCBN", "BIDMC", "UM"),
  n = c(123, 84, 36),
  non_recurrence = c(80, 73, 1),
  recurrence = c(43, 11, 35)
)
cohorts$matched_samples <- cohorts$n * 2

final_cohorts <- data.frame(
  cohort = c("PCBN", "BM"),
  label = c("Discovery", "Validation"),
  n = c(123, 120),
  non_recurrence = c(80, 74),
  recurrence = c(43, 46),
  color = c(colors$discovery, colors$validation)
)
final_cohorts$matched_samples <- final_cohorts$n * 2

# ============================================================================
# FUNCTION: Create cohort box with stacked bar
# ============================================================================
create_cohort_box <- function(name, n, nr, r, matched, border_color = "grey30", 
                              show_label = NULL, label_color = NULL) {
  
  total <- nr + r
  nr_pct <- nr / total
  r_pct <- r / total
  
  # Create data for stacked bar
  bar_data <- data.frame(
    category = factor(c("Non-recurrence", "Recurrence"), 
                      levels = c("Non-recurrence", "Recurrence")),
    count = c(nr, r),
    pct = c(nr_pct, r_pct)
  )
  
  p <- ggplot() +
    # Background rectangle
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "white", color = border_color, linewidth = 1) +
    
    # Cohort label badge (if provided)
    {if (!is.null(show_label)) 
      annotate("label", x = 0.5, y = 0.92, label = toupper(show_label),
               fill = label_color, color = "white", fontface = "bold",
               size = 2.5, label.padding = unit(0.15, "lines"),
               label.r = unit(0.1, "lines"))
    } +
    
    # Cohort name
    annotate("text", x = 0.5, y = ifelse(is.null(show_label), 0.85, 0.78), 
             label = name, fontface = "bold", size = 4.5, color = 
               ifelse(is.null(label_color), "grey20", label_color)) +
    
    # Patient count
    annotate("text", x = 0.5, y = ifelse(is.null(show_label), 0.72, 0.65), 
             label = paste0("n = ", n, " patients"), size = 3, color = "grey40") +
    
    # Matched samples count
    annotate("text", x = 0.5, y = ifelse(is.null(show_label), 0.60, 0.53), 
             label = paste0(matched, " matched samples"), 
             size = 2.8, color = colors$matched, fontface = "bold") +
    
    # Stacked bar - Non-recurrence
    annotate("rect", xmin = 0.1, xmax = 0.1 + 0.8 * nr_pct, 
             ymin = 0.28, ymax = 0.42,
             fill = colors$non_recurrence, color = NA) +
    
    # Stacked bar - Recurrence
    annotate("rect", xmin = 0.1 + 0.8 * nr_pct, xmax = 0.9, 
             ymin = 0.28, ymax = 0.42,
             fill = colors$recurrence, color = NA) +
    
    # Bar labels
    annotate("text", x = 0.1 + 0.8 * nr_pct / 2, y = 0.35, 
             label = nr, color = "white", fontface = "bold", size = 3) +
    annotate("text", x = 0.1 + 0.8 * nr_pct + 0.8 * r_pct / 2, y = 0.35, 
             label = r, color = "white", fontface = "bold", size = 3) +
    
    coord_fixed(ratio = 1.2) +
    xlim(-0.05, 1.05) +
    ylim(0, 1) +
    theme_void()
  
  return(p)
}

# ============================================================================
# CREATE SOURCE COHORT BOXES
# ============================================================================
p_pcbn <- create_cohort_box("PCBN", 123, 80, 43, 246)
p_bidmc <- create_cohort_box("BIDMC", 84, 73, 11, 168)
p_um <- create_cohort_box("UM", 36, 1, 35, 72)

# Source cohorts row
source_row <- p_pcbn + p_bidmc + p_um + 
  plot_layout(ncol = 3, widths = c(1, 1, 1))

# ============================================================================
# CREATE ARROWS AND MERGE BRACKET
# ============================================================================
create_arrow_section <- function() {
  ggplot() +
    # Left arrow (PCBN -> Discovery)
    annotate("segment", x = 0.17, xend = 0.17, y = 0.9, yend = 0.15,
             color = colors$arrow, linewidth = 0.8) +
    annotate("segment", x = 0.14, xend = 0.17, y = 0.2, yend = 0.1,
             color = colors$arrow, linewidth = 0.8) +
    annotate("segment", x = 0.20, xend = 0.17, y = 0.2, yend = 0.1,
             color = colors$arrow, linewidth = 0.8) +
    
    # Merge bracket for BIDMC and UM
    # Left arm (under BIDMC)
    annotate("segment", x = 0.5, xend = 0.5, y = 0.9, yend = 0.7,
             color = colors$arrow, linewidth = 0.8) +
    # Right arm (under UM)
    annotate("segment", x = 0.83, xend = 0.83, y = 0.9, yend = 0.7,
             color = colors$arrow, linewidth = 0.8) +
    # Horizontal connector
    annotate("segment", x = 0.5, xend = 0.83, y = 0.7, yend = 0.7,
             color = colors$arrow, linewidth = 0.8) +
    # Down to merged
    annotate("segment", x = 0.665, xend = 0.665, y = 0.7, yend = 0.15,
             color = colors$arrow, linewidth = 0.8) +
    # Arrowhead
    annotate("segment", x = 0.635, xend = 0.665, y = 0.2, yend = 0.1,
             color = colors$arrow, linewidth = 0.8) +
    annotate("segment", x = 0.695, xend = 0.665, y = 0.2, yend = 0.1,
             color = colors$arrow, linewidth = 0.8) +
    
    # Merge label
    annotate("label", x = 0.665, y = 0.55, label = "Merged",
             fill = "white", color = colors$arrow, size = 2.5,
             label.padding = unit(0.15, "lines")) +
    
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

p_arrows <- create_arrow_section()

# ============================================================================
# CREATE FINAL COHORT BOXES (Discovery & Validation)
# ============================================================================
p_discovery <- create_cohort_box("PCBN", 123, 80, 43, 246, 
                                 border_color = colors$discovery,
                                 show_label = "Discovery", 
                                 label_color = colors$discovery)

p_validation <- create_cohort_box("BM", 120, 74, 46, 240,
                                  border_color = colors$validation,
                                  show_label = "Validation",
                                  label_color = colors$validation)

# Final cohorts row with spacing
final_row <- p_discovery + plot_spacer() + p_validation + 
  plot_layout(ncol = 3, widths = c(1, 0.5, 1))

# ============================================================================
# CREATE MATCHED SAMPLES SECTION
# ============================================================================
create_matched_section <- function() {
  ggplot() +
    # Background box
    annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.1, ymax = 0.9,
             fill = "#e8f6f3", color = colors$matched, linewidth = 1.2) +
    
    # "All Patients" badge
    annotate("label", x = 0.12, y = 0.9, label = "ALL PATIENTS",
             fill = colors$matched, color = "white", fontface = "bold",
             size = 2.5, label.padding = unit(0.2, "lines")) +
    
    # Title
    annotate("text", x = 0.5, y = 0.78, 
             label = "Matched Tumor & Adjacent Benign Tissue Pairs",
             fontface = "bold", size = 4, color = "#0e6655") +
    
    # Patient icon (circle with person symbol approximation)
    annotate("point", x = 0.12, y = 0.5, size = 12, shape = 21,
             fill = "white", color = colors$matched, stroke = 1.5) +
    annotate("text", x = 0.12, y = 0.5, label = "👤", size = 5) +
    annotate("text", x = 0.12, y = 0.32, label = "Each Patient",
             size = 2.5, fontface = "bold", color = "#0e6655") +
    
    # Arrow with "paired" label
    annotate("segment", x = 0.18, xend = 0.28, y = 0.5, yend = 0.5,
             color = colors$matched, linewidth = 1,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    annotate("text", x = 0.23, y = 0.57, label = "PAIRED",
             size = 2, color = colors$matched, fontface = "bold") +
    
    # Tissue pair box
    annotate("rect", xmin = 0.30, xmax = 0.68, ymin = 0.38, ymax = 0.62,
             fill = "white", color = "#a3e4d7", linewidth = 0.8) +
    
    # Tumor tissue
    annotate("point", x = 0.37, y = 0.5, size = 8, shape = 21,
             fill = "#fdeaa8", color = colors$tumor, stroke = 1.5) +
    annotate("text", x = 0.45, y = 0.53, label = "Tumor", 
             fontface = "bold", size = 2.8, hjust = 0) +
    annotate("text", x = 0.45, y = 0.47, label = "tissue", 
             size = 2.5, hjust = 0, color = "grey40") +
    
    # Plus sign
    annotate("text", x = 0.52, y = 0.5, label = "+",
             size = 5, fontface = "bold", color = colors$matched) +
    
    # Adjacent benign tissue
    annotate("point", x = 0.56, y = 0.5, size = 8, shape = 21,
             fill = "#d5dbdb", color = colors$benign, stroke = 1.5) +
    annotate("text", x = 0.595, y = 0.53, label = "Adjacent Benign", 
             fontface = "bold", size = 2.8, hjust = 0) +
    annotate("text", x = 0.595, y = 0.47, label = "tissue", 
             size = 2.5, hjust = 0, color = "grey40") +
    
    # Arrow to RNA-seq
    annotate("segment", x = 0.70, xend = 0.77, y = 0.5, yend = 0.5,
             color = colors$rna_seq, linewidth = 1,
             arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    
    # RNA-seq box
    annotate("rect", xmin = 0.78, xmax = 0.96, ymin = 0.40, ymax = 0.60,
             fill = colors$rna_seq, color = NA) +
    annotate("text", x = 0.87, y = 0.53, label = "RNA-seq",
             color = "white", fontface = "bold", size = 3.5) +
    annotate("text", x = 0.87, y = 0.46, label = "Transcriptome profiling",
             color = "white", size = 2, alpha = 0.9) +
    
    # Total samples summary
    annotate("text", x = 0.35, y = 0.20, 
             label = "Total: 243 patients × 2 samples =",
             size = 3, fontface = "bold", color = "#0e6655") +
    annotate("label", x = 0.62, y = 0.20, label = "486 RNA-seq profiles",
             fill = colors$matched, color = "white", fontface = "bold",
             size = 3, label.padding = unit(0.2, "lines")) +
    
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

p_matched <- create_matched_section()

# ============================================================================
# CREATE LEGEND
# ============================================================================
create_legend <- function() {
  ggplot() +
    # Non-recurrence
    annotate("rect", xmin = 0.30, xmax = 0.34, ymin = 0.4, ymax = 0.6,
             fill = colors$non_recurrence) +
    annotate("text", x = 0.36, y = 0.5, label = "Non-recurrence",
             hjust = 0, size = 3, color = "grey40") +
    
    # Recurrence
    annotate("rect", xmin = 0.55, xmax = 0.59, ymin = 0.4, ymax = 0.6,
             fill = colors$recurrence) +
    annotate("text", x = 0.61, y = 0.5, label = "Recurrence (BCR)",
             hjust = 0, size = 3, color = "grey40") +
    
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void()
}

p_legend <- create_legend()

# ============================================================================
# COMBINE ALL ELEMENTS
# ============================================================================
# Figure title
p_title <- ggplot() +
  annotate("text", x = 0, y = 0.5, label = "A", 
           fontface = "bold", size = 5, hjust = 0) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

# Assemble final figure
final_figure <- p_title / 
  source_row / 
  p_arrows / 
  final_row / 
  p_matched / 
  p_legend +
  plot_layout(heights = c(0.08, 1, 0.4, 1, 0.8, 0.15))

# ============================================================================
# SAVE FIGURE
# ============================================================================
# Save as PDF (vector format, best for journals)
ggsave("Cohorts.pdf", final_figure, path = outdir, 
       width = 8, height = 10, units = "in", dpi = 600)

# Save as PNG (raster format)
ggsave("Cohorts.png", final_figure, path = outdir, 
       width = 8, height = 10, units = "in", dpi = 600)

# Save as TIFF (often required by journals)
ggsave("Cohorts.tiff", final_figure, path = outdir, 
       width = 8, height = 10, units = "in", dpi = 600, compression = "lzw")

# Display figure
print(final_figure)

message("Figure saved as Figure1A.pdf, Figure1A.png, and Figure1A.tiff")
