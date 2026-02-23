## ========== Load packages ==========
set.seed(777)
required_pkgs <- c("patchwork","glmnet", "xgboost", "pROC", "dplyr", "ggplot2", "gtools", "rpart", "nnet", "yardstick", "tibble", "tidyr","xgboost","readxl","reshape","MASS","purrr","ggpubr","RColorBrewer","glue")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, require, character.only = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))
## ========== Input data and output directory ==========
benign_sig_score_file='../data/Signature_score_benign.txt'
PCBN_score_file='../data/PCBN_clinical_ML_score.txt'
BM_score_file='../data/BM_clinical_ML_score.txt'
model_file='../data/ML_model.rds'

outdir='../results/recurency_prediction'

if(!dir.exists(outdir)){
  dir.create(outdir,recursive = T)
}

## ========== Data loading and pre-processing ==========
## Load data 
benign_df=read.delim(benign_sig_score_file)
train_psa_pred=read.delim(PCBN_score_file)
test_psa_pred=read.delim(BM_score_file)

## ------- Learn threshold from Traning data (PCBN) -------
df <- train_psa_pred_sel

# Convert labels to binary
df <- df %>%
  mutate(label = ifelse(obs == "yes", 1, 0))

# reverse B cell score
df$B_cell_33520406 <- -df$B_cell_33520406

# Identify score columns (excluding metadata)
df_plot=df
colnames(df_plot)=gsub('\\-','',colnames(df_plot))
# colnames(df_plot) <- sub("_alpha.*", "", colnames(df_plot))
score_cols <- setdiff(names(df_plot), c("PATIENT", "obs", "label"))

# Prepare a data frame to store cutoff results
cutoff_results <- data.frame(Score_Name = character(), Youden_Cutoff = numeric(), stringsAsFactors = FALSE)

# Loop through each score column
for (col in score_cols) {
  scores <- as.numeric(df_plot[[col]])
  labels <- df_plot$label
  
  # Calculate ROC
  roc_obj <- roc(labels, scores, quiet = TRUE)
  youden_coords <- coords(roc_obj, "best", best.method = "youden", ret = "threshold")
  cutoff <- as.numeric(youden_coords$threshold[1])
  
  # Store result
  cutoff_results <- rbind(cutoff_results, data.frame(Score_Name = col, Youden_Cutoff = as.numeric(cutoff)))
  
  # Generate density plot
  p <- ggplot(df_plot, aes_string(x = col, fill = 'obs')) +
    geom_density(alpha = 0.7) +
    scale_fill_brewer(palette="Dark2")+
    geom_vline(xintercept = as.numeric(cutoff), color = "red", linetype = "dashed") +
    labs(title = paste0(col, " (Youden cutoff = ", round(cutoff, 3), ")"),
         x = "Score", y = "Density") +
    theme_classic()
  p
  ggsave(paste0(col,'.train.pdf'), p, path = outdir,width = 6,height = 4)
}

# Show cutoff results
print(cutoff_results)

# odds ratio function
calc_odds_ratio <- function(cutoff, score, response) {
  pred <- ifelse(score > cutoff, 1, 0)
  tp <- sum(pred == 1 & response == 1, na.rm = TRUE)
  tn <- sum(pred == 0 & response == 0, na.rm = TRUE)
  fp <- sum(pred == 1 & response == 0, na.rm = TRUE)
  fn <- sum(pred == 0 & response == 1, na.rm = TRUE)
  
  denominator <- fp + fn
  if (denominator == 0) return(NA)  # avoid divide-by-zero
  
  odds_ratio <- (tp + tn) / denominator
  return(odds_ratio)
}

# Prepare result dataframe
odds_results <- data.frame(Score_Name = character(), Youden_Cutoff = numeric(), Odds_Ratio = numeric(), stringsAsFactors = FALSE)

# Loop through each score and calculate odds
for (i in 1:nrow(cutoff_results)) {
  col_name <- cutoff_results$Score_Name[i]
  cutoff <- cutoff_results$Youden_Cutoff[i]
  score <- df_plot[[col_name]]
  response <- df_plot$label
  
  # Try to compute odds, safely skip if error
  odds <- tryCatch({
    calc_odds_ratio(cutoff, score, response)
  }, error = function(e) {
    warning(paste("Failed for", col_name))
    return(NA)
  })
  
  # Try to compute odds and CI
  odds_ci <- tryCatch({
    calc_odds_ratio_ci(cutoff, score, response)
  }, error = function(e) {
    warning(paste("Failed for", col_name))
    return(NA)
  })
  odds_results <- rbind(odds_results, data.frame(Score_Name = col_name, Youden_Cutoff = cutoff, Odds_Ratio = odds))
}

## ------- Apply cutoff to Testing cohort -------
# Ensure test_df has same structure and response column as binary
test_df <- test_psa_pred_sel %>%
  mutate(label = ifelse(obs == "yes", 1, 0))

colnames(test_df)=gsub('\\-','',colnames(test_df))

# reverse B cell score
test_df$B_cell_33520406 <- -test_df$B_cell_33520406

# Also ensure df has label (training set)
df <- df %>%
  mutate(label = ifelse(obs == "yes", 1, 0))

# Create results table
all_odds <- data.frame()

# Apply TRAINING cutoffs directly
for (i in 1:nrow(cutoff_results)) {
  score_name <- cutoff_results$Score_Name[i]
  cutoff <- cutoff_results$Youden_Cutoff[i]
  
  # From training set
  train_or <- tryCatch(calc_odds_ratio(cutoff, df_plot[[score_name]], df$label), error = function(e) NA)
  
  # From test set using same cutoff
  tmp_test=test_df[c(score_name,'label')]
  tmp_test=na.omit(tmp_test)
  test_or <- tryCatch(calc_odds_ratio(cutoff, test_df[[score_name]], test_df$label), error = function(e) NA)
  test_or <- tryCatch(calc_odds_ratio(cutoff, tmp_test[[score_name]], tmp_test$label), error = function(e) NA)
  
  # Append to result
  all_odds <- rbind(all_odds,
                    data.frame(Score = score_name, Dataset = "Train", Odds_Ratio = train_or),
                    data.frame(Score = score_name, Dataset = "Test", Odds_Ratio = test_or))
}


# Plot odds ratios
#all_odds$Score=factor(all_odds$Score,levels = c(sel_mls,'PSA'))
all_odds$Dataset=factor(all_odds$Dataset,levels = c('Train','Test'))
all_odds$ML=sapply(all_odds$Score,function(x)sub("_.*", "", x))
all_odds$Combination=sapply(all_odds$Score,function(x)sub("^[A-Za-z0-9]+_", "", x))

sel_comb2=c("PSA","GLEASON_SCORE","CAPRA","DAmico_num","B_cell_33520406","GLMNet_PCBN_alpha0.90_lambda0.04281")

tmp_df_sel=all_odds[which(all_odds$Score%in%sel_comb2),]

tmp_df_sel$Score[which(tmp_df_sel$Score=='GLMNet_PCBN_alpha0.90_lambda0.04281')]='Benign clinical score'
tmp_df_sel$Score[which(tmp_df_sel$Score=='B_cell_33520406')]='B cell sig. score'
tmp_df_sel$Score[which(tmp_df_sel$Score=='PSA')]='PSA'
tmp_df_sel$Score[which(tmp_df_sel$Score=='GLEASON_SCORE')]='Gleason score'
tmp_df_sel$Score[which(tmp_df_sel$Score=='CAPRA')]='CAPRA'
tmp_df_sel$Score[which(tmp_df_sel$Score=='DAmico_num')]='DAmico'

tmp_df_sel$Score=factor(tmp_df_sel$Score,levels = c('Benign clinical score','B cell sig. score','PSA','Gleason score','CAPRA','DAmico'))

p=ggplot(tmp_df_sel, aes(x = Score, y = Odds_Ratio, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(Odds_Ratio, 2)), 
            position = position_dodge(width = 0.7), 
            vjust = -0.4, size = 2.5) +
  labs(title = "",
       x = "",
       y = "Odds Ratio") +
  scale_fill_manual(values = c("Train" = "#1f77b4", "Test" = "#ff7f0e")) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5),
        legend.position = "top")
p
ggsave(paste0('ML','.odds.pdf'), p, path = outdir,width = 6,height = 5)

write.table(all_odds,file.path(outdir,'all_odds.txt'),quote = F,sep = '\t',row.names = F)


## ----- Check parameters and importance -----
# caret train object
model <- readRDS(model_file)

# Use the lambda that caret selected (or set manually if you prefer)
lambda_used <- model$bestTune$lambda
lambda_used
# if you really want to force it:
# lambda_used <- 0.04281

# Extract coefficients from the underlying glmnet model
coef_mat <- as.matrix(coef(model$finalModel, s = lambda_used))

coef_df <- data.frame(
  term = rownames(coef_mat),
  coefficient = as.numeric(coef_mat[, 1]),
  row.names = NULL
)

# Separate intercept and predictors
intercept <- coef_df$coefficient[coef_df$term == "(Intercept)"]
terms_df <- subset(coef_df, term != "(Intercept)" & coefficient != 0)

# Build equation string (logistic example)
term_strings <- paste0(
  sprintf("%.4f", terms_df$coefficient),
  " * ",
  terms_df$term
)

equation_string <- paste0(
  "logit(p) = ",
  sprintf("%.4f", intercept),
  if (nrow(terms_df) > 0) paste0(" + ", paste(term_strings, collapse = " + ")) else ""
)

cat(equation_string, "\n")


## Barplot of importance 
importance_df <- subset(coef_df, term != "(Intercept)" & coefficient != 0)
importance_df$importance <- (importance_df$coefficient)

# Top N features
top_n <- 20
importance_df <- importance_df[order(importance_df$importance, decreasing = TRUE), ]
plot_df <- head(importance_df, top_n)
plot_df$term=c('Gleason score','PSA','B cell sig. score')
plot_df$coef_type <- ifelse(plot_df$coefficient > 0, "Positive", "Negative")

ggplot(plot_df, aes(x = reorder(term, importance), 
                    y = importance, 
                    fill = coef_type)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#e31a1c", "Negative" = "#1f78b4")) +
  xlab("") +
  ylab("Coefficient") +
  ggtitle("") +
  theme_classic() +
  theme(legend.title = element_blank(),legend.position = 'none')

ggsave('GLMNet_PCBN_impoartance.pdf',path = outdir,width = 4.5,height = 2)

## ========== Signature score distribution between recurrence and non-recurrence ==========
df_plot_train=train_psa_pred_sel
df_plot_train$obs=ifelse(df_plot_train$obs=='yes','Rec','Non-rec')
df_plot_train=melt(df_plot_train,id.vars = c('obs'))
df_plot_train$Cohort='PCBN'

df_plot_test=test_psa_pred_sel
df_plot_test$obs=ifelse(df_plot_test$obs=='yes','Rec','Non-rec')
df_plot_test=melt(df_plot_test,id.vars = c('obs'))
df_plot_test$Cohort='BIDMC_UM'

df_merged_sel=rbind(df_plot_train,df_plot_test)

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
    
    p=ggplot(tmp_cohort_var, aes(x = obs, y = value, fill = obs)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, color = 'black') +
      scale_fill_brewer(palette="Dark2")+
      labs(title = i, y = y_lab,x='') +
      theme_classic()+stat_compare_means(label = "p.format",method = 'wilcox',label.x = 1.5)+
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")  # centered
      )
    p
    ggsave(filename = paste0("feature_plot_", i,'_',j, ".pdf"), plot = p,path=outdir, width = 3.5, height = 4)
  }
}

## ========== ROC curve for each sig score ==========
# Function
roc_plot=function(df_score){
  
  # Create a list to store ROC curves
  roc_list <- list(
    `B cell sig. score` = pROC::roc(df_score$obs, df_score$`B cell sig. score`)
  )
  
  # Calculate AUC and update names
  auc_values <- sapply(roc_list, pROC::auc)
  names(roc_list) <- paste0(names(roc_list), " (AUC: ", round(auc_values, 2), ")")
  
  # Define custom colors using a color palette
  color_palette <- brewer.pal(n = length(roc_list), name = "Set1")
  custom_colors <- setNames(color_palette, names(roc_list))
  
  # Plot with ggroc and custom colors
  p = ggroc(roc_list,linewidth = 1) +
    ggtitle("") +
    theme_classic2() +
    labs(color = "Model") +
    scale_color_manual(values = custom_colors) +
    theme(
      legend.position = "right",
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )+
  
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE)+
    geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "black", linewidth = 1)+   # AUC=0.5 baseline
    geom_hline(yintercept = 1, color = "black", linewidth = 1)+         
    geom_vline(xintercept = 0, color = "black", linewidth = 1)         
  
  return(p)
}

# AUC curve for PCBN
df_score=train_psa_pred_sel
colnames(df_score)=c('PSA','Gleason score','CAPRA','DAmico_num' ,'B cell sig. score','Benign clinical score','obs')

p_train=roc_plot(df_score)
ggsave(file.path(outdir,paste0('ROC_AUC_PCBN_B_sig.pdf')),width = 5,height = 3)

# AUC curve for BIDMC_UM
df_score=test_psa_pred_sel
colnames(df_score)=c('PSA','Gleason score','CAPRA','DAmico_num' ,'B cell sig. score','Benign clinical score','obs')

p_test=roc_plot(df_score)
ggsave(file.path(outdir,paste0('ROC_AUC_BIDMC_UM_B_sig.pdf')),width = 5,height = 3)

## ========== AUC of all feature scores ==========
## Function 
auc_prauc=function(df){
  auc_res=data.frame()
  for (i in setdiff(colnames(df),'obs')){
    auc1 <- pROC::auc(df$obs, df[,i], direction = "<")
    prauc1 <- yardstick::pr_auc_vec(as.factor(df$obs), as.numeric(df[,i]), event_level = "second")
    tmp_res=data.frame(AUC=auc1,PRAUC=prauc1,sig_name=i)
    auc_res=rbind(auc_res,tmp_res)
  }
  return(auc_res)
}

## Calculate AUC and PRAUC
df_auc_train=train_psa_pred_sel
colnames(df_auc_train)=c('PSA','Gleason score','CAPRA','DAmico_num' ,'B cell sig. score','Benign clinical score','obs')
# reverse B cell score
df_auc_train$`B cell sig. score` <- -df_auc_train$`B cell sig. score`

df_auc_test=test_psa_pred_sel
colnames(df_auc_test)=c('PSA','Gleason score','CAPRA','DAmico_num' ,'B cell sig. score','Benign clinical score','obs')
# reverse B cell score
df_auc_test$`B cell sig. score` <- -df_auc_test$`B cell sig. score`

## Training dataset
train_auc=auc_prauc(df_auc_train)

## Testing dataset
test_auc=auc_prauc(df_auc_test)

## plot AUC and PRAUC
train_auc$Cohort='PCBN'
test_auc$Cohort='BM'

df_auc=rbind(train_auc,test_auc)
df_auc$sig_name=factor(df_auc$sig_name,levels = c('Benign clinical score','B cell sig. score','PSA','Gleason score','CAPRA','DAmico_num'))
df_auc$Cohort=factor(df_auc$Cohort,levels = c('PCBN','BM'))

# AUC
ggplot(df_auc, aes(x = sig_name, y = AUC)) +
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
ggplot(df_auc, aes(x = sig_name, y = PRAUC)) +
  geom_bar(stat = "identity",fill='#FF7F50',width = 0.75) +
  facet_wrap(~Cohort, scales = "free_y") +
  #geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", PRAUC)),
            position = position_dodge(width = 0.7), vjust = -0.25, size = 3) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))+
  labs(title = "",y = "AUPRC",x='') 
ggsave(file.path(outdir, "PRAUC_merged_cohort.pdf"), width = 8,height = 4.5)

## ========== Correlation between B cell signature score and clinical score ==========
## Function
cor_scatter=function(df,x,y,filename){
  sp <- ggscatter(df, x = x, y = y,
                  add = "reg.line",  # Add regressin line
                  add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                  conf.int = TRUE # Add confidence interval
  )
  # Add correlation coefficient
  sp <- sp + stat_cor(method = "pearson")
  ggsave(filename = filename,plot = sp,path = outdir,width = 4,height = 3)
}

cor_scatter(df=train_psa_pred_sel,x='PSA',y='B_cell_33520406',filename='cor_b_cell_psa_PCBN.pdf')
cor_scatter(df=train_psa_pred_sel,x='GLEASON_SCORE',y='B_cell_33520406',filename='cor_b_cell_gleason_PCBN.pdf')
cor_scatter(df=train_psa_pred_sel,x='PSA',y='GLEASON_SCORE',filename='cor_gleason_psa_PCBN.pdf')
cor_scatter(df=test_psa_pred_sel,x='PSA',y='B_cell_33520406',filename='cor_b_cell_psa_BIDMC_UM.pdf')
cor_scatter(df=test_psa_pred_sel,x='GLEASON_SCORE',y='B_cell_33520406',filename='cor_b_cell_gleason_BIDMC_UM.pdf')
cor_scatter(df=test_psa_pred_sel,x='PSA',y='GLEASON_SCORE',filename='cor_gleason_psa_BIDMC_UM.pdf')

## ========== Plot TP TN FP FN ==========
## Function
confusion_plot=function(df,odds_results,score_id,prefix){
  cutoff=odds_results$Youden_Cutoff[which(odds_results$Score_Name==score_id)]
  colnames(df)[which(colnames(df)==score_id)]='score'
  df2 <- df %>%
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
  
  p = ggplot(conf_df, aes(x = pred_lab, y = truth_lab, fill = n)) +
    geom_tile(color = "black") +
    geom_text(aes(label = label), size = 6) +
    scale_fill_gradient(low = "#deebf7", high = "#08519c") +  # scientific blues
    coord_equal() +
    labs(
      x = "Predicted",
      y = "True",
      fill = "Count",
      title = score_id
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
  ggsave(paste0(prefix,'_',score_id,'_confusion_matrix.pdf'),path = outdir,width = 5,height = 4.5)
}

## PCBN
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='GLMNet_PCBN_alpha0.90_lambda0.04281',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='PSA',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='GLEASON_SCORE',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='B_cell_33520406',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='PSA',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='CAPRA',prefix='PCBN')
confusion_plot(df=train_psa_pred,odds_results=odds_results,score_id='DAmico_num',prefix='PCBN')

## BIDMC_UM
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='GLMNet_PCBN_alpha0.90_lambda0.04281',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='PSA',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='GLEASON_SCORE',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='B_cell_33520406',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='PSA',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='CAPRA',prefix='BIDMC_UM')
confusion_plot(df=test_psa_pred,odds_results=odds_results,score_id='DAmico_num',prefix='BIDMC_UM')

## ========== Multiple variable regression ==========
dat <- benign_df %>%
  transmute(
    RECUR           = as.integer(RECUR),
    B_cell_33520406 = as.numeric(B_cell_33520406),
    PSA             = as.numeric(PSA),
    GLEASON_SCORE   = as.numeric(GLEASON_SCORE),
    SAMPLE_SOURCE   = SAMPLE.SOURCE
  ) %>%
  filter(!is.na(RECUR), !is.na(B_cell_33520406),
         !is.na(PSA), !is.na(GLEASON_SCORE), !is.na(SAMPLE_SOURCE)) %>%
  mutate(
    COHORT = if_else(SAMPLE_SOURCE %in% c("BIDMC","UM"), "BM", "PCBN")
  )

# Helper to fit one cohort --------------------------------------------
fit_one_cohort <- function(df_sub) {
  stopifnot(length(unique(df_sub$RECUR)) == 2)
  fit <- glm(RECUR ~ B_cell_33520406 + PSA + GLEASON_SCORE,
             data = df_sub, family = binomial())
  or_tbl <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::transmute(
      term,
      OR      = estimate,
      ci_low  = conf.low,
      ci_high = conf.high,
      p       = p.value
    )
  auc <- pROC::roc(df_sub$RECUR, fitted(fit))$auc
  list(or = or_tbl, auc = as.numeric(auc), n = nrow(df_sub), rate = mean(df_sub$RECUR))
}

# Fit BM and PCBN separately ------------------------------------------
bm   <- dat %>% filter(COHORT == "BM")
pcbn <- dat %>% filter(COHORT == "PCBN")

res_bm   <- fit_one_cohort(bm)
res_pcbn <- fit_one_cohort(pcbn)

print(res_bm$or)
print(res_pcbn$or)
cat(glue("BM:   n={res_bm$n}, RECUR rate={round(res_bm$rate,3)}, AUC={round(res_bm$auc,3)}\n"))
cat(glue("PCBN: n={res_pcbn$n}, RECUR rate={round(res_pcbn$rate,3)}, AUC={round(res_pcbn$auc,3)}\n"))

# Forest-style figure with two panels ---------------------------------
clean_labels <- function(tbl) {
  tbl %>%
    mutate(term = recode(term,
                         "(Intercept)"      = "Intercept",
                         "B_cell_33520406"  = "B_cell_33520406 (per unit)",
                         "PSA"              = "PSA (per unit)",
                         "GLEASON_SCORE"    = "Gleason score (per unit)"
    )) %>%
    filter(term != "Intercept")
}

bm_plot   <- clean_labels(res_bm$or)
pcbn_plot <- clean_labels(res_pcbn$or)

p1 <- ggplot(bm_plot, aes(x = OR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "grey40") +
  geom_point(size = 2) +
  geom_text(aes(x = ci_high * 1.05,
                label = glue("OR={formatC(OR, digits=2, format='f')} (CI {formatC(ci_low, digits=2, format='f')}-{formatC(ci_high, digits=2, format='f')}), p={formatC(p, digits=3, format='f')}")),
            hjust = 0, vjust = 0.5, size = 3.2) +
  labs(title = glue("BM cohort (n={res_bm$n}, RECUR rate={round(res_bm$rate,2)}; AUC={round(res_bm$auc,3)})"),
       x = "Odds ratio for RECUR", y = NULL) +
  theme_classic(base_size = 12) + theme(panel.grid.minor = element_blank())

p2 <- ggplot(pcbn_plot, aes(x = OR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "grey40") +
  geom_point(size = 2) +
  geom_text(aes(x = ci_high * 1.05,
                label = glue("OR={formatC(OR, digits=2, format='f')} (CI {formatC(ci_low, digits=2, format='f')}-{formatC(ci_high, digits=2, format='f')}), p={formatC(p, digits=3, format='f')}")),
            hjust = 0, vjust = 0.5, size = 3.2) +
  labs(title = glue("PCBN cohort (n={res_pcbn$n}, RECUR rate={round(res_pcbn$rate,2)}; AUC={round(res_pcbn$auc,3)})"),
       x = "Odds ratio for RECUR", y = NULL) +
  theme_classic(base_size = 12) + theme(panel.grid.minor = element_blank())

fig <- p1 + p2
fig
ggsave(file.path(outdir,"RECUR_logit_by_cohort.pdf"), fig, width = 12, height = 3.5)



