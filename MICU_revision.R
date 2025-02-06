setwd("~/Documents/GitHub/Medical_Intensive_Care_Unit/")

load("Data/MICU_Data_Anon.revision.RData")

library(pROC)
library(survminer)
library(survival)
library(Maaslin2)
library(gridExtra)
library(tidyverse)

# maaslin at family, genus and species level between survivor vs. non-survivor --------------------------------------

## maaslin at genus level --------------------------------------------------

maaslin_mat_genus <- t_metaphlan_micu_nocovid %>%
  distinct() %>%
  left_join(taxdmp %>% mutate(taxid = as.character(taxid))) %>%
  mutate(Genus = ifelse(Genus == "", str_extract(Genus, pattern = "([^\\s]+)"), Genus)) %>%
  group_by(shotgunSeq_id, metabolomicsID, db, Genus) %>%
  summarise(pctseqs = sum(pctseqs)) %>%
  pivot_wider(
    id_cols = c(shotgunSeq_id),
    names_from = "Genus",
    values_from = "pctseqs",
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames(var = "shotgunSeq_id")

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_genus,
  input_metadata = data.frame(
    t_metaphlan_micu_nocovid_mat %>%
      rownames_to_column(var = "shotgunSeq_id") %>%
      select(shotgunSeq_id) %>%
      left_join(micu_new_nocovid_oc %>%
                  select(shotgunSeq_id, unique_id)) %>%
      left_join(
        tableone_nocovid_df_filt %>%
          labelled::remove_labels() %>%
          janitor::clean_names() %>%
          mutate(
            race_factor = as.character(race_factor),
            race_factor = ifelse(
              race_factor %in% c("Asian", "More than one race"),
              "Other",
              race_factor
            )
          )
      ) %>%
      column_to_rownames(var = "shotgunSeq_id") %>%
      select(-c(unique_id)) %>%
      mutate(race_factor = as.factor(race_factor))
  ),
  output = "Maaslin_Results/Genus_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("thirtyday_mortality_overall"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

### Run Maaslin2 with Covariates ---------
set.seed(123)

maaslin.meta <- data.frame(
  t_metaphlan_micu_nocovid_mat %>%
    rownames_to_column(var = "shotgunSeq_id") %>%
    select(shotgunSeq_id) %>%
    left_join(micu_new_nocovid_oc %>%
                select(shotgunSeq_id, unique_id)) %>%
    left_join(
      tableone_nocovid_df_filt %>%
        labelled::remove_labels() %>%
        janitor::clean_names() %>%
        mutate(
          race_factor = as.character(race_factor),
          race_factor = case_when(
            race_factor %in% c("Asian", "More than one race") ~ "Other",
            race_factor == "White, non-Hispanic" ~ "White",
            TRUE ~ race_factor
          )
        )
    ) %>%
    column_to_rownames(var = "shotgunSeq_id") %>%
    select(-c(unique_id)) %>%
    mutate(race_factor = as.factor(race_factor))
)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_genus,
  input_metadata = maaslin.meta,
  output = "Maaslin_Results/Genus_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "thirtyday_mortality_overall",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

# Manual p-value adjustment for specific comparsisons
maaslin2_all_results <- maaslin_model$results

maaslin2_results <-
  maaslin2_all_results %>% filter(metadata == "thirtyday_mortality_overall") # Discard covariate associations

maaslin2_results$qval <-
  p.adjust(maaslin2_results$pval, method = "BH") # FDR correction using 'BH'

## maaslin at species level --------------------------------------------------

maaslin_mat_sp <- t_metaphlan_micu_nocovid %>%
  distinct() %>%
  left_join(taxdmp %>% mutate(taxid = as.character(taxid))) %>%
  group_by(shotgunSeq_id, metabolomicsID, db, Species) %>%
  summarise(pctseqs = sum(pctseqs)) %>%
  pivot_wider(
    id_cols = c(shotgunSeq_id),
    names_from = "Species",
    values_from = "pctseqs",
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames(var = "shotgunSeq_id")

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_sp,
  input_metadata = data.frame(
    t_metaphlan_micu_nocovid_mat %>%
      rownames_to_column(var = "shotgunSeq_id") %>%
      select(shotgunSeq_id) %>%
      left_join(micu_new_nocovid_oc %>%
                  select(shotgunSeq_id, unique_id)) %>%
      left_join(
        tableone_nocovid_df_filt %>%
          labelled::remove_labels() %>%
          janitor::clean_names() %>%
          mutate(
            race_factor = as.character(race_factor),
            race_factor = ifelse(
              race_factor %in% c("Asian", "More than one race"),
              "Other",
              race_factor
            )
          )
      ) %>%
      column_to_rownames(var = "shotgunSeq_id") %>%
      select(-c(unique_id)) %>%
      mutate(race_factor = as.factor(race_factor))
  ),
  output = "Maaslin_Results/Species_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("thirtyday_mortality_overall"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

### Run Maaslin2 with Covariates ---------
set.seed(123)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_sp,
  input_metadata = maaslin.meta,
  output = "Maaslin_Results/Species_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "thirtyday_mortality_overall",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

# Manual p-value adjustment for specific comparsisons
maaslin2_all_results <- maaslin_model$results

maaslin2_results <-
  maaslin2_all_results %>% filter(metadata == "thirtyday_mortality_overall") # Discard covariate associations

maaslin2_results$qval <-
  p.adjust(maaslin2_results$pval, method = "BH") # FDR correction using 'BH'

## maaslin at family level: --------------------------------------------------

maaslin_mat_fa <- t_metaphlan_micu_nocovid %>%
  distinct() %>%
  left_join(taxdmp %>% mutate(taxid = as.character(taxid))) %>%
  mutate(Family = ifelse(Family == "", str_extract(Species, pattern = "([^\\s]+)"), Family)) %>%
  group_by(shotgunSeq_id, metabolomicsID, db, Family) %>%
  summarise(pctseqs = sum(pctseqs)) %>%
  pivot_wider(
    id_cols = c(shotgunSeq_id),
    names_from = "Family",
    values_from = "pctseqs",
    values_fill = 0,
    values_fn = sum
  ) %>%
  column_to_rownames(var = "shotgunSeq_id")

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_fa,
  input_metadata = data.frame(
    t_metaphlan_micu_nocovid_mat %>%
      rownames_to_column(var = "shotgunSeq_id") %>%
      select(shotgunSeq_id) %>%
      left_join(micu_new_nocovid_oc %>%
                  select(shotgunSeq_id, unique_id)) %>%
      left_join(
        tableone_nocovid_df_filt %>%
          labelled::remove_labels() %>%
          janitor::clean_names() %>%
          mutate(
            race_factor = as.character(race_factor),
            race_factor = ifelse(
              race_factor %in% c("Asian", "More than one race"),
              "Other",
              race_factor
            )
          )
      ) %>%
      column_to_rownames(var = "shotgunSeq_id") %>%
      select(-c(unique_id)) %>%
      mutate(race_factor = as.factor(race_factor))
  ),
  output = "Maaslin_Results/Family_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("thirtyday_mortality_overall"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

### Run Maaslin2 with Covariates ---------
set.seed(123)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_fa,
  input_metadata = maaslin.meta,
  output = "Maaslin_Results/Family_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "thirtyday_mortality_overall",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("thirtyday_mortality_overall,Survivor")
)

# Manual p-value adjustment for specific comparsisons
maaslin2_all_results <- maaslin_model$results

maaslin2_results <-
  maaslin2_all_results %>% filter(metadata == "thirtyday_mortality_overall") # Discard covariate associations

maaslin2_results$qval <-
  p.adjust(maaslin2_results$pval, method = "BH") # FDR correction using 'BH'

# Samples in training set within 3 days of admission: CoxPH-----------------------------------------------

cox_sub <- cox_df |> 
  filter(`Time to stool sample` <= 3)

## MDS as continuous: CoxPH ---------
sub.mds.cox <- coxph(
  Surv(cox_sub$surv_days, cox_sub$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    `MDS`,
  data = cox_sub
)

base.sub.cox <- sub.mds.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)  %>%
  gtsummary::modify_caption("**Cox Proportional Hazards Regression**")

base.sub.cox

gt::gtsave(gtsummary::as_gt(base.sub.cox), file = "Results/training.subset3days_cont.coxph.png")

## MDS as dichotomous: CoxPH -------- 

cox_sub <- cox_sub |> 
  mutate(
    grouped_md_score = factor(grouped_md_score, levels = c("Low Score","High Score"))     
  )

sub.mds.cox <- coxph(
  Surv(cox_sub$surv_days, cox_sub$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    `grouped_md_score`,
  data = cox_sub
)

base.sub.cox <- sub.mds.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)  %>%
  gtsummary::modify_caption("**Cox Proportional Hazards Regression**")

base.sub.cox

gt::gtsave(gtsummary::as_gt(base.sub.cox), file = "Results/training.subset3days_dich.coxph.png")

## confusion matrix ----

caret::confusionMatrix(table(
  factor(cox_sub$grouped_md_score,
         labels = c("Survivor", "Non-Survivor"),
         levels = c("Low Score", "High Score")),
  factor(
    cox_sub$thirtyday_mortality_overall,
    levels = c("Survivor", "Non-Survivor")
  )
))

#              Survivor Non-Survivor
# Survivor           59            7
# Non-Survivor        8           22
# 
# Accuracy : 0.8438          
# 95% CI : (0.7554, 0.9098)
# No Information Rate : 0.6979          
# P-Value [Acc > NIR] : 0.0007765       
# 
# Kappa : 0.633           
# 
# Mcnemar's Test P-Value : 1.0000000       
#                                           
#             Sensitivity : 0.8806          
#             Specificity : 0.7586          
#          Pos Pred Value : 0.8939          
#          Neg Pred Value : 0.7333          
#              Prevalence : 0.6979          
#          Detection Rate : 0.6146          
#    Detection Prevalence : 0.6875          
#       Balanced Accuracy : 0.8196          
#                                           
#        'Positive' Class : Survivor  

cox_sub |> glimpse()

## KM curve ------

fit_mmp <- survfit(Surv(cox_sub$surv_days, cox_sub$thirtyday_mortality_overall_class)  ~ grouped_md_score, data = cox_sub)

ggs_mmp <- ggsurvplot(
  fit_mmp,
  data = cox_sub,
  size = 1,
  palette = c("#C45258", "#2F4858"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  risk.table.height = 0.4,
  risk.table.y.text.col = TRUE,
  tables.y.text = FALSE,
  risk.table.fontsize = 2.8,
  pval.size = 3.5,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.5, color = "gray90"),
    axis.text.y = et(color = "black", size = 10),
    axis.title.y = et(color = "black")
  ),
  legend.labs = c("High Score", "Low Score")
)

ggs_mmp$table <-
  ggs_mmp$table + labs(x = NULL, y = NULL) + theme(plot.title = eb()) # risk table

ggs_mmp

pdf(
  "cox/training.subset3days_km.pdf",
  height = 6,
  width = 8,
  onefile = FALSE
)
ggs_mmp
invisible(dev.off())

# Samples in test/validation set within 3 days of admission: coxPH + KM curve ---------------

## CoxPH: MDS ---------
cox_vc_sub <- cox_df_vc |> 
  filter(`Time to stool sample` <= 3)

sub.mds.cox.sub <- coxph(
  Surv(cox_vc_sub$surv_days, cox_vc_sub$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    `MDS`,
  data = cox_vc_sub
)

sub.mds.cox.sub.tbl <- sub.mds.cox.sub %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**Cox Proportional Hazards Regression**")

sub.mds.cox.sub.tbl

gt::gtsave(gtsummary::as_gt(sub.mds.cox.sub.tbl), file = "Results/test.subset3days_coxph.png")

## KM curve: MDS ------

fit_mmp_vc <- survfit(Surv(cox_df_vc$surv_days, cox_df_vc$thirtyday_mortality_overall_class)  ~ grouped_md_score, data = cox_df_vc)

ggs_mmp_vc <- ggsurvplot(
  fit_mmp_vc,
  data = cox_df_vc,
  size = 1,
  palette = c("#C45258", "#2F4858"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  risk.table.height = 0.4,
  risk.table.y.text.col = TRUE,
  tables.y.text = FALSE,
  risk.table.fontsize = 2.8,
  pval.size = 3.5,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.5, color = "gray90"),
    axis.text.y = et(color = "black", size = 10),
    axis.title.y = et(color = "black")
  )
)

ggs_mmp_vc$table <-
  ggs_mmp_vc$table + labs(x = NULL, y = NULL) + theme(plot.title = eb()) # risk table

ggs_mmp_vc

pdf(
  "Results/test.subset3days_km.pdf",
  height = 6,
  width = 8,
  onefile = FALSE
)
ggs_mmp_vc
invisible(dev.off())

## validation cohort MMP metrics --------------

mmp_df_vc <- cutpoints_df_vc %>% 
  filter(compound %in% c("deoxycholic acid", 
                         "isodeoxycholic acid", 
                         "lithocholic acid", 
                         "desaminotyrosine")) %>% 
  mutate(cutpoint_prediction = case_when(
    compound == "deoxycholic acid" &
      mvalue__mM >= (89.92/1000) ~ 0,
    compound == "isodeoxycholic acid" &
      mvalue__mM >= (0.97/1000) ~ 0,
    compound == "lithocholic acid" &
      mvalue__mM >= (258.25/1000) ~ 0,
    compound == "desaminotyrosine" &
      mvalue__mM >= (21.31/1000) ~ 0,
    TRUE ~ 1
  )) %>% 
  group_by(metabolomicsID, thirtyday_mortality_overall) %>%
  summarize(mmp_score = sum(cutpoint_prediction)) %>% 
  ungroup() %>% 
  mutate(grouped_mmp_score = ifelse(mmp_score >= 2, "High MMP", "Low MMP"),
         thirtyday_mortality_overall_class = ifelse(thirtyday_mortality_overall == "Survivor", 0, 1))

pROC::roc(mmp_df_vc$thirtyday_mortality_overall_class,
          mmp_df_vc$mmp_score)
# Area under the curve: 0.6627

caret::confusionMatrix(table(
  factor(mmp_df_vc$grouped_mmp_score,
         labels = c("Survivor", "Non-Survivor"),
         levels = c("Low MMP", "High MMP")),
  factor(
    mmp_df_vc$thirtyday_mortality_overall,
    levels = c("Survivor", "Non-Survivor"))
)
)

#              Survivor Non-Survivor
# Survivor           14            3
# Non-Survivor       20           12
# 
# Accuracy : 0.5306          
# 95% CI : (0.3827, 0.6747)
# No Information Rate : 0.6939          
# P-Value [Acc > NIR] : 0.9946247       
# 
# Kappa : 0.1608          
# 
# Mcnemar's Test P-Value : 0.0008492       
#                                           
#             Sensitivity : 0.4118          
#             Specificity : 0.8000          
#          Pos Pred Value : 0.8235          
#          Neg Pred Value : 0.3750          
#              Prevalence : 0.6939          
#          Detection Rate : 0.2857          
#    Detection Prevalence : 0.3469          
#       Balanced Accuracy : 0.6059          
#                                           
#        'Positive' Class : Survivor 


# Shannon CoxPH -------------------------------------------------------
## continuous ---------------------------------------------------------

shannon_cox <-
  coxph(
    Surv(cox_df$surv_days, cox_df$thirtyday_mortality_overall_class) ~
      `Charlson Comorbidity Index` +
      `SOFA Score` +
      # `Time to stool sample` +
      `Shannon Diversity`,
    data = cox_df
  )

coxauc_shannon <- shannon_cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)


coxauc_shannon <- coxauc_shannon %>%
  gtsummary::modify_caption("**Shannon Continuous CoxPH**")

# In case you get an error: "Error in s$close() : attempt to apply non-function", run this code below:
# f <- chromote::default_chromote_object() #get the f object
# f$close()

coxauc_shannon

gt::gtsave(gtsummary::as_gt(coxauc_shannon), file = "Results/shannon_cont.coxPH.training.png")

## dichotomous ---------------------------------------------------------

shannon_cox <-
  coxph(
    Surv(cox_df$surv_days, cox_df$thirtyday_mortality_overall_class) ~
      `Charlson Comorbidity Index` +
      `SOFA Score` +
      # `Time to stool sample` +
      shannon_class,
    data = cox_df
  )

coxauc_shannon <- shannon_cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)


coxauc_shannon <- coxauc_shannon %>%
  gtsummary::modify_caption("**Shannon Group CoxPH**")

coxauc_shannon

gt::gtsave(gtsummary::as_gt(coxauc_shannon), file = "./Results/shannon_class.coxPH.training.png")

# COVID-19 score: CoxPH ---------------------------------------------------

## continuous -------------------------------------------------------------
mmp_df_cox <- mmp_df |> 
  left_join(cox_df |> 
              select(metabolomicsID, surv_days,
                     `Charlson Comorbidity Index`, thirtyday_mortality_overall_class,
                       `SOFA Score`)) |> 
  mutate(grouped_mmp_score = factor(grouped_mmp_score, 
                                    levels = c("Low MMP","High MMP")))

mmp.cox <- coxph(
  Surv(mmp_df_cox$surv_days, mmp_df_cox$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    `mmp_score`,
  data = mmp_df_cox
)

base.mmp.cox <- mmp.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)  %>%
  gtsummary::modify_caption("**Cox Proportional Hazards Regression**")

base.mmp.cox

gt::gtsave(gtsummary::as_gt(base.mmp.cox), file = "Results/training.mmp_cont.coxph.png")

## dichotomous ----------------------------
mmp.cox <- coxph(
  Surv(mmp_df_cox$surv_days, mmp_df_cox$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    `grouped_mmp_score`,
  data = mmp_df_cox
)

base.mmp.cox <- mmp.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE)  %>%
  gtsummary::modify_caption("**Cox Proportional Hazards Regression**")

base.mmp.cox

gt::gtsave(gtsummary::as_gt(base.mmp.cox), file = "Results/training.mmp_class.coxph.png")

## confusion matrix -----

# Area under the curve: 0.6174

pROC::roc(mmp_df_cox$thirtyday_mortality_overall_class,
          mmp_df_cox$mmp_score)

caret::confusionMatrix(table(
  factor(mmp_df_cox$grouped_mmp_score,
         labels = c("Survivor", "Non-Survivor"),
         levels = c("Low MMP", "High MMP")),
  factor(
    mmp_df_cox$thirtyday_mortality_overall,
    levels = c("Survivor", "Non-Survivor"))
  )
)

#             Survivor Non-Survivor
# Survivor           42            8
# Non-Survivor       60           37
# 
# Accuracy : 0.5374          
# 95% CI : (0.4534, 0.6199)
# No Information Rate : 0.6939          
# P-Value [Acc > NIR] : 1               
# 
# Kappa : 0.1769          
# 
# Mcnemar's Test P-Value : 6.224e-10       
#                                           
#             Sensitivity : 0.4118          
#             Specificity : 0.8222          
#          Pos Pred Value : 0.8400          
#          Neg Pred Value : 0.3814          
#              Prevalence : 0.6939          
#          Detection Rate : 0.2857          
#    Detection Prevalence : 0.3401          
#       Balanced Accuracy : 0.6170          
#                                           
#        'Positive' Class : Survivor


# MDS vs alpha diversity --------------------------------------------------

library(ggpmisc)

shan_mds_sp <- cor.test(shannon_mmp_list$md_score,
    shannon_mmp_list$Shannon,
    method = "spearman")
# rho = -0.4869761

shannon_mmp_list |> 
  ggplot(aes(md_score, Shannon)) +
  geom_jitter(width = 0, height = 0.15, alpha = 0.65) +
  stat_poly_line() +
  # stat_poly_eq(label.x = 1, use_label(c("eq", "R2", "p"))) +
  stat_poly_eq(label.x = 1, use_label(c("R2", "p"))) +
  theme_bw() +
  labs(x="MDS", y = "Shannon diversity index")

ggsave("Results/shannon_mds.liner.pdf", height = 6.5, width = 7.5)

## tax + shannon + mds plot ------------------------------------------------

gg_metaphlan

gg_shannon <- shannon_mmp_list |> 
  mutate(thirtyday_mortality_overall = factor(
    thirtyday_mortality_overall,
    levels = c("Non-Survivor", "Survivor"))) |> 
  ggplot(aes(reorder(shotgunSeq_id, Shannon), Shannon)) +
  geom_col(color="black", fill = "grey", width = 0.65) +
  facet_grid(. ~ thirtyday_mortality_overall,
             scales = "free",
             space = "free_x"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = eb(),
    axis.ticks.x = eb(),
    strip.text.x = eb(),
    strip.background = eb(),
    axis.title.y = et(color = "black", size = 14),
    axis.text.y = et(color = "black", size = 12),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = eb(),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 0,
      l = 5
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.065),
    breaks = seq(0, 4, 1)
  ) +
  scale_x_discrete(expand = expansion(add = 1)) +
  ylab("Shannon\n") +
  xlab("") +
  geom_hline(aes(yintercept = 2.16), linetype = 2, color = "red")

gg_mds <- shannon_mmp_list |> 
  mutate(thirtyday_mortality_overall = factor(
    thirtyday_mortality_overall,
    levels = c("Non-Survivor", "Survivor"))) |> 
  ggplot(aes(reorder(shotgunSeq_id, Shannon), md_score)) +
  geom_col(color="black", fill = "grey", width = 0.65) +
  facet_grid(. ~ thirtyday_mortality_overall,
             scales = "free",
             space = "free_x"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = eb(),
    axis.ticks.x = eb(),
    strip.text.x = eb(),
    strip.background = eb(),
    axis.title.y = et(color = "black", size = 14),
    axis.text.y = et(color = "black", size = 12),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = eb(),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 0,
      l = 5
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.065),
    breaks = seq(0, 13, 2.5)
  ) +
  scale_x_discrete(expand = expansion(add = 1)) +
  ylab("MDS\n") +
  xlab("") +
  geom_hline(aes(yintercept = 7.5), linetype = 2, color = "red")

### Enterococcus -----------

library(yingtools2)

entero.gg <- phylo_rel_abd |> 
  filter(organism == "Enterococcus") |> 
  left_join(shannon_mmp_list) |> 
  mutate(thirtyday_mortality_overall = factor(
    thirtyday_mortality_overall,
    levels = c("Non-Survivor", "Survivor"))) |> 
  ggplot(aes(reorder(shotgunSeq_id, Shannon), pctseqs)) +
  geom_col(color="black", fill = "grey", width = 0.65) +
  facet_grid(. ~ thirtyday_mortality_overall,
             scales = "free",
             space = "free_x"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = eb(),
    axis.ticks.x = eb(),
    strip.text.x = eb(),
    strip.background = eb(),
    axis.title.y = et(color = "black", size = 14),
    axis.text.y = et(color = "black", size = 12),
    panel.spacing = unit(0.5, "lines"),
    panel.grid = eb(),
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 0,
      l = 5
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.065),
    breaks = c(0, 0.0012, 0.05, 0.1, 0.5, 1),
    trans = log_epsilon_trans(0.01)
  ) +
  ylab("Enterococcus Abund.\n") +
  xlab("") +
  geom_hline(aes(yintercept = 0.0012), linetype = 2, color = "red")

pdf("Results/shannon_mds.entero.tax.pdf", width = 11, height = 6.5)
gg.stack(gg_metaphlan, gg_shannon, entero.gg, gg_mds, 
         heights = c(1, 0.45, 0.45, 0.45),
         newpage = F)
dev.off()

## color beta diversity with MDS score -------------------------------------

library(vegan)

beta_mds_df <- mds_data %>%
  rownames_to_column(var = "shotgunSeq_id") %>%
  left_join(micu_new_nocovid_oc %>%
              select(shotgunSeq_id, metabolomicsID, thirtyday_mortality_overall)) |> 
  left_join(cox_df |> 
              select(metabolomicsID, grouped_md_score, MDS)) |> 
  as_tibble() |> 
  mutate(grouped_md_score = factor(grouped_md_score,
                                   levels = c("Low Score", "High Score")))

dispersion <-
  permutest(betadisper(beta_dist, beta_mds_df$grouped_md_score)) 

dispersion_pval <- dispersion$tab$`Pr(>F)`[1]

# Stats: PERMANOVA
set.seed(123)
mds_stats <-
  adonis2(
    beta_dist ~ beta_mds_df$grouped_md_score,
    method = "bray-curtis",
    permutations = 999
  )

mds_pval <- mds_stats$`Pr(>F)`[1]

ggplot_mds <-
  ggplot(
    beta_mds_df,
    aes(
      x = MDS1,
      y = MDS2,
      color = grouped_md_score,
      fill = grouped_md_score
    )
  ) +
  stat_ellipse(
    level = 0.1,
    geom = "polygon",
    alpha = 0.35,
    type = "euclid"
  ) +
  geom_point(alpha = 0.65, size = 10) +
  theme_bw() +
  theme(
    axis.title = et(color = "black", size = 72),
    axis.text = et(color = "black", size = 60),
    # plot.subtitle = et(color = "black", size = 79),
    panel.grid.minor = eb(),
    panel.grid.major = eb(),
    # legend.position = "none",
    # plot.margin = margin(
    #   # Top margin
    #   t = 5,
    #   # Right margin
    #   r = 5,
    #   # Bottom margin
    #   b = 5,
    #   # Left margin
    #   l = 5 
    # )
  ) +
  annotate(
    "text",
    x = -0.4,
    y = 1.5,
    hjust = 0,
    size = 18,
    label = paste0(
      "BetaDisper = ", dispersion_pval, "\n",
      "PERMANOVA, p = ", mds_pval
    )
  ) +
  labs(
    y = "NMDS2",
    x = "NMDS1",
    title = ""
  ) +
  ggsci::scale_color_lancet() +
  ggsci::scale_fill_lancet() +
  guides(
    fill = guide_legend("Outcome"),
    color = guide_legend("Outcome")
  )   # +
  # coord_equal(
  #   ylim = c(-1.2, 1.5),
  #   xlim = c(-1.2, 1.5)
  # )

ggplot_mds

ggsave(
  plot = ggplot_mds,
  filename = "Results/Beta_Diversity_BrayCurtis_train.pdf",
  height = 14,
  width = 16.5,
  units = "in"
)

ggplot(
  beta_mds_df,
  aes(
    x = MDS1,
    y = MDS2,
    color = MDS
  )
) +
  geom_point(alpha = 0.85, size = 10) +
  paletteer::scale_color_paletteer_c("ggthemes::Red-Blue Diverging") +
  theme_bw() +
  theme(
    axis.title = et(color = "black", size = 72),
    axis.text = et(color = "black", size = 60),
    # plot.subtitle = et(color = "black", size = 79),
    panel.grid.minor = eb(),
    panel.grid.major = eb(),
    # legend.position = "none",
    # plot.margin = margin(
    #   # Top margin
    #   t = 5,
    #   # Right margin
    #   r = 5,
    #   # Bottom margin
    #   b = 5,
    #   # Left margin
    #   l = 5 
    # )
  ) +
  labs(title = "")

ggsave(
  filename = "Results/Beta_Diversity_BrayCurtis_train_cont.pdf",
  height = 14,
  width = 15,
  units = "in"
)

# MDS vs taxon:  maaslin --------------------------------------------

maaslin.meta.updated <- maaslin.meta |> 
  rownames_to_column(var = "shotgunSeq_id") |> 
  left_join(beta_mds_df |> 
              select(shotgunSeq_id, grouped_md_score, MDS)) |> 
  column_to_rownames(var = "shotgunSeq_id") |> 
  as.data.frame()

## maaslin at genus level --------------------------------------------------

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_genus,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Genus_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("grouped_md_score"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)

### Run Maaslin2 with Covariates ----------
set.seed(123)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_genus,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Genus_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "grouped_md_score",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)

## maaslin at species level --------------------------------------------------

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_sp,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Species_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("grouped_md_score"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)

### Run Maaslin2 with Covariates ------
set.seed(123)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_sp,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Species_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "grouped_md_score",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)

## maaslin at family level --------------------------------------------------

### Run Maaslin without covariates ----
set.seed(123)

maaslin_no_covariates <- Maaslin2(
  input_data = maaslin_mat_fa,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Family_Base",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c("grouped_md_score"),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)

### Run Maaslin2 with Covariates -----------
set.seed(123)

maaslin_model <- Maaslin2(
  input_data = maaslin_mat_fa,
  input_metadata = maaslin.meta.updated,
  output = "mds_metagenome/maaslin/Family_Covariate",
  min_abundance = 0.001,
  # At least 0.1% abundance
  min_prevalence = 0.10,
  # Taxa found in at least 10% of samples
  min_variance = -Inf,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  max_significance = 0.05,
  # p.adj <= 0.05 (qval = padjust)
  random_effects = NULL,
  fixed_effects = c(
    "grouped_md_score",
    "cci_total_sc",
    "sofa_score_total"
  ),
  correction = "BH",
  standardize = TRUE,
  cores = 12,
  plot_heatmap = TRUE,
  plot_scatter = TRUE,
  heatmap_first_n = 50,
  reference = c("grouped_md_score,High Score")
)


# MDS class coxPH  ------------------------------------
## training set -----------------------------------
cox_df <- cox_df |> 
  mutate(grouped_md_score = factor(grouped_md_score,
                                   levels = c("Low Score",
                                              "High Score")))

mds_group.cox <- coxph(
  Surv(cox_df$surv_days, cox_df$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    grouped_md_score,
  data = cox_df
)

mds_group.cox.tbl <- mds_group.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**MDS as group CoxPH**")

mds_group.cox.tbl

gt::gtsave(gtsummary::as_gt(mds_group.cox.tbl), file = "Results/msd_grp.training_coxph.png")

## test set -----------------------------------

cox_df_vc <- cox_df_vc |> 
  mutate(grouped_md_score = factor(grouped_md_score,
                                   levels = c("Low Score",
                                              "High Score")))

mds_group.cox <- coxph(
  Surv(cox_df_vc$surv_days, cox_df_vc$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    grouped_md_score,
  data = cox_df_vc
)

mds_group.cox.tbl <- mds_group.cox %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**MDS as group CoxPH**")

mds_group.cox.tbl

gt::gtsave(gtsummary::as_gt(mds_group.cox.tbl), file = "Results/msd_grp.test_coxph.png")

# new cut point for Enteros ----------------------------------------------------

library(cutpointr)
safe_cutpointr <- possibly(.f = cutpointr, otherwise = "Error")

cutpoints_entero <-
  phylo_rel_abd %>%
  mutate(thirtyday_mortality_overall_class = if_else(thirtyday_mortality_overall == "Survivor",
                                                     0, 1)) |> 
    group_by(organism) |> 
  group_map(
    ~ safe_cutpointr(
      .,
      pctseqs,
      thirtyday_mortality_overall_class,
      method = maximize_metric,
      metric = youden,
      pos_class = 0,
      neg_class = 1,
      boot_runs = 100,
      use_midpoints = TRUE,
      na.rm = T
    )  |> mutate(organism = unique(pull(.x, organism))),
    .keep = TRUE
  ) 

cutpoints_entero_res <- cutpoints_entero |> 
  map_df(as_tibble) |> 
  relocate(organism, optimal_cutpoint)

write_csv(cutpoints_entero_res, "Results/cutpoints_entero_res.csv")

# cutpoints_entero_res <- read_csv("Results/cutpoints_entero_res.csv")

## ROC curve with new cutoffs -------------------

phylo_entero <- phylo_rel_abd %>%
  mutate(thirtyday_mortality_overall_class = if_else(thirtyday_mortality_overall == "Survivor",
                                                     0, 1)) |> 
  spread(organism, pctseqs) 

### Enterococcus -----
pROC_obj_Enterococcus <- pROC::roc(
  phylo_entero$thirtyday_mortality_overall_class,
  phylo_entero$Enterococcus,
  smoothed = TRUE,
  ci = TRUE,
  plot = FALSE,
  auc.polygon = TRUE,
  best.method = TRUE,
  print.auc = TRUE,
  print.auc.col = "black",
  col = "#2F472F",
  auc.polygon.border = "black",
  auc.polygon.col = "gray65",
  print.thres.best.method = "youden"
)

coordinates_Enterococcus <-
  coords(pROC_obj_Enterococcus,
         "best",
         ret = c("acc", "threshold", "sens", "spec", "ppv", "npv")
  )

cairo_pdf(
  "Results/ROC_curve_cutpoint_Enterococcus_train.pdf",
  width = 8,
  height = 6
)

pROC::roc(
  phylo_entero$thirtyday_mortality_overall_class,
  phylo_entero$Enterococcus,
  smoothed = TRUE,
  ci = TRUE,
  plot = TRUE,
  auc.polygon = TRUE,
  best.method = TRUE,
  print.auc = TRUE,
  print.auc.col = "black",
  col = "#2F472F",
  auc.polygon.border = "black",
  auc.polygon.col = "gray65",
  print.thres.best.method = "youden"
)

text(
  paste("ACC:", round(coordinates_Enterococcus$accuracy, 3) * 100, "%"),
  x = 0.5,
  y = 0.45,
  adj = 0
)
text(paste("PPV:", round(coordinates_Enterococcus$ppv, 2)),
     x = 0.5,
     y = 0.41,
     adj = 0
)
text(paste("NPV:", round(coordinates_Enterococcus$npv, 2)),
     x = 0.5,
     y = 0.37,
     adj = 0
)
text(paste("Threshold:", round(coordinates_Enterococcus$threshold, 4)),
     x = 0.5,
     y = 0.33,
     adj = 0
)

roc_plot_Enterococcus <- grDevices::recordPlot()

dev.off()

### Enterobacterales -------

pROC_obj_Enterobacterales <- pROC::roc(
  phylo_entero$thirtyday_mortality_overall_class,
  phylo_entero$Enterobacterales,
  smoothed = TRUE,
  ci = TRUE,
  plot = FALSE,
  auc.polygon = TRUE,
  best.method = TRUE,
  print.auc = TRUE,
  print.auc.col = "black",
  col = "#2F472F",
  auc.polygon.border = "black",
  auc.polygon.col = "gray65",
  print.thres.best.method = "youden"
)

coordinates_Enterobacterales <-
  coords(pROC_obj_Enterobacterales,
         "best",
         ret = c("acc", "threshold", "sens", "spec", "ppv", "npv")
  )

cairo_pdf(
  "Results/ROC_curve_cutpoint_Enterobacterales_train.pdf",
  width = 8,
  height = 6
)

pROC::roc(
  phylo_entero$thirtyday_mortality_overall_class,
  phylo_entero$Enterobacterales,
  smoothed = TRUE,
  ci = TRUE,
  plot = TRUE,
  auc.polygon = TRUE,
  best.method = TRUE,
  print.auc = TRUE,
  print.auc.col = "black",
  col = "#2F472F",
  auc.polygon.border = "black",
  auc.polygon.col = "gray65",
  print.thres.best.method = "youden"
)

text(
  paste("ACC:", round(coordinates_Enterobacterales$accuracy, 3) * 100, "%"),
  x = 0.5,
  y = 0.45,
  adj = 0
)
text(paste("PPV:", round(coordinates_Enterobacterales$ppv, 2)),
     x = 0.5,
     y = 0.41,
     adj = 0
)
text(paste("NPV:", round(coordinates_Enterobacterales$npv, 2)),
     x = 0.5,
     y = 0.37,
     adj = 0
)
text(paste("Threshold:", round(coordinates_Enterobacterales$threshold, 4)),
     x = 0.5,
     y = 0.33,
     adj = 0
)

roc_plot_Enterococcus <- grDevices::recordPlot()

dev.off()

## see distribution with new cutoffs ----

phylo_cutp <- phylo_rel_abd %>%
  mutate(thirtyday_mortality_overall_class = if_else(thirtyday_mortality_overall == "Survivor",
                                                     0, 1)) |> 
  left_join(cutpoints_entero_res |> 
              select(organism, optimal_cutpoint, direction)) |> 
  left_join(cox_df |> 
              select(shotgunSeq_id, MDS, grouped_md_score))

phylo_cutp |> 
  mutate(entero_class = if_else(pctseqs <= optimal_cutpoint,
                                "Low", "High")) |> 
  dplyr::count(organism, thirtyday_mortality_overall, entero_class) |> 
  ggplot(aes(thirtyday_mortality_overall, entero_class)) +
  geom_tile(aes(fill = n)) +
  facet_grid(. ~ organism) +
  geom_text(aes(label = n), color = "white", size = 12) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange") +
  theme_bw()

ggsave("Results/cutpoints_entero_res.hm.pdf", height = 4.5, width = 9.5)

phylo_cutp |> 
  mutate(entero_class = if_else(pctseqs <= optimal_cutpoint,
                                "Low", "High")) |> 
  dplyr::count(organism, grouped_md_score, entero_class) |> 
  ggplot(aes(grouped_md_score, entero_class)) +
  geom_tile(aes(fill = n)) +
  facet_grid(. ~ organism) +
  geom_text(aes(label = n), color = "white", size = 12) +
  paletteer::scale_fill_paletteer_c("ggthemes::Orange") +
  theme_bw()

ggsave("Results/cutpoints_entero_res.MDS.hm.pdf", height = 4.5, width = 9.5)
  
phylo_rel_abd |> 
  ggplot(aes(thirtyday_mortality_overall, pctseqs)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.65) +
  facet_grid(. ~ organism) +
  geom_hline(data = cutpoints_entero_res,
             aes(yintercept = optimal_cutpoint),
             color = "red", linetype = 2 ) +
  theme_bw()

ggsave("Results/cutpoints_entero_res.pdf", height = 6.5, width = 9.5)

## CoxPH + KM curve with new cutoff -----------------------------------------------------

new_cox_df <- phylo_cutp |> 
  mutate(entero_class = if_else(pctseqs <= optimal_cutpoint,
                                paste0("No ",organism, " Domination"),
                                paste0(organism, " Domination"))) |> 
  select(shotgunSeq_id, organism, entero_class) |> 
  spread(organism, entero_class) |> 
  left_join(cox_df |> 
              select(shotgunSeq_id,surv_days,
                     `Charlson Comorbidity Index`,
                       `SOFA Score`,
                     thirtyday_mortality_overall_class)) |> 
  mutate(Enterobacterales = factor(Enterobacterales,
                                   levels = c("No Enterobacterales Domination","Enterobacterales Domination")),
         Enterococcus = factor(Enterococcus,
                                   levels = c("No Enterococcus Domination","Enterococcus Domination")))

### Enterobacterales: as domination or not, dichotomous ---------------
library(survival)

#### CoxPH ------------------
new_enteroba <- coxph(
  Surv(new_cox_df$surv_days, new_cox_df$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    Enterobacterales,
  data = new_cox_df
)

new_enteroba.tbl <- new_enteroba %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**Enterobacterales CoxPH**")

new_enteroba.tbl

gt::gtsave(gtsummary::as_gt(new_enteroba.tbl), file = "Results/Enterobacterales_newCP_coxph.png")

#### KM curve -----------------
surv_object <-
  Surv(new_cox_df$surv_days, new_cox_df$thirtyday_mortality_overall_class)

fit_Enterobac <- survfit(surv_object ~ Enterobacterales, data = new_cox_df)

ggs_Enterobac <- ggsurvplot(
  fit_Enterobac,
  data = new_cox_df,
  size = 1,
  palette = c( "#107b1c","#C45258"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  
  risk.table.height = 0.45,
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 3.25,
  pval.size = 3.5,
  tables.y.text = FALSE,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.25, color = "gray90"),
    axis.text.y = et(color = "black", size = 8),
    axis.title.y = et(color = "black")
  ),
  legend.labs = c("No Domination", "Enterobacterales Domination")
)

# Change table axis labels
ggs_Enterobac$table <-
  ggs_Enterobac$table + labs(x = NULL, y = NULL) + theme(plot.title = eb()) # risk table

ggs_Enterobac

pdf(
  "./Results/kaplan_meier_train_newEnterobac.pdf",
  height = 4,
  width = 6,
  onefile = FALSE
)
ggs_Enterobac
invisible(dev.off())

### Enterobacterales: as continuous variable ---------------

entero_cox_df <- phylo_rel_abd |> 
  select(shotgunSeq_id, organism, pctseqs) |> 
  spread(organism, pctseqs) |> 
  left_join(new_cox_df |> 
              select(-c(Enterobacterales, Enterococcus)))

#### CoxPH -------------
new_enteroba <- coxph(
  Surv(entero_cox_df$surv_days, entero_cox_df$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    Enterobacterales,
  data = entero_cox_df
)

new_enteroba.tbl <- new_enteroba %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**Enterobacterales CoxPH**")

new_enteroba.tbl

gt::gtsave(gtsummary::as_gt(new_enteroba.tbl), file = "Results/Enterobacterales_newCP_cont.coxph.png")

### Enterococcus: as domination or not, dichotomous  ------------------------

#### CoxPH ------------
new_enterocc <- coxph(
  Surv(new_cox_df$surv_days, new_cox_df$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    Enterococcus,
  data = new_cox_df
)

new_enterocc.tbl <- new_enterocc %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**Enterococcus CoxPH**")

new_enterocc.tbl

gt::gtsave(gtsummary::as_gt(new_enterocc.tbl), file = "Results/Enterococcus_newCP_coxph.png")

grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

#### KM curve ---------------
fit_Enteroco <- survfit(surv_object ~ Enterococcus, data = new_cox_df)

ggs_Enteroco <- ggsurvplot(
  fit_Enteroco,
  data = new_cox_df,
  size = 1,
  palette = c( "#107b1c","#C45258"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  tables.y.text = FALSE,
  risk.table.height = 0.45,
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 3.25,
  pval.size = 3.5,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.25, color = "gray90"),
    axis.text.y = et(color = "black", size = 8),
    axis.title.y = et(color = "black")
  ),
  legend.labs = c("No Domination", "Enterococcus Domination")
)

# Change table axis labels
ggs_Enteroco$table <-
  ggs_Enteroco$table + 
  labs(x = NULL, y = NULL) + 
  theme(plot.title = eb()) # risk table

ggs_Enteroco

ggsave(
  "./Results/kaplan_meier_train_newEnteroco.pdf",
  ggs_Enteroco,
  height = 4,
  width = 6,
)

### Enterococcus: as continuous variable ------

#### CoxPH ------------------
new_enterocc <- coxph(
  Surv(entero_cox_df$surv_days, entero_cox_df$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    Enterococcus,
  data = entero_cox_df
)

new_enterocc.tbl <- new_enterocc %>%
  gtsummary::tbl_regression(
    exp = TRUE,
    pvalue_fun = function(x) {
      if_else(is.na(x), NA_character_, if_else(
        x < 0.001,
        format(x,
               digits = 3, scientific = TRUE
        ),
        format(round(x, 3),
               scientific = F
        )
      ))
    }
  ) %>%
  gtsummary::modify_footnote(everything() ~ NA, abbreviation = TRUE) %>%
  gtsummary::modify_caption("**Enterococcus CoxPH**")

new_enterocc.tbl

gt::gtsave(gtsummary::as_gt(new_enterocc.tbl), file = "Results/Enterococcus_newCP_cont.coxph.png")


#### KM curve for MDS + enterococcus: training ------------------------------------

alt_cox_df <- new_cox_df |> 
  left_join(cox_df |> select(shotgunSeq_id, grouped_md_score)) |> 
  mutate(MDS_entero = if_else(
    grouped_md_score == "High Score" | Enterococcus == "Enterococcus Domination",
    "High MDS/Enterococcus Domination", "Neither"
  ))

surv_object <-
  Surv(alt_cox_df$surv_days, alt_cox_df$thirtyday_mortality_overall_class)

fit_Enteroco <- survfit(surv_object ~ MDS_entero, data = alt_cox_df)

ggs_Enteroco <- ggsurvplot(
  fit_Enteroco,
  data = alt_cox_df,
  size = 1,
  palette = c( "#C45258", "#107b1c"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  tables.y.text = FALSE,
  risk.table.height = 0.45,
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 3.25,
  pval.size = 3.5,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.25, color = "gray90"),
    axis.text.y = et(color = "black", size = 8),
    axis.title.y = et(color = "black")
  ) ,
  legend.labs = c("High MDS/Enterococcus Domination", "Neither")
)

# Change table axis labels
ggs_Enteroco$table <-
  ggs_Enteroco$table + 
  labs(x = NULL, y = NULL) + 
  theme(plot.title = eb()) # risk table

ggs_Enteroco

ggsave(
  "./Results/kaplan_meier_train_highMDS_entero.pdf",
  ggs_Enteroco,
  height = 4,
  width = 6,
)

#### KM curve for MDS + enterococcus: test/validation ---------------------------

alt_cox_vc_entero <- cox_df_vc_entero |> 
  left_join(cox_df |> select(shotgunSeq_id, grouped_md_score)) |> 
  mutate(MDS_entero = if_else(
    grouped_md_score == "High Score" | Enterococcus == "Enterococcus Domination",
    "High MDS/Enterococcus Domination", "Neither"
  ))

surv_object <-
  Surv(alt_cox_vc_entero$surv_days, alt_cox_vc_entero$thirtyday_mortality_overall_class)

fit_Enteroco <- survfit(surv_object ~ MDS_entero, data = alt_cox_vc_entero)

ggs_Enteroco <- ggsurvplot(
  fit_Enteroco,
  data = alt_cox_vc_entero,
  size = 1,
  palette = c( "#C45258", "#107b1c"),
  xlab = "Days from Admission",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = "abs_pct",
  legend = "bottom",
  tables.y.text = FALSE,
  risk.table.height = 0.45,
  risk.table.y.text.col = TRUE,
  risk.table.fontsize = 3.25,
  pval.size = 3.5,
  ggtheme = theme_test() + theme(
    panel.grid.major = el(linewidth = 0.25, color = "gray90"),
    axis.text.y = et(color = "black", size = 8),
    axis.title.y = et(color = "black")
  ) ,
  legend.labs = c("High MDS/Enterococcus Domination", "Neither")
)

# Change table axis labels
ggs_Enteroco$table <-
  ggs_Enteroco$table + 
  labs(x = NULL, y = NULL) + 
  theme(plot.title = eb()) # risk table

ggs_Enteroco

ggsave(
  "./Results/kaplan_meier_test_highMDS_entero.pdf",
  ggs_Enteroco,
  height = 4,
  width = 6,
)

### Enteroccus vs mds ------------------------------------------------

entero_pct_mds <- cox_df |> 
  left_join(phylo_cutp |> 
              filter(organism == "Enterococcus") |> 
              select(shotgunSeq_id, pctseqs, organism))

cor.test(entero_pct_mds$MDS,
         entero_pct_mds$pctseqs,
         method = "spearman")
# rho = 0.3136491

library(ggpmisc)

entero_pct_mds |> 
  ggplot(aes(MDS, pctseqs)) +
  geom_point(alpha = 0.65) +
  stat_poly_line() +
  # stat_poly_eq(label.x = 1, use_label(c("eq", "R2", "p"))) +
  stat_poly_eq(label.x = 1, use_label(c("R2", "p"))) +
  theme_bw() +
  labs(x="MDS", y = "Enterococcus abundance")

ggsave("Results/enteroc_mds.liner.pdf", height = 6.5, width = 7.5)

# cox.zph assumption check ------------------------------------------------

library(survminer)
library(survival)
library(gridExtra)

cox_all <- cox_df |> 
  select(-c(`Enterococcus Domination`,`Enterobacterales Domination`)) |> 
  left_join(mmp_df |> 
              select(metabolomicsID, mmp_score)) |> 
  left_join(new_cox_df |> 
              select(shotgunSeq_id,
                     `Enterococcus Domination` = Enterococcus,
                     `Enterobacterales Domination` = Enterobacterales)) |> 
  left_join(phylo_entero |>
              select(shotgunSeq_id, Enterobacterales, Enterococcus))

base.chosen.cox <- coxph(
  Surv(cox_all$surv_days, cox_all$thirtyday_mortality_overall_class) ~
    `Charlson Comorbidity Index` +
    `SOFA Score` +
    # `Time to stool sample` +  # this will violate the ph assumptions
    `Shannon Diversity` +
    mmp_score +
    `Enterococcus Domination` +
    `Enterobacterales Domination` +
    Enterococcus +
    Enterobacterales +
    `MDS`
  ,
  data = cox_all
)

base.chosen.test <- cox.zph(base.chosen.cox)
base.chosen_zph <- ggcoxzph(base.chosen.test)

ggsave("Results/base.chosen.test.not2stool.pdf", 
       arrangeGrob(grobs = base.chosen_zph, ncol = 2),
       height = 18.5, width = 10.5)

ggcoxdiagnostics(base.chosen.cox, type = "dfbeta", linear.predictions = FALSE)
ggsave("Results/base.chosen.test_diag_dfbeta.not2stool.pdf",height = 15.5, width = 14.5)

ggcoxdiagnostics(base.chosen.cox, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
ggsave("Results/base.chosen.test_diag_deviance.not2stool.pdf",height = 6.5, width = 7)

# deltaRMST of KM curve of enterococcus and MDS in the training / validation cohort -----

## training ----------------------------------------------------------------

rmst_mds <-
  survRM2::rmst2(
    cox_df$surv_days,
    cox_df$thirtyday_mortality_overall_class,
    factor(
      cox_df$grouped_md_score,
      levels = c("Low Score", "High Score"),
      labels = c(1, 0)                      # Low Score = 1, High Score = 0, due to area = Low Score - High Score
    ),
    tau = 30
  )

rmst_mds[["unadjusted.result"]][1] # 12.28399

rmst_Enterococcus <-
  survRM2::rmst2(
    new_cox_df$surv_days,
    new_cox_df$thirtyday_mortality_overall_class,
    factor(
      new_cox_df$Enterococcus,
      levels = c("Enterococcus Domination", "No Enterococcus Domination"),
      labels = c(1, 0)  
    ),
    tau = 30
  )

rmst_Enterococcus[["unadjusted.result"]][1] # -5.986439

## testing cohort -------

rmst_mds_vc <-
  survRM2::rmst2(
    cox_df_vc$surv_days,
    cox_df_vc$thirtyday_mortality_overall_class,
    factor(
      cox_df_vc$grouped_md_score,
      levels = c("Low Score", "High Score"),
      labels = c(1, 0)                      # Low Score = 1, High Score = 0, due to area = Low Score - High Score
    ),
    tau = 30
  )

rmst_mds_vc[["unadjusted.result"]][1] # 8.241453

cox_df_vc <- metaphlan |> 
  filter(! shotgunSeq_id %in% cox_df$shotgunSeq_id,
         genLab == "Enterococcus") |> 
  select(shotgun_seq_id = shotgunSeq_id, organism = genLab, pctseqs) |> 
  right_join(cox_df_vc) |> 
  replace_na(list(organism = "Enterococcus",
                  pctseqs = 0))

cox_df_vc_entero <- cox_df_vc |> 
  left_join(cutpoints_entero_res |> 
              select(organism, optimal_cutpoint)) |> 
  mutate(Enterococcus = if_else(pctseqs <= optimal_cutpoint,
                                paste0("No ",organism, " Domination"),
                                paste0(organism, " Domination"))) 
  

rmst_Enterococcus_vc <-
  survRM2::rmst2(
    cox_df_vc_entero$surv_days,
    cox_df_vc_entero$thirtyday_mortality_overall_class,
    factor(
      cox_df_vc_entero$Enterococcus,
      levels = c("Enterococcus Domination", "No Enterococcus Domination"),
      labels = c(1, 0)   
    ),
    tau = 30
  )

rmst_Enterococcus_vc[["unadjusted.result"]][1] # 1.329365


# qual volcano plot with 13 metabolites -----------------------------------

library(EnhancedVolcano)

comps <- tolower( c("Kynurenine",
                    "Glycocholic Acid",
                    "Succinate",
                    "Tryptamine",
                    "Kynurenic Acid",
                    "Tyrosine",
                    "Cholic Acid",
                    "Serotonin",
                    "Indole-3-Propionate",
                    "Isodeoxycholic Acid",
                    "Indole-3-Lactate",
                    "Butyrate",
                    "Indole-3-Carboxaldehyde"))

kval.shape <- c(ifelse(rownames(qual_tot) %in% comps, 17, 16))
names(kval.shape)[kval.shape == 17] <- "MDS 13 metabolites"
names(kval.shape)[kval.shape != 17] <- "Other metabolites"

EnhancedVolcano(
  qual_tot,
  lab = rownames(qual_tot),
  selectLab = comps,
  x = "log2fc_val",
  y = "p.adj",
  title = NULL,
  pCutoff = 0.1,
  FCcutoff = 0.75,
  pointSize = c(ifelse(rownames(qual_tot) %in% comps, 8, 5)),
  labSize = 8,
  axisLabSize = 32,
  # labCol = volcano_labcol$color,
  shapeCustom = kval.shape,
  caption = NULL,
  colAlpha = 0.65,
  col = c("gray85", c("grey40", "grey10", "#F27DFA")),
  legendPosition = "bottom",
  legendLabels = c(
    expression(q > 0.1 * ";" ~ Log[2] ~ FC < "\u00B1" * 0.75),
    expression(q > 0.1 * ";" ~ Log[2] ~ FC >= "\u00B1" *
                 0.75),
    expression(q <= 0.1 * ";" ~ Log[2] ~ FC < "\u00B1" *
                 0.75),
    expression(q <= 0.1 * ";" ~ Log[2] ~ FC >= "\u00B1" *
                 0.75)
  ),
  legendLabSize = 14,
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  arrowheads = TRUE,
  gridlines.minor = FALSE,
  gridlines.major = FALSE,
  max.overlaps = Inf,
  # min.segment.length = 0.5
) +
  theme(
    axis.text = et(color = "black"),
    legend.text = et(hjust = 0, size = 18),
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) +
  labs(subtitle = NULL) +
  annotate(
    "segment",
    x = 0.8,
    xend = 2.5,
    y = 1.5,
    yend = 1.5,
    arrow = arrow(),
    size = 2,
    color = ggsci::pal_lancet(palette = "lanonc")(2)[1]
  ) +
  annotate(
    "text",
    x = 0.8,
    y = 1.6,
    hjust = 0,
    label = "Survivor",
    size = 12,
    color = ggsci::pal_lancet(palette = "lanonc")(2)[1]
  ) +
  annotate(
    "rect",
    xmin = 0.75,
    xmax = Inf,
    ymin = -log(0.1, base = 10),
    ymax = Inf,
    alpha = .1,
    fill = ggsci::pal_lancet(palette = "lanonc")(2)[1]
  ) +
  annotate(
    "segment",
    x = -0.8,
    xend = -2.5,
    y = 1.5,
    yend = 1.5,
    arrow = arrow(),
    size = 2,
    color = ggsci::pal_lancet(palette = "lanonc")(2)[2]
  ) +
  annotate(
    "text",
    x = -1.55,
    y = 1.6,
    hjust = 0.5,
    label = "Non-Survivor",
    size = 12,
    color = ggsci::pal_lancet(palette = "lanonc")(2)[2]
  ) +
  annotate(
    "rect",
    xmin = -0.75,
    xmax = -Inf,
    ymin = -log(0.1, base = 10),
    ymax = Inf,
    alpha = .1,
    fill = ggsci::pal_lancet(palette = "lanonc")(2)[2]
  ) +
  guides(
    color = guide_legend(nrow = 4),
    shape = guide_legend(nrow = 4)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 1.6),
    breaks = seq(0, 2, 0.5)
  )

ggsave(
  filename = "Results/Qual_Metab_Volcano_30_Day_Mortality_train_13MDS.pdf",
  width = 24,
  height = 14
)
