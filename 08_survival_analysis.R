# ============================================================
# 08_survival_analysis.R
# Kaplan-Meier + Cox regression on TCGA-LUAD and TCGA-LUSC
# ============================================================

library(survival); library(survminer); library(dplyr); library(writexl); library(TCGAbiolinks)

# Load TCGA expression data
load_tcga <- function(file_path) {
  expr <- read.table(file_path, header=TRUE, sep="\t",
                     row.names=1, check.names=FALSE)
  colnames(expr) <- substr(colnames(expr), 1, 12)
  expr <- expr[, !duplicated(colnames(expr))]
  return(expr)
}

tcga_luad_expr <- load_tcga("data/TCGA.LUAD.HiSeqV2.gz")
tcga_lusc_expr <- load_tcga("data/TCGA.LUSC.HiSeqV2.gz")

build_surv_df <- function(clinical_df) {
  surv <- data.frame(
    sample_id = clinical_df$submitter_id,
    OS.time   = as.numeric(clinical_df$days_to_death),
    OS        = as.numeric(clinical_df$vital_status == "Dead"))
  surv$OS.time[is.na(surv$OS.time)] <-
    as.numeric(clinical_df$days_to_last_follow_up[is.na(surv$OS.time)])
  surv %>% dplyr::filter(!is.na(OS.time) & OS.time > 0)
}

luad_clin <- GDCquery_clinic("TCGA-LUAD", type="clinical")
lusc_clin <- GDCquery_clinic("TCGA-LUSC", type="clinical")
surv_luad <- build_surv_df(luad_clin)
surv_lusc <- build_surv_df(lusc_clin)

run_km <- function(gene, expr_mat, surv_data, cohort_name) {
  if (!gene %in% rownames(expr_mat)) { cat(gene, "not found\n"); return(NULL) }
  common <- intersect(surv_data$sample_id, colnames(expr_mat))
  surv_sub <- surv_data[surv_data$sample_id %in% common, ]
  expr_vec  <- as.numeric(expr_mat[gene, surv_sub$sample_id])
  surv_sub$Group <- ifelse(expr_vec > median(expr_vec, na.rm=TRUE), "High", "Low")
  surv_sub$time  <- surv_sub$OS.time
  surv_sub$event <- surv_sub$OS
  cox_fit <- coxph(Surv(time, event) ~ Group, data=surv_sub)
  cox_s   <- summary(cox_fit)
  lr_test <- survdiff(Surv(time, event) ~ Group, data=surv_sub)
  p_val   <- round(1 - pchisq(lr_test$chisq, df=1), 4)
  hr      <- round(cox_s$conf.int[1, "exp(coef)"], 2)
  hr_lo   <- round(cox_s$conf.int[1, "lower .95"], 2)
  hr_hi   <- round(cox_s$conf.int[1, "upper .95"], 2)
  cat(gene, cohort_name, "HR:", hr, "(", hr_lo, "-", hr_hi, ") p:", p_val, "\n")
  km_fit <- survfit(Surv(time, event) ~ Group, data=surv_sub)
  p_km   <- ggsurvplot(km_fit, data=surv_sub, pval=TRUE, conf.int=TRUE,
                        risk.table=TRUE, risk.table.height=0.28,
                        palette=c("#E24B4A","#378ADD"), legend.labs=c("High","Low"),
                        title=paste0(gene, " — ", cohort_name,
                                     "\nHR=", hr, " (", hr_lo, "-", hr_hi, ") p=", p_val),
                        xlab="Time (days)", ylab="Overall Survival", ggtheme=theme_bw(12))
  ggsave(paste0("results/survival/KM_", gene, "_", cohort_name, ".png"),
         print(p_km), width=8, height=9, dpi=300)
  return(data.frame(Gene=gene, Cohort=cohort_name, HR=hr,
                    HR_low=hr_lo, HR_high=hr_hi, P_value=p_val))
}

hub_list <- readRDS("data/hub_genes.rds")
lc_hub   <- hub_list$lc

surv_results <- do.call(rbind, lapply(lc_hub, function(g) {
  rbind(run_km(g, tcga_luad_expr, surv_luad, "LUAD"),
        run_km(g, tcga_lusc_expr, surv_lusc, "LUSC"))
}))
surv_results <- surv_results[!is.null(surv_results), ]
write_xlsx(surv_results, "results/tables/Survival_Combined.xlsx")
cat("Step 08 complete. Survival analysis done.\n")
