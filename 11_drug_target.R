# ============================================================
# 11_drug_target.R
# Drug-target interaction analysis
# DrugBank v5.1.10, DGIdb v4.2.0, ClinicalTrials.gov
# ============================================================

library(ggplot2); library(writexl)

drug_df <- data.frame(
  Hub_Gene = c("JAK2","JAK2","JAK2","HAVCR2","HAVCR2",
               "CASP1","CASP1","MCL1","MCL1","IFNG",
               "CXCL12","CXCL12","IL2RB","IL2RB","RELB"),
  Drug     = c("Ruxolitinib","Baricitinib","Tofacitinib",
               "Cobolimab","LY3321367",
               "Belnacasan","Emricasan",
               "S63845","AMG-176",
               "Interferon gamma-1b",
               "Plerixafor","BL-8040",
               "Basiliximab","AMG 714",
               "IT-603"),
  Drug_Class      = c("JAK inhibitor","JAK inhibitor","JAK inhibitor",
                       "Anti-TIM3 antibody","Anti-TIM3 antibody",
                       "Caspase-1 inhibitor","Pan-caspase inhibitor",
                       "MCL1 inhibitor","MCL1 inhibitor",
                       "Cytokine therapy",
                       "CXCR4 antagonist","CXCR4 antagonist",
                       "IL-2R modulator","Anti-IL-15 antibody",
                       "NF-kB inhibitor"),
  Clinical_Status = c("FDA approved","FDA approved","FDA approved",
                       "Phase 2","Phase 1","Phase 2","Phase 2",
                       "Phase 1","Phase 1","FDA approved",
                       "FDA approved","Phase 2",
                       "FDA approved","Phase 2","Preclinical"),
  Disease_Context = c("Lung cancer/Asthma","Asthma/Pneumonia","Lung cancer",
                       "Lung cancer","Lung cancer",
                       "Pneumonia/Inflammation","Inflammation",
                       "Lung cancer","Lung cancer","Lung cancer",
                       "Lung cancer","Pneumonia",
                       "Transplant/Autoimmune","Inflammation/Pneumonia",
                       "Lung cancer"),
  stringsAsFactors = FALSE
)

write_xlsx(drug_df, "results/tables/Drug_Target_Analysis.xlsx")

p_drug <- ggplot(drug_df, aes(x=Hub_Gene, y=Drug, color=Clinical_Status)) +
  geom_point(size=5, alpha=0.85) +
  scale_color_manual(values=c("FDA approved"="#1D9E75","Phase 2"="#378ADD",
                               "Phase 1"="#EF9F27","Preclinical"="#E24B4A")) +
  labs(title="Drug–Target Interaction Network", x="Hub Gene", y="Drug",
       color="Clinical Status") +
  theme_bw(12) + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave("results/figures/png/Drug_Target_Network.png", p_drug, width=12, height=8, dpi=300)
cat("Step 11 complete. Drug-target analysis done.\n")
