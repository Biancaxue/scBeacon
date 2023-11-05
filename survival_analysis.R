# run CIBERSORT
# CIBERSORT R environment needs to be set up before running CIBERSORT function
results <- CIBERSORT("signature matrix file", 
                     "bulk RNA-seq file", perm=100, QN=FALSE, absolute=TRUE, abs_method='no.sumto1')

# surivial analysis
library(mclust)
library(Matrix)
library(survival)
library(survminer)
library(data.table)
subtypes <- read.delim("collected_subtype_annotations.tsv", row.names = 1)
subtypes[,"TCGA_disease"] <- unlist(lapply(as.character(subtypes[,"disease_subtype"]), function(x) unlist(strsplit(x, split=".", fixed=T))[1]))
LGG.samples <- rownames(read.delim("Cibersort_results_EBI-full_TCGA-LGG.txt"))
GBM.samples <- rownames(read.delim("Cibersort_results_EBI-full_TCGA-GBM.txt"))
subtypes[LGG.samples,"TCGA_disease"] <- rep("LGG", length(LGG.samples))
subtypes[GBM.samples,"TCGA_disease"] <- rep("GBM", length(GBM.samples))
survival <- read.delim("Survival_SupplementalTable_S1_20171025_xena_sp", row.names=1)
#load deconvolution results
files <- list.files(path=".", pattern="^Cibersort_results_EBI-full.*txt")
files <- unlist(lapply(names(table(subtypes[,"TCGA_disease"])), function(x) files[grep(x, files)]))
#subtype file does not include DLBC 
bimodel.files <- list.files(path=".", pattern="EBI-full_.*decision-boundary.tsv")
cibersort.files <- list.files(path=".", pattern="Cibersort_results")

# assign patients groups for signatures that passed bimodality model test
patients.group.list <- list()
for (file in files){
  message(file)
  cancer <- gsub("Cibersort_results_EBI-full_TCGA-", "", file)
  cancer <- gsub(".txt", "", cancer)
  cibersort_result <- read.delim(file)
  sigs <- colnames(cibersort_result)[1:as.numeric(ncol(cibersort_result)-4)]
  cibersort_result <- cibersort_result[, sigs]
  bimodel.sig <- read.table(bimodel.files[grep(cancer, cibersort.files)])
  subtype <- subtypes[grep(cancer, subtypes[,"disease_subtype"]),"disease_subtype"]
  intersect.samples <- intersect(names(subtype), rownames(cibersort_result))
  subtype <- subtype[intersect.samples]
  patients.group <- list()
  for (sig in bimodel.sig[,"sigs"]){
    boundary <- bimodel.sig[which(bimodel.sig[,"sigs"] == sig), "boundary"]
    range1.lower <- bimodel.sig[which(bimodel.sig[,"sigs"] == sig), "mean1"]+0*bimodel.sig[which(bimodel.sig[,"sigs"] == sig), "sd1"]
    if (range1.lower > boundary){
      range1.lower = boundary}
    range1 <- c(min(cibersort_result[,paste("X", sig, sep="")]), range1.lower)
    range2.upper <- bimodel.sig[which(bimodel.sig[,"sigs"] == sig), "mean2"]-0*bimodel.sig[which(bimodel.sig[,"sigs"] == sig), "sd2"]
    if (range2.upper < boundary){
      range2.upper = boundary}
    range2 <- c(range2.upper, max(cibersort_result[,paste("X", sig, sep="")]))
    sig <- paste("X", sig, sep="")
    patients.up <- rownames(cibersort_result)[cibersort_result[,sig] %inrange% range2]
    patients.down <- rownames(cibersort_result)[cibersort_result[,sig] %inrange% range1]
    common <- intersect(patients.up, patients.down)
    patients.up <- setdiff(patients.up, common)
    patients.down <- setdiff(patients.down, common)
    patient_group <- matrix(NA, ncol=2, nrow=length(c(patients.up, patients.down)), dimnames=list(c(patients.up, patients.down), c("group", "subtype")))
    patient_group[patients.up, "group"] <- "up"
    patient_group[patients.down, "group"] <- "down"
    patient_group[, "subtype"] <- subtypes[rownames(patient_group), "disease_subtype"]
    patient_group <- patient_group[!is.na(patient_group[,"subtype"]),]
    patients.group[[sig]] <- patient_group
  }
  patients.group.list[[file]] <- patients.group
}
save(patients.group.list, file="EBI-full_TCGA_bimodal_patients_group_assignment_sig217.rda")

# assign patients groups for signatures that didn't pass bimodality model test 
patients.group.list <- list()
for (file in cibersort.files){
  message(file)
  cancer <- gsub("Cibersort_results_EBI-full_TCGA-", "", file)
  cancer <- gsub(".txt", "", cancer)
  cibersort_result <- read.delim(file)
  sigs <- colnames(cibersort_result)[1:as.numeric(ncol(cibersort_result)-4)]
  cibersort_result <- cibersort_result[, sigs]
  bimodel.sig <- read.table(bimodel.files[grep(cancer, cibersort.files)])
  bimodel.sigs <- unlist(lapply(bimodel.sig[,"sigs"], function(x) paste("X", x, sep="")))
  leftover.sig <- sigs[!sigs %in% bimodel.sigs]
  subtype <- subtypes[grep(cancer, subtypes[,"disease_subtype"]),"disease_subtype"]
  intersect.samples <- intersect(names(subtype), rownames(cibersort_result))
  subtype <- subtype[intersect.samples]
  patients.group <- list()
  for (sig in leftover.sig){
    percent.zeros <- length(which(cibersort_result[,sig] == 0))/nrow(cibersort_result)
    if (percent.zeros < 0.5) {
      patients.sorted <- rownames(cibersort_result)[order(cibersort_result[,sig], decreasing = T)]
      patients.up <- head(patients.sorted, floor(length(patients.sorted)*0.5))
      patients.down <- tail(patients.sorted, floor(length(patients.sorted)*0.5))
    } else {
      patients.up <- rownames(cibersort_result)[which(cibersort_result[,sig] > 0)]
      patients.down <- rownames(cibersort_result)[which(cibersort_result[,sig] == 0)]
    }
    patients.group.df <- matrix(NA, ncol=2, nrow=length(c(patients.up, patients.down)), dimnames=list(c(patients.up, patients.down), c("group", "subtype")))
    patients.group.df[patients.up, "group"] <- "up"
    patients.group.df[patients.down, "group"] <- "down"
    patients.group.df[, "subtype"] <- subtypes[rownames(patients.group.df), "disease_subtype"]
    patients.group.df <- patients.group.df[!is.na(patients.group.df[,"subtype"]),]
    patients.group[[sig]] <- patients.group.df
  }
  patients.group.list[[file]] <- patients.group
}
save(patients.group.list, file="EBI-full_TCGA_bimodal_patients_group_assignment_sig217_not_passing_bimodality.rda")

# extract signatures that passed bimodal analysis and have less than 10 samples in up/down group
patients_group <- get(load("EBI-full_TCGA_bimodal_patients_group_assignment_sig217.rda"))
sample_number_list <- list()
for (file in names(patients_group)){
  cancer <- unlist(strsplit(unlist(strsplit(file, split="_"))[4], split=".", fixed=T))[1]
  for (sig in names(patients_group[[file]])){
    n <- table(patients_group[[file]][[sig]][,"group"])
    sample_number_list[[cancer]][[sig]] <- n
  }
}
sample_number_list <- as.data.frame(unlist(sample_number_list))
sample_number_list <- sample_number_list %>% dplyr::mutate(group=unlist(lapply(rownames(sample_number_list), function(x) unlist(strsplit(x, split=".", fixed=T))[3])))
sample_number_list["cancer"] <- unlist(lapply(rownames(sample_number_list), function(i) unlist(strsplit(i, split=".", fixed=T))[1]))
sample_number_list["sig"] <- unlist(lapply(rownames(sample_number_list), function(i) unlist(strsplit(i, split=".", fixed=T))[2]))
sample_number_list <- sample_number_list %>% dplyr::rename(number_of_samples="unlist(sample_number_list)")
sample_number_list %>% dplyr::filter(number_of_samples < 10) %>% dplyr::select(c("cancer", "sig")) %>% unique

# assign patients groups for signatures that passed bimodality model test and have less than 10 samples in up/down group
patients.group.list <- list()
for (cancer_type in unique(sample_number_list %>% dplyr::filter(number_of_samples < 10) %>% dplyr::pull(cancer))){
  message(cancer_type)
  cibersort_result <- read.delim(paste0("Cibersort_results_EBI-full_", cancer_type, ".txt"))
  sigs <- unique(sample_number_list %>% dplyr::filter(number_of_samples < 10) %>% dplyr::filter(cancer == cancer_type) %>% dplyr::pull(sig))
  cibersort_result <- cibersort_result[sigs]
  subtype <- subtypes[grep(unlist(strsplit(cancer_type, split="-"))[2], subtypes[,"disease_subtype"]),"disease_subtype"]
  intersect.samples <- intersect(names(subtype), rownames(cibersort_result))
  subtype <- subtype[intersect.samples]
  patients.group <- list()
  for (sig in sigs){
    percent.zeros <- length(which(cibersort_result[,sig] == 0))/nrow(cibersort_result)
    if (percent.zeros < 0.5) {
      patients.sorted <- rownames(cibersort_result)[order(cibersort_result[,sig], decreasing = T)]
      patients.up <- head(patients.sorted, floor(length(patients.sorted)*0.5))
      patients.down <- tail(patients.sorted, floor(length(patients.sorted)*0.5))
    } else {
      patients.up <- rownames(cibersort_result)[which(cibersort_result[,sig] > 0)]
      patients.down <- rownames(cibersort_result)[which(cibersort_result[,sig] == 0)]
    }
    patients.group.df <- matrix(NA, ncol=2, nrow=length(c(patients.up, patients.down)), dimnames=list(c(patients.up, patients.down), c("group", "subtype")))
    patients.group.df[patients.up, "group"] <- "up"
    patients.group.df[patients.down, "group"] <- "down"
    patients.group.df[, "subtype"] <- subtypes[rownames(patients.group.df), "disease_subtype"]
    patients.group.df <- patients.group.df[!is.na(patients.group.df[,"subtype"]),]
    patients.group[[sig]] <- patients.group.df
  }
  patients.group.list[[paste0("Cibersort_results_EBI-full_", cancer_type, ".txt")]] <- patients.group
}
save(patients.group.list, file="EBI-full_TCGA_bimodal_patients_group_assignment_sig217_pass_bimodality_small_group_number_median_split.rda")


# survival analysis using the three patients groupings
pvalue_list_uni <- list()
pvalue_list_multi <- list()
hazard_ratio_uni <- list()
hazard_ratio_multi <- list()
survival_index_list <- list()
coeff_list_uni <- list()
coeff_list_multi <- list()
for (file in names(patients_group)){
  splots <- list()
  cancer <- gsub("Cibersort_results_EBI-full_TCGA-", "", file)
  cancer <- gsub(".txt", "", cancer)
  message(cancer)
  coeff_list_uni[[cancer]] <- list()
  coeff_list_multi[[cancer]] <- list()
  for (sig in names(patients_group[[file]])){
    message(sig)
    df_for_plot <- as.data.frame(patients_group[[file]][[sig]])
    if (unique(is.na(unique(survival[which(survival[,"cancer.type.abbreviation"] %in% unlist(strsplit(cancer, split="-"))), "PFI"])))){  # if this cancer type does not have PFI data available, use OS survival
      df_for_plot[,"survival.time"] <- survival[unlist(lapply(rownames(df_for_plot), function(x) substr(x, 1, nchar(x)-1))), "OS.time"]
      df_for_plot[,"survival"] <- survival[unlist(lapply(rownames(df_for_plot), function(x) substr(x, 1, nchar(x)-1))), "OS"]
      survival_index <- "OS"
    } else {
      df_for_plot[,"survival.time"] <- survival[unlist(lapply(rownames(df_for_plot), function(x) substr(x, 1, nchar(x)-1))), "PFI.time"] # default use PFI survival
      df_for_plot[,"survival"] <- survival[unlist(lapply(rownames(df_for_plot), function(x) substr(x, 1, nchar(x)-1))), "PFI"]
      survival_index <- "PFI"
    }
    survival_index_list[[cancer]] <- survival_index
    df_for_plot <- df_for_plot[!is.na(df_for_plot[,"survival.time"]),]
    if (length(unique(df_for_plot[,"group"])) == 2 & length(which(df_for_plot[,"group"]=="up")) >= 10 & 
        length(which(df_for_plot[,"group"]=="down")) >= 10){
      # the univariate CoxPH model
      cox_res <- summary(coxph(Surv(survival.time, survival) ~ group, data = df_for_plot))
      # the multivariate CoxPH model, adjusting for subtype
      coxmv_res <- summary(coxph(Surv(survival.time, survival) ~ group + strata(subtype), data = df_for_plot))
      # # make survival plot
      # p2 <- ggsurvplot(
      #   fit = survfit(Surv(survival.time, survival) ~ group, data = df_for_plot),
      #   #pval = T,
      #   #conf.int = T,
      #   title = paste(cancer, sig, sep=" "),
      #   xlab = "Days",
      #   ylab = paste(c("survival probability - survival", " (", survival_index, ")"), collapse = ""))
      # # annotate survival plot with CoxPH results
      # p2$plot <- p2$plot +
      #   ggplot2::annotate(
      #     "text",
      #     x = Inf, y = Inf,
      #     vjust = 1, hjust = 1,
      #     label = paste0("univariate Cox proportional hazards model: \n HR = ",signif(cox_res$conf.int[,"exp(coef)"],digits=3)
      #                    ," (",signif(cox_res$conf.int[,"lower .95"],digits=3),",",signif(cox_res$conf.int[,"upper .95"],digits=3),")"
      #                    ,", p = ",signif(cox_res$sctest["pvalue"],digits=3),
      #                    "\n multivariate CoxPH, adjusting for subtype: \n HR = ",signif(coxmv_res$conf.int[,"exp(coef)"],digits=3)
      #                    ," (",signif(coxmv_res$conf.int[,"lower .95"],digits=3),",",signif(coxmv_res$conf.int[,"upper .95"],digits=3),")"
      #                    ,", p = ",signif(coxmv_res$sctest["pvalue"],digits=3),
      #                    "\n samples in up-group: ",length(which(df_for_plot[,"group"]=="up")),
      #                    "\n samples in down-group: ",length(which(df_for_plot[,"group"]=="down"))
      #     ),
      #     size = 4
      #   )
      # splots[[sig]] <- p2
      pvalue_list_uni[[cancer]][[sig]] <- cox_res$sctest
      pvalue_list_multi[[cancer]][[sig]] <- coxmv_res$sctest
      hazard_ratio_uni[[cancer]][[sig]] <- cox_res$conf.int
      hazard_ratio_multi[[cancer]][[sig]] <- coxmv_res$conf.int
    }
  }
  res <- arrange_ggsurvplots(splots, print = FALSE, ncol = 2, nrow = 2)
  # save the survival plot, naming should reflect the patients grouping method
  ggsave(paste("EBI-full-combined-cancer_survival-CoxPH_", cancer, "bimodal", ".pdf", sep=""), res, width=10, height=10)
  dev.off()
}
survival_index_list <- unlist(survival_index_list)
# save the output of survival analysis, naming should also reflex the patients grouping method
save(survival_index_list, file="survival_index_list_bimodal.rda")
save(pvalue_list_uni, file="pvalue_list_uni_bimodal.rda")
save(pvalue_list_multi, file="pvalue_list_multi_bimodal.rda")
save(hazard_ratio_uni, file="hazard_ratio_uni_bimodal.rda")
save(hazard_ratio_multi, file="hazard_ratio_multi_bimodal.rda")

