# # library(tidyverse)
# # library(survminer)
# # library(survival) 
# source("scripts/utils/pfunc.R")
# 
# # 加载疾病状态
# diseases = read.csv("data/disease_with_pn.csv", row.names = "eid", check.names = FALSE)
# diseases = diseases[, c("Stroke", "Dementia", "Parkinson",
#                         "hypertension", "nonrheumatic aortic valve disorders",
#                         "cardiomyopathy", "atrioventricular and left bundle-branch block",
#                         "Atrial_Fibrillation", "Heart_Failure", "other cardiac arrhythmias",
#                         "Asthma", "COPD", "Chronic kidney disease",
#                         "non-insulin-dependent diabetes mellitus")]
# diseases[diseases == 1] = 2
# diseases[diseases == 0] = 1
# 
# # # 加载疾病时间数据，将其转化为月份
# # diseases_time = read.csv("data/disease_days.csv", row.names = "eid", check.names = FALSE)
# # diseases_time = diseases_time[, c("Stroke", "Dementia", "Parkinson",
# #                                  "hypertension", "nonrheumatic aortic valve disorders",
# #                                  "cardiomyopathy", "atrioventricular and left bundle-branch block",
# #                                  "Atrial_Fibrillation", "Heart_Failure", "other cardiac arrhythmias",
# #                                  "Asthma", "COPD", "Chronic kidney disease",
# #                                  "non-insulin-dependent diabetes mellitus")]
# # colnames(diseases_time) = paste0(colnames(diseases_time), "_time")
# # f = function(x) {ceiling(x/30)}
# # for (i in seq_along(diseases_time)) {
# #   for (j in seq_along(diseases_time[[i]])) {
# #     diseases_time[[i]][j] <- f(diseases_time[[i]][j])
# #   }
# # }
# diseases_time = read.csv("data/disease_months.csv", row.names = "eid", check.names = FALSE)
# 
# 
# # 加载表型状态，按照表型大小对表型进行分组
# pheno = read.csv("data/ukb_new_pheno/whole/ukb_new_pheno.csv", row.names = 1)
# colnames(pheno) = str_remove(colnames(pheno), pattern = "y_res_")
# 
# # 表型分组函数，输入为array
# divide_pheno_by_quantile = function(data = data,
#                                     probs=c(0, 0.2, 0.8, 1),
#                                     labels=c("lower", "middle", "upper")) {
#   # 计算分位数
#   quantiles <- quantile(data, probs = probs, na.rm = TRUE)
# 
#   # 使用cut函数分组
#   groups <- cut(data, breaks = quantiles, labels = labels, include.lowest = TRUE)
#   if (length(data) != length(groups)) {
#     error("输入和输出的array长度不一致")
#   }
#   return(groups)
# }
# 
# df_info = data.frame()
# for (i in seq_along(pheno)) {
#   organ_axis_name = names(pheno)[i]
#   organ_axis = pheno[i]
#   organ_axis = na.omit(organ_axis)
#   organ_axis_dividation = divide_pheno_by_quantile(organ_axis[[1]])
#   organ_axis = data.frame(row.names = rownames(organ_axis), organ_axis_name=organ_axis_dividation)
#   names(organ_axis) = organ_axis_name
# 
#   for (j in seq_along(diseases)) {
#     disease_name = names(diseases)[j]
#     disease = diseases[j]
#     disease = na.omit(disease)
#     disease_time = diseases_time[j]
#     disease_time = na.omit(disease_time)
#     df = inner_join(disease_time %>% rownames_to_column("eid"),
#                     disease %>% rownames_to_column("eid"),
#                     by="eid")
#     df = inner_join(df, organ_axis %>% rownames_to_column("eid"), by="eid")
#     n_disease_minus = nrow(df[(df[[3]] == -1), ])
#     n_disease_plus = nrow(df[(df[[3]] == 2), ])
#     n_nondisease = nrow(df[df[[3]] == 1, ])
# 
#     # 选择after的人进行分析
#     df = df[df[[2]] > 0, ]
# 
#     if (nrow(df) > 0) {
#       # png(paste0("plot/risk_accumulation/", colnames(df)[4], "_", colnames(df)[3], ".png"),
#       #     width = 1600, height = 1000)
# 
#       addbacktick = function(str) {paste0("`", str, "`")}
#       fit_formula = paste0("Surv(",
#                            addbacktick(colnames(df)[2]), ", ",
#                            addbacktick(colnames(df)[3]), ") ~ ",
#                            colnames(df)[4])
#       fit = survfit(as.formula(fit_formula),
#                     data = df)
# 
#       pval = survdiff(as.formula(fit_formula), data = df)$pvalue
#       survival_plot  = ggsurvplot(fit, data = df,
#                                   conf.int = TRUE, # 增加置信区间
#                                   fun = "cumhaz", # 绘制累计风险曲线
#                                   risk.table = TRUE,
#                                   ggtheme = theme_bw(),
#                                   pval = TRUE,
#       )
# 
#       g2 <- ggplotGrob(survival_plot$plot + theme(panel.grid = element_blank(),
#                                                   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#                                                   axis.title.x = element_blank(),
#                                                   axis.title.y = element_text(size = 20, face = "bold", vjust = -10)))
#       g3 <- ggplotGrob(survival_plot$table +
#                          theme(legend.position = "none") +
#                          theme(panel.grid = element_blank(),
#                                axis.text = element_text(size = 15, face = "bold", colour = "black"),
#                                axis.title.y = element_blank(),
#                                axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
#                                )
#                          )
#       min_ncol <- min(ncol(g2), ncol(g3))
#       g = gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
#       ggsave(filename = paste0("plot/risk_accumulation/", disease_name, "/", colnames(df)[4], "_", colnames(df)[3], ".png"),
#              width = 8, height = 6, plot = g)
# 
#       temp_df_info = data.frame(pheno = colnames(df)[4],
#                                 disease = colnames(df)[3],
#                                 n_disease_minus = n_disease_minus,
#                                 n_disease_plus = n_disease_plus,
#                                 n_nondisease = n_nondisease,
#                                 pval = pval)
#       df_info = bind_rows(df_info, temp_df_info)
#     }
#     else {
#       temp_df_info = data.frame(pheno = colnames(df)[4],
#                                 disease = colnames(df)[3],
#                                 n_disease_minus = n_disease_minus,
#                                 n_disease_plus = n_disease_plus,
#                                 n_nondisease = n_nondisease,
#                                 pval = NA)
#       df_info = bind_rows(df_info, temp_df_info)
#     }
# 
#   }
# }
# df_info["p_bonf"] = p.adjust(df_info$pval, method = "bonferroni")
# write.csv(df_info, "data/risk_accumulation.csv")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##organ pca disease accumulation 
# # 加载表型状态，按照表型大小对表型进行分组
# pheno = read.csv("data/organ_pca.csv", row.names = 1)
# 
# 
# df_info = data.frame()
# for (i in seq_along(pheno)) {
#   organ_pca_name = names(pheno)[i]
#   organ_pca = pheno[i]
#   organ_pca = na.omit(organ_pca)
#   organ_pca_dividation = divide_pheno_by_quantile(organ_pca[[1]])
#   organ_pca = data.frame(row.names = rownames(organ_pca), organ_pca_name=organ_pca_dividation)
#   names(organ_pca) = organ_pca_name
# 
#   for (j in seq_along(diseases)) {
#     disease_name = names(diseases)[j]
#     disease = diseases[j]
#     disease = na.omit(disease)
#     disease_time = diseases_time[j]
#     disease_time = na.omit(disease_time)
#     df = inner_join(disease_time %>% rownames_to_column("eid"),
#                     disease %>% rownames_to_column("eid"),
#                     by="eid")
#     df = inner_join(df, organ_pca %>% rownames_to_column("eid"), by="eid")
#     n_disease_minus = nrow(df[(df[[2]] < 0 & df[[3]] == 2), ])
#     n_disease_plus = nrow(df[(df[[2]] > 0 & df[[3]] == 2), ])
#     n_nondisease = nrow(df[df[[3]] == 1, ])
# 
# 
#     df = df[df[[2]] > 0, ]
#     if (nrow(df) > 0) {
#       # png(paste0("plot/risk_accumulation/", colnames(df)[4], "_", colnames(df)[3], ".png"),
#       #     width = 1600, height = 1000)
# 
#       addbacktick = function(str) {paste0("`", str, "`")}
#       fit_formula = paste0("Surv(",
#                            addbacktick(colnames(df)[2]), ", ",
#                            addbacktick(colnames(df)[3]), ") ~ ",
#                            colnames(df)[4])
#       fit = survfit(as.formula(fit_formula),
#                     data = df)
# 
#       pval = survdiff(as.formula(fit_formula), data = df)$pvalue
#       survival_plot  = ggsurvplot(fit, data = df,
#                                   conf.int = TRUE, # 增加置信区间
#                                   fun = "cumhaz", # 绘制累计风险曲线
#                                   risk.table = TRUE,
#                                   ggtheme = theme_bw(),
#       )
# 
#       g2 <- ggplotGrob(survival_plot$plot + theme(panel.grid = element_blank(),
#                                                   axis.text = element_text(size = 15, face = "bold", colour = "black"),
#                                                   axis.title.x = element_blank(),
#                                                   axis.title.y = element_text(size = 20, face = "bold", vjust = -10)))
#       g3 <- ggplotGrob(survival_plot$table +
#                          theme(legend.position = "none") +
#                          theme(panel.grid = element_blank(),
#                                axis.text = element_text(size = 15, face = "bold", colour = "black"),
#                                axis.title.y = element_blank(),
#                                axis.title.x = element_text(size = 20, face = "bold", colour = "black"),
#                          )
#       )
#       min_ncol <- min(ncol(g2), ncol(g3))
#       g = gridExtra::gtable_rbind(g2[, 1:min_ncol], g3[, 1:min_ncol], size="last")
#       ggsave(filename = paste0("test/risk_accumulation/", colnames(df)[4], "_", colnames(df)[3], ".png"),
#              width = 8, height = 6, plot = g)
# 
#       temp_df_info = data.frame(pheno = colnames(df)[4],
#                                 disease = colnames(df)[3],
#                                 n_disease_minus = n_disease_minus,
#                                 n_disease_plus = n_disease_plus,
#                                 n_nondisease = n_nondisease,
#                                 pval = pval)
#       df_info = bind_rows(df_info, temp_df_info)
#     }
#     else {
#       temp_df_info = data.frame(pheno = colnames(df)[4],
#                                 disease = colnames(df)[3],
#                                 n_disease_minus = n_disease_minus,
#                                 n_disease_plus = n_disease_plus,
#                                 n_nondisease = n_nondisease,
#                                 pval = NA)
#       df_info = bind_rows(df_info, temp_df_info)
#     }
# 
#   }
# }
# # df_info["p_bonf"] = p.adjust(df_info$pval, method = "bonferroni")
# # write.csv(df_info, "data/pca_risk_accumulation.csv")
# 








############################### organ axis disease accumulation################
library(tidyverse)
library(survminer)
library(survival) 
source("scripts/utils/pfunc.R")
source("scripts//utils/Func.R")
# 加载疾病状态
diseases = read.csv("data/disease_with_pn.csv", row.names = "eid", check.names = FALSE)
# diseases = diseases[, c("Stroke", "Dementia", "Parkinson", 
#                         "hypertension", "nonrheumatic aortic valve disorders", 
#                         "cardiomyopathy", "atrioventricular and left bundle-branch block", 
#                         "Atrial_Fibrillation", "Heart_Failure", "other cardiac arrhythmias", 
#                         "Asthma", "COPD", "Chronic kidney disease", 
#                         "non-insulin-dependent diabetes mellitus")]
diseases[diseases == 1] = 2
diseases[diseases == 0] = 1


diseases_time = read.csv("data/disease_months.csv", row.names = "eid", check.names = FALSE)
colnames(diseases_time) = paste0(colnames(diseases_time), "_month")
# # 加载疾病时间数据，将其转化为月份
# diseases_time = read.csv("data/disease_days.csv", row.names = "eid", check.names = FALSE)
# diseases_time = diseases_time[, c("Stroke", "Dementia", "Parkinson", 
#                                   "hypertension", "nonrheumatic aortic valve disorders", 
#                                   "cardiomyopathy", "atrioventricular and left bundle-branch block", 
#                                   "Atrial_Fibrillation", "Heart_Failure", "other cardiac arrhythmias", 
#                                   "Asthma", "COPD", "Chronic kidney disease", 
#                                   "non-insulin-dependent diabetes mellitus")]
# colnames(diseases_time) = paste0(colnames(diseases_time), "_time")
# f = function(x) {ceiling(x/30)}
# for (i in seq_along(diseases_time)) {
#   for (j in seq_along(diseases_time[[i]])) {
#     diseases_time[[i]][j] <- f(diseases_time[[i]][j])  
#   }
# }


# 加载表型状态，按照表型大小对表型进行分组
pheno = read.csv("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv", row.names = 1)
colnames(pheno) = str_remove(colnames(pheno), pattern = "y_res_")



df_info_whole = data.frame()
df_info = data.frame()
for (i in seq_along(pheno)) {
  print(i)
  organ_axis_name = names(pheno)[i]
  organ_axis = pheno[i]
  organ_axis = na.omit(organ_axis)
  organ_axis_dividation = divide_pheno_by_quantile(organ_axis[[1]])
  organ_axis = data.frame(row.names = rownames(organ_axis), organ_axis_name=organ_axis_dividation)
  names(organ_axis) = organ_axis_name
  
  for (j in seq_along(diseases)) {
    disease_name = names(diseases)[j]
    disease = diseases[j]
    disease = na.omit(disease)
    disease_time = diseases_time[j]
    disease_time = na.omit(disease_time)
    df = inner_join(disease_time %>% rownames_to_column("eid"), 
                    disease %>% rownames_to_column("eid"), 
                    by="eid")
    df = inner_join(df, organ_axis %>% rownames_to_column("eid"), by="eid")
    n_disease_minus = nrow(df[(df[[3]] == -1), ])
    n_disease_plus = nrow(df[(df[[3]] == 2), ])
    n_nondisease = nrow(df[df[[3]] == 1, ])
    
    # 选择after的人进行分析
    df = df[df[[2]] > 0, ]
    if (sum(df[3] == 2) >= 10) {
      if (nrow(df) > 0) {
        addbacktick = function(str) {paste0("`", str, "`")}
        fit_formula = paste0("Surv(", 
                             addbacktick(colnames(df)[2]), ", ", 
                             addbacktick(colnames(df)[3]), ") ~ ", 
                             colnames(df)[4])
        fit = survfit(as.formula(fit_formula), 
                      data = df)
        
        pval = survdiff(as.formula(fit_formula), data = df)$pvalue
        survival_plot  = ggsurvplot(fit, data = df,
                                    conf.int = TRUE, # 增加置信区间
                                    fun = "cumhaz", # 绘制累计风险曲线
                                    ggtheme = theme_bw(),
                                    title=str_glue("pval={pval}")
        )
        g2 <- ggplotGrob(survival_plot$plot + theme(panel.grid = element_blank(),
                                                    axis.text = element_text(size = 15, face = "bold", colour = "black"),
                                                    axis.title.x = element_text(size = 20, face = "bold"),
                                                    axis.title.y = element_text(size = 20, face = "bold")))


        ggsave(filename = paste0("plot/risk_accumulation/", disease_name, "/", colnames(df)[4], "_", colnames(df)[3], ".png"),
               width = 8, height = 6)
        # 
        temp_df_info = data.frame(pheno = colnames(df)[4], 
                                  disease = colnames(df)[3], 
                                  n_disease_minus = n_disease_minus, 
                                  n_disease_plus = n_disease_plus, 
                                  n_nondisease = n_nondisease, 
                                  pval = pval)
        df_info = bind_rows(df_info, temp_df_info)
      }
      else {
        temp_df_info = data.frame(pheno = colnames(df)[4], 
                                  disease = colnames(df)[3], 
                                  n_disease_minus = n_disease_minus, 
                                  n_disease_plus = n_disease_plus, 
                                  n_nondisease = n_nondisease, 
                                  pval = NA)
        df_info = bind_rows(df_info, temp_df_info)
      }
    }
  }
}
df_info_whole = rbind(df_info_whole, df_info)
df_info["p_bonf"] = p.adjust(df_info$pval, method = "bonferroni")
write.csv(df_info, "data/risk_accumulation.csv")


########################### organ pca disease accumulation #####################
# 加载表型状态，按照表型大小对表型进行分组
pheno = read.csv("data/organ_pca.csv", row.names = 1)

df_info = data.frame()
for (i in seq_along(pheno)) {
  organ_pca_name = names(pheno)[i]
  organ_pca = pheno[i]
  organ_pca = na.omit(organ_pca)
  organ_pca_dividation = divide_pheno_by_quantile(organ_pca[[1]])
  organ_pca = data.frame(row.names = rownames(organ_pca), organ_pca_name=organ_pca_dividation)
  names(organ_pca) = organ_pca_name
  
  for (j in seq_along(diseases)) {
    disease_name = names(diseases)[j]
    disease = diseases[j]
    disease = na.omit(disease)
    disease_time = diseases_time[j]
    disease_time = na.omit(disease_time)
    df = inner_join(disease_time %>% rownames_to_column("eid"), 
                    disease %>% rownames_to_column("eid"), 
                    by="eid")
    df = inner_join(df, organ_pca %>% rownames_to_column("eid"), by="eid")
    
    n_disease_minus = nrow(df[(df[[3]] == -1), ])
    n_disease_plus = nrow(df[(df[[3]] == 2), ])
    n_nondisease = nrow(df[df[[3]] == 1, ])
    
    # 选择after的人进行分析
    df = df[df[[2]] > 0, ]
    if (sum(df[3] == 2) >= 10) {
      if (nrow(df) > 0) {
        # png(paste0("plot/risk_accumulation/", colnames(df)[4], "_", colnames(df)[3], ".png"),
        #     width = 1600, height = 1000)
        
        addbacktick = function(str) {paste0("`", str, "`")}
        fit_formula = paste0("Surv(", 
                             addbacktick(colnames(df)[2]), ", ", 
                             addbacktick(colnames(df)[3]), ") ~ ", 
                             colnames(df)[4])
        fit = survfit(as.formula(fit_formula), 
                      data = df)
        
        pval = survdiff(as.formula(fit_formula), data = df)$pvalue
        survival_plot  = ggsurvplot(fit, data = df,
                                    conf.int = TRUE, # 增加置信区间
                                    fun = "cumhaz", # 绘制累计风险曲线
                                    ggtheme = theme_bw(),
                                    title=str_glue("pval={pval}")
        )
        g2 <- ggplotGrob(survival_plot$plot + theme(panel.grid = element_blank(),
                                                    axis.text = element_text(size = 15, face = "bold", colour = "black"),
                                                    axis.title.x = element_text(size = 20, face = "bold"),
                                                    axis.title.y = element_text(size = 20, face = "bold")))

        ggsave(filename = paste0("plot/risk_accumulation/", disease_name, "/", colnames(df)[4], "_", colnames(df)[3], ".png"),
               width = 8, height = 6)
        
        temp_df_info = data.frame(pheno = colnames(df)[4], 
                                  disease = colnames(df)[3], 
                                  n_disease_minus = n_disease_minus, 
                                  n_disease_plus = n_disease_plus, 
                                  n_nondisease = n_nondisease, 
                                  pval = pval)
        df_info = bind_rows(df_info, temp_df_info)
      }
      else {
        temp_df_info = data.frame(pheno = colnames(df)[4], 
                                  disease = colnames(df)[3], 
                                  n_disease_minus = n_disease_minus, 
                                  n_disease_plus = n_disease_plus, 
                                  n_nondisease = n_nondisease, 
                                  pval = NA)
        df_info = bind_rows(df_info, temp_df_info)
      }
    }
  }
}
df_info_whole = rbind(df_info_whole, df_info)
df_info["p_bonf"] = p.adjust(df_info$pval, method = "bonferroni")
write.csv(df_info, "data/pca_risk_accumulation.csv")


########################### organ_axis_mix disease accumulation #####################
library("glmnet")
# 加载表型状态，按照表型大小对表型进行分组
pheno = read.csv("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv", row.names = 1)
pheno = na.omit(pheno)
colnames(pheno) = str_remove(colnames(pheno), pattern = "y_res_")

# 加载疾病状态
diseases = read.csv("data/disease_with_pn.csv", row.names = "eid", check.names = FALSE)
diseases_time = read.csv("data/disease_months.csv", row.names = "eid", check.names = FALSE)
colnames(diseases_time) = paste0(colnames(diseases_time), "_month")

############## 使用混合的新表型生成新表型，目标是最大化将健康人和疾病分组 ################
df_info = data.frame()
for (j in seq_along(diseases)) {
  disease_name = names(diseases)[j]
  disease = diseases[j]
  disease = na.omit(disease)
  df = inner_join(pheno %>% rownames_to_column("eid"), 
                  disease %>% rownames_to_column("eid"), 
                  by="eid")
  df = column_to_rownames(df, "eid")
  df = na.omit(df)
  n_disease_minus = nrow(df[(df[[ncol(df)]] == -1), ])
  n_disease_plus = nrow(df[(df[[ncol(df)]] == 1), ])
  n_nondisease = nrow(df[df[[ncol(df)]] == 0, ])
  # 选择after的人进行分析
  df = df[df[[ncol(df)]] >= 0, ]
  
  if (sum(df[ncol(df)] == 1) >= 10) {
    formula = my_formula(colnames(df)[ncol(df)], ".")
    # results = glm(formula = formula, data = df, family = "binomial")
    x = df[-ncol(df)]
    x = as.matrix(x)
    y = df[ncol(df)]
    y = as.matrix(y)
    model = cv.glmnet(x, y, family = "binomial", alpha = 0)
    best_lambda = model$lambda.min
    model = glmnet(x, y, lambda = best_lambda, family = "binomial", alpha = 0)
    mixed_pheno = as.data.frame(as.matrix(pheno) * as.matrix(coef(model))[-1, ])
    mixed_pheno = na.omit(mixed_pheno)
    mixed_pheno_dividation = divide_pheno_by_quantile(mixed_pheno[[1]])
    mixed_pheno = data.frame(row.names = rownames(mixed_pheno), mixed_pheno=mixed_pheno_dividation)
    
    
    disease = diseases[j]
    disease[disease == 1] = 2
    disease[disease == 0] = 1
    disease = na.omit(disease)
    disease_time = diseases_time[j]
    disease_time = na.omit(disease_time)
    df = inner_join(disease_time %>% rownames_to_column("eid"), 
                    disease %>% rownames_to_column("eid"), 
                    by="eid")
    df = inner_join(df, mixed_pheno %>% rownames_to_column("eid"), by="eid")
    
    n_disease_minus = nrow(df[(df[[3]] == -1), ])
    n_disease_plus = nrow(df[(df[[3]] == 2), ])
    n_nondisease = nrow(df[df[[3]] == 1, ])
    
    # 选择after的人进行分析
    df = df[df[[2]] > 0, ]
    if (nrow(df) > 0) {
      addbacktick = function(str) {paste0("`", str, "`")}
      fit_formula = paste0("Surv(", 
                           addbacktick(colnames(df)[2]), ", ", 
                           addbacktick(colnames(df)[3]), ") ~ ", 
                           colnames(df)[4])
      fit = survfit(as.formula(fit_formula), 
                    data = df)
      
      pval = survdiff(as.formula(fit_formula), data = df, rho = 0)$pvalue
      survival_plot  = ggsurvplot(fit, data = df,
                                  conf.int = TRUE, # 增加置信区间
                                  fun = "cumhaz", # 绘制累计风险曲线
                                  ggtheme = theme_bw(),pval = TRUE, pval.method = TRUE, 
                                  title=str_glue("pval={pval}")
      )
      g2 <- ggplotGrob(survival_plot$plot + theme(panel.grid = element_blank(),
                                                  axis.text = element_text(size = 15, face = "bold", colour = "black"),
                                                  axis.title.x = element_text(size = 20, face = "bold"),
                                                  axis.title.y = element_text(size = 20, face = "bold")))

      ggsave(filename = paste0("plot/risk_accumulation/", disease_name, "/", colnames(df)[4], "_", colnames(df)[3], ".png"),
             width = 8, height = 6)
      
      temp_df_info = data.frame(pheno = colnames(df)[4], 
                                disease = colnames(df)[3], 
                                n_disease_minus = n_disease_minus, 
                                n_disease_plus = n_disease_plus, 
                                n_nondisease = n_nondisease, 
                                pval = pval)
      df_info = bind_rows(df_info, temp_df_info)
    }
    else {
      temp_df_info = data.frame(pheno = colnames(df)[4], 
                                disease = colnames(df)[3], 
                                n_disease_minus = n_disease_minus, 
                                n_disease_plus = n_disease_plus, 
                                n_nondisease = n_nondisease, 
                                pval = NA)
      df_info = bind_rows(df_info, temp_df_info)
    }
  }
}
df_info_whole = rbind(df_info_whole, df_info)

df_info["p_bonf"] = p.adjust(df_info$pval, method = "bonferroni")
write.csv(df_info, "data/mix_pheno_risk_accumulation.csv")


df_info_split_by_disease = df_info_whole %>% group_split(disease)
names(df_info_split_by_disease) = map_chr(df_info_split_by_disease, ~unique(.x$disease))
.f = function(df) {
  df = df %>% arrange(pval)  
  df$p_bonf = p.adjust(df$pval)
  return(df)
}
df_info_split_by_disease = map(df_info_split_by_disease, .f = .f)
names(df_info_split_by_disease) = substr(names(df_info_split_by_disease), 1, 30)
write_xlsx(df_info_split_by_disease, "data/figure3B_ukb.xlsx")












