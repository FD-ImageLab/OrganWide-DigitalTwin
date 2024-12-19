library(tidyverse)
library(data.table)
library(tools)
library("PMA")
library("mdthemes")
library("ggrepel")
library(writexl)
library("png")
library(grid)
library(gridExtra)

source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")

heart = read.csv("data/heart.csv", check.names = FALSE)
brain = read.csv("data/brain.csv", check.names = FALSE)
kidney = read.csv("data/kidney.csv", check.names = FALSE)
liver = read.csv("data/liver.csv", check.names = FALSE)
lung = read.csv("data/lung.csv", check.names = FALSE)
pancreas = read.csv("data/pancreas.csv", check.names = FALSE)
spleen = read.csv("data/spleen.csv", check.names = FALSE)

heart = heart %>% column_to_rownames(var="eid")
brain = brain %>% column_to_rownames(var="eid")
kidney = kidney %>% column_to_rownames(var="eid")
liver = liver %>% column_to_rownames(var="eid")
lung = lung %>% column_to_rownames(var="eid")
pancreas = pancreas %>% column_to_rownames(var="eid")
spleen = spleen %>% column_to_rownames(var="eid")

organs = list(Brain=brain, Heart=heart, Lung=lung, Liver=liver, 
              Spleen=spleen, Pancreas=pancreas, Kidney=kidney)

# 回归数据BMI和性别
regression_data = read.csv("data/covariates.csv", check.names = FALSE, row.names = "eid")
regression_data$BMI = regression_data$weight / (regression_data$height*0.01 * regression_data$height*0.01)
regression_data = regression_data[c("isMale", "age", "BMI")]
# 回归后的器官数据
regressed_organs = map(organs, get_res_with_regression_data, regression_data = regression_data)


disease = read.csv("data/disease_with_pn.csv", check.names = FALSE, row.names = "eid")
disease_before_df = disease[apply(disease, 1, function(x) any(x == -1)), ]
health_df = disease[apply(disease, 1, function(x) all(x == 0 | x == 1)), ]
health_id = rownames(health_df)
disease_after_df = disease[apply(disease, 1, function(x) any(x == 1)), ]


basic_info = read.csv("data/basic_info.csv", check.names = FALSE, row.names = "eid")
covariates = read.csv("data/covariates.csv", check.names = FALSE, row.names = "eid")


################## 计算器官pc1结果 ###############
.f = function(organ, name_organ) {
  colnames_pc_organ = paste0(name_organ, "_PC1")
  df = prcomp(na.omit(organ), rank.=1)$x
  df = as.data.frame(df)
  colnames(df) = colnames_pc_organ
  return(df)
}
pc_organs_df = map2(organs, names(organs), .f)
pc_organs_df = map(pc_organs_df, rownames_to_column, var="eid")
pc_organs_df = reduce(pc_organs_df, full_join, by="eid")
write.csv(pc_organs_df, "data/organ_pca.csv", row.names = FALSE)



################### 计算ukb新表型 ###################
output_new_pheno_all = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      new_pheno = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(new_pheno) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("eid")
      check_path(paste0("data/ukb_new_pheno_without_regression/ukb_", colnames(new_pheno), ".csv"))
      write.csv(output_new_pheno, paste0("data/ukb_new_pheno_without_regression/ukb_", colnames(new_pheno), ".csv"), row.names = FALSE)
      if (nrow(output_new_pheno_all) == 0) {
        output_new_pheno_all = output_new_pheno
      }
      else {
        output_new_pheno_all = full_join(output_new_pheno_all, output_new_pheno, by="eid")
      }
    }
  }
}
# check_path(paste0("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv"))
# write.csv(output_new_pheno_all, paste0("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv"), row.names = FALSE)



output_new_pheno_all = data.frame()
for (i in seq_along(regressed_organs)) {
  for (j in seq_along(regressed_organs)) {
    if (i < j) {
      organ1 = regressed_organs[[i]]
      organ2 = regressed_organs[[j]]
      organ1_name = names(regressed_organs)[i]
      organ2_name = names(regressed_organs)[j]
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      new_pheno = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(new_pheno) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("eid")
      check_path(paste0("data/ukb_new_pheno_with_regression/ukb_", colnames(new_pheno), ".csv"))
      write.csv(output_new_pheno, paste0("data/ukb_new_pheno_with_regression/ukb_", colnames(new_pheno), ".csv"), row.names = FALSE)
      if (nrow(output_new_pheno_all) == 0) {
        output_new_pheno_all = output_new_pheno
      }
      else {
        output_new_pheno_all = full_join(output_new_pheno_all, output_new_pheno, by="eid")
      }
    }
  }
}
check_path(paste0("data/ukb_new_pheno_with_regression/whole/ukb_new_pheno.csv"))
write.csv(output_new_pheno_all, paste0("data/ukb_new_pheno_with_regression/whole/ukb_new_pheno.csv"), row.names = FALSE)



################生成xlsx文件 ##################
output_cca_edge_all = list()
############## fig 2健康人的器官间热图，使用y_res作为值 #######################
data = read.csv('data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv', row.names = NULL)
data = data[data$eid %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("eid")

data = data.frame(map(data, ~mean(abs(.x), na.rm=TRUE)))


# 提取行名和列名
row_names = color_map$organ[c(1:7)]
col_names = color_map$organ[c(1:7)]

# 创建新的数据框
new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
colnames(new_data) = unique(col_names)
rownames(new_data) = unique(row_names)

# 填充值
for (i in 1:length(row_names)) {
  for (j in 1:length(col_names)) {
    if (i < j) {
      data_col = paste0("y_res_", row_names[i], "_", col_names[j])
      if (!(data_col %in% colnames(data))) {
        new_data[row_names[i], col_names[j]] = 0
      }
      else {
        new_data[row_names[i], col_names[j]] = data[, data_col]

      }
    }
    if (i > j) {
      data_col = paste0("y_res_", col_names[j], "_", row_names[i])
      if (!(data_col %in% colnames(data))) {
        new_data[row_names[i], col_names[j]] = 0
      }
      else {
        new_data[row_names[i], col_names[j]] = data[, data_col]
      }
    }
    if (i == j) {
      data_col = paste0("y_res_", col_names[j], "_", row_names[i])
      new_data[row_names[i], col_names[j]] = 0

    }
  }
}
output_cca_heatmap_matrix = new_data
new_data = new_data * upper.tri(new_data)
new_data[new_data == 0] = NA
new_data = new_data %>% rownames_to_column("organ1")
output_cca_edge = pivot_longer(new_data, cols = c(2:8))
colnames(output_cca_edge) = c("organ1", "organ2", "average_y_res")
output_cca_edge = na.omit(output_cca_edge)
output_cca_edge$organ1 = rename_based_on_df(output_cca_edge$organ1, color_map, "organ", "short_organ")
output_cca_edge$organ2 = rename_based_on_df(output_cca_edge$organ2, color_map, "organ", "short_organ")

output_cca_edge_all[["health"]] = output_cca_edge

draw_organ_complexnet(output_cca_edge,
                      output_name = "plot/fig2/health_ukb_organ_axis_cca_complexnet.png",
                      limits = c(0.5, 1.4))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/fig2/health_ukb_organ_axis_cca_heatmap.png"))
# draw heatmap plot
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"), bias=1.5)(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/fig2/health_ukb_organ_axis_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## fig 2疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease_before_df)) {
  disease_name = colnames(disease_before_df)[k]
  disease_id = rownames(disease_before_df)[disease_before_df[disease_name] == -1]
  data = read.csv('data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv', row.names = NULL)
  data = data[data$eid %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("eid")
  data = data.frame(map(data, ~mean(abs(.x), na.rm=TRUE)))

  # 提取行名和列名
  row_names = color_map$organ[c(1:7)]
  col_names = color_map$organ[c(1:7)]

  # 创建新的数据框
  new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
  colnames(new_data) = unique(col_names)
  rownames(new_data) = unique(row_names)

  # 填充值
  for (i in 1:length(row_names)) {
    for (j in 1:length(col_names)) {
      if (i < j) {
        data_col = paste0("y_res_", row_names[i], "_", col_names[j])
        if (!(data_col %in% colnames(data))) {
          new_data[row_names[i], col_names[j]] = 0
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]

        }
      }
      if (i > j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        if (!(data_col %in% colnames(data))) {
          new_data[row_names[i], col_names[j]] = 0
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]
        }
      }
      if (i == j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        new_data[row_names[i], col_names[j]] = 0

      }
    }
  }
  output_cca_heatmap_matrix = new_data

  new_data = new_data * upper.tri(new_data)
  new_data[new_data == 0] = NA
  new_data = new_data %>% rownames_to_column("organ1")
  output_cca_edge = pivot_longer(new_data, cols = c(2:8))
  colnames(output_cca_edge) = c("organ1", "organ2", "average_y_res")
  output_cca_edge = na.omit(output_cca_edge)

  output_cca_edge$organ1 = rename_based_on_df(output_cca_edge$organ1, color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(output_cca_edge$organ2, color_map, "organ", "short_organ")

  output_cca_edge_all[[disease_name]] = output_cca_edge

  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/fig2/", disease_name, "_ukb_organ_axis_cca_complexnet.png"),
                        limits = c(0.5, 1.4))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/fig2/", disease_name, "_ukb_organ_axis_cca_heatmap.png"))
  
# draw heatmap plot
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"), bias=1.5)(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/fig2/", disease_name, "_ukb_organ_axis_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
write_xlsx(output_cca_edge_all, "data/figure2_ukb.xlsx")






################### axis outlier plot ###################
# for (i in seq_along(regressed_organs)) {
#   for (j in seq_along(regressed_organs)) {
#     if (i < j) {
#       organ1 = regressed_organs[[i]]
#       organ2 = regressed_organs[[j]]
#       
#       organ1_name = names(regressed_organs)[i]
#       organ2_name = names(regressed_organs)[j]
#       
#       enroll_id = intersect(rownames(na.omit(organ1)), rownames(na.omit(organ2)))      enroll_id = intersect(intersect(rownames(na.omit(organ1)), rownames(na.omit(organ2))), health_id)
#       organ1 = organ1[enroll_id, ]
#       organ2 = organ2[enroll_id, ]
#       
#       results = PMA::CCA(organ1, organ2, typex = "standard", typez = "standard")
#       u = results$u
#       v = results$v
#       x_prime = as.matrix(organ1) %*% u
#       y_prime = as.matrix(organ2) %*% v
#       
#       x_prime = scale(x_prime)
#       y_prime = scale(y_prime)
#       
#       merged_data = data.frame(x_prime = x_prime, y_prime = y_prime, row.names = enroll_id)
#       if ((sd(merged_data[[1]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
#         r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
#         p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
#         n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
# 
#         
#         lm_fit = lm(my_formula(y_variable = colnames(merged_data)[2], 
#                                x_variables = colnames(merged_data)[1]), 
#                     data = merged_data)
#         residuals = residuals(lm_fit)
#         intercept = lm_fit$coefficients[1]
#         slope = lm_fit$coefficients[2]
#         top_outliers = names(sort(abs(residuals), decreasing = TRUE)[1:10])
#         outlier_data = merged_data[top_outliers, ]
#         outlier_text = outlier_data
#         outlier_text$label = paste0("eid ", rownames(outlier_data))
#         
#         My_Theme = theme(
#           panel.background = element_blank(), 
#           title = element_text(size = 7), 
#           text = element_text(size = 6), 
#         )
#         
#         risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                       y = get(colnames(merged_data)[2]))) +
#           geom_point(size=0.5, col="#2171b6") + 
#           geom_abline(intercept = intercept, slope = slope, col="red") + 
#           geom_hline(yintercept = 0, col="black") + 
#           theme(axis.line = element_line(color = "black", size = 0.2)) +
#           mdthemes::md_theme_classic() + 
#           labs(x = organ1_name, 
#                y = organ2_name, 
#                title = paste0("scatter plot"), 
#                subtitle = paste0("p = ", p)) + 
#           My_Theme + 
#           geom_point(data = outlier_data, 
#                      aes(x = get(colnames(merged_data)[1]), 
#                          y = get(colnames(merged_data)[2])), 
#                      size = 2, col = "darkorange") + 
#           ggrepel::geom_text_repel(data = outlier_text, 
#                                    aes(x = get(colnames(merged_data)[1]), 
#                                        y = get(colnames(merged_data)[2]), 
#                                        label = label), 
#                                    color = "darkorange", seed = 5678)
#         check_path(paste0("plot/axis_scatter_plot/ukb_scatter_", organ1_name, "_", organ2_name, ".png"))
#         ggsave(paste0("plot/axis_scatter_plot/ukb_scatter_", organ1_name, "_", organ2_name, ".png"))
#       }
#     }
#   }
# }










# ################### sfig 5 axis age, sex, bmi with regression###################
# for (i in seq_along(regressed_organs)) {
#   for (j in seq_along(regressed_organs)) {
#     if (i < j) {
#       organ1 = regressed_organs[[i]]
#       organ2 = regressed_organs[[j]]
#       organ1_name = names(regressed_organs)[i]
#       organ2_name = names(regressed_organs)[j]
#       
#       results = gene_cca_xy_prime(organ1, organ2, organ1_name, organ2_name)
#       x_prime = results$x_prime
#       y_prime = results$y_prime
#       
#       merged_data = data.frame(x_prime = x_prime, y_prime = y_prime)
#       merged_data = merged_data[rownames(merged_data) %in% health_id, ]
#       merged_data = left_join(x = merged_data %>% rownames_to_column("eid"), 
#                               y = basic_info[c("31-0.0", "21001-0.0")] %>% rownames_to_column("eid"), 
#                               by="eid") %>% column_to_rownames("eid")
#       merged_data = left_join(x = merged_data %>% rownames_to_column("eid"), 
#                               y = covariates[c("age")] %>% rownames_to_column("eid"), 
#                               by="eid") %>% column_to_rownames("eid")
#       merged_data = dplyr::rename(merged_data, "sex"="31-0.0", "bmi"="21001-0.0")
#       merged_data$sex = ifelse(merged_data$sex == 1, "male", "female")
#       merged_data$bmi = ifelse(merged_data$bmi > median(merged_data$bmi, na.rm = TRUE), 
#                                "High", "Low")
#       merged_data$age = ifelse(merged_data$age > median(merged_data$age, na.rm = TRUE), 
#                                "High", "Low")
#       merged_data = na.omit(merged_data)
#       if ((sd(merged_data[[1]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
#         r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
#         p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
#         n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
#         
#         
#         lm_fit = lm(my_formula(y_variable = colnames(merged_data)[2], 
#                                x_variables = colnames(merged_data)[1]), 
#                     data = merged_data)
#         residuals = residuals(lm_fit)
#         intercept = lm_fit$coefficients[1]
#         slope = lm_fit$coefficients[2]
#         top_outliers = names(sort(abs(residuals), decreasing = TRUE)[1:10])
#         outlier_data = merged_data[top_outliers, ]
#         outlier_text = outlier_data
#         outlier_text$label = paste0("eid ", rownames(outlier_data))
#         
#         My_Theme = theme(
#           panel.background = element_blank(), 
#           title = element_text(size = 24, face = "bold"), 
#           text = element_text(size = 20), legend.position = "none"
#         )
#         
#         x_range = range(merged_data[[colnames(merged_data)[1]]])
#         y_range = range(merged_data[[colnames(merged_data)[2]]])
#         
#         # 计算 y 轴的 20% 处, x 轴的0.14处
#         x_position <- 0.2 * diff(x_range) + x_range[1]
#         y_position <- 0.14 * diff(y_range) + y_range[1]
#         if (p < 0.01) {
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=sex)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)
#                   ) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#           
#           
#           
#           
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=bmi)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#           
#           
#           
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=age)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#           
#         }
#         if (p >= 0.05) {
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=sex)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p = ", round(p, digits = 2), 
#                                             "\n", 
#                                             "r = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#           
#           
#           
#           
#           
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=bmi)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p = ", round(p, digits = 2), 
#                                             "\n", 
#                                             "r = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#           
# 
#           risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
#                                                         y = get(colnames(merged_data)[2]), 
#                                                         col=age)) +
#             geom_point(size=0.5, alpha = 0.5) + 
#             annotate("text", label = paste0("p = ", round(p, digits = 2), 
#                                             "\n", 
#                                             "r = ", round(r, digits = 2)), 
#                      x=x_position, y=y_position, size=6) +
#             geom_abline(intercept = intercept, slope = slope, col="red") + 
#             geom_hline(yintercept = 0, col="black") + 
#             theme(axis.line = element_line(color = "black", size = 0.2), 
#                   axis.text = element_text(size=20)) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
#           
#         }
#       }
#     }
#   }
# }
# 
# 


################### sfig 5 axis age, sex, bmi without regression###################
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]

      results = gene_cca_xy_prime(organ1, organ2, organ1_name, organ2_name)
      x_prime = results$x_prime
      y_prime = results$y_prime

      merged_data = data.frame(x_prime = x_prime, y_prime = y_prime)
      merged_data = merged_data[rownames(merged_data) %in% health_id, ]
      merged_data = left_join(x = merged_data %>% rownames_to_column("eid"),
                              y = basic_info[c("31-0.0", "21001-0.0")] %>% rownames_to_column("eid"),
                              by="eid") %>% column_to_rownames("eid")
      merged_data = left_join(x = merged_data %>% rownames_to_column("eid"),
                              y = covariates[c("age")] %>% rownames_to_column("eid"),
                              by="eid") %>% column_to_rownames("eid")
      merged_data = dplyr::rename(merged_data, "sex"="31-0.0", "bmi"="21001-0.0")
      merged_data$sex = ifelse(merged_data$sex == 1, "male", "female")
      merged_data$bmi = ifelse(merged_data$bmi > median(merged_data$bmi, na.rm = TRUE),
                               "High", "Low")
      merged_data$age = ifelse(merged_data$age > median(merged_data$age, na.rm = TRUE),
                               "High", "Low")
      merged_data = na.omit(merged_data)
      merged_data = merged_data[sample(rownames(merged_data)), ]
      if ((sd(merged_data[[1]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
        r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
        p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
        n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1


        lm_fit = lm(my_formula(y_variable = colnames(merged_data)[2],
                               x_variables = colnames(merged_data)[1]),
                    data = merged_data)
        residuals = residuals(lm_fit)
        intercept = lm_fit$coefficients[1]
        slope = lm_fit$coefficients[2]
        top_outliers = names(sort(abs(residuals), decreasing = TRUE)[1:10])
        outlier_data = merged_data[top_outliers, ]
        outlier_text = outlier_data
        outlier_text$label = paste0("eid ", rownames(outlier_data))

        My_Theme = theme(
          panel.background = element_blank(),
          title = element_text(size = 24, face = "bold"),
          text = element_text(size = 20), legend.position = "none", 
          plot.margin = margin(0.1, 0.1, 0, 0.1, "cm")
        )

        x_range = range(merged_data[[colnames(merged_data)[1]]])
        y_range = range(merged_data[[colnames(merged_data)[2]]])

        # 计算 y 轴的 20% 处, x 轴的0.14处
        x_position <- 0.2 * diff(x_range) + x_range[1]
        y_position <- 0.14 * diff(y_range) + y_range[1]
        if (p < 0.01) {
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=sex)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2), "\n p < 0.01"),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width=4, height=3, dpi = 1200)





          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=bmi)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2), "\n p < 0.01"),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)



          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=age)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2), "\n p < 0.01"),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)


        }
        if (p >= 0.01) {
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=sex)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2),
                                            "\n",
                                            "p = ", round(p, digits = 2)),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width=4, height=3, dpi = 1200)






          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=bmi)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2),
                                            "\n",
                                            "p = ", round(p, digits = 2)),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width=4, height=3, dpi=1200)







          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]),
                                                        y = get(colnames(merged_data)[2]),
                                                        col=age)) +
            geom_point(size=0.5, alpha = 0.5) +
            annotate("text", label = paste0("r = ", round(r, digits = 2),
                                            "\n",
                                            "p = ", round(p, digits = 2)),
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") +
            geom_hline(yintercept = 0, col="black") +
            theme(axis.line = element_line(color = "black", size = 0.2),
                  axis.text = element_text(size=20, colour="black")) +
            mdthemes::md_theme_classic() +
            labs(x = organ1_name,
                 y = organ2_name,
            ) +
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/ukb_scatter_age_", organ1_name, "_", organ2_name, ".png"), width=4, height=3)
        }
      }
    }
  }
}

################## sfig5 combine png #################
sort_path = function(path) {
  organ1 = map_chr(path, ~unlist(strsplit(str_remove(basename(.x), pattern = ".png"), "_"))[4])
  organ2 = map_chr(path, ~unlist(strsplit(str_remove(basename(.x), pattern = ".png"), "_"))[5])
  df = data.frame(path = path, organ1 = organ1, organ2 = organ2)
  df$organ1 = factor(df$organ1, levels = color_map$organ)
  df$organ2 = factor(df$organ2, levels = color_map$organ)
  df = df %>% arrange(organ1, organ2)
  path = df$path
  return(path)
}
combine_png_tri_21 = function(image_files, output_path = NULL) {
  # 创建一个空的Grobs列表来保存图片
  image_grobs <- vector("list", length(image_files))

  # 从所有文件中读取图片并添加到Grobs列表
  for (i in seq_along(image_files)) {
    img <- readPNG(image_files[i])
    image_grobs[[i]] <- rasterGrob(img, interpolate=TRUE)
  }

  # 创建向左对齐的倒三角形布局
  row_counts <- 1:6  # 各行图片数从6递减至1
  grid_list <- list()
  max_count <- max(row_counts)
  for (i in seq_along(row_counts)) {
    # 计算这一行的起始和结束索引
    start_index <- sum(row_counts[row_counts > row_counts[i]]) + 1
    end_index <- start_index + row_counts[i] - 1
    # 获取这一行的图片
    row_images <- image_grobs[start_index:end_index]
    
    
    # 使用`nullGrob()`填充行的剩余部分
    row_layout <- c(rep(list(nullGrob()), times = max(row_counts) - row_counts[i]), row_images)
    grid_list[[length(row_counts) - i + 1]] <- do.call(gridExtra::arrangeGrob, c(row_layout, ncol = max(row_counts)))
  }

  # 将所有行整合成最终的向左对齐的倒三角形布局
  final_layout <- do.call(gridExtra::arrangeGrob, c(grid_list, ncol = 1))

  # 绘制结果
  grid.newpage()
  grid.draw(final_layout)

  # 将整个布局导出到一个新的PNG文件
  ggsave(output_path, final_layout, width = 10, height = 8, units = "in")

}

# image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
# image_files = image_paths[grep("bmi", image_paths)]
# image_files = sort_path(image_files)
# combine_png_tri_21(image_files, "plot/axis_scatter_plot/with_regression_bmi_combined.png")
# 
# 
# image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
# image_files = image_paths[grep("sex", image_paths)]
# image_files = sort_path(image_files)
# combine_png_tri_21(image_files, "plot/axis_scatter_plot/with_regression_sex_combined.png")
# 
# 
# image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
# image_files = image_paths[grep("age", image_paths)]
# image_files = sort_path(image_files)
# combine_png_tri_21(image_files, "plot/axis_scatter_plot/with_regression_age_combined.png")




image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("bmi", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/without_regression_bmi_combined.png")


image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("sex", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/without_regression_sex_combined.png")


image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("age", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/without_regression_age_combined.png")





# ################### sfig 5 variance #################
# data = read.csv("data/ukb_new_pheno_with_regression/whole/ukb_new_pheno.csv", row.names = "eid")
# data = data[rownames(data) %in% health_id, ]
# df = data.frame(organ1=NA, organ2=NA, pheno = colnames(data), var=NA)
# df$organ1 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][3])
# df$organ2 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][4])
# 
# df$var = map_dbl(data, sd, na.rm=TRUE)
# 
# df = pivot_longer(df, cols = c("organ1", "organ2"), values_to = "organ")
# 
# df = df %>%
#   group_by(organ) %>%
#   summarise(
#     mean = mean(var, na.rm = TRUE),
#     sd = sd(var, na.rm = TRUE),
#     lower = mean - sd,
#     upper = mean + sd
#   )
# df$organ_color = rename_based_on_df(df$organ, nmapdf = color_map, from = "organ", to = "color")
# 
# p <- ggplot(df, aes(x = organ, y = mean, fill = organ_color)) +
#   geom_bar(stat = 'identity') +  # 条形图
#   scale_fill_identity() +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # 错误条
#   coord_flip() +  # 翻转坐标轴，使条形图水平
#   theme_minimal() +  # 使用简洁的主题
#   theme(legend.position = "none")  # 不显示图例
# ggsave("plot/sfig5/variance_with_regression.png", width = 10, height = 10)
# 
# 
# 
# data = read.csv("data/ukb_new_pheno_without_regression/whole/ukb_new_pheno.csv", row.names = "eid")
# data = data[rownames(data) %in% health_id, ]
# df = data.frame(organ1=NA, organ2=NA, pheno = colnames(data), var=NA)
# df$organ1 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][3])
# df$organ2 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][4])
# 
# df$var = map_dbl(data, sd, na.rm=TRUE)
# 
# df = pivot_longer(df, cols = c("organ1", "organ2"), values_to = "organ")
# 
# df = df %>%
#   group_by(organ) %>%
#   summarise(
#     mean = mean(var, na.rm = TRUE),
#     sd = sd(var, na.rm = TRUE),
#     lower = mean - sd,
#     upper = mean + sd
#   )
# df$organ_color = rename_based_on_df(df$organ, nmapdf = color_map, from = "organ", to = "color")
# 
# p <- ggplot(df, aes(x = organ, y = mean, fill = organ_color)) +
#   geom_bar(stat = 'identity') +  # 条形图
#   scale_fill_identity() +
#   geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # 错误条
#   coord_flip() +  # 翻转坐标轴，使条形图水平
#   theme_minimal() +  # 使用简洁的主题
#   theme(legend.position = "none")  # 不显示图例
# ggsave("plot/sfig5/variance_without_regression.png", width = 10, height = 10)
# 
# 



##################### sfig 5 r value ######################
r_value_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      results = gene_cca_xy_prime(organ1, organ2, organ1_name, organ2_name)
      x_prime = results$x_prime
      y_prime = results$y_prime
      
      merged_data = data.frame(x_prime = x_prime, y_prime = y_prime)
      merged_data = merged_data[rownames(merged_data) %in% health_id, ]
      r = cor.test(merged_data$x_prime, merged_data$y_prime)$estimate
      tmp_df = data.frame(organ1 = organ1_name, organ2 = organ2_name, r=r)
      r_value_df= rbind(r_value_df, tmp_df)
    }
  }
}
r_value_df = pivot_longer(r_value_df, cols = c("organ1", "organ2"), values_to = "organ")
r_value_df = r_value_df %>%
  group_by(organ) %>%
  summarise(r_mean = mean(r), 
            r_sd = sd(r), 
            lower = r_mean - r_sd,
            upper = r_mean + r_sd)


r_value_df$organ_color = rename_based_on_df(r_value_df$organ, nmapdf = color_map, 
                                            from = "organ", to = "color")
r_value_df = r_value_df[order(match(r_value_df$organ, color_map$organ)), ]
r_value_df$organ = factor(r_value_df$organ)
check_path("data/sfig5/ukb_r.csv")
write.csv(r_value_df, "data/sfig5/ukb_r.csv")

p <- ggplot(r_value_df, aes(x = organ, y = r_mean, fill = organ_color)) +
  geom_bar(stat = 'identity', width=0.5) +  # 条形图
  scale_fill_identity() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +  # 错误条
  coord_flip() +  # 翻转坐标轴，使条形图水平
  theme_minimal() +  # 使用简洁的主题
  theme(legend.position = "none", 
        axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30), 
        axis.title = element_blank())  # 不显示图例
ggsave("plot/sfig5/r_without_regression.png", width = 6, height = 8)








################## figure 4A UKB ##########################
feature_dict = read.csv("Sandbox/ukb_category/feature.csv")

features = read.csv("data/feature.csv", row.names = "eid", check.names = FALSE)
# features = features[colnames(features) %in% feature_dict$UDI]
feature_sd = map_dbl(features, ~sd(.x, na.rm = TRUE))
features = features[names(na.omit(feature_sd[feature_sd != 0]))]


# 小分子表型
moleculars = read.csv("data/molecular.csv", row.names = "eid", check.names = FALSE)


df_all = data.frame()
for (new_pheno_file in list.files("data/ukb_new_pheno_without_regression/", full.names = TRUE, pattern = ".csv")) {
  print(paste0("calculating ", basename(new_pheno_file)))
  male_df = read.csv(new_pheno_file, row.names = "eid")
  df = data.frame()
  
  n = 0
  for (i in seq_along(features)) {
    for (j in seq_along(male_df)) {
      feature = features[i]
      feature_name = colnames(features)[i]
      if (class(features[[i]]) == "numeric" | class(features[[i]]) == "integer") {
        axis = male_df[j]
        axis_name = str_remove(colnames(male_df)[j], "y_res_")
        
        merged_data = inner_join(feature %>% rownames_to_column("eid"), 
                                 axis %>% rownames_to_column("eid"), 
                                 by="eid") %>% 
          column_to_rownames("eid") %>% 
          na.omit()
        if ((sd(merged_data[[1]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
          r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
          p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
          n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
          temp_df = data.frame(pheno = colnames(male_df)[j], axis = axis_name, 
                               feature = feature_name, r = r, p = p, n = n)
          df = bind_rows(df, temp_df)
          n = n + 1
        } 
      }
    }
  }
  df["organ1"] = map_chr(str_split(df$axis, "_"), ~.x[1])
  df["organ2"] = map_chr(str_split(df$axis, "_"), ~.x[2])
  df["p_bonf"] = p.adjust(df$p, method = "bonferroni", n = nrow(df))
  write.csv(df, paste0("data/5/ukb_", colnames(male_df), "_data.csv"), row.names = FALSE)
  df_all = bind_rows(df_all, df)
}
df_all["p_bonf"] = p.adjust(df_all$p, method = "bonferroni", n = nrow(df_all))
write.csv(df_all, paste0("data/5/0_ukb_all_data.csv"), row.names = FALSE)






data = read.csv("data/5/0_ukb_all_data.csv")
feature_dict = read.csv("Sandbox/ukb_category/feature.csv")
data = data[data$feature %in% feature_dict$UDI, ]
data = data[data$p_bonf < 0.05, ]
new_data = as.data.frame(matrix(NA, nrow=length(unique(data$feature)), ncol=7, dimnames = list(c(unique(data$feature)), color_map$organ[1:7])))

start = Sys.time()
for (i in c(1:nrow(new_data))) {
  for (j in c(1:ncol(new_data))) {
    new_data[i, j] = mean(abs(data$r[((data$feature == rownames(new_data)[i]) & grepl(colnames(new_data)[j], data$axis, ignore.case = TRUE))]))
  }
}
end = Sys.time()

copy_new_data = new_data


new_data = copy_new_data
new_data$Sum = rowMeans(new_data)
df_sorted = new_data[order(new_data$Sum, decreasing = TRUE), ]
df_sorted$Sum = NULL
temp = df_sorted
temp$feature = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'names')
temp  = temp %>% dplyr::select(feature, everything())
temp$category = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'category_name')
write.csv(temp, file = "data/feature_net_sort_by_wholenet.csv", row.names = TRUE)

a = read.csv("Sandbox/ukb_category/feature.csv")
temp = read.csv("data/feature_net_sort_by_wholenet.csv", row.names = 1)
temp$name = rename_based_on_df(rownames(temp), nmapdf = a, from = "UDI", to = "names")
temp$category = rename_based_on_df(rownames(temp), nmapdf = a, from = "UDI", to = "category_name")
write.csv(temp, 'a.csv')


df_sorted = read.csv("data/feature_net_sort_by_wholenet.csv", row.names = 1)
df_sorted = df_sorted[!duplicated(df_sorted$feature), ]
df_sorted = df_sorted[1:30, ]
# for (i in c(1:nrow(df_sorted))) {
#   df = data[data$feature == rownames(df_sorted)[i], ]
#   df$feature = rename_based_on_df(df$feature, nmapdf = feature_dict, from = 'UDI', to = 'names')
#   feature_name = unique(df$feature)
#   df = df %>% dplyr::select(organ1, organ2, r)
#   df$r = abs(df$r)
#   df$organ1 = str_to_title(df$organ1)
#   df$organ2 = str_to_title(df$organ2)
#   draw_organ_complexnet(df, paste0("plot/complexnet/feature/whole_net/", gsub("[\"'/:*?,<>|%]", "_", feature_name), ".png"), 
#                         limit=c(0, 0.35))
# }
df_sorted_with_UDI = df_sorted
rownames(df_sorted) = NULL
df_sorted = df_sorted %>% column_to_rownames("feature")
df_sorted$category = NULL
hM <- format(round(as.matrix(df_sorted), 2))#对数据保留2位小数


pheatmap::pheatmap(as.matrix(df_sorted), cluster_rows = FALSE, cluster_cols = FALSE, 
                   number_color = "black",
                   fontsize_number = 14,
                   fontsize = 20,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50), 
                   cellwidth=30,
                   cellheight=30,
                   filename = paste0("plot/heatmap/whole_net.png"),
                   display_numbers = hM, 
                   width = 16,
                   height = 16,
                   units = "in", border_color = "black", 
                   res = 600)

# ################### Fig 4B UKB ###################
# data = read.csv("data/mediation/top10_whole.csv", row.names = 1)
# data = data[data$x %in% unique(data$x)[1:10], ]
# feature_dict = read.csv("Sandbox/ukb_category/feature.csv")
# molecular_dict = read.csv("Sandbox/ukb_category/molecular.csv")
# pheno_dict = read.csv("scripts/pheno_dict.csv")
# data$y = rename_based_on_df(data$y, pheno_dict, from = "names", to = "description")
# data$molecular_category = rename_based_on_df(data$mediator, molecular_dict, from = "UDI", "Category_name")
# 
# data_list = data %>%
#   filter(prop.mediated.p <= 0.05) %>%
#   group_split(x) 
# 
# .f = function(df) {
#   df = df %>%
#     arrange(y, mediator) %>%
#     dplyr::select(x, y, mediator, prop.mediated, prop.mediated.p, molecular_category, N)
#   df$BP = c(1:nrow(df))
#   df$mht_col = rename_based_on_df(df$molecular_category, nmapdf = ukb_molecular_map, 
#                                   from = "molecular", to = "color")
#   return(df)
# }
# df = map(data_list, .f)
# names(df) = map_chr(data_list, ~unique(.x$x))
# writexl::write_xlsx(df, "data/mediation/mht.xlsx")
# 
# 
# # plot complexnet
# for (i in seq_along(df)) {
#   output_png_name = paste0("plot/fig4B/complexnet/", names(df)[i], ".png")
#   output_csv_name = paste0("data/fig4B/complexnet/", names(df)[i], ".csv")
#   
#   png_df = df[[i]]
#   
#   png_df$organ1 = map_chr(png_df$y, ~str_split(.x, "_")[[1]][1])
#   png_df$organ2 = map_chr(png_df$y, ~str_split(.x, "_")[[1]][2])
#   
#   results = png_df %>%
#     group_by(organ1, organ2) %>%
#     dplyr::summarise(mean_prop=mean(prop.mediated, na.rm=TRUE)) %>%
#     ungroup()
#   
#   check_path(output_csv_name)
#   check_path(output_png_name)
#   write.csv(results, file = output_csv_name)
#   draw_organ_complexnet(results, output_name = output_png_name)
# }
# ################## micro phenotype scatter plot ##########################
# 
# micro_dict = read.csv("Sandbox/ukb_category/molecular.csv")
# 
# micros = read.csv("data/molecular.csv", row.names = "eid", check.names = FALSE)
# micro_sd = map_dbl(micros, ~sd(.x, na.rm = TRUE))
# micros = micros[names(na.omit(micro_sd[micro_sd != 0]))]
# 
# 
# df_all = data.frame()
# for (new_pheno_file in list.files("data/ukb_new_pheno_without_regression", full.names = TRUE, pattern = ".csv")) {
#   axises = read.csv(new_pheno_file) %>% column_to_rownames(var = "eid")
#   df = data.frame()
#   n = 0
#   for (i in seq_along(micros)) {
#     for (j in seq_along(axises)) {
#       micro = micros[i]
#       micro_name = colnames(micros)[i]
#       if (class(micros[[i]]) == "numeric" | class(micros[[i]]) == "integer") {
#         axis = axises[j]
#         axis_name = str_remove(colnames(axises)[j], "y_res_")
#         
#         merged_data = inner_join(micro %>% rownames_to_column("eid"), axis %>% rownames_to_column("eid")) %>% na.omit()
#         rownames(merged_data) = NULL
#         merged_data = merged_data %>% column_to_rownames("eid")
# 
#         if ((sd(merged_data[[1]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
#           r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
#           p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
#           n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
#           temp_df = data.frame(pheno = colnames(axises)[j], axis = axis_name, 
#                                micro_udi = micro_name, 
#                                micro_name = gsub("[\"'/:*?,<>|%]", "_", micro_dict$Label[micro_dict$UDI %in% colnames(micro)]), 
#                                r = r, p = p, n = n)
#           df = bind_rows(df, temp_df)
#     
#           n = n + 1
#         } 
#       }
#     }
#   }
#   df["organ1"] = map_chr(str_split(df$axis, "_"), ~.x[1])
#   df["organ2"] = map_chr(str_split(df$axis, "_"), ~.x[2])
#   df["p_bonf"] = p.adjust(df$p, method = "bonferroni", n = nrow(df))
#   write.csv(df, paste0("data/5/ukb_micro_", colnames(axises), "_data.csv"), row.names = FALSE)
#   df_all = bind_rows(df_all, df)
# }
# df_all["p_bonf"] = p.adjust(df_all$p, method = "bonferroni", n = nrow(df_all))
# write.csv(df_all, paste0("data/5/0_ukb_micro_all_data.csv"), row.names = FALSE)
# 
# 
# 
# ukb_field = read.delim("data/ukb_field_added.txt", header = FALSE, check.names = TRUE)
# unique(ukb_field$V3)
# 
# 
# 
# 
# 
# 
# 

##################### 修改后的fig2 #################
install.packages("UpSetR")
library("UpSetR")

setbarcolor <- c("#2e409a","#942d8d","#d75427","#006b7b","#4da0a0","#9b3a74")#设置色条颜色
upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
      order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE))


data = read.csv("data/disease_odd_before.csv")
dict = read.csv("scripts/pheno_dict.csv")
data = data[data$p <= 0.05, ]
data$pheno = rename_based_on_df(data$pheno, dict, from = "names", to = "description")
data = data %>% dplyr::select(pheno, disease)
data = data %>%
  separate_rows(pheno, sep = "_")
# 需求好像不一样


####################### figure 2 柱状图 disease before ############
data = read.csv("data/disease_odd_before.csv")
dict = read.csv("scripts/pheno_dict.csv")
data$disease_organ = rename_based_on_df(data$disease, nmapdf = dict, from = "description", to = "organ")
data$disease_color = rename_based_on_df(data$disease_organ, nmapdf = color_map, from = "organ", to = "color")
# 按疾病分类，计算odds均值及p值显著的pheno个数（定义p < 0.01 为显著）
result = data %>%
  group_by(disease, disease_color, disease_organ) %>%
  summarise(
    mean_odds = mean(odds),                      # odds的均值
    significant_pheno_count = sum(p < 0.01)      # p值显著的pheno个数
  )
write.csv(result, "data/figure2_ukb_bar.csv", row.names = FALSE)

# result$disease = rename_based_on_df(result$disease, dict, from = "description", to = "abbreviation")
# result = data[data$disease %in% c("T2D", "Chronic Renal Failure", "Asthma", "AVLBB", "Myocardial Infarction", 
#                                 "FIB/CIRR", "Stroke", "Inflammatory Prostate Diseases", 
#                                 "Inflammatory Liver Diseases","Heart Failure", 
#                                 "Cardiomyopathy", "Viral Pneumonia"), ]
# 绘制柱状图
ggplot(result, aes(x = reorder(disease, significant_pheno_count),
                   y = significant_pheno_count,
                   label = significant_pheno_count,
                   fill = disease_color)) +
  geom_col() +
  geom_text(aes(label = disease),  # 添加疾病名称的标签
            position = position_stack(vjust = 0),  # 调整标签的位置
            size = 3, hjust=0) +
  scale_fill_identity() +
  labs(x = "Disease", y = "Significant axis count") +
  # labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20, colour = "black"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 24, face = "bold"),
        axis.title.y = element_text(size = 24, face = "bold"),
        panel.grid = element_blank()) +
  coord_flip()
ggsave("plot/figure2_ukb_bar.png", width = 6, height = 6)


####################### figure 2 柱状图2 disease before ############
data = read.csv("data/disease_odd_before.csv")
dict = read.csv("scripts/pheno_dict.csv")


significant_organs = data %>%
  separate(col = pheno, into = c("y", "res", "Organ1", "Organ2"), sep = "_") %>%  # 分割pheno列
  filter(p < 0.01) %>%  # 筛选显著的疾病（p值小于0.01）
  select(Organ1, Organ2, p, disease) %>%  # 选择感兴趣的列
  gather(key = "OrganKey", value = "Organ", Organ1, Organ2) %>%  # 把器官列转换为长格式
  group_by(Organ) %>%  # 按器官分组
  summarise(Count = n())  # 计算每个器官关联的显著疾病的数量
significant_organs$organ_color = rename_based_on_df(significant_organs$Organ,
                                                    nmapdf = color_map, from = "organ", to = "color")
significant_organs$Organ = factor(significant_organs$Organ, levels = rev(color_map$organ))
write.csv(significant_organs, "data/figure2_ukb_bar2.csv", row.names = FALSE)
# 绘制柱状图
ggplot(significant_organs, aes(x = Organ, y = Count,
                               label = Count, fill=organ_color)) +
  geom_col() +
  coord_flip() + 
  scale_x_discrete(position = "top") + 
  scale_fill_identity() + 
  scale_y_reverse() + 
  geom_text(size = 8, hjust=-1) +  # 在柱子上方标示显著的pheno个数
  labs(x = "Organ", y = "Significant Count") +
  # geom_segment(aes(x = -Inf, xend = Inf, y = -Inf, yend = -Inf), 
  #              color = "black", size = 1) + 
  geom_segment(aes(x = -Inf, xend = -Inf, y = -Inf, yend = Inf), 
               color = "black", size = 1) + 
  # ylim(0, max(significant_organs$Count) * 1.1) +
  theme(axis.text.x = element_text(size = 30, colour = "black"),
        axis.title.x = element_text(size = 36, face = "bold"),  # 设置 X 轴标题文字大小
        axis.ticks.length.y = unit(0, "cm"), 
        
        axis.text.y = element_text(size = 30, colour = "black", hjust=0),
        axis.title.y = element_blank(), 
        plot.margin = unit(c(1.1, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white")) 
ggsave("plot/figure2_ukb_bar2.png", width = 6, height = 6)


# ####################### figure 3 气泡图 disease after ############
# 
# data = read.csv("data/disease_odd_after.csv")
# dict = read.csv("scripts/pheno_dict.csv")
# data = data %>%
#   separate(col = pheno, into = c("y", "res", "Organ1", "Organ2"), sep = "_")  # 分割pheno列
# data$sig = case_when(
#   data$p > 0.05 ~ " ",
#   (data$p > 0.01 & data$p <= 0.05) ~ "*",
#   (data$p > 0.001 & data$p <= 0.01) ~ "**",
#   (data$p <= 0.001) ~ "***"
# )
# data$sig = as.factor(data$sig)
# data$logp = -log10(data$p)
# 
# 
# data = data %>%
#   dplyr::select(Organ1, Organ2, disease, odds, logp, sig)
# data$disease = factor(data$disease)
# data = data %>%
#   gather(key = "OrganKey", value = "Organ", Organ1, Organ2)
# data$y_shift = ifelse(data$OrganKey == "Organ1", -0.05, 0.05)
# data$y_shift = data$y_shift * as.numeric(data$sig)
# data$organ_color = rename_based_on_df(data$Organ, color_map, from = "organ", to = "color")
# write.csv(data, "data/fig3A_ukb.csv", row.names = FALSE)
# 
# ggplot(data, aes(x = odds, y = disease, color=organ_color, size=sig)) +
#   geom_point(alpha = 1, position = position_nudge(y=data$y_shift)) + # 设置透明度
#   scale_size_discrete(range = c(1, 4)) + # 调整气泡的大小范围
#   scale_color_identity() +
#   labs(x = "OR",
#        y = "Diseases") +
#   theme_minimal() +
#   geom_vline(xintercept = 1)
# ggsave("plot/figure3A_ukb.png", width = 12, height = 6)
# 





#################毕设####################

data = read.csv("data/5/0_ukb_all_data.csv")
feature_dict = read.csv("Sandbox/ukb_category/feature.csv")
data = data[data$feature %in% feature_dict$UDI, ]
# data = data[data$p_bonf < 0.05, ]
data = data[grepl("Heart", data$axis), ]


new_data = as.data.frame(matrix(NA, nrow=length(unique(data$feature)), ncol=6, dimnames = list(c(unique(data$feature)), unique(data$axis))))

for (i in c(1:nrow(new_data))) {
  for (j in c(1:ncol(new_data))) {
    if (!identical(data$r[((data$feature == rownames(new_data)[i]) & grepl(colnames(new_data)[j], data$axis, ignore.case = TRUE))], numeric(0))) {
      new_data[i, j] = data$r[((data$feature == rownames(new_data)[i]) & grepl(colnames(new_data)[j], data$axis, ignore.case = TRUE))]
      
    }
  }
}



new_data_p = as.data.frame(matrix(NA, nrow=length(unique(data$feature)), ncol=6, dimnames = list(c(unique(data$feature)), unique(data$axis))))

for (i in c(1:nrow(new_data_p))) {
  for (j in c(1:ncol(new_data_p))) {
    if (!identical(data$p[((data$feature == rownames(new_data_p)[i]) & grepl(colnames(new_data_p)[j], data$axis, ignore.case = TRUE))], numeric(0))) {
      new_data_p[i, j] = data$p[((data$feature == rownames(new_data_p)[i]) & grepl(colnames(new_data_p)[j], data$axis, ignore.case = TRUE))]
      
    }
  }
}
p_bonf = 0.05 / dim(new_data_p)[1]*dim(new_data_p)[2]
new_data[new_data_p > p_bonf] = 0



new_data$sum = rowSums(abs(new_data))
new_data = new_data %>% slice_max(sum, n = 30)
new_data$feature = rename_based_on_df(rownames(new_data), nmapdf = feature_dict, from = 'UDI', to = 'names')









new_data = copy_new_data
new_data$Sum = rowMeans(new_data)
df_sorted = new_data[order(new_data$Sum, decreasing = TRUE), ]
df_sorted$Sum = NULL
temp = df_sorted
temp$feature = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'names')
temp  = temp %>% dplyr::select(feature, everything())
temp$category = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'category_name')
write.csv(temp, file = "data/feature_net_sort_by_wholenet.csv", row.names = TRUE)

a = read.csv("Sandbox/ukb_category/feature.csv")
temp = read.csv("data/feature_net_sort_by_wholenet.csv", row.names = 1)
temp$name = rename_based_on_df(rownames(temp), nmapdf = a, from = "UDI", to = "names")
temp$category = rename_based_on_df(rownames(temp), nmapdf = a, from = "UDI", to = "category_name")
write.csv(temp, 'a.csv')


df_sorted = read.csv("data/feature_net_sort_by_wholenet.csv", row.names = 1)
df_sorted = df_sorted[!duplicated(df_sorted$feature), ]
df_sorted = df_sorted[1:30, ]
# for (i in c(1:nrow(df_sorted))) {
#   df = data[data$feature == rownames(df_sorted)[i], ]
#   df$feature = rename_based_on_df(df$feature, nmapdf = feature_dict, from = 'UDI', to = 'names')
#   feature_name = unique(df$feature)
#   df = df %>% dplyr::select(organ1, organ2, r)
#   df$r = abs(df$r)
#   df$organ1 = str_to_title(df$organ1)
#   df$organ2 = str_to_title(df$organ2)
#   draw_organ_complexnet(df, paste0("plot/complexnet/feature/whole_net/", gsub("[\"'/:*?,<>|%]", "_", feature_name), ".png"), 
#                         limit=c(0, 0.35))
# }
df_sorted_with_UDI = df_sorted
rownames(df_sorted) = NULL
df_sorted = df_sorted %>% column_to_rownames("feature")
df_sorted$category = NULL
hM <- format(round(as.matrix(df_sorted), 2))#对数据保留2位小数


pheatmap::pheatmap(as.matrix(df_sorted), cluster_rows = FALSE, cluster_cols = FALSE, 
                   number_color = "black",
                   fontsize_number = 14,
                   fontsize = 20,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50), 
                   cellwidth=30,
                   cellheight=30,
                   filename = paste0("plot/heatmap/whole_net.png"),
                   display_numbers = hM, 
                   width = 16,
                   height = 16,
                   units = "in", border_color = "black", 
                   res = 600)

