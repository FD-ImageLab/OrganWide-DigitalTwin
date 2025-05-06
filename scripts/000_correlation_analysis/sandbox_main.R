source("scripts/utils/Func.R")
source("scripts/utils/read_sandbox_data.R", encoding = "UTF-8")
source("scripts/utils/pfunc.R", encoding = "UTF-8")

library(tidyverse)
library(data.table)
library(tools)
library("PMA")
library(writexl)

organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]


full_organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
full_organs = full_organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")


male_organs = subset(full_organs, names(full_organs) != "Uterus")
female_organs = subset(full_organs, names(full_organs) != "Prostate")
organs = subset(full_organs, ((names(full_organs) != "Prostate") & names(full_organs) != "Uterus"))

# regress data
regression_data = read_regression_data()
regressed_organs = map(organs, get_res_with_regression_data, regression_data = regression_data)


disease = read.csv("Sandbox/disease_sandbox.csv", check.names = FALSE, row.names = "buId")
health_df = disease[apply(disease, 1, function(x) all(x == 0)), ]
health_id = rownames(health_df)




################## 计算器官pc1结果 ###############
.f = function(organ, name_organ) {
  colnames_pc_organ = paste0(name_organ, "_PC1")
  df = prcomp(na.omit(organ), rank.=1)$x
  df = as.data.frame(df)
  colnames(df) = colnames_pc_organ
  return(df)
}
pc_organs_df = map2(organs, names(organs), .f)
pc_organs_df = map(pc_organs_df, rownames_to_column, var="buId")
pc_organs_df = reduce(pc_organs_df, full_join, by="buId")
write.csv(pc_organs_df, "data/sandbox_organ_pca.csv", row.names = FALSE)



############ male ind ###########
male_df = data.frame()
for (i in seq_along(male_organs)) {
  for (j in seq_along(male_organs)) {
    if (i < j) {
      x = male_organs[[i]]
      y = male_organs[[j]]

      x_name = names(male_organs)[i]
      y_name = names(male_organs)[j]
      results = gene_new_pheno(x, y, x_name, y_name)
      results = as.data.frame(results) %>% rownames_to_column("buId")
      if (nrow(male_df) == 0) {
        male_df = results
      }
      else {
        male_df = full_join(male_df, results,by="buId")
      }
    }
  }
}
male_df = male_df %>% column_to_rownames("buId")
write.csv(x = male_df, file = "data/male_indnet_data.csv")








############ female ind ###########
female_df = data.frame()
for (i in seq_along(female_organs)) {
  for (j in seq_along(female_organs)) {
    if (i < j) {
      x = female_organs[[i]]
      y = female_organs[[j]]
      x_name = names(female_organs)[i]
      y_name = names(female_organs)[j]
      results = gene_new_pheno(x, y, x_name, y_name)
      results = as.data.frame(results) %>% rownames_to_column("buId")
      if (nrow(female_df) == 0) {
        female_df = results
      }
      else {
        female_df = full_join(female_df, results,by="buId")
      }
    }
  }
}
female_df = female_df %>% column_to_rownames("buId")
write.csv(x = female_df, file = "data/female_indnet_data.csv")




############ without regression #############
without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]

      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      check_path(paste0("data/sandbox_new_pheno_without_regression/sandbox_", colnames(results), ".csv"))
      write.csv(output_new_pheno, paste0("data/sandbox_new_pheno_without_regression/sandbox_", colnames(results), ".csv"), row.names = FALSE, quote=FALSE)
      if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}
# without_sex_df = without_sex_df %>% column_to_rownames("buId")
check_path("data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv")
write.csv(x = without_sex_df, file = "data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv", row.names = FALSE, quote=FALSE)


############ with regression #############
without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]

      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results)
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      check_path(paste0("data/sandbox_new_pheno_with_regression/sandbox_", colnames(results), ".csv"))
      write.csv(output_new_pheno, paste0("data/sandbox_new_pheno_with_regression/sandbox_", colnames(results), ".csv"), row.names = FALSE, quote=FALSE)
      if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}
# without_sex_df = without_sex_df %>% column_to_rownames("buId")
check_path("data/sandbox_new_pheno_with_regression/whole/sandbox_new_pheno.csv")
write.csv(x = without_sex_df, file = "data/sandbox_new_pheno_with_regression/whole/sandbox_new_pheno.csv", row.names = FALSE, quote=FALSE)











############### 健康人的器官间热图，使用y_res作为值 ###########
output_cca_edge_all = list()

data = read.csv("data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv")
data = data[data$buId %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("buId")

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
                      output_name = "plot/fig2/health_sandbox_organ_axis_cca_complexnet.png",
                      limits = c(0.5, 1.5))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/fig2/health_sandbox_organ_axis_cca_heatmap.png"))
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/fig2/health_sandbox_organ_axis_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## 疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease)) {
  disease_name = colnames(disease)[k]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  
  data = read.csv("data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv")
  data = data[data$buId %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
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
                        output_name = paste0("plot/fig2/", disease_name, "_sandbox_organ_axis_cca_complexnet.png"),
                        limits = c(0.5, 1.5))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/fig2/", disease_name, "_sandbox_organ_axis_cca_heatmap.png"))
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/fig2/", disease_name, "_sandbox_organ_axis_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
write_xlsx(output_cca_edge_all, "data/figure2/sandbox.xlsx")













# ################### axis age with regression###################
# for (i in seq_along(regressed_organs)) {
#   for (j in seq_along(regressed_organs)) {
#     if (i < j) {
#       organ1 = regressed_organs[[i]]
#       organ2 = regressed_organs[[j]]
#       organ1_name = names(regressed_organs)[i]
#       organ2_name = names(regressed_organs)[j]
#       
#       results = gene_cca_xy_prime(organ1, organ2, organ1_name, organ2_name)
#       
#       x_prime = results$x_prime
#       y_prime = results$y_prime
#       
#       merged_data = data.frame(x_prime = x_prime, y_prime = y_prime)
#       merged_data = left_join(x = merged_data %>% rownames_to_column("buId"), 
#                               y = regression_data %>% rownames_to_column("buId"), 
#                               by="buId") %>% column_to_rownames("buId")
#       merged_data = merged_data[rownames(merged_data) %in% health_id, ]
#       merged_data = dplyr::rename(merged_data, "sex"="性别", "bmi"="BMI", "age"="年龄")
#       
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
#             ) +
#             mdthemes::md_theme_classic() + 
#             labs(x = organ1_name, 
#                  y = organ2_name, 
#             ) + 
#             My_Theme
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
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
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
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
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
#           
#           
#         }
#         if (p >= 0.01) {
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
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
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
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
#           
#           
#           
#           
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
#           check_path(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"))
#           ggsave(paste0("plot/axis_scatter_plot/with_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
#         }
#       }
#     }
#   }
# }



################### axis age, sex, bmi without regression###################
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
      
      merged_data = left_join(x = merged_data %>% rownames_to_column("buId"), 
                              y = regression_data %>% rownames_to_column("buId"), 
                              by="buId") %>% column_to_rownames("buId")
      merged_data = merged_data[rownames(merged_data) %in% health_id, ]
      merged_data = dplyr::rename(merged_data, "sex"="性别", "bmi"="BMI", "age"="年龄")
      
      merged_data$sex = ifelse(merged_data$sex == 1, "male", "female")
      merged_data$bmi = ifelse(merged_data$bmi > median(merged_data$bmi, na.rm = TRUE), 
                               "High", "Low")
      merged_data$age = ifelse(merged_data$age > median(merged_data$age, na.rm = TRUE), 
                               "High", "Low")
      merged_data = na.omit(merged_data)
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
          text = element_text(size = 20), legend.position = "none"
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
            annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
          
          
          
          
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
                                                        y = get(colnames(merged_data)[2]), 
                                                        col=bmi)) +
            geom_point(size=0.5, alpha = 0.5) + 
            annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
          
          
          
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
                                                        y = get(colnames(merged_data)[2]), 
                                                        col=age)) +
            geom_point(size=0.5, alpha = 0.5) + 
            annotate("text", label = paste0("p < 0.01", "\nr = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
          
        }
        if (p >= 0.01) {
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
                                                        y = get(colnames(merged_data)[2]), 
                                                        col=sex)) +
            geom_point(size=0.5, alpha = 0.5) + 
            annotate("text", label = paste0("p = ", round(p, digits = 2), 
                                            "\n", 
                                            "r = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_sex_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
          
          
          
          
          
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
                                                        y = get(colnames(merged_data)[2]), 
                                                        col=bmi)) +
            geom_point(size=0.5, alpha = 0.5) + 
            annotate("text", label = paste0("p = ", round(p, digits = 2), 
                                            "\n", 
                                            "r = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_bmi_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
          
          
          
          
          
          
          risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), 
                                                        y = get(colnames(merged_data)[2]), 
                                                        col=age)) +
            geom_point(size=0.5, alpha = 0.5) + 
            annotate("text", label = paste0("p = ", round(p, digits = 2), 
                                            "\n", 
                                            "r = ", round(r, digits = 2)), 
                     x=x_position, y=y_position, size=6) +
            geom_abline(intercept = intercept, slope = slope, col="red") + 
            geom_hline(yintercept = 0, col="black") + 
            theme(axis.line = element_line(color = "black", size = 0.2), 
                  axis.text = element_text(size=20)) +
            mdthemes::md_theme_classic() + 
            labs(x = organ1_name, 
                 y = organ2_name, 
            ) + 
            My_Theme
          check_path(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"))
          ggsave(paste0("plot/axis_scatter_plot/without_regression/sandbox_scatter_age_", organ1_name, "_", organ2_name, ".png"), width = 4, height = 3)
          
        }
        
        
      }
    }
  }
}


################### combine png #################
library("png")
library(grid)
library(gridExtra)

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
  ggsave(output_path, final_layout, width = 10, height = 10, units = "in")
  
}

image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
image_files = image_paths[grep("bmi", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_with_regression_bmi_combined.png")


image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
image_files = image_paths[grep("sex", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_with_regression_sex_combined.png")


image_paths = list.files("plot/axis_scatter_plot/with_regression/", full.names = TRUE)
image_files = image_paths[grep("age", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_with_regression_age_combined.png")




image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("bmi", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_without_regression_bmi_combined.png")


image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("sex", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_without_regression_sex_combined.png")


image_paths = list.files("plot/axis_scatter_plot/without_regression/", full.names = TRUE)
image_files = image_paths[grep("age", image_paths)]
image_files = sort_path(image_files)
combine_png_tri_21(image_files, "plot/axis_scatter_plot/sandbox_without_regression_age_combined.png")





################### variance #################
data = read.csv("data/sandbox_new_pheno_with_regression/whole/sandbox_new_pheno.csv", row.names = "buId")
data = data[rownames(data) %in% health_id, ]
df = data.frame(organ1=NA, organ2=NA, pheno = colnames(data), var=NA)
df$organ1 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][3])
df$organ2 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][4])

df$var = map_dbl(data, sd, na.rm=TRUE)

df = pivot_longer(df, cols = c("organ1", "organ2"), values_to = "organ")

df = df %>%
  group_by(organ) %>%
  summarise(
    mean = mean(var, na.rm = TRUE),
    sd = sd(var, na.rm = TRUE),
    lower = mean - sd,
    upper = mean + sd
  )
df$organ_color = rename_based_on_df(df$organ, nmapdf = color_map, from = "organ", to = "color")
check_path("data/sfig4/sandbox_with_regression_variance.csv")
write.csv(df, "data/sfig4/sandbox_with_regression_variance.csv")

p <- ggplot(df, aes(x = organ, y = mean, fill = organ_color)) +
  geom_bar(stat = 'identity') +  # 条形图
  scale_fill_identity() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # 错误条
  coord_flip() +  # 翻转坐标轴，使条形图水平
  theme_minimal() +  # 使用简洁的主题
  theme(legend.position = "none")  # 不显示图例
check_path("plot/sfig4/variance_with_regression.png")
ggsave("plot/sfig4/variance_with_regression.png", width = 10, height = 10)




data = read.csv("data/sandbox_new_pheno_without_regression/whole/sandbox_new_pheno.csv", row.names = "buId")
data = data[rownames(data) %in% health_id, ]
df = data.frame(organ1=NA, organ2=NA, pheno = colnames(data), var=NA)
df$organ1 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][3])
df$organ2 = map_chr(df$pheno, ~str_split(.x, pattern = "_")[[1]][4])

df$var = map_dbl(data, sd, na.rm=TRUE)

df = pivot_longer(df, cols = c("organ1", "organ2"), values_to = "organ")

df = df %>%
  group_by(organ) %>%
  summarise(
    mean = mean(var, na.rm = TRUE),
    sd = sd(var, na.rm = TRUE),
    lower = mean - sd,
    upper = mean + sd
  )
df$organ_color = rename_based_on_df(df$organ, nmapdf = color_map, from = "organ", to = "color")
write.csv(df, "data/sfig4/sandbox_without_regression_variance.csv")


p <- ggplot(df, aes(x = organ, y = mean, fill = organ_color)) +
  geom_bar(stat = 'identity') +  # 条形图
  scale_fill_identity() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +  # 错误条
  coord_flip() +  # 翻转坐标轴，使条形图水平
  theme_minimal() +  # 使用简洁的主题
  theme(legend.position = "none")  # 不显示图例
ggsave("plot/sfig4/variance_without_regression.png", width = 10, height = 10)










####################### figure 2 柱状图 disease before ############
data = read.csv("data/sandbox_disease_odd.csv")
data$disease_organ = map_chr(data$disease, ~str_split(.x, "-")[[1]][2])
data$disease_organ = tools::toTitleCase(data$disease_organ)
data$disease_organ = map_chr(data$disease_organ, ~str_replace(.x, "Head", replacement = "Brain"))
data$disease_color = rename_based_on_df(data$disease_organ, nmapdf = color_map, from = "organ", to = "color")
# 按疾病分类，计算odds均值及p值显著的pheno个数（定义p < 0.05 为显著）
result = data %>%
  group_by(disease, disease_color, disease_organ) %>%
  summarise(
    mean_odds = mean(odds),                      # odds的均值
    significant_pheno_count = sum(p < 0.05)      # p值显著的pheno个数
  )
write.csv(result, "data/figure2_sandbox_bar.csv", row.names = FALSE)
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
ggsave("plot/figure2_sandbox_bar.png", width = 6, height = 6)



####################### figure 2 柱状图2 disease before ############
data = read.csv("data/sandbox_disease_odd.csv")


significant_organs = data %>%
  separate(col = pheno, into = c("y", "res", "Organ1", "Organ2"), sep = "_") %>%  # 分割pheno列
  filter(p < 0.05) %>%  # 筛选显著的疾病（p值小于0.05）
  select(Organ1, Organ2, p, disease) %>%  # 选择感兴趣的列
  gather(key = "OrganKey", value = "Organ", Organ1, Organ2) %>%  # 把器官列转换为长格式
  group_by(Organ) %>%  # 按器官分组
  summarise(Count = n())  # 计算每个器官关联的显著疾病的数量
significant_organs$organ_color = rename_based_on_df(significant_organs$Organ, 
                                                    nmapdf = color_map, from = "organ", to = "color")
significant_organs$Organ = factor(significant_organs$Organ, levels = color_map$organ)
write.csv(significant_organs, "data/figure2_sandbox_bar2.csv", row.names = FALSE)
# 绘制柱状图
ggplot(significant_organs, aes(x = Organ, y = Count, 
                               label = Count, fill=organ_color)) +
  geom_col() +
  scale_fill_identity() +
  geom_text(vjust = -0.3, size = 8) +  # 在柱子上方标示显著的pheno个数
  labs(x = "Organ", y = "Significant Count") +
  theme_minimal() +
  ylim(0, max(significant_organs$Count) * 1.1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 30, colour = "black"),
        axis.text.y = element_text(size = 30, colour = "black"),
        axis.title.x = element_text(size = 36, face = "bold"),  # 设置 X 轴标题文字大小
        axis.title.y = element_text(size = 36, face = "bold"), 
        plot.margin = unit(c(1.1, 0.5, 0.5, 0.5), "cm"), 
        panel.grid.major = element_blank(),  # 去除主要网格线
        panel.grid.minor = element_blank())  # 设置 Y 轴标题文字大小
ggsave("plot/figure2_sandbox_bar2.png", width = 6, height = 6)





####################### figure 3 气泡图 disease after ############
data = read.csv("data/sandbox_disease_odd.csv")
data = data %>%
  separate(col = pheno, into = c("y", "res", "Organ1", "Organ2"), sep = "_")  # 分割pheno列
data$sig = case_when(
  data$p > 0.05 ~ " ", 
  (data$p > 0.01 & data$p <= 0.05) ~ "*",
  (data$p > 0.001 & data$p <= 0.01) ~ "**", 
  (data$p <= 0.001) ~ "***"
)
data$sig = as.factor(data$sig)
data$logp = -log10(data$p)


data = data %>%
  dplyr::select(Organ1, Organ2, disease, odds, logp, sig)
data$disease = factor(data$disease)
data = data %>%
  gather(key = "OrganKey", value = "Organ", Organ1, Organ2)
data$y_shift = ifelse(data$OrganKey == "Organ1", -0.05, 0.05)
data$y_shift = data$y_shift * as.numeric(data$sig)
data$organ_color = rename_based_on_df(data$Organ, color_map, from = "organ", to = "color")
write.csv("data/fig3A_sandbox.csv", row.names = FALSE)

ggplot(data, aes(x = odds, y = disease, color=organ_color, size=sig)) +
  geom_point(alpha = 1, position = position_nudge(y=data$y_shift)) + # 设置透明度
  scale_size_discrete(range = c(1, 4)) + # 调整气泡的大小范围
  scale_color_identity() +
  labs(x = "OR",
       y = "Diseases") +
  theme_minimal() +
  geom_vline(xintercept = 1)
ggsave("plot/figure3A_sandbox.png", width = 12, height = 6)










# 
# read_feature_data = function() {
#   BMI <- read.csv(file = "Sandbox/Demography/DXA体成分-体成分高分辨率数字成像化系统.csv",
#                   header = TRUE, row.names = 1, check.names = FALSE)
#   BMI <- BMI[, 1, drop=FALSE]
#   sex_age <- read.csv(file = "Sandbox/Demography/基本信息.csv",
#                       header = TRUE, row.names = 1, check.names = FALSE)
#   sexual <- sex_age[, 1, drop=FALSE]
#   sexual <- map_df(sexual, ~ ifelse(. == "男", 1, 0))
#   rownames(sexual) <- rownames(sex_age)
#   age <- sex_age[, 2, drop = FALSE]
#   partial <- cbind(BMI, sexual, age)
# 
#   nc = read.csv(file = "Sandbox/Demography/人体外观测量-三维人体扫描分析系统.csv",
#                   header = TRUE, row.names = 1, check.names = FALSE)
# 
#   tjdata <- read.csv("Sandbox/Serum/志愿者招募-常规体检.csv",
#                      header = TRUE, check.names = FALSE, row.names = 1)
#   tjdata <- tjdata[c("血压(mm/Hg)收缩压（mm/Hg）", "血压(mm/Hg)舒张压（mm/Hg）", "葡萄糖（mmol/L）",
#                      "心率")]
#   features <- bind_cols(partial, nc, tjdata)
#   mapcnames <- c("BMI", "Gender", "Age", "NC", "SBP", "DBP", "Glucose", "HR")
#   colnames(features) <- mapcnames
#   return(features)
# }
# features = read_feature_data()
# 
# 
# for (i in seq_along(features)) {
#   df = data.frame()
#   for (j in seq_along(male_df)) {
#     feature = features[i]
#     feature_name = colnames(features)[i]
# 
#     axis = male_df[j]
#     axis_name = str_remove(colnames(male_df)[j], "y_res_")
# 
#     merged_data = merge_and_divide(feature, axis)
#     feature = merged_data[[1]]
#     axis = merged_data[[2]]
#     r = cor.test(feature[[1]], axis[[1]])$estimate
#     p = cor.test(feature[[1]], axis[[1]])$p.value
#     n = cor.test(feature[[1]], axis[[1]])$parameter + 1
#     temp_df = data.frame(axis = axis_name, feature = feature_name, r = r, p = p, n = n)
#     df = bind_rows(df, temp_df)
#   }
#   df["organ1"] = map_chr(str_split(df$axis, "_"), ~.x[1])
#   df["organ2"] = map_chr(str_split(df$axis, "_"), ~.x[2])
#   df["p_bonf"] = p.adjust(df$p, method = "bonferroni", n = length(features) * length(male_df))
# 
#   write.csv(df, paste0("data/5/y_res_male_", feature_name, "_data.csv"), row.names = FALSE)
# }











################# figure 4A sandbox ####################
features = read.csv("data/sandbox_feature.csv", row.names = "buId")
feature_dict = read.csv("Sandbox/Category/feature.csv", row.names = NULL)
feature_sd = map_dbl(features, ~sd(.x, na.rm = TRUE))
features = features[names(na.omit(feature_sd[feature_sd != 0]))]


# 小分子表型
moleculars = ReadMicrophenomics(dataPath = "Sandbox/")
moleculars = reduce(moleculars, my_join)


df_all = data.frame()
for (new_pheno_file in list.files("data/sandbox_new_pheno_without_regression/", full.names = TRUE, pattern = ".csv")) {
  print(paste0("calculating ", basename(new_pheno_file)))
  male_df = read.csv(new_pheno_file, row.names = "buId")
  df = data.frame()
  
  n = 0
  for (i in seq_along(features)) {
    for (j in seq_along(male_df)) {
      feature = features[i]
      feature_name = colnames(features)[i]
      if (class(features[[i]]) == "numeric" | class(features[[i]]) == "integer") {
        axis = male_df[j]
        axis_name = str_remove(colnames(male_df)[j], "y_res_")
        
        merged_data = inner_join(feature %>% rownames_to_column("buId"), 
                                 axis %>% rownames_to_column("buId"), 
                                 by="buId") %>% 
          column_to_rownames("buId") %>% 
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
  write.csv(df, paste0("data/5/sandbox_", colnames(male_df), "_data.csv"), row.names = FALSE)
  df_all = bind_rows(df_all, df)
}
df_all["p_bonf"] = p.adjust(df_all$p, method = "bonferroni", n = nrow(df_all))
write.csv(df_all, paste0("data/5/0_sandbox_all_data.csv"), row.names = FALSE)


data = read.csv("data/5/0_sandbox_all_data.csv")

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
new_data$Sum = rowSums(new_data)
df_sorted = new_data[order(new_data$Sum, decreasing = TRUE), ]
df_sorted$Sum = NULL
temp = df_sorted
temp$feature = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'names')
temp  = temp %>% dplyr::select(feature, everything())
temp$category = rename_based_on_df(rownames(temp), nmapdf = feature_dict, from = 'UDI', to = 'category_name')
write.csv(temp, file = "data/sandbox_feature_net_sort_by_wholenet.csv", row.names = TRUE)


df_sorted = read.csv("data/sandbox_feature_net_sort_by_wholenet.csv", row.names = 1)
df_sorted = df_sorted[!duplicated(df_sorted$feature), ]
df_sorted = df_sorted[1:30, ]

df_sorted_with_UDI = df_sorted
rownames(df_sorted) = NULL
df_sorted = df_sorted %>% column_to_rownames("feature")
df_sorted$category = NULL
hM <- format(round(as.matrix(df_sorted), 2))#对数据保留2位小数


check_path("plot/sandbox_heatmap/whole_net.png")
pheatmap::pheatmap(as.matrix(df_sorted), cluster_rows = FALSE, cluster_cols = FALSE, 
                   number_color = "black",
                   fontsize_number = 14,
                   fontsize = 20,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50), 
                   cellwidth=30,
                   cellheight=30,
                   filename = paste0("plot/sandbox_heatmap/whole_net.png"),
                   display_numbers = hM, 
                   width = 16,
                   height = 16,
                   units = "in", border_color = "black", 
                   res = 600)


# ############### 个体的图 #######################
############ 所有人的图 ############

for (i_buId in seq_along(rownames(disease))) {
  buId = rownames(disease)[i_buId]
  if ((identical(regression_data$性别[rownames(regression_data) %in% buId], 0)) & (buId %in% rownames(female_df))) {
    data = female_df[buId, ]
    data = abs(data)
    # 提取行名和列名
    row_names = color_map$organ[c(1:7, 9)]
    col_names = color_map$organ[c(1:7, 9)]
  }
  else if ((identical(regression_data$性别[rownames(regression_data) %in% buId], 1)) & (buId %in% rownames(male_df))) {
    data = male_df[buId, ]
    data = abs(data)
    # 提取行名和列名
    row_names = color_map$organ[c(1:8)]
    col_names = color_map$organ[c(1:8)]
  }
  else {next}
  
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
          new_data[row_names[i], col_names[j]] = NA
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]
          
        }
      }
      if (i > j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        if (!(data_col %in% colnames(data))) {
          new_data[row_names[i], col_names[j]] = NA
        }
        else {
          new_data[row_names[i], col_names[j]] = data[, data_col]
        }
      }
      if (i == j) {
        data_col = paste0("y_res_", col_names[j], "_", row_names[i])
        new_data[row_names[i], col_names[j]] = NA
        
      }
    }
  }
  
  # 打印结果
  print(new_data)
  
  
  new_data[is.na(new_data)] = 0
  check_path(paste0("plot/sfig8/complexnet/", buId, ".png"))
  check_path(paste0("plot/sfig8/heatmap/", buId, ".png"))
  draw_organ_complexnet(new_data, paste0("plot/sfig8/complexnet/", buId, ".png"))
  
  pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     main = buId,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     filename = paste0("plot/sfig8/heatmap/", buId, ".png"))
}


for (path in list.files("data/sfig8/", recursive = TRUE, full.names = TRUE)) {
  try(
    {
      data = read.csv(path, row.names = 1)
      data = data[1: (nrow(data) - 1), 1: (ncol(data) - 1)]
      new_data = data %>% rownames_to_column("organ1")
      new_data = pivot_longer(new_data, cols = -1, names_to = "organ2")
      new_data$organ1 = factor(new_data$organ1, levels = color_map$organ, ordered = TRUE)
      new_data$organ2 = factor(new_data$organ2, levels = color_map$organ, ordered = TRUE)
      new_data = new_data[new_data$organ1 < new_data$organ2, ]
      new_data$organ1 = as.character(new_data$organ1)
      new_data$organ2 = as.character(new_data$organ2)
      
      
      disease_name = basename(dirname(path))
      buId = file_path_sans_ext(basename(path))
      
      new_data[is.na(new_data)] = 0
      check_path(paste0("plot/sfig8/complexnet/", disease_name, "/", buId, ".png"))
      check_path(paste0("plot/sfig8/heatmap/", disease_name, "/", buId, ".png"))
      draw_organ_complexnet(new_data, paste0("plot/sfig8/complexnet/",disease_name, "/",  buId, ".png"))
      
      
      bk = seq(0, 2, length.out=50)
      pheatmap::pheatmap(data, cellwidth = 15, cellheight = 15,
                         cluster_rows = FALSE, cluster_cols = FALSE,
                         main = buId,
                         color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                         filename = paste0("plot/sfig8/heatmap/", disease_name, "/", buId, ".png"), 
                         breaks = bk,
                         legend_breaks = c(0, 1, 2)
                         )
      
    }
  )

}








# ##### sample 1 Kidney01-囊肿，lung02结节 #####
# 
# data = female_df["01210708V440", ]
# data = abs(data)
# 
# # 提取行名和列名
# row_names = color_map$organ[c(1:7, 9)]
# col_names = color_map$organ[c(1:7, 9)]
# 
# # 创建新的数据框
# new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
# colnames(new_data) = unique(col_names)
# rownames(new_data) = unique(row_names)
# 
# # 填充值
# for (i in 1:length(row_names)) {
#   for (j in 1:length(col_names)) {
#     if (i < j) {
#       data_col = paste0("y_res_", row_names[i], "_", col_names[j])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#         
#       }
#     }
#     if (i > j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#       }
#     }
#     if (i == j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       new_data[row_names[i], col_names[j]] = NA
#       
#     }
#   }
# }
# 
# # 打印结果
# print(new_data)
# 
# 
# new_data[is.na(new_data)] = 0
# pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15, 
#                    cluster_rows = FALSE, cluster_cols = FALSE, 
#                    main = "01210708V440", 
#                    color = colorRampPalette(c("#FEFEBD","firebrick3"))(50), 
#                    filename = "E:/OneDrive/桌面/sample1.png")
# 
# 
# 
# 
# 
# ########## sample 2 liver02 囊肿#########
# data = female_df["01210708V303", ]
# data = abs(data)
# 
# # 提取行名和列名
# row_names = color_map$organ[c(1:7, 9)]
# col_names = color_map$organ[c(1:7, 9)]
# 
# # 创建新的数据框
# new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
# colnames(new_data) = unique(col_names)
# rownames(new_data) = unique(row_names)
# 
# # 填充值
# for (i in 1:length(row_names)) {
#   for (j in 1:length(col_names)) {
#     if (i < j) {
#       data_col = paste0("y_res_", row_names[i], "_", col_names[j])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#         
#       }
#     }
#     if (i > j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#       }
#     }
#     if (i == j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       new_data[row_names[i], col_names[j]] = NA
#       
#     }
#   }
# }
# 
# # 打印结果
# print(new_data)
# 
# 
# new_data[is.na(new_data)] = 0
# pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15, 
#                    cluster_rows = FALSE, cluster_cols = FALSE, 
#                    main = "01210708V303", 
#                    color = colorRampPalette(c("#FEFEBD","firebrick3"))(50), 
#                    filename = "E:/OneDrive/桌面/sample2.png")
# 
# 
# 
# 
# 
# 
# 
# 
# ########## sample 3 liver05 胆结石#########
# data = male_df["01210625V006", ]
# data = abs(data)
# 
# # 提取行名和列名
# row_names = color_map$organ[c(1:8)]
# col_names = color_map$organ[c(1:8)]
# 
# # 创建新的数据框
# new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
# colnames(new_data) = unique(col_names)
# rownames(new_data) = unique(row_names)
# 
# # 填充值
# for (i in 1:length(row_names)) {
#   for (j in 1:length(col_names)) {
#     if (i < j) {
#       data_col = paste0("y_res_", row_names[i], "_", col_names[j])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#         
#       }
#     }
#     if (i > j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#       }
#     }
#     if (i == j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       new_data[row_names[i], col_names[j]] = NA
#       
#     }
#   }
# }
# 
# # 打印结果
# print(new_data)
# 
# 
# new_data[is.na(new_data)] = 0
# pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15, 
#                    cluster_rows = FALSE, cluster_cols = FALSE, 
#                    main = "01210625V006", 
#                    color = colorRampPalette(c("#FEFEBD","firebrick3"))(50), 
#                    filename = "E:/OneDrive/桌面/sample3.png")
# 
# 
# 
# 
# 
# 
# ########## sample 4 liver02 肝囊肿 ovary01 囊肿#########
# data = female_df["01210806V002", ]
# data = abs(data)
# 
# # 提取行名和列名
# row_names = color_map$organ[c(1:7, 9)]
# col_names = color_map$organ[c(1:7, 9)]
# 
# # 创建新的数据框
# new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
# colnames(new_data) = unique(col_names)
# rownames(new_data) = unique(row_names)
# 
# # 填充值
# for (i in 1:length(row_names)) {
#   for (j in 1:length(col_names)) {
#     if (i < j) {
#       data_col = paste0("y_res_", row_names[i], "_", col_names[j])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#         
#       }
#     }
#     if (i > j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#       }
#     }
#     if (i == j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       new_data[row_names[i], col_names[j]] = NA
#       
#     }
#   }
# }
# 
# # 打印结果
# print(new_data)
# 
# 
# new_data[is.na(new_data)] = 0
# pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15, 
#                    cluster_rows = FALSE, cluster_cols = FALSE, 
#                    main = "01210806V002", 
#                    color = colorRampPalette(c("#FEFEBD","firebrick3"))(50), 
#                    filename = "E:/OneDrive/桌面/sample4.png")
# 
# 
# 
# 
# 
# 
# 
# ########### sample 5 uterus02 囊肿#########
# data = female_df["01210602V002", ]
# data = abs(data)
# 
# # 提取行名和列名
# row_names = color_map$organ[c(1:7, 9)]
# col_names = color_map$organ[c(1:7, 9)]
# 
# # 创建新的数据框
# new_data <- data.frame(matrix(NA, nrow = length(unique(row_names)), ncol = length(unique(col_names))))
# colnames(new_data) = unique(col_names)
# rownames(new_data) = unique(row_names)
# 
# # 填充值
# for (i in 1:length(row_names)) {
#   for (j in 1:length(col_names)) {
#     if (i < j) {
#       data_col = paste0("y_res_", row_names[i], "_", col_names[j])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#         
#       }
#     }
#     if (i > j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       if (!(data_col %in% colnames(data))) {
#         new_data[row_names[i], col_names[j]] = NA
#       }
#       else {
#         new_data[row_names[i], col_names[j]] = data[, data_col]
#       }
#     }
#     if (i == j) {
#       data_col = paste0("y_res_", col_names[j], "_", row_names[i])
#       new_data[row_names[i], col_names[j]] = NA
#       
#     }
#   }
# }
# 
# # 打印结果
# print(new_data)
# 
# 
# new_data[is.na(new_data)] = 0
# pheatmap::pheatmap(new_data, cellwidth = 15, cellheight = 15, 
#                    cluster_rows = FALSE, cluster_cols = FALSE, 
#                    main = "01210602V002", 
#                    color = colorRampPalette(c("#FEFEBD","firebrick3"))(50), 
#                    filename = "E:/OneDrive/桌面/sample5.png")
# 










# 
# 
# 
# ################## popNet plot ###################
# output_cca_heatmap_matrix = matrix(NA, nrow =7, ncol = 7, dimnames = list(c("Heart", "Brain", "Kidney", "Liver", 
#                                                                             "Lung", "Pancreas", "Spleen"), 
#                                                                           c("Heart", "Brain", "Kidney", "Liver", 
#                                                                             "Lung", "Pancreas", "Spleen")))
# output_cca_heatmap = as.data.frame(output_cca_heatmap_matrix)
# output_cca_edge = data.frame()
# organs = list(Heart = without_sex_organs$Heart, Brain = without_sex_organs$Brain, Kidney = without_sex_organs$Kidney, 
#               Liver = without_sex_organs$Liver, Lung = without_sex_organs$Lung, Pancreas = without_sex_organs$Pancreas, 
#               Spleen = without_sex_organs$Spleen)
# for (i in seq_along(organs)) {
#   for (j in seq_along(organs)) {
#     if (i < j) {
#       organ1 = organs[[i]]
#       organ2 = organs[[j]]
#       organ1_name = names(organs)[i]
#       organ2_name = names(organs)[j]
#       results = get_cca_result(organ1, organ2)
#       output_cca_heatmap_matrix[i, j] = results$cca
#       output_cca_heatmap_matrix[j, i] = results$cca
#       
#       temp_output_cca_edge = data.frame(organ1 = organ1_name, organ2 = organ2_name, 
#                                         r = results$cca, p = results$cca.p, n=results$N)
#       output_cca_edge = bind_rows(output_cca_edge, temp_output_cca_edge)
#     }
#     else if (i == j) {
#       output_cca_heatmap_matrix[i, j] = 1
#       
#     }
#   }
# }
# 
# draw_organ_complexnet(output_cca_edge, "plot/complexnet/sandbox_organ_axis_cca_complexnet.png")
# pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE, 
#                    color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50), 
#                    cellwidth=12,
#                    cellheight=12,
#                    filename = paste0("plot/heatmap/sandbox_organ_axis_cca_heatmap.png"),
#                    width = 3,
#                    height = 3,
#                    units = "in",
#                    res = 600)
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
# 
# 
# 
# ################### risk factor plot #################
# features = read.csv("data/sandbox_feature.csv", row.names = "buId", check.names = FALSE)
# feature_sd = map_dbl(features, ~sd(.x, na.rm = TRUE))
# features = features[names(na.omit(feature_sd[feature_sd != 0]))]
# 
# 
# df_all = data.frame()
# for (new_pheno_file in list.files("data/sandbox_new_pheno", full.names = TRUE, pattern = ".csv")) {
#   male_df = read.csv(new_pheno_file) %>% column_to_rownames(var = "buId")
#   df = data.frame()
#   n = 0
#   for (i in seq_along(features)) {
#     for (j in seq_along(male_df)) {
#       feature = features[i]
#       feature_name = colnames(features)[i]
#       if (class(features[[i]]) == "numeric" | class(features[[i]]) == "integer") {
#         axis = male_df[j]
#         axis_name = str_remove(colnames(male_df)[j], "y_res_")
#         # axis_name = str_remove(axis_name, "x_res_")
#         # axis_name = str_remove(axis_name, "dis_")
#         
#         merged_data = inner_join(feature %>% rownames_to_column("buId"), axis %>% rownames_to_column("buId")) %>% na.omit() 
#         rownames(merged_data) = NULL
#         merged_data = merged_data %>% column_to_rownames("buId")
#         if ((sd(merged_data[[2]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
#             r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
#             p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
#             n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
#             temp_df = data.frame(pheno = colnames(male_df)[j], axis = axis_name, feature = feature_name, r = r, p = p, n = n)
#             df = bind_rows(df, temp_df)
#             
#             
#             lm_fit = lm(my_formula(y_variable = colnames(merged_data)[2], 
#                                    x_variables = colnames(merged_data)[1]),
#                         data = merged_data)
#             residuals = residuals(lm_fit)
#             intercept = lm_fit$coefficients[1]
#             slope = lm_fit$coefficients[2]
#             top_outliers <- names(sort(abs(residuals), decreasing = TRUE)[1:10])
#             outlier_data = merged_data[top_outliers, ]
#             outlier_text = outlier_data
#             outlier_text$label = paste0("buId ", rownames(outlier_data))
#             
#             My_Theme = theme(
#               panel.background = element_blank(), 
#               title = element_text(size = 7),
#               text = element_text(size = 6), 
#             )
#             
#             
#             risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2]))) + 
#               geom_point(size=0.5, col="#2171b5") +
#               geom_abline(intercept = intercept, slope = slope, col="red") +
#               geom_hline(yintercept = 0, col="black") +
#               theme(axis.line = element_line(color="black", size = 0.2)) +
#               mdthemes::md_theme_classic() +
#               labs(x = colnames(merged_data)[1], 
#                    y = colnames(merged_data)[2],
#                    title=paste0("scatter plot"), 
#                    subtitle = paste0("p = ", p)) +
#               My_Theme + 
#               geom_point(data=outlier_data, 
#                          aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2])),
#                          size=2, col="darkorange")+
#               ggrepel::geom_text_repel(data=outlier_text,
#                               aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2]), 
#                                   label=label),color="darkorange",seed=5678)
#             ggsave(paste0("plot/", gsub("[\"'/:*?,<>|%]", "_", colnames(merged_data)[1]), "_", colnames(merged_data)[2], ".png"))
#             n = n + 1
#         }
#       }
#     }
#   }
#   df["organ1"] = map_chr(str_split(df$axis, "_"), ~.x[1])
#   df["organ2"] = map_chr(str_split(df$axis, "_"), ~.x[2])
#   df["p_bonf"] = p.adjust(df$p, method = "bonferroni", n = nrow(df))
#   write.csv(df, paste0("data/5/sandbox_risk_", colnames(male_df), "_data.csv"), row.names = FALSE)
#   df_all = bind_rows(df_all, df)
# }
# df_all["p_bonf"] = p.adjust(df_all$p, method = "bonferroni", n = nrow(df_all))
# check_path(paste0("data/5/0_sandbox_risk_data.csv"))
# write.csv(df_all, paste0("data/5/0_sandbox_risk_data.csv"), row.names = FALSE)
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
# 
# 
# 
# 
# 
# 
# 
# 
# ################### micro phenotype plot #################
# micro_files = list.files("data/sandbox/", full.names = TRUE)
# for (micro_file in micro_files) {
#   micro_cate = file_path_sans_ext(basename(micro_file))
#   
#   micros = read.csv(micro_file, row.names = "buId", check.names = FALSE)
#   micro_sd = map_dbl(micros, ~sd(.x, na.rm = TRUE))
#   micros = micros[names(na.omit(micro_sd[micro_sd != 0]))]
# 
#   
#   df_all = data.frame()
#   for (new_pheno_file in list.files("data/sandbox_new_pheno", full.names = TRUE, pattern = ".csv")) {
#     male_df = read.csv(new_pheno_file) %>% column_to_rownames(var = "buId")
#     df = data.frame()
#     n = 0
#     for (i in seq_along(micros)) {
#       for (j in seq_along(male_df)) {
#         micro = micros[i]
#         micro_name = colnames(micros)[i]
#         if (class(micros[[i]]) == "numeric" | class(micros[[i]]) == "integer") {
#           axis = male_df[j]
#           axis_name = str_remove(colnames(male_df)[j], "y_res_")
#           # axis_name = str_remove(axis_name, "x_res_")
#           # axis_name = str_remove(axis_name, "dis_")
#           
#           merged_data = inner_join(micro %>% rownames_to_column("buId"), axis %>% rownames_to_column("buId")) %>% na.omit() 
#           rownames(merged_data) = NULL
#           merged_data = merged_data %>% column_to_rownames("buId")
#           if ((sd(merged_data[[2]], na.rm = TRUE) != 0) & (nrow(merged_data) > 100)) {
#             r = cor.test(merged_data[[1]], merged_data[[2]])$estimate
#             p = cor.test(merged_data[[1]], merged_data[[2]])$p.value
#             n = cor.test(merged_data[[1]], merged_data[[2]])$parameter + 1
#             temp_df = data.frame(pheno = colnames(male_df)[j], axis = axis_name, micro = micro_name, micro_cate = micro_cate, 
#                                  r = r, p = p, n = n)
#             df = bind_rows(df, temp_df)
#             # lm_fit = lm(my_formula(y_variable = colnames(merged_data)[2], 
#             #                        x_variables = colnames(merged_data)[1]),
#             #             data = merged_data)
#             # residuals = residuals(lm_fit)
#             # intercept = lm_fit$coefficients[1]
#             # slope = lm_fit$coefficients[2]
#             # top_outliers <- names(sort(abs(residuals), decreasing = TRUE)[1:10])
#             # outlier_data = merged_data[top_outliers, ]
#             # outlier_text = outlier_data
#             # outlier_text$label = paste0("buId ", rownames(outlier_data))
#             # 
#             # My_Theme = theme(
#             #   panel.background = element_blank(), 
#             #   title = element_text(size = 7),
#             #   text = element_text(size = 6), 
#             # )
#             # 
#             # 
#             # risk_scatter = ggplot(data = merged_data, aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2]))) + 
#             #   geom_point(size=0.5, col="#2171b5") +
#             #   geom_abline(intercept = intercept, slope = slope, col="red") +
#             #   geom_hline(yintercept = 0, col="black") +
#             #   theme(axis.line = element_line(color="black", size = 0.2)) +
#             #   mdthemes::md_theme_classic() +
#             #   labs(x = colnames(merged_data)[1], 
#             #        y = colnames(merged_data)[2],
#             #        title=paste0("scatter plot"), 
#             #        subtitle = paste0("p = ", p)) +
#             #   My_Theme + 
#             #   geom_point(data=outlier_data, 
#             #              aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2])),
#             #              size=2, col="darkorange")+
#             #   ggrepel::geom_text_repel(data=outlier_text,
#             #                            aes(x = get(colnames(merged_data)[1]), y = get(colnames(merged_data)[2]), 
#             #                                label=label),color="darkorange",seed=5678)
#             # ggsave(paste0("plot/sandbox_micro_", micro_cate, "/", gsub("[\"'/:*?,<>|%]", "_", colnames(merged_data)[1]), "_", colnames(merged_data)[2], ".png"))
#             n = n + 1
#           }
#         }
#       }
#     }
#     df["organ1"] = map_chr(str_split(df$axis, "_"), ~.x[1])
#     df["organ2"] = map_chr(str_split(df$axis, "_"), ~.x[2])
#     df["p_bonf"] = p.adjust(df$p, method = "bonferroni", n = nrow(df))
#     write.csv(df, paste0("data/5/sandbox_micro_", colnames(male_df)[1], "_", micro_cate, "_data.csv"), row.names = FALSE)
#     df_all = bind_rows(df_all, df)
#   }
#   df_all["p_bonf"] = p.adjust(df_all$p, method = "bonferroni", n = nrow(df_all))
#   check_path(paste0("data/5/0_sandbox_micro_", micro_cate, "_data.csv"))
#   write.csv(df_all, paste0("data/5/0_sandbox_micro_", micro_cate, "_data.csv"), row.names = FALSE)
# }
# 









################## cytoscope data modality 1#######################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


prefix = rep(names(organs), as.numeric(map(organs, length)))
organs = flatten(organs)
names(organs) = paste0(prefix, "_", names(organs))


category_organ = names(organs)
category_organ[grep("Morphology", names(organs))] = "Morphology"
category_organ[grep("Texture", names(organs))] = "Texture"
category_organ[grep("Density", names(organs))] = "Water"
category_organ[grep("Water", names(organs))] = "Water"
category_organ[grep("Fat", names(organs))] = "Fat"
category_organ[grep("Iron", names(organs))] = "Iron"
category_organ[grep("Diffusion", names(organs))] = "Diffusion"
category_organ[grep("Function", names(organs))] = "Function"
category_organ[grep("Vessel", names(organs))] = "Vessel"
category_organ[grep("ECG", names(organs))] = "ECG"
category_organ[grep("EEG", names(organs))] = "EEG"

df = data.frame()
for (i in seq_along(organs)) {
  for(j in seq_along(organs)) {
    if (i < j) {
      category_x = category_organ[i]
      category_y = category_organ[j]
      x = organs[[i]]
      y = organs[[j]]
      x_name = names(organs)[i]
      y_name = names(organs)[j]
      x_organ = str_split(x_name, pattern = "_")[[1]][1]
      y_organ = str_split(y_name, pattern = "_")[[1]][1]
      
      results = get_cca_result(x, y)
      cca = results$cca
      cca.p = results$cca.p
      temp_df = data.frame(name1 = x_name, 
                           name1_category = category_x, 
                           organ1 = x_organ,
                           name2 = y_name, 
                           name2_category = category_y,
                           organ2 = y_organ, 
                           r = cca, p = cca.p
      )
      df = rbind(df, temp_df)
    }
  }
}

node = data.frame(node=c(df$name1, df$name2), 
                  node_category=c(df$name1_category, df$name2_category), 
                  node_organ=c(df$organ1, df$organ2))
node = node[!duplicated(node$node), ]
edge = df

# 使用col_map来使用两个色系
edge$col_map = ifelse(edge$organ1 == edge$organ2, edge$r, edge$r + 10)

node_names = node$node

for (node_name in node_names){
  node$size[which(node$node %in% node_name)]=mean(edge$r[which(edge$name1 %in% node_name|edge$name2 %in% node_name)])
}

node$node_organ = factor(node$node_organ, levels = color_map$organ[1:7])


# 增加pos列
pos = get_position(category = node$node_organ, r_n2 = 0.3, coef = 100)
pos[, 1] = pos[, 1] - min(pos[, 1])
pos[, 2] = pos[, 2] - min(pos[, 2])
node = cbind(node, as.data.frame(pos))

network = list(node = node, edge = edge)
check_path(str_glue("data/circle_circle/modality/cca.xlsx"))
write_xlsx(network, path = str_glue("data/circle_circle/modality/cca.xlsx"))


################## cytoscope data modality 2#######################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


prefix = rep(names(organs), as.numeric(map(organs, length)))
organs = flatten(organs)
names(organs) = paste0(prefix, "_", names(organs))
category_organ = names(organs)


all_data = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name, col_seq = "__")
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      if (nrow(all_data) == 0) {
        all_data = output_new_pheno
      }
      else {
        all_data = full_join(all_data, output_new_pheno,by="buId")
      }
    }
  }
}

health_and_disease_id_list = list()
health_and_disease_id_list$health = health_id
for (i in seq_along(disease)) {
  disease_name = colnames(disease)[i]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  health_and_disease_id_list[[disease_name]] = disease_id
}

for (i in seq_along(health_and_disease_id_list)) {
  id_name = names(health_and_disease_id_list)[i]
  id = health_and_disease_id_list[[i]]
  data = all_data[all_data$buId %in% id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
  
  data = data.frame(map(data, ~mean(abs(.x), na.rm=TRUE)))
  colnames(data) = str_remove(colnames(data), "y_res_")
  data = pivot_longer(data, cols = everything())
  data = data %>%
    separate(col = name, into = c("name1", "name2"), sep = "__")
  data$name1_category = map_chr(data$name1, ~str_split(.x, pattern = "_")[[1]][2])
  data$name2_category = map_chr(data$name2, ~str_split(.x, pattern = "_")[[1]][2])
  data$organ1 = map_chr(data$name1, ~str_split(.x, pattern = "_")[[1]][1])
  data$organ2 = map_chr(data$name2, ~str_split(.x, pattern = "_")[[1]][1])
  
  
  node = data.frame(node=c(data$name1, data$name2), 
                    node_category=c(data$name1_category, data$name2_category), 
                    node_organ=c(data$organ1, data$organ2))
  node = node[!duplicated(node$node), ]
  edge = data
  # 使用col_map来使用两个色系
  edge$col_map = ifelse(edge$organ1 == edge$organ2, edge$value, edge$value + 10)
  node_names = node$node
  for (node_name in node_names){
    node$size[which(node$node %in% node_name)]=mean(edge$value[which(edge$name1 %in% node_name|edge$name2 %in% node_name)])
  }
  node$node_organ = factor(node$node_organ, levels = color_map$organ[1:7])
  
  # 增加pos列
  pos = get_position(category = node$node_organ, r_n2 = 0.3, coef = 100)
  pos[, 1] = pos[, 1] - min(pos[, 1])
  pos[, 2] = pos[, 2] - min(pos[, 2])
  node = cbind(node, as.data.frame(pos))
  
  network = list(node = node, edge = edge)
  
  check_path(str_glue("data/circle_circle/modality/{id_name}.xlsx"))
  write_xlsx(network, path = str_glue("data/circle_circle/modality/{id_name}.xlsx"))
}


#################### cytoscope data subarea #######################
organs = read_data("scripts/category/subarea/", class = "SubArea", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


prefix = rep(names(organs), as.numeric(map(organs, length)))
organs = flatten(organs)
names(organs) = paste0(prefix, "_", names(organs))


category_organ = names(organs)
category_organ = sub("^[^_]*_", "", category_organ)


df = data.frame()
for (i in seq_along(organs)) {
  for(j in seq_along(organs)) {
    if (i < j) {
      category_x = category_organ[i]
      category_y = category_organ[j]
      x = organs[[i]]
      y = organs[[j]]
      x_name = names(organs)[i]
      y_name = names(organs)[j]
      x_organ = str_split(x_name, pattern = "_")[[1]][1]
      y_organ = str_split(y_name, pattern = "_")[[1]][1]
      
      results = get_cca_result(x, y)
      cca = results$cca
      cca.p = results$cca.p
      temp_df = data.frame(name1 = x_name, 
                           name1_category = category_x, 
                           organ1 = x_organ,
                           name2 = y_name, 
                           name2_category = category_y,
                           organ2 = y_organ, 
                           r = cca, p = cca.p
      )
      df = rbind(df, temp_df)
    }
  }
}

node = data.frame(node=c(df$name1, df$name2), 
                  node_category=c(df$name1_category, df$name2_category), 
                  node_organ=c(df$organ1, df$organ2))
node = node[!duplicated(node$node), ]
edge = df

# 使用col_map来使用两个色系
edge$col_map = ifelse(edge$organ1 == edge$organ2, edge$r, edge$r + 10)


node_names = node$node

for (node_name in node_names){
  node$size[which(node$node %in% node_name)]=mean(edge$r[which(edge$name1 %in% node_name|edge$name2 %in% node_name)])
}

node$node_organ = factor(node$node_organ, levels = color_map$organ[1:7])


# 增加pos列
pos = get_position(category = node$node_organ, r_n2 = 0.3, coef = 100)
pos[, 1] = pos[, 1] - min(pos[, 1])
pos[, 2] = pos[, 2] - min(pos[, 2])
node = cbind(node, as.data.frame(pos))

check_path("data/circle_circle/subarea/cca.xlsx")
network = list(node = node, edge = edge)
write_xlsx(network, path = "data/circle_circle/subarea/cca.xlsx")




################## cytoscope data subarea 2#######################
organs = read_data("scripts/category/subarea/", class = "SubArea", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


prefix = rep(names(organs), as.numeric(map(organs, length)))
organs = flatten(organs)
names(organs) = paste0(prefix, "_", names(organs))
category_organ = names(organs)


all_data = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name, col_seq = "__")
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      if (nrow(all_data) == 0) {
        all_data = output_new_pheno
      }
      else {
        all_data = full_join(all_data, output_new_pheno,by="buId")
      }
    }
  }
}

health_and_disease_id_list = list()
health_and_disease_id_list$health = health_id
for (i in seq_along(disease)) {
  disease_name = colnames(disease)[i]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  health_and_disease_id_list[[disease_name]] = disease_id
}

for (i in seq_along(health_and_disease_id_list)) {
  id_name = names(health_and_disease_id_list)[i]
  id = health_and_disease_id_list[[i]]
  data = all_data[all_data$buId %in% id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
  
  data = data.frame(map(data, ~mean(abs(.x), na.rm=TRUE)))
  colnames(data) = str_remove(colnames(data), "y_res_")
  data = pivot_longer(data, cols = everything())
  data = data %>%
    separate(col = name, into = c("name1", "name2"), sep = "__")
  
  data$name1_category = sub("^[^_]*_", "", data$name1)
  data$name2_category = sub("^[^_]*_", "", data$name2)
  data$organ1 = map_chr(data$name1, ~str_split(.x, pattern = "_")[[1]][1])
  data$organ2 = map_chr(data$name2, ~str_split(.x, pattern = "_")[[1]][1])
  
  
  node = data.frame(node=c(data$name1, data$name2), 
                    node_category=c(data$name1_category, data$name2_category), 
                    node_organ=c(data$organ1, data$organ2))
  node = node[!duplicated(node$node), ]
  edge = data
  # 使用col_map来使用两个色系
  edge$col_map = ifelse(edge$organ1 == edge$organ2, edge$value, edge$value + 10)
  node_names = node$node
  for (node_name in node_names){
    node$size[which(node$node %in% node_name)]=mean(edge$value[which(edge$name1 %in% node_name|edge$name2 %in% node_name)])
  }
  node$node_organ = factor(node$node_organ, levels = color_map$organ[1:7])
  
  # 增加pos列
  pos = get_position(category = node$node_organ, r_n2 = 0.3, coef = 100)
  pos[, 1] = pos[, 1] - min(pos[, 1])
  pos[, 2] = pos[, 2] - min(pos[, 2])
  node = cbind(node, as.data.frame(pos))
  
  network = list(node = node, edge = edge)
  
  check_path(str_glue("data/circle_circle/subarea/{id_name}.xlsx"))
  write_xlsx(network, path = str_glue("data/circle_circle/subarea/{id_name}.xlsx"))
}





########################### sfig 7 Morphology ############################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = map(organs, ~.x$Morphology)

organs = subset(organs, ((names(organs) != "Prostate") & names(organs) != "Uterus"))



without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}



############### 健康人的器官间热图，使用y_res作为值 ###########
output_cca_edge_all = list()

data = without_sex_df
data = data[data$buId %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("buId")

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

output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
# output_cca_edge = na.omit(output_cca_edge)

output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")

output_cca_edge_all[["health"]] = output_cca_edge

draw_organ_complexnet(output_cca_edge,
                      output_name = "plot/sfig7/morphology/health_sandbox_organ_cca_complexnet.png",
                      limits = c(0.5, 1.5))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/sfig7/morphology/health_sandbox_organ_cca_heatmap.png"))
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/sfig7/morphology/health_sandbox_organ_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## 疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease)) {
  disease_name = colnames(disease)[k]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  
  data = without_sex_df
  data = data[data$buId %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
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
  output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
  output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
  output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
  # output_cca_edge = na.omit(output_cca_edge)
  
  output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")
  
  output_cca_edge_all[[disease_name]] = output_cca_edge
  
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/morphology/", disease_name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/sfig7/morphology/", disease_name, "_sandbox_organ_cca_heatmap.png"))
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/sfig7/morphology/", disease_name, "_sandbox_organ_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
check_path("data/sfig7/morphology/sandbox.xlsx")
write_xlsx(output_cca_edge_all, "data/sfig7/morphology/sandbox.xlsx")









########################### sfig 7 Function ############################
# Composition/Function的类别包括ECG/EEG/LungFunction/Water
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Prostate") & names(organs) != "Uterus"))

organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
organs = organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")
f = function(df) {
  colnames(df) = gsub(pattern = "ECG", replacement = "Function", x = colnames(df))
  colnames(df) = gsub(pattern = "EEG", replacement = "Function", x = colnames(df))
  colnames(df) = gsub(pattern = "Water", replacement = "Function", x = colnames(df))
  df = df[, grepl(pattern = "Function", colnames(df))]
  return(df)
}
organs = map(organs, f)
organs = organs[!map_vec(organs, is.null)]


without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}



############### 健康人的器官间热图，使用y_res作为值 ###########
output_cca_edge_all = list()

data = without_sex_df
data = data[data$buId %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("buId")

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

output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
# output_cca_edge = na.omit(output_cca_edge)

output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")

output_cca_edge_all[["health"]] = output_cca_edge

draw_organ_complexnet(output_cca_edge,
                      output_name = "plot/sfig7/function/health_sandbox_organ_cca_complexnet.png",
                      limits = c(0.5, 1.5))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/sfig7/function/health_sandbox_organ_cca_heatmap.png"))
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/sfig7/function/health_sandbox_organ_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## 疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease)) {
  disease_name = colnames(disease)[k]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  
  data = without_sex_df
  data = data[data$buId %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
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
  output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
  output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
  output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
  # output_cca_edge = na.omit(output_cca_edge)
  
  output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")
  
  output_cca_edge_all[[disease_name]] = output_cca_edge
  
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/function/", disease_name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/sfig7/function/", disease_name, "_sandbox_organ_cca_heatmap.png"))
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/sfig7/function/", disease_name, "_sandbox_organ_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
check_path("data/sfig7/function/sandbox.xlsx")
write_xlsx(output_cca_edge_all, "data/sfig7/function/sandbox.xlsx")








########################### sfig 7 Vessel ############################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = map(organs, ~.x$Vessel)

organs = subset(organs, ((names(organs) != "Prostate") & names(organs) != "Uterus"))



without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
       if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}



############### 健康人的器官间热图，使用y_res作为值 ###########
output_cca_edge_all = list()

data = without_sex_df
data = data[data$buId %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("buId")

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

output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
# output_cca_edge = na.omit(output_cca_edge)

output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")

output_cca_edge_all[["health"]] = output_cca_edge

draw_organ_complexnet(output_cca_edge,
                      output_name = "plot/sfig7/vessel/health_sandbox_organ_cca_complexnet.png",
                      limits = c(0.5, 1.5))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/sfig7/vessel/health_sandbox_organ_cca_heatmap.png"))
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/sfig7/vessel/health_sandbox_organ_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## 疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease)) {
  disease_name = colnames(disease)[k]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  
  data = without_sex_df
  data = data[data$buId %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
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
  output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
  output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
  output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
  # output_cca_edge = na.omit(output_cca_edge)
  
  output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")
  
  output_cca_edge_all[[disease_name]] = output_cca_edge
  
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/vessel/", disease_name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/sfig7/vessel/", disease_name, "_sandbox_organ_cca_heatmap.png"))
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/sfig7/vessel/", disease_name, "_sandbox_organ_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
check_path("data/sfig7/vessel/sandbox.xlsx")
write_xlsx(output_cca_edge_all, "data/sfig7/vessel/sandbox.xlsx")






########################### sfig 7 Texture ############################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = map(organs, ~.x$Texture)

organs = subset(organs, ((names(organs) != "Prostate") & names(organs) != "Uterus"))
organs = organs[!map_vec(organs, is.null)]


without_sex_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      health_id = rownames(health_df)
      # 使用健康的人计算u,v得到健康和疾病人的y_res
      results = gene_new_pheno(organ1, organ2, organ1_name, organ2_name)
      output_new_pheno = as.data.frame(results) 
      output_new_pheno = output_new_pheno %>% rownames_to_column("buId")
      if (nrow(without_sex_df) == 0) {
        without_sex_df = output_new_pheno
      }
      else {
        without_sex_df = full_join(without_sex_df, output_new_pheno,by="buId")
      }
    }
  }
}



############### 健康人的器官间热图，使用y_res作为值 ###########
output_cca_edge_all = list()

data = without_sex_df
data = data[data$buId %in% health_id, ]
rownames(data) = NULL
data = data %>% column_to_rownames("buId")

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

output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
# output_cca_edge = na.omit(output_cca_edge)

output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")

output_cca_edge_all[["health"]] = output_cca_edge

draw_organ_complexnet(output_cca_edge,
                      output_name = "plot/sfig7/texture/health_sandbox_organ_cca_complexnet.png",
                      limits = c(0.5, 1.5))
bk = seq(0.49, 1.51, length.out=50)
check_path(paste0("plot/sfig7/texture/health_sandbox_organ_cca_heatmap.png"))
pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                   color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                   cellwidth=12,
                   cellheight=12,
                   filename = paste0("plot/sfig7/texture/health_sandbox_organ_cca_heatmap.png"),
                   width = 3,
                   height = 3,
                   units = "in",
                   res = 600,
                   breaks = bk,
                   legend_breaks = c(0.5, 1, 1.5))


########## 疾病人的器官间热图，使用average y_res作为值 ###########
for (k in seq_along(disease)) {
  disease_name = colnames(disease)[k]
  disease_id = rownames(disease)[disease[disease_name] == 1]
  
  data = without_sex_df
  data = data[data$buId %in% disease_id, ]
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
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
  output_cca_edge$organ1 = factor(output_cca_edge$organ1, levels = color_map$organ, ordered = TRUE)
  output_cca_edge$organ2 = factor(output_cca_edge$organ2, levels = color_map$organ, ordered = TRUE)
  output_cca_edge = output_cca_edge[output_cca_edge$organ1 < output_cca_edge$organ2, ]
  # output_cca_edge = na.omit(output_cca_edge)
  
  output_cca_edge$organ1 = rename_based_on_df(as.character(output_cca_edge$organ1), color_map, "organ", "short_organ")
  output_cca_edge$organ2 = rename_based_on_df(as.character(output_cca_edge$organ2), color_map, "organ", "short_organ")
  
  output_cca_edge_all[[disease_name]] = output_cca_edge
  
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/texture/", disease_name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5))
  bk = seq(0.49, 1.51, length.out=50)
  check_path(paste0("plot/sfig7/texture/", disease_name, "_sandbox_organ_cca_heatmap.png"))
  pheatmap::pheatmap(output_cca_heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("#FFFF33","orange","firebrick3"))(50),
                     cellwidth=12,
                     cellheight=12,
                     filename = paste0("plot/sfig7/texture/", disease_name, "_sandbox_organ_cca_heatmap.png"),
                     width = 3,
                     height = 3,
                     units = "in",
                     res = 600,
                     breaks = bk,
                     legend_breaks = c(0.5, 1, 1.5))
}
names(output_cca_edge_all) = substr(names(output_cca_edge_all), 1, 30)
check_path("data/sfig7/texture/sandbox.xlsx")
write_xlsx(output_cca_edge_all, "data/sfig7/texture/sandbox.xlsx")























file_path  = "liumeng_240114/data/sfig7/morphology/sandbox.xlsx"
sheet_names <- excel_sheets(file_path)
output_cca_edge_all = lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})
names(output_cca_edge_all) = sheet_names
for (i in seq_along(output_cca_edge_all)) {
  output_cca_edge = output_cca_edge_all[[i]]
  name = names(output_cca_edge_all)[i]
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/morphology/", name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5), node_size = c(0.5, 2), col = "#7FCDC0")
}






file_path  = "liumeng_240114/data/sfig7/vessel/sandbox.xlsx"
sheet_names <- excel_sheets(file_path)
output_cca_edge_all = lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})
names(output_cca_edge_all) = sheet_names
for (i in seq_along(output_cca_edge_all)) {
  output_cca_edge = output_cca_edge_all[[i]]
  name = names(output_cca_edge_all)[i]
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/vessel/", name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5), node_size = c(0.5, 2), col = "#E37F83")
}







file_path  = "liumeng_240114/data/sfig7/function/sandbox.xlsx"
sheet_names <- excel_sheets(file_path)
output_cca_edge_all = lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})
names(output_cca_edge_all) = sheet_names
for (i in seq_along(output_cca_edge_all)) {
  output_cca_edge = output_cca_edge_all[[i]]
  name = names(output_cca_edge_all)[i]
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/function/", name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5), node_size = c(0.5, 2), col = "#9DAAC4")
}






file_path  = "liumeng_240114/data/sfig7/texture/sandbox.xlsx"
sheet_names <- excel_sheets(file_path)
output_cca_edge_all = lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})
names(output_cca_edge_all) = sheet_names
for (i in seq_along(output_cca_edge_all)) {
  output_cca_edge = output_cca_edge_all[[i]]
  name = names(output_cca_edge_all)[i]
  draw_organ_complexnet(output_cca_edge,
                        output_name = paste0("plot/sfig7/texture/", name, "_sandbox_organ_cca_complexnet.png"),
                        limits = c(0.5, 1.5), node_size = c(0.5, 2), col = "#C59C95")
}

