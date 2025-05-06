library("ggplot2")
library("tidyverse")

force_all = FALSE
############# ukb #######################
######### 计算ukb的表型的关联
if (force_all) {
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
  
  full_organs = organs
  
  
  # 设置用于并行计算的核心数
  no_cores <- 120
  
  # 使用mclapply进行并行计算
  results_list <- mclapply(seq_along(full_organs), function(i) {
    i_organ = full_organs[[i]]
    i_organ_name = names(full_organs)[i]
    
    mclapply(seq_along(full_organs), function(j) {
      j_organ = full_organs[[j]]
      j_organ_name = names(full_organs)[j]
      
      temp_df_list <- list()
      
      for (k in seq_along(i_organ)) {
        for (l in seq_along(j_organ)) {
          ik_name = colnames(i_organ)[k]
          jl_name = colnames(j_organ)[l]
          ik = i_organ[, k, drop=FALSE]
          jl = j_organ[, l, drop=FALSE]
          
          intersect_rows = intersect(rownames(ik), rownames(jl))
          results = cor.test(ik[intersect_rows, ], jl[intersect_rows, ])
          
          temp_df_list[[length(temp_df_list) + 1]] <- data.frame(organ1=i_organ_name, organ2=j_organ_name, 
                                                                 organ1_pheno=ik_name, organ2_pheno=jl_name, 
                                                                 r=results$estimate, p=results$p.value)
        }
      }
      do.call(rbind, temp_df_list)
    }, mc.cores = no_cores)
  }, mc.cores = no_cores)
  
  # 组合结果
  df <- do.call(rbind, do.call(rbind, results_list))
  
  
  check_path("data/organ_pheno_organ_pheno/ukb_r.csv")
  write.csv(df, "data/organ_pheno_organ_pheno/ukb_r.csv", row.names = FALSE)
}







data = read.csv("data/organ_pheno_organ_pheno/sandbox_r.csv")
data$p = p.adjust(data$p, method = "bonferroni")
data = data[data$p < 0.01, ]


a = sum((data$organ1 == data$organ2) & (data$p <= 0.01))/dim(data)[1]
b = sum((data$organ1 != data$organ2) & (data$p <= 0.01))/dim(data)[1]
c = sum(data$p > 0.01)/dim(data)[1]

df = data.frame(label = c("Significant intra-organ",
                          "Significant organ-to-organ",
                          "Not significantly"),
                percent = c(a, b, c))



check_path("plot/pie/sandbox_pie.png")
png("plot/pie/sandbox_pie.png", width = 1600, height = 1600)
ggpubr::ggdonutchart(df, "percent", fill="label", label = df$label, font.family = "sans",
                     lab.font = c(5, "bold", "black"), palette = c("grey", "#619CFF", "#F8766D"), 
                     color="white")
dev.off()






data_filter = data %>% filter(p <= 0.01)

interorgan_filter = as.data.frame(data_filter[data_filter$organ1 == data_filter$organ2, ])
crossorgan_filter = as.data.frame(data_filter[data_filter$organ1 != data_filter$organ2, ])

interorgan_filter = interorgan_filter %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())
crossorgan_filter = crossorgan_filter %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())


interorgan = as.data.frame(data[data$organ1 == data$organ2, ])
crossorgan = as.data.frame(data[data$organ1 != data$organ2, ])

interorgan = interorgan %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())
crossorgan = crossorgan %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())



interorgan_filter$percent = interorgan_filter$n/interorgan$n
crossorgan_filter$percent = crossorgan_filter$n/crossorgan$n


interorgan_filter$col = rename_based_on_df(interorgan_filter$value, color_map, from = "organ", to = "color")
crossorgan_filter$col = rename_based_on_df(crossorgan_filter$value, color_map, from = "organ", to = "color")


interorgan_filter$str_per = sprintf("%.2f%%", interorgan_filter$percent * 100)
crossorgan_filter$str_per = sprintf("%.2f%%", crossorgan_filter$percent * 100)

  
crossorgan_filter = crossorgan_filter %>% dplyr::rename(category="value", value="percent", label="str_per", color="col")
crossorgan_filter$category = factor(crossorgan_filter$category, levels = rev(color_map$organ))
interorgan_filter = interorgan_filter %>% dplyr::rename(category="value", value="percent", label="str_per", color="col")
interorgan_filter$category = factor(interorgan_filter$category, levels = rev(color_map$organ))


crossorgan_filter_sandbox = crossorgan_filter
interorgan_filter_sandbox = interorgan_filter

check_path("plot/pie/sandbox_pie_interorgan.png")
png("plot/pie/sandbox_pie_interorgan.png", width = 1600, height = 1600)
draw_pie(interorgan_filter_sandbox)
dev.off()


check_path("plot/pie/sandbox_pie_crossorgan.png")
png("plot/pie/sandbox_pie_crossorgan.png", width = 1600, height = 1600)
draw_pie(crossorgan_filter_sandbox)
dev.off()




  
# check_path("plot/pie/sandbox_pie_interorgan.png")
# png("plot/pie/sandbox_pie_interorgan.png", width = 1600, height = 1600)
# ggpubr::ggdonutchart(interorgan_filter, "percent", fill="value", font.family = "sans",palette = interorgan_filter$col, label = "str_per",
#                      lab.font = c(20, "bold", "black"),
#                      color="white")
# dev.off()
# 
# 
# check_path("plot/pie/sandbox_pie_crossorgan.png")
# png("plot/pie/sandbox_pie_crossorgan.png", width = 1600, height = 1600)
# ggpubr::ggdonutchart(crossorgan_filter, "percent", fill="value", font.family = "sans", palette = crossorgan_filter$col, label = "str_per",
#                      lab.font = c(20, "bold", "black"),
#                      color="white")
# dev.off()








##################################################


data = read.csv("data/organ_pheno_organ_pheno/ukb_r.csv")
data$p = p.adjust(data$p, method = "bonferroni")

a = sum((data$organ1 == data$organ2) & (data$p <= 0.01))/dim(data)[1]
b = sum((data$organ1 != data$organ2) & (data$p <= 0.01))/dim(data)[1]
c = sum(data$p > 0.01)/dim(data)[1]

df = data.frame(label = c("Significant intra-organ",
                          "Significant organ-to-organ",
                          "Not significantly"),
                percent = c(a, b, c))


check_path("plot/pie/ukb_pie.png")
png("plot/pie/ukb_pie.png", width = 1600, height = 1600)
ggpubr::ggdonutchart(df, "percent", fill="label", label = df$label, font.family = "sans",
                     lab.font = c(5, "bold", "black"),palette = c("grey", "#619CFF", "#F8766D"), 
                     color="white")
dev.off()



data_filter = data %>% filter(p <= 0.01)

interorgan_filter = as.data.frame(data_filter[data_filter$organ1 == data_filter$organ2, ])
crossorgan_filter = as.data.frame(data_filter[data_filter$organ1 != data_filter$organ2, ])

interorgan_filter = interorgan_filter %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())
crossorgan_filter = crossorgan_filter %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())


interorgan = as.data.frame(data[data$organ1 == data$organ2, ])
crossorgan = as.data.frame(data[data$organ1 != data$organ2, ])

interorgan = interorgan %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())
crossorgan = crossorgan %>% pivot_longer(cols = c(organ1, organ2)) %>% group_by(value) %>% summarise(n=n())



interorgan_filter$percent = interorgan_filter$n/interorgan$n
crossorgan_filter$percent = crossorgan_filter$n/crossorgan$n


interorgan_filter$col = rename_based_on_df(interorgan_filter$value, color_map, from = "organ", to = "color")
crossorgan_filter$col = rename_based_on_df(crossorgan_filter$value, color_map, from = "organ", to = "color")

interorgan_filter$str_per = sprintf("%.2f%%", interorgan_filter$percent * 100)
crossorgan_filter$str_per = sprintf("%.2f%%", crossorgan_filter$percent * 100)




crossorgan_filter = crossorgan_filter %>% dplyr::rename(category="value", value="percent", label="str_per", color="col")
crossorgan_filter$category = factor(crossorgan_filter$category, levels = rev(color_map$organ))
interorgan_filter = interorgan_filter %>% dplyr::rename(category="value", value="percent", label="str_per", color="col")
interorgan_filter$category = factor(interorgan_filter$category, levels = rev(color_map$organ))


crossorgan_filter_ukb = crossorgan_filter
interorgan_filter_ukb = interorgan_filter


check_path("plot/pie/ukb_pie_interorgan.png")
png("plot/pie/ukb_pie_interorgan.png", width = 1600, height = 1600)
draw_pie(interorgan_filter_ukb)
dev.off()


check_path("plot/pie/ukb_pie_crossorgan.png")
png("plot/pie/ukb_pie_crossorgan.png", width = 1600, height = 1600)
draw_pie(crossorgan_filter_ukb)
dev.off()


circlize::circos.heatmap()

# check_path("plot/pie/ukb_pie_interorgan.png")
# png("plot/pie/ukb_pie_interorgan.png", width = 1600, height = 1600)
# ggpubr::ggdonutchart(interorgan_filter, "percent", fill="value", font.family = "sans",palette = interorgan_filter$col, label = "str_per", 
#                      lab.font = c(20, "bold", "black"),
#                      color="white")
# dev.off()
# 
# 
# check_path("plot/pie/ukb_pie_crossorgan.png")
# png("plot/pie/ukb_pie_crossorgan.png", width = 1600, height = 1600)
# ggpubr::ggdonutchart(crossorgan_filter, "percent", fill="value", font.family = "sans", palette = crossorgan_filter$col, label = "str_per", 
#                      lab.font = c(20, "bold", "black"),
#                      color="white")
# dev.off()






col_fun1 = colorRamp2(c(-0.6, 0, 0.6), c("#EEE57E", "#EB943D", "#E60A13"))



