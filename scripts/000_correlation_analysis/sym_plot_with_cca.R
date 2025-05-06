# Plot network figure
# 2022.05.02
rm(list=ls())

library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(scales)


source("scripts/utils/pfunc.R", encoding = "UTF-8")
source("scripts/utils/Func.R")


Left_Heart = read.csv("data/left_heart.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Right_Heart = read.csv("data/right_heart.csv", header = TRUE, row.names = "eid", check.names = FALSE)

Left_Brain = read.csv("data/left_brain.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Right_Brain = read.csv("data/right_brain.csv", header = TRUE, row.names = "eid", check.names = FALSE)

Left_Kidney = read.csv("data/left_kidney.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Right_Kidney = read.csv("data/right_kidney.csv", header = TRUE, row.names = "eid", check.names = FALSE)

Heart = read.csv("data/heart.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Brain = read.csv("data/brain.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Kidney = read.csv("data/kidney.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Liver = read.csv("data/liver.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Lung = read.csv("data/lung.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Pancreas = read.csv("data/pancreas.csv", header = TRUE, row.names = "eid", check.names = FALSE)
Spleen = read.csv("data/spleen.csv", header = TRUE, row.names = "eid", check.names = FALSE)






draw_center_network_plot = function(df) {
  data = data %>% 
    # filter(CC1.p <= p_threshold) %>%
    rename(r_mean=CC1)
  
  data$organ1 = str_replace(data$organ1, "_", "\n")
  data$organ2 = str_replace(data$organ2, "_", "\n")
  data$organ1 = str_to_title(data$organ1)
  data$organ2 = str_to_title(data$organ2)
  
  
  data$organ1 = rename_based_on_df(data$organ1, color_map, "organ", "short_organ")
  data$organ2 = rename_based_on_df(data$organ2, color_map, "organ", "short_organ")
  
  
  # build network
  exports <- data %>% distinct(organ1) %>% rename(label = organ1)
  imports <- data %>% distinct(organ2) %>% rename(label = organ2)
  nAttr <- nrow(full_join(exports, imports, by = "label"))
  
  nodes <- full_join(exports, imports, by = "label")
  nodes$label = nodes$label[order(match(nodes$label, color_map$short_organ))]
  nodes = nodes %>% 
    mutate(id = 1:nAttr) %>%
    dplyr::select(id, everything())
  
  # change node color
  nodes$label = factor(nodes$label, levels=nodes$label)
  node_actual_name = as.character(nodes$label)
  
  node_color = rename_based_on_df(node_actual_name, color_map, 
                                  from = "short_organ", to = "color")
  node_color = rename_based_on_df(node_color, lr_organ_color_map, 
                                  from = "organ", to = "color")
  
  
  edges<-data%>%
    left_join(nodes, by = c("organ1" = "label")) %>% 
    rename(from = id)
  
  edges <- edges %>% 
    left_join(nodes, by = c("organ2" = "label")) %>% 
    rename(to = id)
  
  edges <- dplyr::select(edges, from, to, r_mean)
  center_left_index=grep("Left", nodes$label) 
  center_right_index=grep("Right", nodes$label) 
  edges$r_mean[which(edges$from %in% c(center_left_index, center_right_index) | edges$to %in% c(center_left_index, center_right_index))] = 1.5 * edges$r_mean[which(edges$from %in% c(center_left_index, center_right_index) | edges$to %in% c(center_left_index, center_right_index))]
  network <- tbl_graph(
    nodes = nodes, edges = edges, directed = TRUE
  )
  
  
  x0 = seq(0,nAttr-3)
  x1=6*sin(x0*360/(nAttr-2)*pi/180)
  y1=6*cos(x0*360/(nAttr-2)*pi/180)
  x<-c(x1,-2.6,2.6)
  y<-c(y1,0,0)
  
  center_left_index=grep("Left", nodes$label) 
  center_right_index=grep("Right", nodes$label) 
  tem<-c(x[center_left_index])
  x[center_left_index]<-c(x[nAttr-1])
  x[nAttr-1]<-c(tem)
  tem<-c(x[center_right_index])
  x[center_right_index]<-c(x[nAttr])
  x[nAttr]<-c(tem)
  tem<-c(y[center_left_index])
  y[center_left_index]<-c(y[nAttr-1])
  y[nAttr-1]<-c(tem)
  tem<-c(y[center_right_index])
  y[center_right_index]<-c(y[nAttr])
  y[nAttr]<-c(tem)
  
  coords1 <- data.frame(
    x,y
  )
  
  edge_color<-c()
  edge_color[which(edges$from!=center_left_index&edges$to!=center_left_index&edges$from!=center_right_index&edges$to!=center_right_index)]<- c('#1f78b4')
  edge_color[which(edges$from==center_left_index|edges$to==center_left_index|edges$from==center_right_index|edges$to==center_right_index)]<-c('#d3311f')
  r_alpha<-c()
  r_alpha[which(edges$from!=center_left_index&edges$to!=center_left_index&edges$from!=center_right_index&edges$to!=center_right_index)]<-c(edges$r_mean[which(edges$from!=center_left_index&edges$to!=center_left_index&edges$from!=center_right_index&edges$to!=center_right_index)])
  r_alpha[which(edges$from==center_left_index|edges$to==center_left_index|edges$from==center_right_index|edges$to==center_right_index)]<-c(1)
  point_size<-c()
  for (i in 1:nAttr){
    point_size[i]=c(sum(edges$r[which(edges$from==i|edges$to==i)]))/4
  }
  # point_size = rescale(point_size, c(1, 3))
  
  ggraph(network,layout=coords1) + 
    geom_edge_fan(aes(edge_width=(r_mean),edge_alpha=r_alpha,edge_color=edge_color
    ),show.legend=F) + scale_edge_color_identity()+
    #  geom_edge_fan(aes(edge_width=(rmean),filter=to==8
    #  ),color="#990000",show.legend=F) +  
    geom_node_circle(aes(color=label,fill=label,r=point_size),
                     show.legend = F) +
    coord_fixed() +
    # scale_edge_width(range = c(0.1,4)) +
    scale_edge_alpha(range = c(0.5, 1)) +
    geom_node_text(aes(label = label),fontface="bold",repel = F,size=7, family=FONTFAMILY) +
    #scale_color_identity(str_replace_all(, "_", "\n"),aesthetics =",label")+
    scale_color_manual(values = node_color)+
    scale_fill_manual(values = node_color)+
    theme_graph()
}



################## ukb ###########################
################### left right heart #######################
organs = list(`Left_Heart` = Left_Heart, `Right_Heart` = Right_Heart, 
              Brain = Brain, Kidney = Kidney, Liver = Liver, Lung = Lung, 
              Pancreas = Pancreas, Spleen = Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/ukb_", "Heart", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/ukb_", "Heart", "Center.png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)





################### left right brain #######################
organs = list(`Left_Brain` = Left_Brain, `Right_Brain` = Right_Brain, 
              Heart = Heart, Kidney = Kidney, Liver = Liver, Lung = Lung, 
              Pancreas = Pancreas, Spleen = Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/ukb_", "Brain", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/ukb_", "Brain", "Center.png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)









################### left right kidney #######################
organs = list(`Left_Kidney` = Left_Kidney, `Right_Kidney` = Right_Kidney, 
              Heart = Heart, Brain = Brain, Liver = Liver, Lung = Lung, 
              Pancreas = Pancreas, Spleen = Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/ukb_", "Kidney", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/ukb_", "Kidney", "Center.png"), 
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)




###################### sandbox ############################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]


organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
organs = organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")

male_organs = subset(organs, names(organs) != "Uterus")
female_organs = subset(organs, names(organs) != "Prostate")
without_sex_organs = subset(organs, ((names(organs) != "Uterus") & 
                                            (names(organs) != "Prostate")))

# 左右器官做y轴
read_left_right_data = function() {
  left_right_organs = list()
  Brain = organs$Brain
  left_brain = Brain[c(grep("_l", colnames(Brain)), grep("LMCA", colnames(Brain)), 
                       grep("LIC", colnames(Brain)), grep("LPCA", colnames(Brain)), 
                       grep("LSCA", colnames(Brain)), grep("LVA", colnames(Brain)))]
  right_brain = Brain[c(grep("_r", colnames(Brain)), grep("RMCA", colnames(Brain)), 
                        grep("RIC", colnames(Brain)), grep("RPCA", colnames(Brain)), 
                        grep("RSCA", colnames(Brain)), grep("RVA", colnames(Brain)))]
  left_right_organs$Brain = list(left_brain = left_brain, right_brain = right_brain)
  
  Heart = organs$Heart
  left_heart = Heart[c(grep("LA", colnames(Heart)), grep("LV", colnames(Heart)), 
                       grep("LPA", colnames(Heart)), grep("LPV", colnames(Heart)))]
  right_heart = Heart[c(grep("RA", colnames(Heart)), grep("RV", colnames(Heart)), 
                        grep("RPA", colnames(Heart)), grep("RPV", colnames(Heart)))]
  left_right_organs$Heart = list(left_heart = left_heart, right_heart = right_heart)
  
  Kidney = organs$Kidney
  left_kidney = Kidney[c(grep("L_", colnames(Kidney)), grep("LRA", colnames(Kidney)), 
                         grep("LRV", colnames(Kidney)))]
  right_kidney = Kidney[c(grep("R_", colnames(Kidney)), grep("RRA", colnames(Kidney)), 
                          grep("RRV", colnames(Kidney)))]
  left_right_organs$Kidney = list(left_kidney = left_kidney, right_kidney = right_kidney)
  
  Lung = organs$Lung
  left_lung = Lung[grep("L_", colnames(Lung))]
  right_lung = Lung[grep("R_", colnames(Lung))]
  left_right_organs$Lung = list(left_lung = left_lung, right_lung = right_lung)
  return(left_right_organs)
}
left_right_organs = read_left_right_data()



################### left right heart #######################
organs = list(`Left_Heart` = left_right_organs$Heart$left_heart, 
              `Right_Heart` = left_right_organs$Heart$left_heart, 
              Brain = without_sex_organs$Brain, 
              Kidney = without_sex_organs$Kidney, 
              Liver = without_sex_organs$Liver, 
              Lung = without_sex_organs$Lung, 
              Pancreas = without_sex_organs$Pancreas, 
              Spleen = without_sex_organs$Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/sandbox_", "Heart", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/sandbox_", "Heart", "Center.png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)





################### left right brain #######################
organs = list(`Left_Brain` = left_right_organs$Brain$left_brain, 
              `Right_Brain` = left_right_organs$Brain$right_brain, 
              Heart = without_sex_organs$Heart, 
              Kidney = without_sex_organs$Kidney, 
              Liver = without_sex_organs$Liver, 
              Lung = without_sex_organs$Lung, 
              Pancreas = without_sex_organs$Pancreas, 
              Spleen = without_sex_organs$Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/sandbox_", "Brain", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/sandbox_", "Brain", "Center.png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)









################### left right kidney #######################
organs = list(`Left_Kidney` = left_right_organs$Kidney$left_kidney, 
              `Right_Kidney` = left_right_organs$Kidney$right_kidney, 
              Heart = without_sex_organs$Heart, 
              Brain = without_sex_organs$Brain, 
              Liver = without_sex_organs$Liver, 
              Lung = without_sex_organs$Lung, 
              Pancreas = without_sex_organs$Pancreas, 
              Spleen = without_sex_organs$Spleen)


output_cancor_coef_df = data.frame()
for (i in seq_along(organs)) {
  for (j in seq_along(organs)) {
    if (i < j) {
      organ1 = organs[[i]]
      organ2 = organs[[j]]
      organ1_name = names(organs)[i]
      organ2_name = names(organs)[j]
      
      temp = get_cca_result(organ1, organ2)
      cor_main_coef = temp$cca
      cor_main_coef.p = temp$cca.p
      N = nrow(temp$X_prime[[1]])
      temp_df = tibble(organ1=organ1_name, organ2=organ2_name,
                       CC1=cor_main_coef, CC1.p=cor_main_coef.p, N=N)
      output_cancor_coef_df <- bind_rows(output_cancor_coef_df, temp_df)
    }
  }
}

data = output_cancor_coef_df
draw_center_network_plot(data)
check_path(paste0("plot/sym_network/sandbox_", "Kidney", "Center.png"))
ggsave(
  filename = paste0("plot/sym_network/sandbox_", "Kidney", "Center.png"), 
  width = 7,
  height = 7,
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)