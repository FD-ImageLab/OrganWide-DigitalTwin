# Plot network figure
# 2022.05.02
rm(list=ls())

# install.packages("ggraph")
# install.packages("tidyverse")
# install.packages("tidygraph")
# install.packages("RColorBrewer")
# install.packages("extrafont")

library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(scales)

# set path and load data
source("scripts/utils/Func.R")

Path = paste0("data/sfig13/organ_cell_organ_mediation_data.csv")
data<-read.csv(Path,header = T, fileEncoding = "UTF-8-BOM")

exports <- data %>% distinct(organ1) %>% rename(label = organ1)
imports <- data %>% distinct(organ2) %>% rename(label = organ2)

n1_label = exports$label[exports$label %in% color_map$organ]
n1_label = n1_label[order(match(n1_label, color_map$organ))]

n2_label = exports$label[!(exports$label %in% color_map$organ)]
#may cause some problems bacause the order of n2_label
n2_label = n2_label[order(match(n2_label, micro_cell_color_map$type))]


n3_label = setdiff(imports$label, exports$label)
n3_label = str_remove(n3_label, "_layer3")
n3_label = n3_label[order(match(n3_label, color_map$organ))]
n3_label = paste0(n3_label, "_layer3")

nAttr <- length(n1_label) + length(n2_label) + length(n3_label)

#################################
# 计算圈大小
temp1 = data %>% filter(organ1 %in% n1_label) %>%
  group_by(organ1) %>%
  summarise(node_weight=sum(r_mean)) %>%
  rename(organ=organ1)
temp1 = temp1[order(match(temp1$organ, color_map$organ)), ]

temp_data = data %>% rename(organ1=organ2, organ2=organ1)
temp2 = bind_rows(data, temp_data) %>% filter(organ1 %in% c("DCs", "Mono", "NK", 
                                                            "Treg", "T-cell", "Th", "B-cell")) %>% 
  group_by(organ1) %>% 
  summarise(node_weight=sum(r_mean)/2) %>%
  rename(organ=organ1)
temp2 = temp2[order(match(temp2$organ, micro_cell_color_map$type)), ]

temp3 = data %>% filter(organ2 %in% n3_label) %>%
  group_by(organ2) %>%
  summarise(node_weight=sum(r_mean)) %>%
  rename(organ=organ2)
temp3$organ = str_remove(temp3$organ, "_layer3")
temp3 = temp3[order(match(temp3$organ, color_map$organ)), ]
temp3$organ = paste0(temp3$organ, "_layer3")


node_weight = bind_rows(temp1, temp2, temp3)

node_weight = c(temp1$node_weight, 
                temp2$node_weight, 
                temp3$node_weight)
node_weight = rescale(node_weight, c(0.5, 1.7))
#################################

nodes <- data.frame(id=1:nAttr, label=c(n1_label, n2_label, n3_label))
nodes$label = factor(nodes$label, levels=nodes$label)
node_actual_name = str_remove(nodes$label, pattern = "_layer3")

node_actual_name = rename_based_on_df(node_actual_name, color_map, "organ", "short_organ")
node_color = rename_based_on_df(node_actual_name, color_map, 
                                from = "short_organ", to = "color")
node_color = rename_based_on_df(node_color, micro_cell_color_map, 
                                from = "type", to = "color")

edges<- data %>%
  left_join(nodes, by = c("organ1" = "label")) %>% 
  rename(from = id)

edges <- edges %>% 
  left_join(nodes, by = c("organ2" = "label")) %>% 
  rename(to = id)

edges <- dplyr::select(edges, from, to, r_mean)

network <- tbl_graph(
  nodes = nodes, edges = edges, directed = TRUE
)


n1 = length(n1_label)
n2 = length(n2_label)
n3 = length(n3_label)
# 调整文字位置
str_slip = c(rep(max(c(temp1$node_weight, temp3$node_weight)) + 0.7, n1), 
             rep(0, n2), 
             rep(-(max(c(temp1$node_weight, temp3$node_weight)) + 0.7), n3))
x_wide=16
x2_wide=18
x1 = seq(0, x_wide, length.out=n1)
x2 = seq(-2, x2_wide, length.out=n2 + 2)
x3 = seq(0, x_wide, length.out=n3)

y0 <- seq(0,nAttr)*0
y1 <- y0
y2 <- y0+5
y3 <- y0+10
x <- c(x1, x2[2:(length(x2)-1)], x3)
y <- c(y3[0:n1],y2[0:n2],y1[0:n3])
coords1 <- data.frame(x,y)


ggraph(network,layout=coords1) + 
  geom_edge_fan(aes(edge_width=(r_mean), edge_alpha=(r_mean)
  ),color="#1f78b4",show.legend=F) + 
  #  geom_edge_fan(aes(edge_width=(rmean),filter=to==8
  #  ),color="#990000",show.legend=F) +  
  geom_node_circle(aes(color=label,fill=label,r=node_weight),
                   show.legend = F) +
  coord_fixed() +
  scale_edge_width(range = c(0, 2)) +
  geom_node_text(aes(label = node_actual_name),fontface="bold",
                 family=FONTFAMILY, repel = F,size=8) + # , nudge_y = str_slip
  #scale_color_identity(str_replace_all(, "_", "\n"),aesthetics ="label")+
  scale_color_manual(values = node_color)+
  scale_fill_manual(values = node_color)+
  theme_graph()

ggsave(
  filename = paste0("plot/sfig13/ConnectionGraph.png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
  width = 8,             # ???
  height = 8,            # ???
  units = "in",          # 单位
  dpi = 600              # 分辨率DPI
)









############################## feature #############################

files = list.files("data/sfig13/feature/", full.names = TRUE, include.dirs = FALSE, 
                   recursive = TRUE, pattern = "organ_cell_organ_mediation_data.csv")
for (i in seq_along(files)) {
  feature = basename(dirname(files[i]))
  
  data<-read.csv(files[i],header = T, fileEncoding = "UTF-8-BOM")
  
  exports <- data %>% distinct(organ1) %>% rename(label = organ1)
  imports <- data %>% distinct(organ2) %>% rename(label = organ2)
  
  n1_label = exports$label[exports$label %in% color_map$organ]
  n1_label = n1_label[order(match(n1_label, color_map$organ))]
  
  n2_label = exports$label[!(exports$label %in% color_map$organ)]
  #may cause some problems bacause the order of n2_label
  n2_label = n2_label[order(match(n2_label, micro_cell_color_map$type))]
  
  
  n3_label = setdiff(imports$label, exports$label)
  n3_label = str_remove(n3_label, "_layer3")
  n3_label = n3_label[order(match(n3_label, color_map$organ))]
  n3_label = paste0(n3_label, "_layer3")
  
  nAttr <- length(n1_label) + length(n2_label) + length(n3_label)
  
  #################################
  # 计算圈大小
  temp1 = data %>% filter(organ1 %in% n1_label) %>%
    group_by(organ1) %>%
    summarise(node_weight=sum(r_mean)) %>%
    rename(organ=organ1)
  temp1 = temp1[order(match(temp1$organ, color_map$organ)), ]
  
  temp_data = data %>% rename(organ1=organ2, organ2=organ1)
  temp2 = bind_rows(data, temp_data) %>% filter(organ1 %in% c("DCs", "Mono", "NK", 
                                                              "Treg", "T-cell", "Th", "B-cell")) %>% 
    group_by(organ1) %>% 
    summarise(node_weight=sum(r_mean)/2) %>%
    rename(organ=organ1)
  temp2 = temp2[order(match(temp2$organ, micro_cell_color_map$type)), ]
  
  temp3 = data %>% filter(organ2 %in% n3_label) %>%
    group_by(organ2) %>%
    summarise(node_weight=sum(r_mean)) %>%
    rename(organ=organ2)
  temp3$organ = str_remove(temp3$organ, "_layer3")
  temp3 = temp3[order(match(temp3$organ, color_map$organ)), ]
  temp3$organ = paste0(temp3$organ, "_layer3")
  
  
  node_weight = bind_rows(temp1, temp2, temp3)
  
  node_weight = c(temp1$node_weight, 
                  temp2$node_weight, 
                  temp3$node_weight)
  node_weight = rescale(node_weight, c(0.5, 1.7))
  #################################
  
  nodes <- data.frame(id=1:nAttr, label=c(n1_label, n2_label, n3_label))
  nodes$label = factor(nodes$label, levels=nodes$label)
  node_actual_name = str_remove(nodes$label, pattern = "_layer3")
  
  node_actual_name = rename_based_on_df(node_actual_name, color_map, "organ", "short_organ")
  node_color = rename_based_on_df(node_actual_name, color_map, 
                                  from = "short_organ", to = "color")
  node_color = rename_based_on_df(node_color, micro_cell_color_map, 
                                  from = "type", to = "color")
  
  edges<- data %>%
    left_join(nodes, by = c("organ1" = "label")) %>% 
    rename(from = id)
  
  edges <- edges %>% 
    left_join(nodes, by = c("organ2" = "label")) %>% 
    rename(to = id)
  
  edges <- dplyr::select(edges, from, to, r_mean)
  
  network <- tbl_graph(
    nodes = nodes, edges = edges, directed = TRUE
  )
  
  
  n1 = length(n1_label)
  n2 = length(n2_label)
  n3 = length(n3_label)
  # 调整文字位置
  str_slip = c(rep(max(c(temp1$node_weight, temp3$node_weight)) + 0.7, n1), 
               rep(0, n2), 
               rep(-(max(c(temp1$node_weight, temp3$node_weight)) + 0.7), n3))
  x_wide=16
  x2_wide=18
  x1 = seq(0, x_wide, length.out=n1)
  x2 = seq(-2, x2_wide, length.out=n2 + 2)
  x3 = seq(0, x_wide, length.out=n3)
  
  y0 <- seq(0,nAttr)*0
  y1 <- y0
  y2 <- y0+5
  y3 <- y0+10
  x <- c(x1, x2[2:(length(x2)-1)], x3)
  y <- c(y3[0:n1],y2[0:n2],y1[0:n3])
  coords1 <- data.frame(x,y)
  
  
  ggraph(network,layout=coords1) + 
    geom_edge_fan(aes(edge_width=(r_mean), edge_alpha=(r_mean)
    ),color="#1f78b4",show.legend=F) + 
    #  geom_edge_fan(aes(edge_width=(rmean),filter=to==8
    #  ),color="#990000",show.legend=F) +  
    geom_node_circle(aes(color=label,fill=label,r=node_weight),
                     show.legend = F) +
    coord_fixed() +
    scale_edge_width(range = c(0, 2)) +
    geom_node_text(aes(label = node_actual_name),fontface="bold",
                   family=FONTFAMILY, repel = F,size=8) + # , nudge_y = str_slip
    #scale_color_identity(str_replace_all(, "_", "\n"),aesthetics ="label")+
    scale_color_manual(values = node_color)+
    scale_fill_manual(values = node_color)+
    theme_graph()
  
  check_path(paste0("plot/sfig13/feature/", feature, "/ConnectionGraph.png"))
  ggsave(
    filename = paste0("plot/sfig13/feature/", feature, "/ConnectionGraph.png"), 
    width = 8,            
    height = 8,           
    units = "in",         
    dpi = 600             
  )
}

