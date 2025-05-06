library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(scales)

source('scripts/utils/Func.R')
Paths = list.files(paste0("data/5"), full.names = TRUE)

for (Path in Paths) {
  data = read.csv(Path, header = T)
  feature = unique(data$feature)
  # process data
  data = data %>%
    filter(p_bonf <= 0.05)
  if (nrow(data) > 0) {
    data$organ1 = rename_based_on_df(data$organ1, color_map, "organ", "short_organ")
    data$organ2 = rename_based_on_df(data$organ2, color_map, "organ", "short_organ")
    
    # build network
    exports <- data %>% distinct(organ1) %>% rename(label = organ1)
    imports <- data %>% distinct(organ2) %>% rename(label = organ2)
    nAttr <- nrow(full_join(exports, imports, by = "label"))
    
    nodes <- full_join(exports, imports, by = "label")
    nodes$label = nodes$label[order(match(nodes$label, color_map$organ))]
    nodes = nodes %>% 
      mutate(id = 1:nAttr) %>%
      dplyr::select(id, everything())
    
    # change node color
    nodes$label = factor(nodes$label, levels=nodes$label)
    node_actual_name = as.character(nodes$label)
    node_color = rename_based_on_df(node_actual_name, color_map, 
                                    from = "short_organ", to = "color")
    
    
    edges<-data%>%
      left_join(nodes, by = c("organ1" = "label")) %>% 
      rename(from = id)
    
    edges <- edges %>% 
      left_join(nodes, by = c("organ2" = "label")) %>% 
      rename(to = id)
    
    edges <- dplyr::select(edges, from, to, r)
    
    network <- tbl_graph(
      nodes = nodes, edges = edges, directed = TRUE
    )
    
    x0 = seq(0,nAttr-1)
    coords1 <- data.frame(
      x<-2*sin(x0*360/nAttr*pi/180),
      y<-2*cos(x0*360/nAttr*pi/180)
    )
    
    r_alpha<-c(edges$r)
    point_size<-c()
    for (i in 1:nAttr){
      point_size[i]<-(c(sum(edges$r[which(edges$from==i|edges$to==i)])) + 2) * 0.03
    }
    # point_size = rescale(point_size, c(0.45, 0.7))
    
    ggraph(network,layout=coords1) + 
      geom_edge_fan(aes(edge_width=r,edge_color="#334A8D"
      ),show.legend=F) + scale_edge_color_identity() +
      #  geom_edge_fan(aes(edge_width=(rmean),filter=to==8
      #  ),color="#990000",show.legend=F) +   
      geom_node_circle(aes(color=label,fill=label, r=point_size),
                       show.legend = F) +
      coord_fixed() +
      scale_edge_width(range = c(0.1,5)) +
      scale_edge_alpha(range = c(0.1, 1)) +
      geom_node_text(aes(label = label),fontface="bold",repel = F,size=15, family=FONTFAMILY) +
      #scale_color_identity(str_replace_all(, "_", "\n"),aesthetics ="label")+
      scale_color_manual(values = node_color)+
      scale_fill_manual(values = node_color)+
      theme_graph()
    
    check_path(paste0("plot/", str_remove(basename(Path), "_data.csv"), ".png"))
    ggsave(
      filename = paste0("plot/", str_remove(basename(Path), "_data.csv"), ".png"), # 保存的文件名称。通过后缀来决定生成什么格式的图片
      width = 7,             
      height = 7,            
      units = "in",          
      dpi = 300              
    )
  }
}
