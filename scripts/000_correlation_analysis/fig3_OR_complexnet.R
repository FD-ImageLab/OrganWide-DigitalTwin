source("scripts/utils/pfunc.R", encoding = "UTF-8")
library(tidyverse)

disease_odd = read.csv("data/disease_odd_before.csv", header = TRUE)
disease_odd["organ1"] = map_chr(str_split(disease_odd$pheno, pattern = "_"), ~.x[3])
disease_odd["organ2"] = map_chr(str_split(disease_odd$pheno, pattern = "_"), ~.x[4])
diseases_odd = disease_odd %>% group_split(disease)


draw_complexnet = function(data, output_name = "plot/complexnet/ukb_organ_axis_cca_complexnet.png", 
                                 limits = c(0, 4)) {

  data$organ1 = rename_based_on_df(data$organ1, color_map, "organ", "short_organ")
  data$organ2 = rename_based_on_df(data$organ2, color_map, "organ", "short_organ")
  
  edge_color = ifelse(data$sign == "neg", yes = "#D68485", no = "#D68485")
  edge_color = ifelse(data$p < 0.01, yes=edge_color, no="grey")
  
  # build network
  exports = data %>% distinct(organ1) %>% rename(label = organ1)
  imports = data %>% distinct(organ2) %>% rename(label = organ2)
  nAttr = nrow(full_join(exports, imports, by = "label"))
  
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
  
  
  edges = data%>%
    left_join(nodes, by = c("organ1" = "label")) %>% 
    rename(from = id)
  
  edges = edges %>% 
    left_join(nodes, by = c("organ2" = "label")) %>% 
    rename(to = id)
  
  edges = dplyr::select(edges, from, to, r)
  
  network = tbl_graph(
    nodes = nodes, edges = edges, directed = TRUE
  )
  
  x0 = seq(0,nAttr-1)
  coords1 = data.frame(
    x=2*sin(x0*360/nAttr*pi/180),
    y=2*cos(x0*360/nAttr*pi/180)
  )
  
  r_alpha=c(edges$r)
  point_size=c()
  for (i in 1:nAttr){
    point_size[i]=(c(mean(edges$r[which(edges$from==i|edges$to==i)])) + 0.7) / 2
  }
  print(point_size)
  g = ggraph(network,layout=coords1) + 
    geom_edge_fan(aes(edge_width=(r),edge_color=edge_color, edge_alpha=r_alpha), 
                  show.legend=F) + scale_edge_color_identity() +
    geom_node_circle(aes(color=label,fill=label,r=point_size),
                     show.legend = F) +
    coord_fixed() +
    scale_edge_width(range = c(1,10), limits = limits) +
    
    scale_edge_alpha(range = c(0.1, 1)) +
    geom_node_text(aes(label = label),fontface="bold",repel = F,size=15, family=FONTFAMILY) +
    #scale_color_identity(str_replace_all(, "_", "\n"),aesthetics ="label")+
    scale_color_manual(values = node_color)+
    scale_fill_manual(values = node_color)+
    theme_graph()
  print(g)
  check_path(paste0(output_name))
  ggsave(plot = g, 
         filename = paste0(output_name), # 保存的文件名称。通过后缀来决定生成什么格式的图片
         width = 7,             # 宽
         height = 7,            # 高
         units = "in",          # 单位
         dpi = 300              # 分辨率DPI
  )
}





for (i in seq_along(diseases_odd)) {
  disease_odd = diseases_odd[[i]]
  disease_name = unique(disease_odd$disease)
  disease_odd["distance2one"] = abs(disease_odd$odds - 1)
  disease_odd["sign"] = ifelse(((disease_odd$odds -1) >= 0), "pos", "neg")
  data = disease_odd %>% dplyr::select(organ1, organ2, distance2one, p, sign)
  data = dplyr::rename(data, r=distance2one)
  draw_complexnet(data, output_name = str_glue("plot/fig2/OR_005/{disease_name}.png"), limits = NULL)
}
