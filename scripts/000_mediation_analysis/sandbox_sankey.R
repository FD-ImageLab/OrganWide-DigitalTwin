library("networkD3")
library("htmlwidgets")
library("tidyverse")


data = read.csv("data/mediation/sandbox_med_top50_sig.csv")
data = data %>% filter(prop.mediated.p < 0.05)
data$y = map_chr(data$y, ~str_remove(.x, pattern = "y_res_"))

feature_dict = read.csv("Sandbox/Category/feature.csv")
data = data[!(data$x %in% feature_dict$Short_name[feature_dict$Category == "Physical_Measurements"]), ]

top_mo = data %>% 
  group_by(mediator) %>%
  dplyr::summarise(mean=mean(abs(prop.mediated))) %>%
  ungroup() %>%
  dplyr::arrange(desc(mean))
  
data = data[data$mediator %in% top_mo$mediator[1:10], ]



xy = data %>%
  group_by(x, mediator) %>%
  dplyr::summarise(mean = mean(abs(prop.mediated))) %>%
  ungroup() %>%
  dplyr::select(x, mediator, mean)
xy = dplyr::rename(xy, source=x, target=mediator)

yz = data %>%
  group_by(y, mediator) %>%
  dplyr::summarise(mean = mean(abs(prop.mediated))) %>%
  ungroup() %>%
  dplyr::select(mediator, y, mean)
yz = dplyr::rename(yz, source=mediator, target=y)

node = c(unique(xy$source), unique(xy$target), unique(yz$target))
node = data.frame(node=node, ID=c(0:(length(node)-1)))

df = rbind(xy, yz)
df = as.data.frame(df)



draw_sankey(df, source = "source", target = "target", value = "mean", 
            node_color_map = node_color_map, from = "key", to = "value",
            output_path = paste0("plot/sankey/sandbox/subarea/", organ_name, ".html"))
webshot::webshot(paste0("plot/sankey/sandbox/subarea/", organ_name, ".html"), 
                 paste0("plot/sankey/sandbox/subarea/", organ_name, ".png"), vwidth = 3200, vheight = 800)




# 
# 
# df$source = rename_based_on_df(df$source, node, from = "node", to = "ID")
# df$target = rename_based_on_df(df$target, node, from = "node", to = "ID")
# 
# df$source = as.double(df$source)
# df$target = as.double(df$target)
# 
# 
# sankey = sankeyNetwork(df, node, Source = 'source',
#               Target = 'target', Value = 'mean', NodeID = 'node', 
#               units = 'TWh', fontSize = 12, nodeWidth = 30, sinksRight = FALSE)
# saveWidget(sankey, "plot/sankey/sandbox/sankey.html")
