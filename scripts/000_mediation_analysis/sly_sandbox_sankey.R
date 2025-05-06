color = read_excel("data/color_map.xlsx")

data = read.csv("data/mediation/sandbox_med_top50_sig_with_serum.csv")
data = data[data$cell_category == "Serum", ]
data = data %>% filter(prop.mediated.p < 0.05)
data$y = map_chr(data$y, ~str_remove(.x, pattern = "y_res_"))

feature_dict = read.csv("Sandbox/Category/feature_sly.csv")
data = data[!(data$x %in% feature_dict$Short_name[feature_dict$Category == "Physical_Measurements"]), ]

top_mo = data %>% 
  group_by(mediator) %>%
  dplyr::summarise(mean=mean(abs(prop.mediated))) %>%
  ungroup() %>%
  dplyr::arrange(desc(mean))

data = data[data$mediator %in% top_mo$mediator[1:10], ]

data  = data[data$x %in% feature_dict$Short_name[feature_dict$Select == "TRUE"], ]
data$x = rename_based_on_df(data$x, feature_dict, from = "Short_name", to = "DispName")

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


draw_sankey(df, source = "source",target = "target",value = "mean",
            node_color_map = color_map, from = "organ", to = "color", 
            output_path = "plot/sankey/sandbox/sankey.html") # 改颜色 参考变量color_map
