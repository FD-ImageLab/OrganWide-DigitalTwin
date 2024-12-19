library("networkD3")
library("htmlwidgets")
library("tidyverse")


data = read.csv("data/mediation/top50_whole.csv", row.names = 1)
feature_dict = read.csv("scripts/pheno_dict.csv")

data$y = rename_based_on_df(data$y, feature_dict, from = "names", to = "description")
data = data %>% filter(prop.mediated.p < 0.05)

top_mo = data %>% 
  group_by(mediator) %>%
  dplyr::summarise(mean=mean(abs(prop.mediated))) %>%
  ungroup() %>%
  dplyr::arrange(desc(mean))

data = data[data$mediator %in% top_mo$mediator[1:10], ]

dict = read.csv("Sandbox/ukb_category/feature.csv")
mo_dict = read.csv("Sandbox/ukb_category/molecular.csv")
data$x = rename_based_on_df(data$x, dict, from = "UDI", to = "names")
data$y = rename_based_on_df(data$y, dict, from = "UDI", to = "names")

data$mediator = rename_based_on_df(data$mediator, mo_dict, from = "UDI", to = "Short_name")

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
df$source = rename_based_on_df(df$source, node, from = "node", to = "ID")
df$target = rename_based_on_df(df$target, node, from = "node", to = "ID")

df$source = as.double(df$source)
df$target = as.double(df$target)


sankey = sankeyNetwork(df, node, Source = 'source',
              Target = 'target', Value = 'mean', NodeID = 'node', 
              units = 'TWh', fontSize = 12, nodeWidth = 30, sinksRight = FALSE, )
check_path("plot/sankey/ukb/sankey.html")
saveWidget(sankey, "plot/sankey/ukb/sankey.html")
