library(tidyverse)
source("scripts/utils/Func.R")


pheno_dict = read.csv("scripts/pheno_dict2.csv")
feature_dict = read.csv("Sandbox/Category/feature_sly.csv")
molecular_dict = read.csv("Sandbox/Category/molecular.csv")
data = read.csv("data/mediation/sandbox_med_top50_sig.csv")
data = data[data$cell_category == "IF&MF", ]

data = data[(data$x %in% feature_dict$Short_name[feature_dict$Select == "TRUE"]), ]
data$x = rename_based_on_df(data$x, nmapdf = feature_dict, from = "Short_name", to = "DispName")

data$y = rename_based_on_df(data$y, nmapdf = pheno_dict, from = "names", to = "description")
data$cell_category = rename_based_on_df(data$mediator, nmapdf = molecular_dict, from = "micro", to = "category2")
data = data %>%
  separate(y, into=c("organ1", "organ2"), sep = "_")
data = data %>%
  pivot_longer(c("organ1", "organ2"), values_to = "organ")


f = function(df) {
  df %>%
    mutate(prop.mediated=abs(prop.mediated), significant = ifelse(prop.mediated.p <= 0.05, TRUE, FALSE)) %>%
    arrange(desc(significant), desc(prop.mediated)) %>%
    dplyr::distinct(mediator, .keep_all = TRUE) %>%
    slice(1:5) %>%
    arrange(desc(prop.mediated)) %>%
    select(organ, mediator, prop.mediated, cell_category, significant) %>%
    rename(value=prop.mediated, individual=mediator, group=cell_category) %>%
    rownames_to_column(var = "id") %>%
    right_join(data.frame(id=as.character(c(1:5))))
}
med_barplot_data = data %>%
  group_split(organ) %>%
  map(f)
func = function(df) {
  df$value[is.na(df$value)] = 0
  df$individual[is.na(df$individual)] = " "
  return(df)
}
med_barplot_data = map(med_barplot_data, func)


tf = function(df) {
  check_path(paste0("data/sfig13/100/", "/Organ/", unique(na.omit(df$organ)), ".csv"))
  write.csv(df, paste0("data/sfig13/100", "/Organ/", unique(na.omit(df$organ)), ".csv"), 
            fileEncoding = "UTF-8")
}
walk(med_barplot_data, tf)
