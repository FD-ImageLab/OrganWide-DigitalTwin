library("networkD3")
library("htmlwidgets")
library("tidyverse")
install.packages("webshot")
webshot::install_phantomjs()



source("scripts/utils/Func.R")
source("scripts/utils/read_sandbox_data.R", encoding = "UTF-8")
source("scripts/utils/pfunc.R", encoding = "UTF-8")


organs = read_data("scripts/category/subarea/", class = "SubArea", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


full_organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
full_organs = full_organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")


full_organs = subset(full_organs, ((names(full_organs) != "Prostate") & names(full_organs) != "Uterus"))


disease = read.csv("Sandbox/disease_sandbox.csv", check.names = FALSE, row.names = "buId")
health_df = disease[apply(disease, 1, function(x) all(x == 0)), ]
health_id = rownames(health_df)


features = read.csv("data/sandbox_feature.csv", row.names = "buId", check.names = FALSE)
################################ subarea #############################
whole_df = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = get_cca_result(i_part, full_organ)
        cca = results$cca
        cca.p = results$cca.p
        temp_df = data.frame(organ_subarea = i_name, 
                             subarea = i_part_name, 
                             organ2 = j_full_organ_name, 
                             r = cca, p = cca.p)
        output_df = rbind(output_df, temp_df)
      }
    }
    whole_df = rbind(whole_df, output_df)
  }
  check_path(paste0("data/circos/cca/", i_name, ".csv"))
  write.csv(output_df, paste0("data/circos/cca/", i_name, ".csv"), row.names = FALSE)
}

whole_df = whole_df %>%
  group_split(organ_subarea)

for (i in seq_along(whole_df)) {
  df = whole_df[[i]]
  organ_name = unique(df$organ_subarea)
  
  node_color_map = data.frame(key = c(color_map$organ, radiomics_color_map$radiomics), 
                              value = c(color_map$color, radiomics_color_map$color))
  draw_sankey(df, source = "subarea", target = "organ2", value = "r", 
              node_color_map = node_color_map, from = "key", to = "value",
              output_path = paste0("plot/sankey/sandbox/subarea/", organ_name, ".html"))
  webshot::webshot(paste0("plot/sankey/sandbox/subarea/", organ_name, ".html"), 
                   paste0("plot/sankey/sandbox/subarea/", organ_name, ".png"), vwidth = 3200, vheight = 800)
}


################################ subarea with disease##############################
all_data = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    i_part_name_with_organ = paste0(i_name, "_", i_part_name)
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = gene_new_pheno(i_part, full_organ, 
                                 i_part_name_with_organ, j_full_organ_name, 
                                 col_seq = "__")
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
  
  data$organ_subarea = map_chr(data$name1, ~str_split(.x, pattern = "_")[[1]][1])
  data$subarea = sub("^[^_]*_", "", data$name1)
  data$organ2 = data$name2
  data = data %>% dplyr::select(organ_subarea, subarea, organ2, value)
  
  check_path(str_glue("data/circos/subarea/disease/{id_name}.csv"))
  write.csv(data, str_glue("data/circos/subarea/disease/{id_name}.csv"))
}




################################ subarea with risk_factor##############################
all_data = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    i_part_name_with_organ = paste0(i_name, "_", i_part_name)
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = gene_new_pheno(i_part, full_organ, 
                                 i_part_name_with_organ, j_full_organ_name, 
                                 col_seq = "__")
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
}



for (i in seq_along(features)) {
  df = data.frame()
  feature = features[, i, drop=FALSE]
  feature_name = colnames(features)[i]
  data = all_data
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
  colnames(data) = str_remove(colnames(data), "y_res_")
  
  for (j in seq_along(data)) {
    new_pheno = data[, j, drop=FALSE]
    results = cor.test(feature[intersect(rownames(feature), rownames(new_pheno)), ], 
                       new_pheno[intersect(rownames(feature), rownames(new_pheno)), ])
    temp_df = data.frame(organ_subarea = str_split(colnames(new_pheno), "_")[[1]][1], 
                         subarea = sub("^[^_]*_", "", str_split(colnames(new_pheno), "__")[[1]][1]), 
                         organ2 = str_split(colnames(new_pheno), "__")[[1]][2], 
                         r = results$estimate, 
                         p = results$p.value)
    df = rbind(df, temp_df)
  }
  check_path(str_glue("data/circos/subarea/feature/{feature_name}.csv"))
  write.csv(df, str_glue("data/circos/subarea/feature/{feature_name}.csv"))
}













############################# modality ##########################
organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


full_organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))


full_organs = full_organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")


full_organs = subset(full_organs, ((names(full_organs) != "Prostate") & names(full_organs) != "Uterus"))


################################ modality #############################
whole_df = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = get_cca_result(i_part, full_organ)
        cca = results$cca
        cca.p = results$cca.p
        temp_df = data.frame(organ_modality = i_name, 
                             modality = i_part_name, 
                             organ2 = j_full_organ_name, 
                             r = cca, p = cca.p)
        output_df = rbind(output_df, temp_df)
      }
    }
    whole_df = rbind(whole_df, output_df)
  }
  check_path(paste0("data/circos/modality/cca/", i_name, ".csv"))
  write.csv(output_df, paste0("data/circos/modality/cca/", i_name, ".csv"), row.names = FALSE)
}


whole_df = whole_df %>%
  group_split(organ_modality)

for (i in seq_along(whole_df)) {
  df = whole_df[[i]]
  organ_name = unique(df$organ_modality)
  
  node_color_map = data.frame(key = c(color_map$organ, radiomics_color_map$radiomics), 
                              value = c(color_map$color, radiomics_color_map$color))
  draw_sankey(df, source = "modality", target = "organ2", value = "r", 
              node_color_map = node_color_map, from = "key", to = "value",
              output_path = paste0("plot/sankey/sandbox/modality/", organ_name, ".html"))
  webshot::webshot(paste0("plot/sankey/sandbox/modality/", organ_name, ".html"), 
                   paste0("plot/sankey/sandbox/modality/", organ_name, ".png"), vwidth = 3200, vheight = 800)
  
}


################################ modality with disease##############################
all_data = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    i_part_name_with_organ = paste0(i_name, "_", i_part_name)
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = gene_new_pheno(i_part, full_organ, 
                                 i_part_name_with_organ, j_full_organ_name, 
                                 col_seq = "__")
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
  
  data$organ_modality = map_chr(data$name1, ~str_split(.x, pattern = "_")[[1]][1])
  data$modality = sub("^[^_]*_", "", data$name1)
  data$organ2 = data$name2
  data = data %>% dplyr::select(organ_modality, modality, organ2, value)
  
  check_path(str_glue("data/circos/modality/disease/{id_name}.csv"))
  write.csv(data, str_glue("data/circos/modality/disease/{id_name}.csv"))
}


################################ modality with risk_factor##############################
all_data = data.frame()
for (i in seq_along(organs)) {
  output_df = data.frame()
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  for (k in seq_along(i_organ)) {
    i_part = i_organ[[k]]
    i_part_name = names(i_organ)[k]
    i_part_name_with_organ = paste0(i_name, "_", i_part_name)
    for (j in seq_along(full_organs)) {
      if (i != j) {
        full_organ = full_organs[[j]]
        j_full_organ_name = names(full_organs)[j]
        results = gene_new_pheno(i_part, full_organ, 
                                 i_part_name_with_organ, j_full_organ_name, 
                                 col_seq = "__")
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
}



for (i in seq_along(features)) {
  df = data.frame()
  feature = features[, i, drop=FALSE]
  feature_name = colnames(features)[i]
  data = all_data
  rownames(data) = NULL
  data = data %>% column_to_rownames("buId")
  colnames(data) = str_remove(colnames(data), "y_res_")
  
  for (j in seq_along(data)) {
    new_pheno = data[, j, drop=FALSE]
    results = cor.test(feature[intersect(rownames(feature), rownames(new_pheno)), ], 
                       new_pheno[intersect(rownames(feature), rownames(new_pheno)), ])
    temp_df = data.frame(organ_modality = str_split(colnames(new_pheno), "_")[[1]][1], 
                         modality = sub("^[^_]*_", "", str_split(colnames(new_pheno), "__")[[1]][1]), 
                         organ2 = str_split(colnames(new_pheno), "__")[[1]][2], 
                         r = results$estimate, 
                         p = results$p.value)
    df = rbind(df, temp_df)
  }
  check_path(str_glue("data/circos/modality/feature/{feature_name}.csv"))
  write.csv(df, str_glue("data/circos/modality/feature/{feature_name}.csv"))
}





# 
# 
# 
# 
# ############################ 计算减去健康人的结果 ##############################
# for (file in list.files("data/circos/modality/disease/", full.names = TRUE)) {
#   df = read.csv(file)
#   health_df = read.csv("data/circos/modality/disease/health.csv")
#   df$value = df$value - health_df$value
#   write.csv(df, file, row.names = FALSE)
# }
# 
# 
# 
# for (file in list.files("data/circos/subarea/disease/", full.names = TRUE)) {
#   df = read.csv(file)
#   health_df = read.csv("data/circos/subarea/disease/health.csv")
#   df$value = df$value - health_df$value
#   write.csv(df, file, row.names = FALSE)
# }
# 
# 
