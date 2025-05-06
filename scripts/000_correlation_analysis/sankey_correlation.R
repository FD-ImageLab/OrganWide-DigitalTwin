library("networkD3")
library("htmlwidgets")
library("tidyverse")


source("scripts/utils/Func.R")
source("scripts/utils/read_sandbox_data.R", encoding = "UTF-8")
source("scripts/utils/pfunc.R", encoding = "UTF-8")


organs = read_data("scripts/category/modality/", class = "Modality", part_name="Prefix")
organs = organs[order(match(names(organs), color_map$organ))]
organs = subset(organs, ((names(organs) != "Uterus") & (names(organs) != "Prostate")))


full_organs = organs %>%
  map(~map(.x, rownames_to_column, var="buId"))

join_by_buId = function(x, y) {inner_join(x, y, by="buId")}
reduce_join = function(l) {l %>% reduce(join_by_buId)}
full_organs = full_organs %>%
  map(reduce_join) %>%
  map(column_to_rownames, var="buId")


# male_organs = subset(full_organs, names(full_organs) != "Uterus")
# female_organs = subset(full_organs, names(full_organs) != "Prostate")
full_organs = subset(full_organs, ((names(full_organs) != "Prostate") & names(full_organs) != "Uterus"))



prefix = rep(names(organs), as.numeric(map(organs, length)))
organs = flatten(organs)
names(organs) = paste0(prefix, "_", names(organs))


category_organ = names(organs)
category_organ[grep("Morphology", names(organs))] = "Morphology"
category_organ[grep("Texture", names(organs))] = "Texture"
category_organ[grep("Density", names(organs))] = "Water"
category_organ[grep("Water", names(organs))] = "Water"
category_organ[grep("Fat", names(organs))] = "Fat"
category_organ[grep("Iron", names(organs))] = "Iron"
category_organ[grep("Diffusion", names(organs))] = "Diffusion"
category_organ[grep("Function", names(organs))] = "Function"
category_organ[grep("Vessel", names(organs))] = "Vessel"
category_organ[grep("ECG", names(organs))] = "ECG"
category_organ[grep("EEG", names(organs))] = "EEG"

output_df = data.frame()
for (i in seq_along(organs)) {
  i_organ = organs[[i]]
  i_name = names(organs)[i]
  i_organ_name = str_split(names(organs)[i], "_")[[1]][1]
  i_category_name = str_split(names(organs)[i], "_")[[1]][2]
  for (j in seq_along(full_organs)) {
    full_organ = full_organs[[j]]
    j_full_organ_name = names(full_organs)[j]
    results = get_cca_result(i_organ, full_organ)
    cca = results$cca
    cca.p = results$cca.p
    temp_df = data.frame(name1 = i_name,
                         organ1 = i_organ_name,
                         category1 = i_category_name, 
                         name2 = j_full_organ_name,
                         r = cca, p = cca.p)
    output_df = rbind(output_df, temp_df)
  }
}

output_df = output_df[!(output_df$organ1 == output_df$name2), ]
whole_df = output_df %>%
  group_split(organ1)

for (i in seq_along(whole_df)) {
  df = whole_df[[i]]
  organ_name = unique(df$organ1)
  
  node_color_map = data.frame(key = c(color_map$organ, radiomics_color_map$radiomics), 
                              value = c(color_map$color, radiomics_color_map$color))
  draw_sankey(df, source = "category1", target = "name2", value = "r", 
              node_color_map = node_color_map, from = "key", to = "value",
              output_path = paste0("plot/sankey/sandbox/", organ_name, ".html"))
}


####################### pie plot #######################
# 设置用于并行计算的核心数
library("parallel")
no_cores <- 64

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



check_path("data/organ_pheno_organ_pheno_r.csv")
write.csv(df, "data/organ_pheno_organ_pheno_r.csv", row.names = FALSE)
