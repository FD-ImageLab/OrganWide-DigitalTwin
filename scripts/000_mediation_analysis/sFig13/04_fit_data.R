library("tidyverse")

data = read.csv(paste0("data/sfig13/organ_cell_organ_mean.csv"), 
                header = TRUE)
data1 = data[data$y_fix == "", ]
data1 = rename(data1, organ1=x_fix, organ2=cell_category, r_mean=mean_prop.mediated)
data1$y_fix = NULL
data2 = data[data$x_fix == "", ]
data2 = rename(data2, organ1=cell_category, organ2=y_fix, r_mean=mean_prop.mediated)
data2$x_fix = NULL
data2$organ2 = paste0(data2$organ2,"_layer3")

data = bind_rows(data1, data2)
write.csv(data, paste0("data/sfig13/organ_cell_organ_mediation_data.csv"), row.names = FALSE)







library("tidyverse")

files = list.files("data/sfig13/feature/", full.names = TRUE, include.dirs = FALSE, 
                   recursive = TRUE, pattern = "organ_cell_organ_mean.csv")
for (i in seq_along(files)) {
  feature = basename(dirname(files[i]))
  data = read.csv(files[i], 
                  header = TRUE)
  data1 = data[data$y_fix == "", ]
  data1 = rename(data1, organ1=x_fix, organ2=cell_category, r_mean=mean_prop.mediated)
  data1$y_fix = NULL
  data2 = data[data$x_fix == "", ]
  data2 = rename(data2, organ1=cell_category, organ2=y_fix, r_mean=mean_prop.mediated)
  data2$x_fix = NULL
  data2$organ2 = paste0(data2$organ2,"_layer3")
  
  data = bind_rows(data1, data2)
  write.csv(data, paste0("data/sfig13/feature/", feature, "/organ_cell_organ_mediation_data.csv"), row.names = FALSE)
  
  
}




