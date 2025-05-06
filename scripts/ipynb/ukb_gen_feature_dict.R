library(tidyverse)
source("scripts/utils/Func.R")

dict_all = read_tsv("data/Data_Dictionary_Showcase.tsv")


dict = data.frame()
fid = dict_all$FieldID[((dict_all$ValueType == "Integer" | 
                           dict_all$ValueType == "Continuous" |
                           dict_all$ValueType == "Categorical single" |
                           dict_all$ValueType == "Categorical multiple") & 
                          (dict_all$ItemType == "Data"))]

names = dict_all$Field[((dict_all$ValueType == "Integer" | 
                           dict_all$ValueType == "Continuous" |
                           dict_all$ValueType == "Categorical single" |
                           dict_all$ValueType == "Categorical multiple") & 
                          (dict_all$ItemType == "Data"))]

description = dict_all$Field[((dict_all$ValueType == "Integer" | 
                                 dict_all$ValueType == "Continuous" |
                                 dict_all$ValueType == "Categorical single" |
                                 dict_all$ValueType == "Categorical multiple") & 
                                (dict_all$ItemType == "Data"))]

category = dict_all$Category[((dict_all$ValueType == "Integer" | 
                                 dict_all$ValueType == "Continuous" |
                                 dict_all$ValueType == "Categorical single" |
                                 dict_all$ValueType == "Categorical multiple") & 
                                (dict_all$ItemType == "Data"))]
dict = data.frame(UDI=NA, fid=fid, names=names, description=description, 
                  Category=category, Type="feature")


need_fids = dict$fid


udi = names(read.csv("data/ukb670788_all.csv", nrows = 1, check.names = FALSE, header = TRUE))

fids = unlist(map(str_split(udi, "-"), ~.x[1]))

dict$UDI = udi[match(need_fids, fids)]
# dict$UDI = udi[(!duplicated(fids)) & (fids %in% need_fids)]
dict = dict[!is.na(dict$UDI), ]

write.csv(dict, "Sandbox/ukb_category/feature.csv", row.names = FALSE)

# bai_dict = read.delim("data/ukb_field_added.txt", header = FALSE)
# feature_dict$category_name = map_chr(feature_dict$fid, ~ifelse(.x %in% bai_dict$V2, bai_dict[which(bai_dict$V2 %in% .x), "V3"], NA))


