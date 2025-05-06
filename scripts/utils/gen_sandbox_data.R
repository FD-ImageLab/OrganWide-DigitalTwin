# user-defined path
mainPath <- "E:/seafile/Seafile/SandBox/"  # local
# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
setwd(paste0(mainPath, "/Sandbox/"))
library(tidyverse)

positions = list.files("Organ/", full.names = TRUE, pattern = "csv")
f = function(position) {tryCatch({
        temp = read.csv(position, row.names=1, header=TRUE)
        return(rownames(temp))}, 
        error = function(e) {})
}
a = map(positions, f)
b = a %>%
        map(is.null) %>%
        as.logical()
a[!b] %>% reduce(union) %>% sample(450) -> buId


temp = readxl::read_xlsx(paste0(mainPath, "/Sandbox/Demography/scanlist_20220713.xlsx"), sheet = "张江扫描名单")
temp = temp[!duplicated(temp$编号),]
temp = temp[!is.na(temp$编号),]
temp_df = temp %>% 
        dplyr::select(编号, 性别, 年龄) %>%
        rename(buId=编号)
write.csv(temp_df, 'Demography/基本信息.csv', row.names = FALSE)

buId = temp$编号[1: 450]
sex = temp$性别[1: 450]
age = as.numeric(temp$年龄[1: 450])




temp = readxl::read_xlsx(paste0(mainPath, "/Sandbox/Demography/zhangjiang.xlsx"), sheet = "Sheet1")
temp$`Height（m）` = as.numeric(temp$`Height（m）`)
temp$`Weight（kg）` = as.numeric(temp$`Weight（kg）`)
a = temp %>%
        mutate(BMI = `Weight（kg）`/(`Height（m）`*`Height（m）`)) %>%
        dplyr::select(c("编号", "BMI")) %>%
        rename(buId=编号)
dxa = temp_df %>%
        left_join(a, by = "buId") %>%
        dplyr::select('buId', 'BMI')
# BMI = runif(nrow(temp_df), 15, 30)
# temp_df = data.frame(temp_df$buId, BMI)
write.csv(dxa, "Demography/DXA体成分-体成分高分辨率数字成像化系统.csv", 
            row.names=FALSE)


temp_df = data.frame(buId = temp_df$buId, matrix(runif(nrow(temp_df) * 20, 10, 50), nrow=nrow(temp_df), ncol=20))
colnames(temp_df) = c("buId", "血压(mm/Hg)收缩压（mm/Hg）", "血压(mm/Hg)舒张压（mm/Hg）", "葡萄糖（mmol/L）", 
                      "总胆固醇（mmol/L）", "高密度脂蛋白胆固醇（mmol/L）", "低密度脂蛋白胆固醇（mmol/L）", 
                      "高敏感C反应蛋白（mg/L）", "淋巴细胞数（*10^9/L）", 
                      "总胆红素（μmol/L）", "尿酸（μmol/L）", "载脂蛋白A-I（g/L）", 
                      "载脂蛋白B（g/L）", "载脂蛋白E（mg/L）", "肌酐（μmol/L）", 
                      "癌胚抗原（ng/ml）", "胰岛素空腹（uU/ml）", "甲胎蛋白（ng/ml）", 
                      "丙氨酸氨基转移酶（U/L）", "门冬氨酸氨基转移酶（U/L）", "心率")
write.csv(temp_df, "Serum/志愿者招募-常规体检.csv", row.names=FALSE)


temp_df = data.frame(buId = temp_df$buId, matrix(runif(nrow(temp_df) * 19, 10, 50), nrow=nrow(temp_df), ncol=1))
colnames(temp_df) = c("buId", "中间的脖子围度（cm）")
write.csv(temp_df, "Demography/人体外观测量-三维人体扫描分析系统.csv", row.names = FALSE)


temp_df = data.frame(buId = buId, matrix(runif(2700, 10, 50), nrow=450, ncol=6))
colnames(temp_df) = c("buId", "左颈总动脉", "右颈总动脉", "左颈外动脉", "右颈外动脉", 
                      "左颈内动脉", "右颈内动脉")
write.csv(temp_df, "Vessel/超声-多普勒超声诊断仪.csv", row.names=FALSE)


temp_df = data.frame(ID = buId, matrix(runif(3600, 10, 50), nrow=450, ncol=8))
colnames(temp_df) = c("ID", "Active_GLP-1", "C-peptide", "FGF-21", "FGF-23", 
                      "Glucagon", "Insulin", "Leptin", "human_GIP-active")
write.csv(temp_df, "MicroPhenomics/metabolic_factor/Metabolic_factors.csv", 
          row.names = FALSE)


temp_df = data.frame(matrix(runif(4500, 10, 50), nrow=450, ncol=10), 
                     row.names = buId)
colnames(temp_df) = c("IFN-γ", "IL-10", "IL-12p70", "IL-13", "IL-1β", 
                      "IL-2", "IL-4", "IL-6", "IL-8", "TNF-α")
write.csv(temp_df, "MicroPhenomics/inflammatory_factor/Inflammatory_factors.csv", 
          fileEncoding = "UTF-8")



temp_df = data.frame(ID = buId, matrix(runif(900, 10, 50), nrow=450, ncol=2))
colnames(temp_df) = c("ID", "IL-6_control", "IL-6_Pam")
write.csv(temp_df, "MicroPhenomics/ELISA/ELISA_Cytokines.csv", 
          row.names = FALSE)

temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "WB_Events_P1_CD45+", "WB_Events_P1_Mono", 
                      "WB_Events_P2_CD45+", "WB_Events_P3_CD45", 
                      "WB_Events_P4_CD45", "WB_Events_P5_CD45", 
                      "WB_Events_P1_CD16+_mDCs", "WB_Events_P1_pDCs", 
                      "WB_Events_P1_Classical_monocytes", "WB_Events_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/Events.csv", 
          row.names = FALSE)




temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "CD64_MFI_in_P1_CD45+", "CD64_MFI_in_P1_Mono", 
                      "CD64_MFI_in_P2_CD45+", "CD64_MFI_in_P3_CD45", 
                      "CD64_MFI_in_P4_CD45", "CD64_MFI_in_P5_CD45", 
                      "CD64_MFI_in_P1_CD16+_mDCs", "CD64_MFI_in_P1_pDCs", 
                      "CD64_MFI_in_P1_Classical_monocytes", "CD64_MFI_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/MFI.csv", 
          row.names = FALSE)


temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "FSC_MFI_in_P1_CD45+", "FSC_MFI_in_P1_Mono", 
                      "FSC_MFI_in_P2_CD45+", "SSC_MFI_in_P3_CD45", 
                      "SSC_MFI_in_P4_CD45", "FSC_MFI_in_P5_CD45", 
                      "FSC_MFI_in_P1_CD16+_mDCs", "SSC_MFI_in_P1_pDCs", 
                      "SSC_MFI_in_P1_Classical_monocytes", "FSC_MFI_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/Morphology.csv", 
          row.names = FALSE)


temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "WB_level2_in_P1_CD45+", "WB_level2_in_P1_Mono", 
                      "WB_level2_in_P2_CD45+", "WB_level2_in_P3_CD45", 
                      "WB_level2_in_P4_CD45", "WB_level2_in_P5_CD45", 
                      "WB_level2_in_P1_CD16+_mDCs", "WB_level2_in_P1_pDCs", 
                      "WB_level2_in_P1_Classical_monocytes", "WB_level2_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/Percentage_grandParent.csv", 
          row.names = FALSE)



temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "WB_Percentage_in_P1_CD45+", "WB_Percentage_in_P1_Mono", 
                      "WB_Percentage_in_P2_CD45+", "WB_Percentage_in_P3_CD45", 
                      "WB_Percentage_in_P4_CD45", "WB_Percentage_in_P5_CD45", 
                      "WB_Percentage_in_P1_CD16+_mDCs", "WB_Percentage_in_P1_pDCs", 
                      "WB_Percentage_in_P1_Classical_monocytes", "WB_Percentage_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/Percentage_in_CD45_cells.csv", 
          row.names = FALSE)



temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "WB_level1_in_P1_CD45+", "WB_level1_in_P1_Mono", 
                      "WB_level1_in_P2_CD45+", "WB_level1_in_P3_CD45", 
                      "WB_level1_in_P4_CD45", "WB_level1_in_P5_CD45", 
                      "WB_level1_in_P1_CD16+_mDCs", "WB_level1_in_P1_pDCs", 
                      "WB_level1_in_P1_Classical_monocytes", "WB_level1_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/Percentage_parent.csv", 
          row.names = FALSE)


temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "P1_CD45+/Lym", "P1_Mono/Lym", 
                      "P2_CD45+/Lym", "P3_CD45/Lym", 
                      "P4_CD45/Lym", "P5_CD45/Lym", 
                      "P1_CD16+_mDCs/Lym", "P1_pDCs/Lym", 
                      "P1_Classical_monocytes/Lym", "P6_CD45/Lym")
write.csv(temp_df, "MicroPhenomics/FACS/Ratio.csv", 
          row.names = FALSE)



temp_df = data.frame(ID = buId, matrix(runif(4500, 10, 50), nrow=450, ncol = 10))
colnames(temp_df) = c("ID", "CD64_MFI_in_P1_CD45+", "CD64_MFI_in_P1_Mono", 
                      "CD64_MFI_in_P2_CD45+", "CD64_MFI_in_P3_CD45", 
                      "CD64_MFI_in_P4_CD45", "CD64_MFI_in_P5_CD45", 
                      "CD64_MFI_in_P1_CD16+_mDCs", "CD64_MFI_in_P1_pDCs", 
                      "CD64_MFI_in_P1_Classical_monocytes", "CD64_MFI_in_P6_CD45")
write.csv(temp_df, "MicroPhenomics/FACS/MFI.csv", 
          row.names = FALSE)









temp_df = data.frame(ID = buId, `是否吸食香烟` = sample(c("是", "否", NA),   size=450, replace=TRUE))
write.csv(temp_df, "Demography/健康问卷调查-问卷-1.csv", row.names = FALSE)

