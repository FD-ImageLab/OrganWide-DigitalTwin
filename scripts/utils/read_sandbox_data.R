read_carotidartery = function(dataPath=dataPath) {
  CarotidArtery = read.csv(paste0(dataPath, "Vessel/超声-多普勒超声诊断仪.csv"), 
                           header = TRUE, check.names = FALSE)
  right_common_carotid = CarotidArtery[intersect(grep("右", colnames(CarotidArtery)), grep("颈总动脉", colnames(CarotidArtery)))]
  right_external_carotid_artery = CarotidArtery[intersect(grep("右", colnames(CarotidArtery)), grep("颈外动脉", colnames(CarotidArtery)))]
  right_internal_carotid = CarotidArtery[intersect(grep("右", colnames(CarotidArtery)), grep("颈内动脉", colnames(CarotidArtery)))]
  left_common_carotid = CarotidArtery[intersect(grep("左", colnames(CarotidArtery)), grep("颈总动脉", colnames(CarotidArtery)))]
  left_external_carotid_artery = CarotidArtery[intersect(grep("左", colnames(CarotidArtery)), grep("颈外动脉", colnames(CarotidArtery)))]
  left_internal_carotid = CarotidArtery[intersect(grep("左", colnames(CarotidArtery)), grep("颈内动脉", colnames(CarotidArtery)))]
  CarotidArtery = list(RCC = right_common_carotid, 
                       # RECA = right_external_carotid_artery, 
                       RICA = right_internal_carotid, 
                       LCC = left_common_carotid, 
                       # LECA = left_external_carotid_artery, 
                       LICA = left_internal_carotid)
  removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y))), drop=FALSE]}
  CarotidArtery = map(CarotidArtery, removeColsAllNa) %>% map(as.data.frame)
  return(CarotidArtery)
}