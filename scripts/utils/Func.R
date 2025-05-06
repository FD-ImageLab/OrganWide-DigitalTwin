library("tools")
library("tidyverse")
library("PMA")
# mainPath<-"D:/seafile/Seafile/SandBox/"
# mainPath <- "E:/seafile/Seafile/SandBox/"  # local
# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox




# 调整饱和度
desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}
# for (i in seq(1, 0, -0.01)) {
#   barplot(1: 9, col = desat(color_map$color, i))
#   title(str_glue("Saturation is {i}"))
#   Sys.sleep(0.05)
# }


saturation = 0.7

ukb_molecular_map = data.frame(molecular=c("NMR metabolomics", "Urine assays", "Blood chemistry"), 
                               color=c("red", "blue", "orange"))

color_map = data.frame(organ = c("Brain", "Heart", "Lung", "Liver", "Spleen", 
                                 "Pancreas", "Kidney", "Prostate", "Uterus"), 
                       short_organ = c("B", "H", "Lu", "Li", "S", "Pa", "K", "Pr", "U"), 
                       color = c("#BDB0A5", "#EB8677", "#BDDD78", 
                                 "#F2B670", "#7DBFA6", "#BDBBD7", 
                                 "#EE924F", "#7AADD2", "#DA8FC0"))
# color_map$color = desat(color_map$color, sat = saturation)
lr_organ_color_map = data.frame(organ = c("Left\nBrain", "Right\nBrain", 
                                          "Left\nHeart", "Right\nHeart", 
                                          "Left\nLung", "Right\nLung",
                                          "Left\nKidney", "Right\nKidney"
                                          ), 
                                color = c("#BDB0A5", "#BDB0A5", 
                                          "#EB8677", "#EB8677", 
                                          "#BDDD78", "#BDDD78", 
                                          "#EE924F", "#EE924F"))
# lr_organ_color_map$color = desat(lr_organ_color_map$color, sat = saturation)
vessel_color_map = data.frame(vessel = c("CA", "CAA", "AA", "SVC", 
                                         "CV", "PV", "HV", "RV", "IA"), 
                              color = c("#BDB0A5", "#AEDAE5", "#E8A79C", 
                                        "#D4A1D4", "#EB8677", "#BDDD78", 
                                        "#F2B670", "#EE924F", "#7AADD2"))
# vessel_color_map$color = desat(vessel_color_map$color, sat = saturation)
cell_color_map = data.frame(cell = c("Mono", "NK", "Treg", "T-cell", "Th", "B-cell"), 
                            cate = c("P1", "P2", "P3", "P4", "P5", "P6"), 
                            color = c("#B3CBDD", "#BF9D97", "#C9D2A5", 
                                      "#E2A0A2", "#EAC39C", "#C2B5D3"))

micro_cell_color_map = data.frame(type = c("DCs", 
                                           "Mono", "NK", "Treg", "T-cell",
                                           "Th", "B-cell", 
                                           "MF", "IF"), 
                                  color = c("#EE924F", 
                                            "#B3CBDD", "#BF9D97", "#C9D2A5", 
                                            "#E2A0A2", "#EAC39C", "#C2B5D3", 
                                            "#DBDA96", "#A5D4DD"))


micro_color_map = data.frame(micro = c("Demography", 
                                       "Serum Routine", 
                                       "Inflammatory Factor", 
                                       "Metabolic Factor", 
                                       "Immunocyte FACS",
                                       "Immunocyte ELISA"), 
                             before_name = c("Demography", 
                                             "physical_examination", 
                                             "inflammatory_factor", 
                                             "metabolic_factor", 
                                             "facs", 
                                             "elisa"), 
                             sn = c("DG", 
                                    "Serum", 
                                    "IF", 
                                    "MF", 
                                    "FACS", 
                                    "ELISA"), 
                             # color = c("#F27970", "#BB9727",
                             #           "#54B345", "#32B897", "#05B9E2","#8983BF"),
                             color = c("#F27970", 
                                       "#ED926B",
                                       "#A5D4DD", 
                                       "#DBDA96", 
                                       "#C1B1D2",
                                       "#E2CB9B"))


                             # color = c("#F6931C", "#FEDD72", "#109474", "#822F28", "#8EC63F", "#202F4D"))

# radiomics_color_map = data.frame(radiomics = c("Shape", "Function", 
#                                                "Density", "Texture", "Vessel"), 
#                                  color = c("#91CBC0", "#9FAAC2",  
#                                            "#F3CCC1", "#BF9E97", "#D68485"))
radiomics_color_map = data.frame(radiomics = c("Morphology", "Function",
                                               "Water", "Diffusion", 
                                               "Fat", "Iron", 
                                               "Texture","Vessel"), 
                                 color = c("#7FCDC0", "#9DAAC4",  
                                           "#FACBBF", "#C55A11", 
                                           "#C3BF3F", "#FBCA73", 
                                           "#C59C95", "#E37F83"))


# check_path(paste0(mainPath, "/Data/ColorMap/organ_color_map.csv"))
# write.csv(color_map, paste0(mainPath, "/Data/ColorMap/organ_color_map.csv"), row.names = FALSE)
# write.csv(vessel_color_map, paste0(mainPath, "/Data/ColorMap/vessel_color_map.csv"), row.names = FALSE)
# write.csv(lr_organ_color_map, paste0(mainPath, "/Data/ColorMap/lr_organ_color_map.csv"), row.names = FALSE)
# write.csv(cell_color_map, paste0(mainPath, "/Data/ColorMap/cell_color_map.csv"), row.names = FALSE)
# write.csv(micro_color_map, paste0(mainPath, "/Data/ColorMap/micro_color_map.csv"), row.names = FALSE)
# write.csv(radiomics_color_map, paste0(mainPath, "/Data/ColorMap/radiomics_color_map.csv"), row.names = FALSE)
# write.csv(micro_cell_color_map, paste0(mainPath, "/Data/ColorMap/micro_cell_color_map.csv"), row.names = FALSE)
# 默认用行名进行merge
merge_and_divide <- function(X, Y, by=0) {
  if (by == 0) {
    temp <- merge(X, Y, by=0)
    rownames(temp) <- temp$Row.names
    temp$Row.names <- NULL
  }
  else {
    temp <- merge(X, Y, by=by)
  }
  temp <- na.omit(temp)
  xcnames <- colnames(X)
  ycnames <- colnames(Y)
  X <- temp[xcnames]
  Y <- temp[ycnames]
  return(list(X, Y))
}

my_formula <- function(y_variable, x_variables, showEnv=FALSE) {
  addbacktick = function(str) {paste0("`", str, "`")}
  fo <- paste0(addbacktick(y_variable), " ~ ", paste(addbacktick(x_variables), collapse = " + "))
  fmla <- as.formula(fo)
  return(fmla)
}

# 根据t2对t1的每一列进行回归;t1 and t2是df
predict_lm_tibble <- function(t1, t2) {
  if (!(is.data.frame(t1) & is.data.frame(t2))) {
    stop("t1 or t2 is not dataframe")
  }
  if (identical(intersect(rownames(t1), rownames(t2)), character(0))) {
    stop("t1 and t2 have no intersection")
  }
  ct2 <- colnames(t2)
  rt1 <- rownames(t1)
  rt1_df = tibble(id = rt1)
  t2 = rownames_to_column(t2, var = "id")
  tf = function(temp_t1, ct1) {
    ttemp_t1 = tibble(id = rt1, !!enquo(ct1) := temp_t1)
    temp = left_join(ttemp_t1, t2, by="id") %>%
      column_to_rownames(var = "id") %>%
      lm(my_formula(ct1, ct2), data = .) %>%
      predict
    temp = tibble(id = names(temp), !!enquo(ct1) := temp) %>%
      left_join(x = rt1_df, y = ., by = "id") %>%
      column_to_rownames(var = "id")
  }
  pt1 = imap_dfc(t1, tf)
  return(pt1)
}

# 按照category_df的一列对表型进行分类，var参数必须传入字符串，select_category为
# var这一列中选择的分类，且传入时必须使用select_category命名
# 字符串使用sym function使str变成symbol
# category_df必须包含Name列，值是data_df的列名，如果需要对data_df的列名进行映射， 
# 加入fnmap/feature name map col其值为category_df中另一列的列名，如Short_name,
# 默认为NULL，即不需要对data_df的列名进行映射
divide_data_by_var_df = function(data_df, category_df, var, select_category=NULL, fnmap = NULL) {
  if (class(var) == "character") {
    var = sym(var)
    cate_df = category_df %>% 
      group_split(!!var) %>%
      lapply(function(x) {
        if (all(x$Name %in% colnames(data_df))) {
          if (is.null(fnmap)) {
            tp = data_df[x$Name]
            return(tp)
          }
          else if (fnmap %in% colnames(category_df)) {
            tp = data_df[x$Name]
            tp = rename_based_on_df(tp, category_df, "Name", fnmap)
            return(tp)
          }
          else if ((!(fnmap %in% colnames(category_df))) & !is.null(fnmap)) {
            warning(str_glue("{fnmap} is not in the colnames(category_df),\\
                             please check the variable \"fnmap\"."))
            tp = data_df[x$Name]
            return(tp)
          }
        }
        else {
          warning(str_glue("\"Name\" in \"{unique(x[[var]])}\" is not all\\
                           contained in the colnames of data_df, please pay attention."))
          if (is.null(fnmap)) {
            tp = data_df[x$Name[x$Name %in% colnames(data_df)]]
            return(tp)
          }
          else if (fnmap %in% colnames(category_df)) {
            tp = data_df[x$Name[x$Name %in% colnames(data_df)]]
            tp = rename_based_on_df(tp, category_df, "Name", fnmap)
            return(tp)
          }
          else if ((!(fnmap %in% colnames(category_df))) & !is.null(fnmap)) {
            warning(str_glue("{fnmap} is not in the colnames(category_df),\\
                             please check the variable \"fnmap\"."))
            tp = data_df[x$Name]
            return(tp)
          }
        }
      })
    names(cate_df) = category_df %>%
      group_split(!!var) %>%
      lapply(function(x) {
        if (all(x$Name %in% colnames(data_df))) {
          return(unique(x[[var]]))
        }
        else {
          return(unique(x[[var]]))
        }
      })
    if (!is.null(select_category)) {
      cate_df = cate_df[select_category[select_category %in% names(cate_df)]]
    }
    return(cate_df)
  }
}

# 按照Category和Region进行分类，如果select_category为空，则将所有分类都提取出来
divide_data_by_category_and_region_df = function(data_df, category_df, select_category=NULL) {
  l1 = divide_data_by_var_df(data_df, category_df, "Category")
  l2 = divide_data_by_var_df(data_df, category_df, "Region")
  ll = c(l1, l2)
  if (!is.null(select_category)) {
    ll = ll[names(ll) %in% select_category]
  }
  return(ll)
}

# 根据nmapdf的from和to列进行df/list/chr的colnames/names的替换
# from/to 可以是变量
rename_based_on_df = function(df, nmapdf, from, to) {
  if (class(df) == "data.frame") {
    for (i in seq_along(df)) {
      if (colnames(df)[i] %in% nmapdf[[from]]) {
        colnames(df)[i] = nmapdf[[to]][which(nmapdf[[from]] == colnames(df)[i])]
      } 
      else {warning(str_glue("{colnames(df)[i]} is not in nmapdf"))}
    }
  }
  if (class(df) == "list") {
    oneton_len = lapply(unique(nmapdf[[from]]), function(x) {
      length(nmapdf[[to]][which(nmapdf[[from]] == x)]) > 1
    })
    if (sum(unlist(oneton_len)) > 1) {
      stop("1 to n is over 1 pair, please transform the nmapdf")
    }
    
    j = 1
    for (i in seq_along(df)) {
      if (names(df)[i] %in% nmapdf[[from]]) {
        if (length(nmapdf[[to]][which(nmapdf[[from]] == names(df)[i])]) == 1) {
          names(df)[i] = nmapdf[[to]][which(nmapdf[[from]] == names(df)[i])]
        }
        else {
          names(df)[i] = nmapdf[[to]][which(nmapdf[[from]] == names(df)[i])][j]
          j = j + 1
        }
      }
      else {warning(str_glue("{colnames(df)[i]} is not in nmapdf"))}
    }
  }
  if (class(df) == "character") {
    for (i in seq_along(df)) {
      if (df[i] %in% nmapdf[[from]]) {
        df[i] = nmapdf[[to]][which(nmapdf[[from]] == df[i])]
      } 
    }
  }
  return(df)
}

# 输入为selected表格
# part_name是df中的列名，表示生成部分的名称，这个字段会被加到生成的df的列名中，
# 主要目的是为了防止列名重复，所以df的part_name必须不能重复。
# class是df中的列名，表示生成部分是如何分类的，gen_part返回的是list，list的元素
# 是不同的df，df的名称命名为class中的类别。
# df的列名必须包括File、RefFile、SelectedColumn、SelectedName， 
gen_part = function(df, File_encoding = "", File_fileEncoding = "", 
                    RefFile_encoding = "", RefFile_fileEncoding = "", 
                    part_name = "Prefix", class = "Modality", 
                    fnmap = NULL) {
  file_positions = df$File
  data = lapply(c(1: nrow(df)), function(r) {
    file_position = df$File[r]
    if (!file.exists(file_position)) {
      warning(str_glue("{file_position} not exist!"))
      part_data = list(NULL)
      names(part_data) = df[[part_name]][r]
      return(part_data)
    }
    else {
      if (file_ext(file_position) == "csv") {
        File = read.csv(file_position, header = TRUE, fileEncoding=File_fileEncoding, 
                        encoding=File_encoding, check.names = FALSE)
        File = File[!duplicated(File[[1]]),]
        rownames(File) = File[[1]]
        File = File[-1]
      } 

      reffile_position = df$RefFile[r]
      RefFile = read.csv(reffile_position, header = TRUE, 
                         encoding = RefFile_encoding, 
                         fileEncoding = RefFile_fileEncoding)
      SelectedColumn = df$SelectedColumn[r]
      SelectedName = df$SelectedName[r]
      gen_part_data = function(File, RefFile, SelectedColumn, SelectedName, fnmap) {
        temp = df[r,,drop=FALSE]
        temp_result = divide_data_by_var_df(File, RefFile, SelectedColumn, SelectedName, fnmap = fnmap)
        temp_result = rename_based_on_df(temp_result, temp, "SelectedName", part_name)
        add_preffix = function(df, preffix) {
          colnames(df) = paste0(preffix, "_", colnames(df))
          return(df)
        }
        preffix = df[[part_name]][r] # 这边的Prefix可以改成part_name
        temp_result = map(temp_result, add_preffix, preffix = preffix)
        return(temp_result)
      }
      part_data = gen_part_data(File, RefFile, SelectedColumn, SelectedName, fnmap)
      return(part_data)
    }
  })
  data = flatten(data)
  #增加df合并模块，将再class中相同的df合并
  f = function(data, new_data) {
    class_idxs = get_index_of_cate_from_array(df[[class]])
    new_data = vector("list", length(class_idxs))
    names(new_data) = names(class_idxs)
    for (i in seq_along(class_idxs)) {
      class_idx = class_idxs[[i]]
      tmp_data = data[class_idx]
      not_null = !unlist(map(tmp_data, is.null))
      
      tmp_data = tmp_data[not_null]
      if (length(tmp_data) == 0) {
        new_data[[i]] = NULL
      }
      else {
        tmp_data = map(tmp_data, rownames_to_column, var="buId")
        tmp_data = reduce(tmp_data, inner_join)
        new_data[[i]] = column_to_rownames(tmp_data, var = "buId")
      }
    }
    return(new_data)
  }
  new_data = f(data)
  return(new_data)
}
# 输入为两张df，默认使用行名进行匹配，即输入的两个df的行名需要有相同的部分，
# 如果相同的行少于20个，报warning，返回的结果为NA
# error 为啥输出为cca=1，可能是X和Y有列相同导致报错？
get_cca_result = function(X, Y) {
  colnames(X) = paste0("i_", colnames(X))
  colnames(Y) = paste0("j_", colnames(Y))
  temp_xy_data = merge_and_divide(X, Y)
  X  = temp_xy_data[[1]]
  Y = temp_xy_data[[2]]
  if (nrow(X) >= 20 & (!is.null(X)) & (!is.null(Y))) {
    X <- X[c(map_dfc(X, var, na.rm=TRUE) != 0)]
    Y <- Y[c(map_dfc(Y, var, na.rm=TRUE) != 0)]
    X_len = length(X)
    Y_len = length(Y)
    if (X_len == 1) {
      X["icopy_col"] = X[[1]]
    }
    if (Y_len == 1) {
      Y["jcopy_col"] = Y[[1]]
    }
    tryCatch({cc_result <- PMA::CCA(X, Y, typex = "standard", typez = "standard",
                                    penaltyx = 1, penaltyz = 1, trace = FALSE, 
                                    standardize = TRUE)
    u <- cc_result$u
    v <- cc_result$v
    cor_main_coef <- cc_result$cors
    cc1_organ <- as.matrix(X) %*% u
    cc2_organ <- as.matrix(Y) %*% v
    cor_main_coef.p = cor.test(cc1_organ, cc2_organ)$p.value
    results = list(cca = cor_main_coef, cca.p = cor_main_coef.p, 
                   X_prime = list(cc1_organ), Y_prime = list(cc2_organ), N = nrow(X))
    }, 
    error = function(e) {
      results = list(cca = 1, cca.p = 0, X_prime = list(cc1_organ), 
                     Y_prime = list(cc2_organ), N = nrow(X))
    })
  }
  else {
    warning("observations are not enough")
    results = list(cca = NA, cca.p = NA, X_prime = NA, 
                   Y_prime = NA, N = nrow(X))}
}

# 删除全是na的列(x is dataframe)
removeColsAllNa  <- function(df=df, na_rate=1){
  df[, apply(df, 2, function(y) {(sum(is.na(y))/length(y)) <= na_rate})]
  return(df)
}

removeRowsAllNa  <- function(df=df, na_rate=1){
  r.df = rownames(df)[apply(df, 1, function(y) {(sum(is.na(y))/length(y)) <= na_rate})]
  df = df[apply(df, 1, function(y) {(sum(is.na(y))/length(y)) <= na_rate}) ,]
  rownames(df) = r.df
  return(df)
}

removeVar0 = function(df=df) {
  df = df[, which(apply(df, 2, var, na.rm=TRUE) != 0)]
  return(df)
}

# 根据regression_data回归data数据
get_res_with_regression_data = function(data=data, regression_data = regression_data) {
  if (is.null(data)) {
    res = data.frame()
  }
  else if (identical(intersect(rownames(data), rownames(regression_data)), character(0))) {
    res = data
    res[res==res] = 0
  }
  else {
    res = data - predict_lm_tibble(data, regression_data)
  }
  return(res)
}

# 修改merge函数，在x或y长度为0时返回为另一个
my_merge <- function(x, y, by='') {
  if (dim(x)[1] == 0) {
    z <- y
  } else if (dim(y)[1] == 0) {
    z <- x
  } else {
    z <- merge(x, y, by=by)
  }
  return(z)
}


# check file_path是否存在，不存在就warning并创建文件夹
check_path = function(s, mute=FALSE, creat_dir=FALSE) {
  if (grepl("\\.", basename(s))) {
    if (!dir.exists(dirname(s))) {
      dir.create(dirname(s), recursive = TRUE)
      if (!mute) {
        warning(str_glue("path: {dirname(s)} has been created!"))
      }
    }
    else {
      sprintf("path: %s exists", dirname(s))
    }
  }
  else {
    if (!dir.exists(s)) {
      dir.create(s, recursive = TRUE)
      if (!mute) {
        warning(str_glue("path: {s} has been created!"))
      }
    }
    sprintf("path: %s exists", s)
  }
}

ReadMicrophenomics = function(dataPath=paste0(mainPath, "/Sandbox/")) {
  inflammatory_factor = read.csv(file = paste0(dataPath, "/MicroPhenomics/inflammatory_factor/Inflammatory_factors.csv"),
                                   header = TRUE, row.names = 1,
                                   check.names = FALSE, encoding = "UTF-8")
  
  matabolic_factor = read.csv(file = paste0(dataPath, "/MicroPhenomics/metabolic_factor/Metabolic_factors.csv"),
                              header = TRUE, row.names = 1, check.names = FALSE)
  
  elisa = read.csv(file = paste0(dataPath, "/MicroPhenomics/ELISA/ELISA_Cytokines.csv"), 
                   header = TRUE, row.names = 1, check.names = FALSE)
  facs_files = list.files(paste0(dataPath, "/MicroPhenomics/FACS/"), pattern = ".csv", 
                          full.names = TRUE)
  facs = map_dfc(facs_files, read.csv, row.names = 1, header = TRUE, 
                 check.names = FALSE)
  cnames <- c("高敏感C反应蛋白（mg/L）", "淋巴细胞数（*10^9/L）", "总胆红素（μmol/L）", "尿酸（μmol/L）", "载脂蛋白A-I（g/L）", 
              "载脂蛋白B（g/L）", "载脂蛋白E（mg/L）", "肌酐（μmol/L）", "癌胚抗原（ng/ml）", 
              "胰岛素空腹（uU/ml）", "葡萄糖（mmol/L）", "高密度脂蛋白胆固醇（mmol/L）", 
              "低密度脂蛋白胆固醇（mmol/L）", "甲胎蛋白（ng/ml）", 
              "丙氨酸氨基转移酶（U/L）", "门冬氨酸氨基转移酶（U/L）")
  mapcnames <- c("CRP", "lymph", "TB", "UA", "ApoA", "ApoB", "ApoE", "CRE", "CEA", "insulin", "glucose", "HDL", "LDL", 
                 "ALT", "AST", "AFP")
  physical_examination <- read.csv(paste0(dataPath, "/Serum/志愿者招募-常规体检.csv"), 
                                     header = TRUE, row.names = 1, check.names = FALSE)
  physical_examination <- physical_examination %>% 
    dplyr::select(dplyr::all_of(cnames))
  physical_examination <- as.data.frame(lapply(physical_examination, as.numeric), row.names = rownames(physical_examination))
  colnames(physical_examination) <- mapcnames  
  
  microphenomics <- list(`Serum` = physical_examination, 
                         `IF` = inflammatory_factor,
                         `MF` = matabolic_factor,
                         `FACS` = facs, 
                         `ELISA` = elisa
                         )
  microphenomics = map(microphenomics, removeColsAllNa)
  microphenomics = map(microphenomics, removeRowsAllNa)
  return(microphenomics)
}



# 返回与key列的value值最接近的n个数据的id号
select_id = function(df, key="BMI", value=20, otrm = TRUE, n=100) {
  key = sym(key)
  value = as.numeric(value)
  if (otrm == TRUE) {
    df = df %>% 
      filter(!df[[key]] %in% boxplot.stats(df[[key]])$out)
  }
  if (nrow(df) <= n) {
    warning(str_glue("n is larger than nrow of df"))
    temp_df = df %>%
      dplyr::mutate(px = abs(!!key - value)) %>%
      dplyr::arrange(px) %>%
      slice_min(px, n = nrow(df), with_ties = FALSE) %>%
      select(-px)
    return(rownames(temp_df))
  }
  temp_df = df %>%
    dplyr::mutate(px = abs(!!key - value)) %>%
    dplyr::arrange(px) %>%
    slice_min(px, n = n, with_ties = FALSE) %>%
    dplyr::select(-px)
  return(rownames(temp_df))
}
# select_id(regression_data, "年龄", value=25, n = 100)



# 读取category文件，path为category文件夹路径，其下包括多个器官的selected.csv文件，
# 返回为list，元素是多个器官的list，器官的list包括器官的不同部位的df，part_name
# 是在selected.csv文件的一列，是在返回的df中，给列名中添加的字段
read_data = function(path, part_name="Prefix", class="Modality", fnmap = NULL) {
  ll = list.files(path, full.names = TRUE, pattern = "csv")
  data = lapply(ll, function(l) {
    temp_df = read.csv(l, fileEncoding="UTF-8-BOM", header=TRUE)
    organ = gen_part(df = temp_df, File_encoding = "GBK", File_fileEncoding = "UTF-8-BOM",
                     RefFile_fileEncoding = "UTF-8-BOM", part_name = part_name, 
                     class = class, fnmap = fnmap)
    return(organ)
  })
  data_names = lapply(ll, function(l) {
    temp_df = read.csv(l, fileEncoding="UTF-8-BOM", header=TRUE)
    organ_name = names(table(temp_df$Organ[temp_df$Organ > 0]))[which.max(table(temp_df$Organ[temp_df$Organ > 0]))]
    return(organ_name)
  })
  names(data) = data_names
  return(data)
}


read_regression_data = function(SEX_=TRUE, BMI_=TRUE, AGE_=TRUE) {
  #回归数据
  BMI <- read.csv(file = "Sandbox/Demography/DXA体成分-体成分高分辨率数字成像化系统.csv",
                  header = TRUE, row.names = 1, check.names = FALSE)
  BMI <- BMI[, 1, drop=FALSE]
  sex_age <- read.csv(file = "Sandbox/Demography/基本信息.csv",
                      header = TRUE, row.names = 1, check.names = FALSE)
  sexual <- sex_age[, 1, drop=FALSE]
  sexual <- map_df(sexual, ~ ifelse(. == "男", 1, 0))
  rownames(sexual) <- rownames(sex_age)
  age <- sex_age[, 2, drop = FALSE]
  partial = NULL
  if (SEX_) {
    partial = bind_cols(partial, sexual)
  }
  if (BMI_) {
    partial = bind_cols(partial, BMI)
  }
  if (AGE_) {
    partial = bind_cols(partial, age)
  }
  return(partial)
}

read_feature_data = function() {
  BMI <- read.csv(file = paste0(mainPath, "/Sandbox/Demography/DXA体成分-体成分高分辨率数字成像化系统.csv"),
                  header = TRUE, row.names = 1, check.names = FALSE)
  BMI <- BMI[, 1, drop=FALSE]
  sex_age <- read.csv(file = paste0(mainPath, "/Sandbox/Demography/基本信息.csv"),
                      header = TRUE, row.names = 1, check.names = FALSE)
  sexual <- sex_age[, 1, drop=FALSE]
  sexual <- map_df(sexual, ~ ifelse(. == "男", 1, 0))
  rownames(sexual) <- rownames(sex_age)
  age <- sex_age[, 2, drop = FALSE]
  partial <- cbind(BMI, sexual, age)
  
  nc = read.csv(file = paste0(mainPath, "/Sandbox/Demography/人体外观测量-三维人体扫描分析系统.csv"), 
                header = TRUE, row.names = 1, check.names = FALSE)
  
  tjdata <- read.csv(paste0(mainPath, "/Sandbox/Serum/志愿者招募-常规体检.csv"), 
                     header = TRUE, check.names = FALSE, row.names = 1)
  tjdata <- tjdata[c("血压(mm/Hg)收缩压（mm/Hg）", "血压(mm/Hg)舒张压（mm/Hg）", "葡萄糖（mmol/L）", 
                     "心率")]
  features <- bind_cols(partial, nc, tjdata)
  mapcnames <- c("BMI", "Gender", "Age", "NC", "SBP", "DBP", "Glucose", "HR")
  colnames(features) <- mapcnames
  return(features)
}

# ggplot theme
# newtheme <- theme_bw() + theme(plot.title = element_text(color = "darkred"))




# 实现matrix的两种索引方式互换，索引方式1是将matrix当做array，索引是matrix的第几个元素，
# 索引方式2是将matrix当做dataframe，索引方式是第几行第几列
indexA2indexB = function(mat, single_index=NULL, multi_index=NULL) {
  if (!is.null(single_index)) {
    if (!is.numeric(single_index)) {
      stop("single index is not a index")
    }
    whichrow = ifelse(single_index %% dim(mat)[1] == 0, dim(mat)[1], single_index %% dim(mat))
    whichcol = ifelse(single_index %% dim(mat)[1] == 0, single_index %/% dim(mat)[1], single_index %/% dim(mat)[1] + 1)
    
    indexB = c(whichrow, whichcol)
    return(indexB)
  }
  if (!is.null(multi_index)) {
    whichrow = multi_index[1]
    whichcol = multi_index[2]
    indexA = whichcol * (dim(mat)[1] - 1) + whichrow
    return(indexA)
  }
}


# 如array=c("A", "B", "A", "B", "A", "C")，
# 将array传入函数后返回的结果为list(A=c(1, 3, 5), B=c(2, 4), C=c(6))
get_index_of_cate_from_array = function(array) {
  cates = unique(array)
  return_list = vector("list", length = length(cates))
  names(return_list) = cates
  for (cate in cates) {
    idx = which(array==cate)
    return_list[[cate]] = idx
  }
  return(return_list)
}





data_preprocess = function(df, col_na_rate=1, row_na_rate=1) {
  df = df %>% 
    rownames_to_column(var="buId") %>% 
    removeColsAllNa(na_rate = col_na_rate) %>%
    column_to_rownames(var="buId") %>%
    removeRowsAllNa(na_rate = row_na_rate) %>%
    removeVar0()
  return(df)
}





FONTFAMILY="sans"






