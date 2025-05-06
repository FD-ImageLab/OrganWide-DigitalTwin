# 加载必要的库
library(data.table)

paths = list.files("output_231227/gwas_results/", full.names = TRUE)
files = list.files("output_231227/gwas_results/")


for (i in seq_along(paths)) {
  path = paths[i]
  file = files[i]
  output_path = paste0("fuma/sandbox/", file)
  # 读取gz文件
  gwas_data <- fread(path)
  
  # 修改chr列，移除"chr"文本
  gwas_data[, chr := as.numeric(gsub("chr", "", chr))]
  
  # 将修改后的数据框写回到新的.gz文件
  fwrite(gwas_data, output_path, sep = "\t")
}


