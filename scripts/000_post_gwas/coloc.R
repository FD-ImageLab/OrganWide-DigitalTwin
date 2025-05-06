files = list.files("fuma/ukb/coloc/", full.names = TRUE)

pheno_dict = read.csv("scripts/pheno_dict.csv")
f = function(file) {
  df = read.csv(file, row.names = 1)
  if (nrow(df) == 0) {
    df = data.frame()
  }
  else {
    df$disease = str_split(file_path_sans_ext(basename(file)), "_")[[1]][1]
    df$pheno = str_split(file_path_sans_ext(basename(file)), "_")[[1]][2]
    return(df)
  }
}

dfs = map(files, f)
df = reduce(dfs, rbind)
df$disease_name = rename_based_on_df(df$disease, nmapdf = pheno_dict, from = "fid", to = "names")
df$pheno_name = rename_based_on_df(df$pheno, nmapdf = pheno_dict, from = "fid", to = "description")
df = df %>%
  dplyr::select(snp, disease, disease_name, pheno, pheno_name, SNP.PP.H4)
write.csv(df, "C:/Users/liumeng/OneDrive/Desktop/coloc.csv")


library("biomaRt")
# 选择ENSEMBL的数据库
ensembl <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp", host = "grch37.ensembl.org")

# 根据rsid查询染色体位置
rsids <- "rs34872471"  # 替换为你的rsid列表

attributes <- c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end')

results <- getBM(attributes = attributes, values = rsids, mart = ensembl)

