library(tidyverse)
source("scripts/utils/Func.R")

files = list.files("scripts/ldsc/h2_log", full.names=TRUE)
pheno_dict = read.csv("scripts/pheno_dict.csv")

h2_df = data.frame()
for (file in files) {
    pheno = str_split(basename(file), "_")[[1]][1]
    d = readLines(file)
    d = d[grep("Total Observed scale h2:", d)]
    h2 = str_extract(d, "(?<=Total Observed scale h2: ).*?(?= \\()")
    se = str_extract(d, "(?<=\\().*(?=\\))")
    temp_h2_df = data.frame(pheno = pheno, h2 = h2, se=se)
    h2_df = rbind(h2_df, temp_h2_df)
    print(h2_df)    
}
h2_df$pheno = rename_based_on_df(h2_df$pheno, nmapdf=pheno_dict, from="fid", to="description")
h2_df$organ1 = h2_df$pheno %>% str_split("_") %>% map_chr(1)
h2_df$organ2 = h2_df$pheno %>% str_split("_") %>% map_chr(2)
h2_df$h2 = as.numeric(h2_df$h2)
h2_df$se = as.numeric(h2_df$se)


rev_h2_df = h2_df
rev_h2_df = rev_h2_df %>% dplyr::rename(organ1=organ2, organ2=organ1) %>% dplyr::select(pheno, h2, se, organ1, organ2)
rev_h2_df$pheno = paste(rev_h2_df$organ1, rev_h2_df$organ2, sep="_")

h2_df = rbind(h2_df, rev_h2_df)
h2_df$organ1 = factor(h2_df$organ1, levels=color_map$organ)
h2_df$organ2 = factor(h2_df$organ2, levels=color_map$organ)
h2_df = h2_df %>% arrange(organ1, organ2)
h2_df$organ1 = as.character(h2_df$organ1)
h2_df$organ2 = as.character(h2_df$organ2)
h2_df$col = rename_based_on_df(h2_df$organ1, nmapdf=color_map, from="organ", to="color")
h2_df$pheno = factor(h2_df$pheno, levels = h2_df$pheno)


png(filename = "h2.png", width = 6000, height = 1200)
ggplot(h2_df, aes(x = pheno, y = h2)) +
  geom_bar(stat = "identity", aes(fill=col), position = position_dodge(width = 0.7), width = 0.7) +
  geom_errorbar(aes(ymin = h2-se, ymax = h2+se), 
                width = .2,                    # Width of the error bars
                position = position_dodge(.7)) +
  scale_fill_identity() +
  theme_minimal() +
  labs(x = "Measure", y = "Heritability") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 50, colour = "black"), 
        axis.text.y = element_text(size=50, colour="black"), 
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=60)) 
dev.off()
