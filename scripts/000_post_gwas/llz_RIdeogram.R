
# 安装RIdeogram（如果你还没有安装的话）
if (!requireNamespace("RIdeogram", quietly = TRUE)) {
  install.packages("RIdeogram")
}

library(RIdeogram)

# 创建虚拟数据
genomic_data <- data.frame(
  chr = c("1", "1", "2", "2"),
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  color = c("red", "blue", "green", "purple")
)

# 绘制染色体图
# 您可能需要根据您的染色体大小调整 ylim 参数
ideogram(
  df = genomic_data,  # 数据框
  chr = "chr",        # 染色体列名
  start = "start",    # 区域起始位置列名
  end = "end",        # 区域结束位置列名
  name = "name",      # 区域名称列名
  RColor = "color"    # 区域颜色列名
)


# 安装RIdeogram（如果你还没有安装的话）
if (!requireNamespace("RIdeogram", quietly = TRUE)) {
  install.packages("RIdeogram")
}

library(RIdeogram)

# 创建虚拟数据
genomic_data <- data.frame(
  chr = c("1", "1", "2", "2"),
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  RColor = c("red", "blue", "green", "purple")
)

# 准备ideogram的输入数据
ideo <- data.frame(
  chr = as.factor(genomic_data$chr),
  start = genomic_data$start,
  end = genomic_data$end,
  color = genomic_data$RColor,
  name = genomic_data$name
)

# 绘制ideogram
ideogram(
  cytoband = FALSE,  # 如果没有染色体带数据，设置为FALSE
  chrIndex = TRUE,   # 是否显示染色体索引
  ideoData = ideo    # 提供整理好的数据框
)


# 安装RIdeogram（如果还没有安装的话）
if (!requireNamespace("RIdeogram", quietly = TRUE)) {
  install.packages("RIdeogram")
}

library(RIdeogram)

# 创建虚拟数据
data <- data.frame(
  chr = c("1", "1", "2", "2"),
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  color = c("red", "blue", "green", "purple")
)

# 将颜色和标签列转换为因子
data$color <- as.factor(data$color)
data$name <- as.factor(data$name)

# 绘制染色体图
# 检查RIdeogram文档以确保参数正确
RIdeogram::fread.cytoband(
  data = data,
  species = "hg19",            # 选择物种，如果您的数据与hg19不同，请相应调整
  cytoband = FALSE,            # 设置为FALSE，因为我们不使用标准的cytoband数据
  centromere = FALSE,          # 设置为FALSE，除非您有着丝粒数据
  outfile = "ideogram_plot"    # 输出文件的名称
)


library(RIdeogram)

# 创建虚拟数据
data <- data.frame(
  chr = c("1", "1", "2", "2"),
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  color = c("red", "blue", "green", "purple")
)

# 绘制染色体图
ideoView(
  dfChr = data,
  centromere = NULL,
  chrIndex = TRUE,
  color = "color"
)



# 安装和加载ggplot2包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# 创建虚拟数据
data <- data.frame(
  chr = factor(c("1", "1", "2", "2"), levels = c("1", "2")), # 将染色体号转换为因子以确保排序
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  color = c("red", "blue", "green", "purple")
)

# 绘制基因组注释图
ggplot(data) +
  geom_segment(aes(x = start, xend = end, y = chr, yend = chr, colour = color), size = 5) +
  theme_minimal() +
  labs(x = "Position", y = "Chromosome", title = "Genomic Annotation") +
  theme(axis.text.y = element_text(angle = 0)) # 调整Y轴文字方向


# 安装和加载ggplot2包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(ggplot2)

# 创建虚拟数据
data <- data.frame(
  chr = factor(c("1", "1", "2", "2"), levels = c("1", "2")), # 将染色体号转换为因子以确保排序
  start = c(3000000, 10000000, 7000000, 25000000),
  end = c(5000000, 12000000, 9000000, 26000000),
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),
  color = c("red", "blue", "green", "purple")
)

# 绘制基因组注释图
ggplot(data) +
  geom_segment(aes(x = start, xend = end, y = chr, yend = chr, colour = color), linewidth = 5) +
  theme_minimal() +
  labs(x = "Position", y = "Chromosome", title = "Genomic Annotation") +
  theme(axis.text.y = element_text(angle = 0)) # 调整Y轴文字方向




# 加载RIdeogram包
library(RIdeogram)

# 创建虚拟数据
# 假设我们有两个染色体（1和2）和几个特定区域
data <- data.frame(
  chr = c("1", "1", "2", "2"),
  start = c(50000000, 100000000, 50000000, 150000000),  # 起始位置
  end = c(55000000, 105000000, 55000000, 155000000),    # 结束位置
  name = c("Region1", "Region2", "Region3", "Region4")  # 区域名称
)

# 绘制染色体和标记特定区域
RIdeogram::plot.ideogram(
  df = data, 
  species = "hg38",  # 选择物种的基因组版本，这里使用人类hg38版本
  chr.index = TRUE   # 显示染色体编号
)


install.packages("RIdeogram")
library(RIdeogram)

# 创建虚拟数据
data <- data.frame(
  chr = c("1", "1", "2", "2"),        # 染色体号
  start = c(3000000, 10000000, 7000000, 25000000),  # 起始位置
  end = c(5000000, 12000000, 9000000, 26000000),    # 结束位置
  name = c("RegionA", "RegionB", "RegionC", "RegionD")  # 区域名称
)

# 准备 ideogram 输入数据
ideo <- data.frame(
  chr = as.factor(data$chr),
  start = data$start,
  end = data$end,
  name = data$name
)

# 绘制染色体图
plot_ideogram(
  dfChr = ideo, 
  centromere = NULL, 
  chrIndex = TRUE, 
  is.add.band = FALSE, 
  is.add.label = TRUE
)


# 安装RIdeogram（如果还没有安装的话）
if (!requireNamespace("RIdeogram", quietly = TRUE)) {
  install.packages("RIdeogram")
}

library(RIdeogram)

# 创建虚拟数据
data <- data.frame(
  chr = c("1", "1", "2", "2"),        # 染色体号
  start = c(3000000, 10000000, 7000000, 25000000),  # 起始位置
  end = c(5000000, 12000000, 9000000, 26000000),    # 结束位置
  name = c("RegionA", "RegionB", "RegionC", "RegionD"),  # 区域名称
  color = c("red", "blue", "green", "purple")  # 区域颜色
)

# 绘制染色体图
ideogram(dfChr = data, species = "hg19")



#加载 RIdeogram
require(RIdeogram)
#加载包中的测试数据
data(human_karyotype, package="RIdeogram")
data(gene_density, package="RIdeogram")
data(Random_RNAs_500, package="RIdeogram")

#当然，我们也可以通过read.table()导入数据绘图。
install.packages("readxl")
library(readxl)

human_karyotype <- read.table("karyotype.txt", sep = "\t", header = T, stringsAsFactors = F)
gene_density <- read.table("data_1.txt", sep = "\t", header = T, stringsAsFactors = F)

SNP <- read_excel("/Users/lanlizhen/Documents/20FDU/phd/CrosstalkPaper/gwas/data\ /SNP_job.xlsx")

SNP$Shape <- as.character(SNP$Shape)
SNP$color <- as.character(SNP$color)

ordered_indices <- order(SNP$Shape, SNP$color)

#在计算gene_density时，包中贴心的给出了计算函数。只需要下载gff3文件，利用下边的代码句即可得到。
gene_density <- GFFex(input = "gencode.v32.annotation.gff3.gz", karyotype = "human_karyotype.txt", feature = "gene", window = 1000000)

setwd("/Users/lanlizhen/Documents/20FDU/phd/CrosstalkPaper/indNet")

#绘制染色体
ideogram(karyotype = human_karyotype)
convertSVG("chromosome.svg", device = "png")

#绘制gene密度
ideogram(karyotype = human_karyotype, overlaid = gene_density)
convertSVG("chromosome.svg", device = "png")

#绘制marker，可以是SNP、QTL或gene等marker
ideogram(karyotype = human_karyotype, label = SNP, label_type = "marker")
convertSVG("chromosome.svg", device = "png")

#同时绘制gene密度和marker
ideogram(karyotype = human_karyotype[1:22,], overlaid = gene_density, label = SNP, label_type = "marker")
convertSVG("chromosome.svg", device = "png")


#通过colorset1更改heatmap颜色
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = SNP, label_type = "marker", colorset1 = c("#fc8d59", "#ffffbf", "#91bfdb"))
convertSVG("chromosome.svg", device = "png")


#如果不知道centromere的位置，可绘制不包含centromere位置的染色体
human_karyotype <- human_karyotype[,1:3]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = SNP, label_type = "marker")
convertSVG("chromosome.svg", device = "png")

#当绘制的染色体比较少时，在画板固定width=170情况下，染色体的宽度就会增加（左下图）。
#通过调整width=100，染色体宽度变窄，视觉上比较好看（右下图）。
human_karyotype <- human_karyotype[1:10,]
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker")
convertSVG("chromosome.svg", device = "png")


#调整width后发现legend太靠右了，所以我们也要调整这个参数Lx和Ly，设置Lx = 80, Ly = 25。
ideogram(karyotype = human_karyotype, overlaid = gene_density, label = Random_RNAs_500, label_type = "marker", width = 100, Lx = 80, Ly = 25)
convertSVG("chromosome.svg", device = "png")

