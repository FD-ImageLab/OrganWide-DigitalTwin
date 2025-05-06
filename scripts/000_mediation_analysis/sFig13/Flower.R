# install.packages('plotrix')
library(plotrix) 

col_alpha = 0.85

# set path and load data
# flower_dat <- read.csv(dataPath, header = T, sep = ',')#read data
# flower_dat = flower_dat[order(match(flower_dat$organ_names, 
#                                     c("Brain", "Uterus", "Prostate", "Kidney", 
#                                       "Pancreas", "Spleen", "Liver", "Lung", "Heart"))), ]
# sample_id <- c(flower_dat[1:9,1])
# otu_num<-c(flower_dat[1:9,2])
# 
# center_data = read.csv(paste(mainPath,"/Data/Flower/flower_center.csv",sep=""), 
#                        header = T, sep = ',')
# center_data  = center_data %>%
#   filter(total!=0) %>%
#   arrange(desc(inter_organ_num))

# core_9 <- center_data$total[center_data$inter_organ_num == 9]
# core_8 <- center_data$total[center_data$inter_organ_num == 8]

# ellipse_col <- c('#bea589CC','#d3311fCC','#5f803fCC','#b8802bCC','#327487CC','#d9712aCC','#5360a0CC','#387e24CC','#dd7094CC')#set color

# ellipse_col = map_chr(flower_dat$organ_names, ~ color_map$color[which(color_map$organ == .x)])
# ellipse_col = alpha(ellipse_col, col_alpha)

flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(0,0,0,0))
  plot(c(2,8),c(2,8),type='n')
  n<-length(sample)
  deg <- 360/n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg*(t - 1))*pi/180), #plot
                 y = 5 + sin((start + deg*(t - 1))*pi/180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
  })

  res <- lapply(1:n, function(t){
    text(x = 5 + 2.65 *cos((start + deg * (t - 1)) * pi / 180),#number
         y = 5 + 2.65 *sin((start + deg * (t - 1)) * pi / 180),
         # bquote(frac(.(str_split(otu_num[t], "/")[[1]][1]), .(str_split(otu_num[t], "/")[[1]][2]))),
         otu_num[t], 
         cex = 1.9,
    )
    # text(x = 5 + 2 *cos((start + deg * (t - 1)) * pi / 180),#number
    #      y = 5 + 2 *sin((start + deg * (t - 1)) * pi / 180),
    #      sample[t],
    #      cex=0.7,
    #      col="white")
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 1.8 * cos((start + deg * (t - 1)) * pi / 180),#name
           y = 5 + 1.8 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           # srt = deg * (1 - t) ,
           adj = 0.5,
           cex = 1.9, font = 2.2
      )
    } else {
      text(x = 5 + 1.8 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 1.8 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           # srt = deg * (1 - t) ,
           adj = 0.5,
           cex = 1.9, font = 2.2
      )
    }
  })
  
  
  
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)#circle
  text(x = 5, y = 5.6, str_glue('Shared by {center_data$inter_organ_num[1]} '),cex = 1.6)#circle name
  text(x = 4.8, y = 5.2, 'organs:',cex = 1.6)#circle name
  text(x = 5.7, y = 5.25, center_data$total[1],cex = 1.6,col='red')#circle name
  text(x = 5, y = 4.8, str_glue('Shared by {center_data$inter_organ_num[2]} '),cex = 1.6 )#circle name
  text(x = 4.8, y = 4.4, 'organs:',cex = 1.6)#circle name
  text(x = 5.7, y = 4.45, center_data$total[2],cex = 1.6,col='red')#circle name
}

# check_path(paste0(resultPath, "/flower.png"))
# png(paste0(resultPath, "/flower.png"), width = 1500, height = 1500, res = 300)
# flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
#             start = 90, a = 1, b = 2.1, r = 1.3, ellipse_col = ellipse_col, circle_col = 'white')
# 
# dev.off()



######################################Organ#####################################

#####################################FACS#######################################

col_alpha = 0.85


flower_dat <- read.csv("data/sfig13/FACS_flower.csv", header = T, sep = ',')#read data
flower_dat$organ_names = rename_based_on_df(flower_dat$organ_names, 
                                            nmapdf = color_map, 
                                            from = "organ", to = "short_organ")
flower_dat = flower_dat[order(match(flower_dat$organ_names, 
                                    c("B", "U", "Pr", "K", 
                                      "Pa", "S", "Li", "Lu", "H"))), ]
sample_id <- c(flower_dat[1:7,1])
otu_num<-c(flower_dat[1:7,2])
otu_num = map_chr(str_split(otu_num, "/"), ~.x[2])

center_data = read.csv("data/sfig13/FACS_flower_center.csv", 
                       header = T, sep = ',')
center_data  = center_data %>%
  filter(total!=0) %>%
  arrange(desc(inter_organ_num))

ellipse_col = map_chr(flower_dat$organ_names, ~ color_map$color[which(color_map$short_organ == .x)])
ellipse_col = alpha(ellipse_col, col_alpha)

check_path("plot/sfig13/FACS_flower.png")
png("plot/sfig13/FACS_flower.png", width = 1500, height = 1500, res = 300)
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 1, b = 2.1, r = 1.3, ellipse_col = ellipse_col, circle_col = 'white')

dev.off()





#####################################IFMF#######################################
col_alpha = 0.85

# set path and load data
flower_dat <- read.csv("data/sfig13/IFMF_flower.csv", header = T, sep = ',')#read data
flower_dat$organ_names = rename_based_on_df(flower_dat$organ_names, 
                                            nmapdf = color_map, 
                                            from = "organ", to = "short_organ")
flower_dat = flower_dat[order(match(flower_dat$organ_names, 
                                    c("B", "U", "Pr", "K", 
                                      "Pa", "S", "Li", "Lu", "H"))), ]
sample_id <- c(flower_dat[1:7,1])
otu_num<-c(flower_dat[1:7,2])
otu_num = map_chr(str_split(otu_num, "/"), ~.x[2])


center_data = read.csv("data/sfig13/IFMF_flower_center.csv", 
                       header = T, sep = ',')
center_data  = center_data %>%
  filter(total!=0) %>%
  arrange(desc(inter_organ_num))

ellipse_col = map_chr(flower_dat$organ_names, ~ color_map$color[which(color_map$short_organ == .x)])
ellipse_col = alpha(ellipse_col, col_alpha)

check_path("plot/sfig13/IFMF_flower.png")
png("plot/sfig13/IFMF_flower.png", width = 1500, height = 1500, res = 300)
flower_plot(sample = sample_id, otu_num = otu_num, core_otu = core_num,
            start = 90, a = 1, b = 2.1, r = 1.3, ellipse_col = ellipse_col, circle_col = 'white')

dev.off()







