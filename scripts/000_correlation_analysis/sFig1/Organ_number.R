# install.packages('plotrix')
library(plotrix) 
# user-defined path
mainPath <- "/mnt/d/University/fdurop/remote_ws/indNet/"  # local
# mainPath <- "E:/seafile/Seafile/SandBox/"  # local
# mainPath <- "/Users/apple/Seafile/Research/SandBox"  # wangcy
# mainPath <- "/public/sandbox/workdir/wangchy/SharedFolder"  # sandbox
# mainPath <- "/public/sandbox/workdir/liumeng/SandBox"  # sandbox
resultPath = paste0(mainPath, "/plot/sfig1/")
source(paste(mainPath,"scripts/utils/Func.R",sep=""), encoding = "UTF-8")

















col_alpha = 0.85
color_map = data.frame(organ = c("Brain", "Heart", "Lung", "Liver", "Spleen", 
                                 "Pancreas", "Kidney", "Prostate", "Uterus"), 
                       color = c("#BDB0A5", "#EB8677", "#BDDD78", 
                                 "#F2B670", "#7DBFA6", "#BDBBD7", 
                                 "#EE924F", "#7AADD2", "#DA8FC0"))
# set path and load data
sample_id <- c("Brain", "Kidney", "Pancreas","Spleen","Liver","Lung","Heart")

# ellipse_col <- c('#bea589CC','#d3311fCC','#5f803fCC','#b8802bCC','#327487CC','#d9712aCC','#5360a0CC','#387e24CC','#dd7094CC')#set color


ellipse_col =  c("#BDB0A5", "#EE924F","#BDBBD7", "#7DBFA6", "#F2B670","#BDDD78","#EB8677")

ellipse_col = alpha(ellipse_col, col_alpha)
flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col,center) {
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
    text(x = 5 + 2.8 *cos((start + deg * (t - 1)) * pi / 180),#number
         y = 5 + 2.8*sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t],
         cex = 1.3,
         
    )
    # text(x = 5 + 2 *cos((start + deg * (t - 1)) * pi / 180),#number
    #      y = 5 + 2 *sin((start + deg * (t - 1)) * pi / 180),
    #      sample[t],
    #      cex=0.7,
    #      col="white")
    
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 1.9 * cos((start + deg * (t - 1)) * pi / 180),#name
           y = 5 + 1.9 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           # srt = deg * (1 - t) ,
           adj = 0.5,
           cex = 1.3,
           
      )
    } else {
      text(x = 5 + 1.9 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 1.9 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           # srt = deg * (1 - t) ,
           adj = 0.5,
           cex = 1.3,
      )
    }
  })
  
  
  
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)#circle
  text(x = 5, y = 5.2, str_glue('Shared by 7 '),cex = 1.4)#circle name
  text(x = 4.8, y = 4.8, 'organs:',cex = 1.4)#circle name
  text(x = 5.7, y = 4.85, center,cex = 1.4,col='red')#circle name
}


# check_path(paste0(resultPath, "/flower.png"))
png(paste0(resultPath, "/flower_sandbox.png"), width =1500, height = 1500, res = 300)
flower_plot(sample = sample_id, otu_num = c(1065,1075,1032,1032,1075,1075,1091) - 693, core_otu = core_num,
            start = 90, a = 1, b = 2.2, r = 1.25, ellipse_col = ellipse_col, circle_col = 'white',center = c(693))
dev.off()


png(paste0(resultPath, "/flower_ukb.png"), width =1500, height = 1500, res = 300)
flower_plot(sample = sample_id, otu_num = c(46393,36137,29464,28475,27418,23154,32462) - 10060, core_otu = core_num,
            start = 90, a = 1, b = 2.2, r = 1.25, ellipse_col = ellipse_col, circle_col = 'white',center = c(10060))
dev.off()


# png(paste0(resultPath, "/flower_male.png"), width = 1500, height = 1500, res = 300)
# flower_plot(sample = sample_id[-which(sample_id=="Uterus")], otu_num = c(118,119,116,117,117,117,115,158), core_otu = core_num,
#             start = 90, a = 1, b = 2.1, r = 1.3, ellipse_col = ellipse_col[-2], circle_col = 'white',center = c(188))
# dev.off()
# png(paste0(resultPath, "/flower_female.png"), width = 1500, height = 1500, res = 300)
# flower_plot(sample = sample_id[-which(sample_id=="Prostate")], otu_num = c(192,191,192,192,192,192,199,216), core_otu = core_num,
#             start = 90, a = 1, b = 2.1, r = 1.3, ellipse_col = ellipse_col[-3], circle_col = 'white',center = c(281))
# dev.off()
