if (!requireNamespace("RColorBrewer", quietly = TRUE)) 
  install.packages("RColorBrewer",dependencies = TRUE)
library(RColorBrewer)

??RColorBrewer

scale_color_manual(values = c("part1" = rgb(0, 0, 1, alpha = 0.5),   # 50% 透明度的蓝色
                              "part2" = rgb(70, 130, 180, alpha = 0.5),  # 50% 透明度的steelblue
                              "part3" = rgb(0, 1, 0, alpha = 0.5),  # 50% 透明度的绿色
                              "part4" = rgb(1, 0, 0, alpha = 0.5),  # 50% 透明度的红色
                              "other" = rgb(0.5, 0.5, 0.5, alpha = 0.5)))  # 50% 透明度的灰色




## create a sequential palette for usage and show colors(这相当于一个向量)
mypalette<-brewer.pal(9,"Greens")
image(1:9,1,as.matrix(1:9),col=mypalette,xlab="Greens (sequential)",
      ylab="",xaxt="n",yaxt="n",bty="n")

mypalette <- brewer.pal(9, "Greens")
selected_colors <- mypalette[1:2]


## display a divergent palette
display.brewer.pal(7,"BrBG")
devAskNewPage(ask=TRUE)

## display a qualitative palette
display.brewer.pal(7,"Accent")
devAskNewPage(ask=TRUE)

## display a palettes simultanoeusly
display.brewer.all(n=10)
devAskNewPage(ask=TRUE)

display.brewer.all()
devAskNewPage(ask=TRUE)

display.brewer.all(type="div")
devAskNewPage(ask=TRUE)

display.brewer.all(type="seq")
devAskNewPage(ask=TRUE)

display.brewer.all(type="qual") 
devAskNewPage(ask=TRUE)

display.brewer.all(n=5,type="div",exact.n=TRUE)
devAskNewPage(ask=TRUE)

display.brewer.all(colorblindFriendly=TRUE)
devAskNewPage(ask=TRUE)

