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
mypalette <- colorRampPalette(c("red","white", "blue"))(9)

scale_colour_gradient2(low = mypalette[1], mid = adjustcolor("grey", alpha.f = 0.5), high = mypalette[9], midpoint = 0.00010035, guide = "colourbar")











# 调用方式
colors()
par(mar=c(1,1,1,1))  # 设置边距为1
n <- length(colors())
par(mfrow=c(25, 10))  # 将图形窗口分为25行和10列
for (i in 1:n) {
  plot(0, 0, col = colors()[i], pch = 19, cex = 3, xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.5, labels = colors()[i], col = "black", cex = 0.8)
}


# 透明度可以改变两个渐变色的对比度一个0.9，一个0.
scale_fill_manual(values = c("#EDB749","#3CB2EC","#9C8D58"))
c("palevioletred", "peachpuff", "lightgreen","lightblue", "plum")
values = c("turquoise", "salmon", "mediumpurple1"
mypalette=c(rgb(1, 0, 0, 0.9), rgb(1, 1, 1, 0))
mypalette <- colorRampPalette(c("red","white", "blue"))(9)

color = colorRampPalette(c("steelblue", "white", "salmon"))(100),

mypalette <- brewer.pal(4,"BrBG")[1:4]
# 使用 alpha() 函数调节颜色的透明度
my_color <- alpha(brewer.pal(4, "Set1"), 0.6)
my_color <- alpha(brewer.pal(4, "Set2"),0.8)
my_color <- alpha(brewer.pal(4, "Set3"), 1)
my_color <- alpha(brewer.pal(4, "Pastel2"), 1)
my_color <- alpha(brewer.pal(4, "Pastel1"), 1)

my_color <- alpha(brewer.pal(4, "Paired"), 0.8)
my_color <- alpha(brewer.pal(4, "Accent"),0.8)
my_color <- alpha(brewer.pal(4, "Dark2"), 0.7)



## 展示
image(1:9,100,as.matrix(1:9),col=my_color,xlab="Greens (sequential)",
      ylab="",xaxt="n",yaxt="n",bty="n")




## display a divergent palette
display.brewer.pal(4,"BrBG")
devAskNewPage(ask=TRUE)

## display a qualitative palette
display.brewer.pal(4,"Accent")
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

