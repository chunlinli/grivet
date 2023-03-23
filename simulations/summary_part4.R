viv <- rep(c(1:6),3)
group <- rep(1:3,each=6)
library(ggplot2)


### Discovery Part
path_used <- "./primary_results/part4/discovery/"
stats_hub1 <- read.csv(file = paste(path_used,"stats_hub1.csv",sep=""),header=TRUE)
stats_hub2 <- read.csv(file = paste(path_used,"stats_hub2.csv",sep=""),header=TRUE)
stats_hub3 <- read.csv(file = paste(path_used,"stats_hub3.csv",sep=""),header=TRUE)
stats_hub4 <- read.csv(file = paste(path_used,"stats_hub4.csv",sep=""),header=TRUE)
stats_hub5 <- read.csv(file = paste(path_used,"stats_hub5.csv",sep=""),header=TRUE)
stats_hub6 <- read.csv(file = paste(path_used,"stats_hub6.csv",sep=""),header=TRUE)

result.mat.discovery <- rbind(colMeans(stats_hub1),colMeans(stats_hub2),colMeans(stats_hub3),colMeans(stats_hub4),colMeans(stats_hub5),colMeans(stats_hub6))
rownames(result.mat.discovery) <- c("(H,C,n=1000)","(H,C,n=500)","(H,C,n=200)","(H,D,n=1000)","(H,D,n=500)","(H,D,n=200)")
colnames(result.mat.discovery) <- rep(c("FDR","JCI","TPR","SHD"),6)

## Plot for TPR, continuous
mat1 <- result.mat.discovery[1:3,(1:6)*4-1]
TPR <- as.numeric(t(mat1))
mat1 <- cbind(viv,TPR,group)
mat1 <- as.data.frame(mat1)
mat1$group <- as.factor(mat1$group)
p1<-ggplot(mat1, aes(x=viv, y=TPR, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  ylim(0,1)+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average true positive rate")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/ctpr.pdf",
  width = 6, 
  height = 5
)
p1
dev.off()

## Plot for TPR, discrete
mat1 <- result.mat.discovery[4:6,(1:6)*4-1]
TPR <- as.numeric(t(mat1))
mat1 <- cbind(viv,TPR,group)
mat1 <- as.data.frame(mat1)
mat1$group <- as.factor(mat1$group)
p2<-ggplot(mat1, aes(x=viv, y=TPR, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  ylim(0,1)+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average true positive rate")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/dtpr.pdf",
  width = 6, 
  height = 5
)
p2
dev.off()

## Plot for JCI, continuous
mat1 <- result.mat.discovery[1:3,(1:6)*4-2]
JCI <- as.numeric(t(mat1))
mat1 <- cbind(viv,JCI,group)
mat1 <- as.data.frame(mat1)
mat1$group <- as.factor(mat1$group)
p3<-ggplot(mat1, aes(x=viv, y=JCI, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  ylim(0,1)+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average Jaccard index")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/cjci.pdf",
  width = 6, 
  height = 5
)
p3
dev.off()

## Plot for JCI, discrete
mat1 <- result.mat.discovery[4:6,(1:6)*4-2]
JCI <- as.numeric(t(mat1))
mat1 <- cbind(viv,JCI,group)
mat1 <- as.data.frame(mat1)
mat1$group <- as.factor(mat1$group)
p4<-ggplot(mat1, aes(x=viv, y=JCI, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  ylim(0,1)+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average Jaccard index")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/djci.pdf",
  width = 6, 
  height = 5
)
p4
dev.off()


pdf("./summary_results/part4/structure_learning.pdf", width = 10, height = 10)
p2 <- p2 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p3 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p4 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
legend_p1 <- cowplot::get_legend(p1)
y_axis_p1 <- ggplot() + labs(y = "average true positive rate") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(r = 0, l = 5, unit = "pt")
  )
y_axis_p3 <- ggplot() + labs(y = "average Jaccard index") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(r = 0, l = 5, unit = "pt")
  )
x_axis_p1 <- ggplot() + labs(x = "the number of valid IVs") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 0, b = 5, unit = "pt")
  )
p1 <- p1 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())

grid_layout1 <- gridExtra::arrangeGrob(
  y_axis_p1, p1, p2,
  y_axis_p3, p3, p4,
  ncol = 3,
  nrow = 2,
  widths = c(0.2,2.7,2.7), heights = c(2.7,2.7)
)
grid_layout2 <- gridExtra::arrangeGrob(
  grid_layout1,
  x_axis_p1,
  ncol = 1,
  nrow = 2,
  widths = c(2.7), heights = c(2.7,0.1)
) 
gridExtra::grid.arrange(grid_layout2, legend_p1, nrow = 1, widths = c(2.7, 0.3))
dev.off()



### Estimation Part
path_used <- "./primary_results/part4/estimation/"
stats_hub1 <- read.csv(file = paste(path_used,"stats_hub1.csv",sep=""),header=TRUE)
stats_hub2 <- read.csv(file = paste(path_used,"stats_hub2.csv",sep=""),header=TRUE)
stats_hub3 <- read.csv(file = paste(path_used,"stats_hub3.csv",sep=""),header=TRUE)
stats_hub4 <- read.csv(file = paste(path_used,"stats_hub4.csv",sep=""),header=TRUE)
stats_hub5 <- read.csv(file = paste(path_used,"stats_hub5.csv",sep=""),header=TRUE)
stats_hub6 <- read.csv(file = paste(path_used,"stats_hub6.csv",sep=""),header=TRUE)

result.mat.estimation <- rbind(colMeans(stats_hub1),colMeans(stats_hub2),colMeans(stats_hub3),colMeans(stats_hub4),colMeans(stats_hub5),colMeans(stats_hub6))
rownames(result.mat.estimation) <- c("(H,C,n=1000)","(H,C,n=500)","(H,C,n=200)","(H,D,n=1000)","(H,D,n=500)","(H,D,n=200)")
colnames(result.mat.estimation) <- rep(c("MaaD","MeaD","MesD"),6)

## Plot for maad, continuous
mat1 <- result.mat.estimation[1:3,(1:6)*3-2]
MaaD <- as.numeric(t(mat1))
mat1 <- cbind(viv,MaaD,group)
mat1 <- as.data.frame(mat1)
mat1$group <- as.factor(mat1$group)
p1<-ggplot(mat1, aes(x=viv, y=MaaD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average maximum absolute deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/cmaad.pdf",
  width = 6, 
  height = 5
)
p1
dev.off()

## Plot for maad, discrete
mat4 <- result.mat.estimation[4:6,(1:6)*3-2]
MaaD <- as.numeric(t(mat4))
mat4 <- cbind(viv,MaaD,group)
mat4 <- as.data.frame(mat4)
mat4$group <- as.factor(mat4$group)
p2<-ggplot(mat4, aes(x=viv, y=MaaD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average maximum absolute deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/dmaad.pdf",
  width = 6, 
  height = 5
)
p2
dev.off()

## Plot for mesd, continuous
mat3 <- result.mat.estimation[1:3,(1:6)*3]
MeSD <- as.numeric(t(mat3))
mat3 <- cbind(viv,MeSD,group)
mat3 <- as.data.frame(mat3)
mat3$group <- as.factor(mat3$group)
p3<-ggplot(mat3, aes(x=viv, y=MeSD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average mean squared deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/cmesd.pdf",
  width = 6, 
  height = 5
)
p3
dev.off()

## Plot for mesd, discrete
mat6 <- result.mat.estimation[4:6,(1:6)*3]
MeSD <- as.numeric(t(mat6))
mat6<- cbind(viv,MeSD,group)
mat6 <- as.data.frame(mat6)
mat6$group <- as.factor(mat6$group)
p4<-ggplot(mat6, aes(x=viv, y=MeSD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average mean squared deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 16))
pdf(
  file ="./summary_results/part4/dmesd.pdf",
  width = 6, 
  height = 5
)
p4
dev.off()

pdf("./summary_results/part4/coefficient_estimation.pdf", width = 10, height = 10)
p2 <- p2 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
p3 <- p3 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
p4 <- p4 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())
legend_p1 <- cowplot::get_legend(p1)
y_axis_p1 <- ggplot() + labs(y = "average maximum absolute deviation") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(r = 0, l = 5, unit = "pt")
  )
y_axis_p3 <- ggplot() + labs(y = "average mean squared deviation") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    plot.margin = margin(r = 0, l = 5, unit = "pt")
  )
x_axis_p1 <- ggplot() + labs(x = "the number of valid IVs") +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 20),
    plot.background = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(t = 0, b = 5, unit = "pt")
  )
p1 <- p1 + theme(legend.position = "none",axis.title.x = element_blank(), axis.title.y = element_blank())

grid_layout1 <- gridExtra::arrangeGrob(
  y_axis_p1, p1, p2,
  y_axis_p3, p3, p4,
  ncol = 3,
  nrow = 2,
  widths = c(0.2,2.7,2.7), heights = c(2.7,2.7)
)
grid_layout2 <- gridExtra::arrangeGrob(
  grid_layout1,
  x_axis_p1,
  ncol = 1,
  nrow = 2,
  widths = c(2.7), heights = c(2.7,0.1)
) 
gridExtra::grid.arrange(grid_layout2, legend_p1, nrow = 1, widths = c(2.7, 0.3))
dev.off()
