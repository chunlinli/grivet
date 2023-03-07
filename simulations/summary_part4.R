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

print(result.mat.discovery)

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
viv <- rep(c(1:6),3)
group <- rep(1:3,each=6)
library(ggplot2)

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
  theme_bw()
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
p4<-ggplot(mat4, aes(x=viv, y=MaaD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average maximum absolute deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()
pdf(
  file ="./summary_results/part4/dmaad.pdf",
  width = 6, 
  height = 5
)
p4
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
  theme_bw()
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
p6<-ggplot(mat6, aes(x=viv, y=MeSD, group=group)) +
  geom_line(aes(color=group))+
  geom_point(aes(color=group))+
  scale_color_manual(labels = c("n=1000", "n=500", "n=200"),values = c("#E64B35B2","#4DBBD5B2","#00A087B2"))+
  xlab("the number of valid IVs")+
  ylab("average mean squared deviation")+
  scale_x_continuous(breaks=1:6)+
  theme_bw()
pdf(
  file ="./summary_results/part4/dmesd.pdf",
  width = 6, 
  height = 5
)
p6
dev.off()
