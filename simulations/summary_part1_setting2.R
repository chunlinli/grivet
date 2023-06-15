## This file reads the preliminary simulation results and summarize the evaluations of performance in table 1 of the supplementary materials.
path_used <- "./primary_results/part1/setting2/"

stats_hub1 <- read.csv(file = paste(path_used,"stats_hub1.csv",sep=""),header=TRUE)
stats_hub2 <- read.csv(file = paste(path_used,"stats_hub2.csv",sep=""),header=TRUE)
stats_hub3 <- read.csv(file = paste(path_used,"stats_hub3.csv",sep=""),header=TRUE)
stats_hub4 <- read.csv(file = paste(path_used,"stats_hub4.csv",sep=""),header=TRUE)
stats_hub5 <- read.csv(file = paste(path_used,"stats_hub5.csv",sep=""),header=TRUE)
stats_hub6 <- read.csv(file = paste(path_used,"stats_hub6.csv",sep=""),header=TRUE)
stats_random1 <- read.csv(file = paste(path_used,"stats_random1.csv",sep=""),header=TRUE)
stats_random2 <- read.csv(file = paste(path_used,"stats_random2.csv",sep=""),header=TRUE)
stats_random3 <- read.csv(file = paste(path_used,"stats_random3.csv",sep=""),header=TRUE)
stats_random4 <- read.csv(file = paste(path_used,"stats_random4.csv",sep=""),header=TRUE)
stats_random5 <- read.csv(file = paste(path_used,"stats_random5.csv",sep=""),header=TRUE)
stats_random6 <- read.csv(file = paste(path_used,"stats_random6.csv",sep=""),header=TRUE)

result.mat <- rbind(colMeans(stats_hub1),colMeans(stats_hub2),colMeans(stats_hub3),colMeans(stats_hub4),colMeans(stats_hub5),colMeans(stats_hub6),
                    colMeans(stats_random1),colMeans(stats_random2),colMeans(stats_random3),colMeans(stats_random4),colMeans(stats_random5),colMeans(stats_random6))
rownames(result.mat) <- c("(H,C,n=1000)","(H,C,n=500)","(H,C,n=200)","(H,D,n=1000)","(H,D,n=500)","(H,D,n=200)","(R,C,n=1000)","(R,C,n=500)","(R,C,n=200)","(R,D,n=1000)","(R,D,n=500)","(R,D,n=200)")
colnames(result.mat) <- c("FDR(RFCI)","JI(RFCI)","TPR(RFCI)","SHD(RFCI)","FDR(LrPS+GES)","JI(LrPS+GES)","TPR(LrPS+GES)","SHD(LrPS+GES)","FDR(proposed)","JI(proposed)","TPR(proposed)","SHD(proposed)")
result.mat1 <- result.mat[,1:4]
result.mat2 <- result.mat[,5:8]
result.mat3 <- result.mat[,9:12]

print(result.mat1)
write.csv(result.mat1,file.path("./summary_results/part1/setting2/","RFCI.csv"),row.names = TRUE)

print(result.mat2)
write.csv(result.mat2,file.path("./summary_results/part1/setting2/","lrps_ges.csv"),row.names = TRUE)

print(result.mat3)
write.csv(result.mat3,file.path("./summary_results/part1/setting2/","grivet.csv"),row.names = TRUE)
