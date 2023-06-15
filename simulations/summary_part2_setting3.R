## This file reads the preliminary simulation results and summarize the evaluations of parameter estimation in table 3 of the supplementary materials.
path_used <- "./primary_results/part2/setting3/"

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
colnames(result.mat) <- c("Max abs diff(Direct)","Mean abs diff(Direct)","Mean sq diff(Direct)","Max abs diff(Proposed)","Mean abs diff(Proposed)","Mean sq diff(Proposed)")

print(result.mat)
write.csv(result.mat,file.path("./summary_results/part2/setting3/","result.csv"),row.names = TRUE)