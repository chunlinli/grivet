path_used <- "./primary_results/part5/"
stats_random.list <- vector("list",6)
for (k in 1:6){
  stats_random.list[[k]] <- read.csv(file = paste(path_used,"stats_random",k,".csv",sep=""),header=TRUE)
}

library(ggplot2)

## Process empirical rejection probabilities under H0
emp <- numeric(300)
for(k in 1:6){
  mat <- stats_random.list[[k]]
  for (j in 1:50){
    stats <- mat[mat[,j]==0,j+50]
    emp[50*(k-1)+j] <- sum(stats>qchisq(0.95,df=1))/length(stats) 
  }
}
n <- rep(rep(c("500","400","300"),50),2)
intervention <- rep(c("continuous","discrete"),each=150)
emp.h0 <-data.frame(emp,intervention,n)

cols <- c("#E64B35B2","#4DBBD5B2","#00A087B2")
p1 <- ggplot(emp.h0, aes(x=n, y=emp, fill = n)) + 
  geom_boxplot()+
  facet_wrap(~intervention)+
  scale_fill_manual(values = cols)+
  ylab("empirical rejection probabilities")+
  theme_bw()
pdf(
  file ="./summary_results/part5/rejh0.pdf",
  width = 6, 
  height = 5
)
p1
dev.off()

## Process empirical rejection probabilities under Ha
emp <- numeric(300)
for(k in 1:6){
  mat <- stats_random.list[[k]]
  for (j in 1:50){
    stats <- mat[mat[,j]==1,j+50]
    emp[50*(k-1)+j] <- sum(stats>qchisq(0.95,df=1))/length(stats) 
  }
}
n <- rep(rep(c("500","400","300"),50),2)
intervention <- rep(c("continuous","discrete"),each=150)
emp.ha <-data.frame(emp,intervention,n)

cols <- c("#E64B35B2","#4DBBD5B2","#00A087B2")
p2 <- ggplot(emp.ha, aes(x=n, y=emp, fill = n)) + 
  geom_boxplot()+
  facet_wrap(~intervention)+
  scale_fill_manual(values = cols)+
  ylab("empirical rejection probabilities")+
  theme_bw()
pdf(
  file ="./summary_results/part5/rejha.pdf",
  width = 6, 
  height = 5
)
p2
dev.off()

