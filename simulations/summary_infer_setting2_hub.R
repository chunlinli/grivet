path_used <- "./primary_results/part3/setting2/"

stats.list <- vector("list",6)
for (k in 1:6){
  stats.list[[k]] <- read.csv(file = paste(path_used,"stats_hub",k,".csv",sep=""),header=TRUE)
}

library(ggplot2)

dfs <- c(1,3,5,1,3,5)
ns <- rep(c(1000,500,200),2)
result.mat <- matrix(0,6,6)
for (i in 1:6){
  cate <- "hub graph,|D|="
  if (i <= 3){
    inter <- ",continuous X"
  }else{
    inter <- ",discrete X"
  }
  for (j in 1:6){
    mat <- stats.list[[i]]
    result.mat[i,j] <- sum(mat[,j]>qchisq(0.95,df=dfs[j]))/sum(mat[,j]>=0)
    if (j < 4){
      pdf(
        file =paste("./summary_results/part3/setting2/figures/",i,"_",j,".pdf",sep=""), # The file name you want to save the plot in
        width = 5.5, # The width & height of the plot in inches
        height = 6.5
      )
      plt <- ggplot(data = data.frame(mat[,j]),aes(mat[,j])) +
        xlim(0,20) +
        geom_histogram(aes(y = after_stat(density)),binwidth = 0.2,boundary = 0,color="#e9ecef")+
        stat_function(fun = dchisq,args = list(df=dfs[j]),colour = "deeppink")+
        labs(title = paste(cate,dfs[j],",m=200,n=",ns[i],inter,sep=""),x = "test statistic")
      print(plt)
      dev.off()
    }
  }
}

colnames(result.mat) <- c("size(|D|=1)","size(|D|=3)","size(|D|=5)","power(|D|=1)","power(|D|=3)","power(|D|=5)")
rownames(result.mat) <- c("(H,C,m=200,n=1000)","(H,C,m=200,n=500)","(H,C,m=200,n=200)","(H,D,m=200,n=1000)","(H,D,m=200,n=500)","(H,D,m=200,n=200)")
print(result.mat)

write.csv(result.mat,file.path("./summary_results/part3/setting2/","hub.csv"),row.names = TRUE)

