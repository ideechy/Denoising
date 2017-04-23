library(Denoising)
set.seed(123)
#test with rank 1 matrix
X=matrix(1,10,10)
X=GenerateNoise(X, 0.2, 1)

par(mfrow=c(3,1))
contour(X-1, main="origin")
contour(RPCA(X)$E, main="RPCA")
contour(X-POPCAUnif(X,1,0.2), main="POPCA")


# figures initialize
library(EBImage)
im <-  readImage("NCState.png")
im <- channel(im, "gray")

## time-image size
size <- c(16,32,64,128,256,512)
maxMSE <- numeric(length(size))
RPCAfig <- POPCAfig1 <- POPCAfig2 <- list()
for (i in 1:length(size)){
    O <- resize(im, size[i], size[i])
    X <- GenerateNoise(O, 0.3, 1)
    maxMSE[i] <- norm(X-O, "F")
    #RPCA
    tol <- c(2^(7:0)/1000,10^-(4:8))
    MSE <- numeric(length(tol))
    time <- numeric(length(tol))
    for (j in 1:length(tol)){
        time[j] <- system.time(A <- RPCA(X, tol=tol[j])$A)[1]
        MSE[j] <- norm(A-O, "F")
    }
    RPCAfig[[i]] <- list(MSE=MSE, time=time, tol=tol, maxMSE=maxMSE)
    #POPCA rank
    MSE <- time <- numeric(size[i])
    for (j in (1:size[i])){
        time[j] <- system.time(A <- POPCAUnif(X, j, 0.3))[1]
        MSE[j] <- norm(A-O, "F")
    }
    POPCAfig1[[i]] <- list(MSE=MSE, time=time, rank=1:size[i], maxMSE=maxMSE)
    optrank <- which.min(MSE)
    #POPCA p
    p <- c((1:10)*0.02,(1:10)*0.05+0.2)
    MSE <- time <- numeric(length(p))
    for (j in 1:length(p)){
        time[j] <- system.time(A <- POPCAUnif(X, optrank, p[j]))[1]
        MSE[j] <- norm(A-O, "F")
    }
    POPCAfig2[[i]] <- list(MSE=MSE, time=time, p=p, maxMSE=maxMSE, rank=optrank)
}

O <- resize(im, 64, 64)
X <- GenerateNoise(O, 0.3, 1)
p <- c((1:10)*0.02,(1:10)*0.05+0.2)
rank <- 1:64
MSE <- matrix(nrow=length(p), ncol=64)
for (i in 1:length(p)){
    for (j in 1:64){
        A <- POPCAUnif(X, j, p[i])
        MSE[i,j] <- norm(A-O, "F")
    }
}

#plot
library(ggplot2)
library(reshape2)
library(RColorBrewer)
dp1 <- c()
for (i in 2:6){
    d <- data.frame(rank=POPCAfig1[[i]]$rank,
                    MSE=POPCAfig1[[i]]$MSE,
                    noiseMSE=POPCAfig1[[i]]$maxMSE[i],
                    size=size[i])
    dp1 <- rbind(dp1,d)
}
dp1$size <- as.factor(dp1$size)

qplot(rank, MSE/noiseMSE, data=dp1, geom="line", colour=size,
      main="Choice of rank in POPCA (p=optimal)")
ggsave("figure/popca_mse_to_rank.eps",w=4,h=3)
ggsave("figure/popca_mse_to_png.eps",w=4,h=3)

dp2 <- c()
for (i in 2:6){
    d <- data.frame(p=POPCAfig2[[i]]$p,
                    MSE=POPCAfig2[[i]]$MSE,
                    noiseMSE=POPCAfig2[[i]]$maxMSE[i],
                    size=size[i])
    dp2 <- rbind(dp2,d)
}
dp2$size <- as.factor(dp2$size)

qplot(p, MSE, data=dp2, geom="line", colour=size,
      main="Choice of p in POPCA (rank=optimal)")
ggsave("figure/popca_mse_to_p.eps",w=4,h=3)
ggsave("figure/popca_mse_to_p.png",w=4,h=3)

dp3 <- melt(MSE, value.name="MSE")
dp3$Var1=p[dp3$Var1]
ggplot(dp3,aes(x=Var1,y=Var2, fill=MSE))+geom_tile()+
    ggtitle("Simultaneous choice of p and rank (size=64)")+
    xlab("p")+ylab("rank")
ggsave("figure/popca_rank_to_p.eps",w=4,h=3)
ggsave("figure/popca_rank_to_p.png",w=4,h=3)

dp4 <- c()
for (i in 2:6){
    d <- data.frame(time=RPCAfig[[i]]$time,
                    MSE=RPCAfig[[i]]$MSE,
                    noiseMSE=RPCAfig[[i]]$maxMSE[i],
                    size=size[i])
    dp4 <- rbind(dp4,d)
}
dp4$size <- as.factor(dp4$size)

qplot(log(time), MSE/noiseMSE, data=dp4, geom="line", colour=size,
      main="Convergence speed of RPCA")
ggsave("figure/rpca_mse_to_time.eps",w=4,h=3)
ggsave("figure/rpca_mse_to_time.png",w=4,h=3)

dp5 <- data.frame(size, RMSE=0, POMSE=0, Rtime=0, POtime=0,
                  noiseMSE=RPCAfig[[6]]$maxMSE)
for (i in 1:6){
    dp5$RMSE[i] <- min(RPCAfig[[i]]$MSE)
    dp5$Rtime[i] <- max(RPCAfig[[i]]$time)
    #dp5$Rtime[i] <- RPCAfig[[i]]$time[min(which(RPCAfig[[i]]$MSE<dp5$RMSE[i]*1.01))]
    #dp5$RMSE[i] <- RPCAfig[[i]]$MSE[min(which(RPCAfig[[i]]$time>=dp5$Rtime[i]))]
    dp5$POMSE[i] <- min(POPCAfig2[[i]]$MSE)
    dp5$POtime[i] <- POPCAfig2[[i]]$time[which.min(POPCAfig2[[i]]$MSE)]
}
dp5=data.frame(size, noiseMSE=dp5$noiseMSE,
            MSE=c(dp5$RMSE, dp5$POMSE),
            time=c(dp5$Rtime, dp5$POtime),
            method=as.factor(c(rep("RPCA",6),rep("POPCA",6))))

qplot(size, MSE, data=dp5, geom=c("point", "line"), colour=method)
ggsave("figure/compare_mse_to_size.eps",w=4,h=3)
ggsave("figure/compare_mse_to_size.png",w=4,h=3)
qplot(size, time, data=dp5, geom=c("point", "line"), colour=method)
ggsave("figure/compare_time_to_size.eps",w=4,h=3)
ggsave("figure/compare_time_to_size.png",w=4,h=3)
