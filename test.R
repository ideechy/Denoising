library(Denoising)

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
## POPCA MSE-rank (time-rank)
O <- resize(im, 64, 64)
rank <- qr(O)$rank
X <- GenerateNoise(O, 0.3, 1)
MSE <- numeric(64); time <- numeric(64)
for (r in (1:64)){
    time[r] <- system.time(A <- POPCAUnif(X, r, 0.3))[1]
    MSE[r] <- norm(A-O, "F")
}
figure1 <- list(MSE=MSE, time=time, rank=1:64)
plot(figure1$rank, figure1$MSE); plot(figure1$rank, figure1$time)
## POPCA MSE-p
O <- resize(im, 64, 64)
X <- GenerateNoise(O, 0.3, 1)
maxMSE <- norm(X-O, "F")
MSE <- numeric(19); time <- numeric(19)
p <- c((1:10)*0.02,(1:10)*0.05+0.2)
for (i in 1:20){
    time[i] <- system.time(A <- POPCAUnif(X, 11, p[i]))[1]
    MSE[i] <- norm(A-O, "F")
}
figure2 <- list(MSE=MSE, time=time, p=p)
plot(figure2$p, figure2$MSE, type='o'); plot(figure2$p, figure2$time)

