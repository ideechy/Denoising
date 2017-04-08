rm(list=ls())
library(EBImage)
im <-  readImage("NCState.png")
display(im)
imBlackNoise <-  GenerateNoise(im, p=0.2)
display(imBlackNoise)
display(Denoise(imBlackNoise, 10, p=0.2))
display(Denoise(imBlackNoise, 20, p=0.2))
display(Denoise(imBlackNoise, 40, p=0.2))
display(Denoise(imBlackNoise, 20, p=0.5))
display(Denoise(imBlackNoise, 20, p=0.9))
imColorNoise <-  GenerateNoise(im, p=0.4, type=1)
display(imColorNoise)
display(Denoise(imColorNoise, 10, p=0.2))
display(Denoise(imColorNoise, 20, p=0.2))
display(Denoise(imColorNoise, 20, p=0.5))
display(Denoise(imColorNoise, 20, p=0.9))
display(DenoiseUnif(imColorNoise,40, p=0.4))
X1=(Denoise(imColorNoise, 40, p=0.4))
X2=(DenoiseUnif(imColorNoise, 40, p=0.4))
X3=(Denoise(imColorNoise, 40, p=0))
res=c()
for (r in c(c(1:18),5*c(4:10))){
    X = DenoiseUnif(imColorNoise, r, p=0.4)
    res = c(res,sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X[,,i],"F")^2})))
}

sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X[,,i],"F")^2}))

sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X2[,,i],"F")}))
sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X3[,,i],"F")}))
dash = sqrt(sum(sapply(1:4,FUN=function(i){norm(im[,,i]-imColorNoise[,,i],"F")^2})))
plot(sqrt(res[11:35]), type='l', ylim=c(200,320), xlab='rank', ylab='Frobenius norm')
abline(h=309,lty=2,col='red')

DetectNoise <- function(X){
    # X is an input m*n matrix
    m <- nrow(X)
    n <- ncol(X)
    LeftRight <- cbind(0,0,X) + cbind(X,0,0)
    LeftRight <- LeftRight[,2:(n+1)]
    UpDown <- rbind(0,0,X) + rbind(X,0,0)
    UpDown <- UpDown[2:(m+1),]
    denominator <- matrix(rep(4,m*n),m,n)
    denominator[c(1,m),] <- 3; denominator[,c(1,n)] <- 3;
    denominator[c(1,m),c(1,n)] <- 2
    PixelAvg <- (LeftRight + UpDown)/denominator
    Diff <- (X - PixelAvg)^2
}

ReduceRank <- function(X, rank, p=0){
    # X is an input m*n matrix
    # rank is the desired rank of the output matrix
    # p is the corruption probability of every pixel
    #m <- nrow(X)
    n <- ncol(X)
    # o is the probability of observe a uncorrupted pixel
    o <- 1-p
    # Calculate the covariance matrix
    XtX = X %*% t(X)
    C <- 1/n * (1/o^2 * XtX + (1/o-1/o^2) * diag(diag(XtX)))
    # Get eigenvectors correspond to the biggest #=rank eigenvalues
    V <- eigen(C)$vector[,1:rank]
    # Project X to the column space of V
    P <- V %*% t(V) %*% X
    return(P)
}
ReduceRankUnif <- function(X, rank, p=0){
    # X is an input m*n matrix
    # rank is the desired rank of the output matrix
    # p is the corruption probability of every pixel
    #m <- nrow(X)
    n <- ncol(X)
    # o is the probability of observe a uncorrupted pixel
    o <- 1-p
    # Calculate the covariance matrix
    Cdiag = (diag(X %*% t(X)) - n*p/3)/o
    Coffdiag = 1/o^2 * (X-p/2) %*% t(X-p/2)
    C <- 1/n * (Coffdiag - diag(diag(Coffdiag)) + diag(Cdiag))
    # Get eigenvectors correspond to the biggest #=rank eigenvalues
    V <- eigen(C)$vector[,1:rank]
    # Project X to the column space of V
    P <- V %*% t(V) %*% X
    return(P)
}
Denoise <- function(X, rank, p=0){
    # X is an input m*n*d array, preferably of class "Image"
    # rank is the desired rank of the output array
    # p is the corruption probability of every pixel
    if (colorMode(X) == 0){
        return(ReduceRank(X, rank, p))
    } else {
        for (d in 1:dim(X)[3]){
            X[,,d] <- ReduceRank(X[,,d], rank, p)
        }
        return(X)
    }
}
DenoiseUnif <- function(X, rank, p=0){
    # X is an input m*n*d array, preferably of class "Image"
    # rank is the desired rank of the output array
    # p is the corruption probability of every pixel
    if (colorMode(X) == 0){
        return(ReduceRankUnif(X, rank, p))
    } else {
        for (d in 1:dim(X)[3]){
            X[,,d] <- ReduceRankUnif(X[,,d], rank, p)
        }
        return(X)
    }
}
GenerateNoise <- function(X, p, type=0){
    # X is an input m*n matrix
    # p is the corruption probability of every pixel
    # type is the type of noise, 0=white and 1=color
    m <- nrow(X)
    n <- ncol(X)
    NoiseLocation <- which(rbinom(m*n, 1, p)==1)
    NoiseNumber <- length(NoiseLocation)
    if (colorMode(X) == 0){
        X[NoiseLocation] <- runif(n=NoiseNumber, max=type)
    } else {
        for (d in 1:dim(X)[3]){
            X[NoiseLocation+(d-1)*m*n] <- runif(n=NoiseNumber, max=type)
        }
    }
    return(X)
}
