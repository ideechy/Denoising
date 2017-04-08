library(EBImage)
setwd('C:/Users/Billy/Desktop/Courses/ST790-Advanced Computing/Project/Haoyu Chen')
im <-  readImage("NCState.png")
im <- readImage('parrot.jpg')
display(im)
# imBlackNoise <-  GenerateNoise(im, p=0.2)
# display(imBlackNoise)
# display(Denoise(imBlackNoise, 10, p=0.2))
# display(Denoise(imBlackNoise, 20, p=0.2))
# display(Denoise(imBlackNoise, 40, p=0.2))
# display(Denoise(imBlackNoise, 20, p=0.5))
# display(Denoise(imBlackNoise, 20, p=0.9))
imColorNoise <-  GenerateNoise(im, p=0.2, type=1)
display(imColorNoise)
# display(Denoise(imColorNoise, 10, p=0.2))
# display(Denoise(imColorNoise, 20, p=0.2))
# display(Denoise(imColorNoise, 20, p=0.5))
# display(Denoise(imColorNoise, 20, p=0.9))
# display(DenoiseUnif(imColorNoise,40, p=0.4))
X1=(Denoise(imColorNoise, 40, p=0.2))
X2=(DenoiseUnif(imColorNoise, 40, p=0.2))
# X3=(Denoise(imColorNoise, 40, p=0))
# res=c()
# 
# for (r in 5*c(1:10)){
#   X = DenoiseUnif(imColorNoise, r, p=0.4)
#   res = c(res,sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X[,,i],"F")})))
# }

# sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X2[,,i],"F")}))
# sum(sapply(1:4,FUN=function(i){norm(im[,,i]-X3[,,i],"F")}))
# sum(sapply(1:4,FUN=function(i){norm(im[,,i]-imColorNoise[,,i],"F")}))


Dist <- function(diff){
  # Compute the Forbenious distance between images
  # diff is the difference of two images
  sum(sapply(1:dim(diff)[3],FUN=function(i){norm(diff[ , , 3],"F")^2}))
}

ExpandRow <- function(image){
  # line up the layers vertically
  expand <- image[ , , 1]
  for ( i in 2:dim(image)[3]){
    expand <- rbind(expand, image[, , i])
  }
  return(expand)
}

ExpandCol <- function(image){
  # line up the layers horizontally
  expand <- image[ , , 1]
  for ( i in 2:dim(image)[3]){
    expand <- cbind(expand, image[, , i])
  }
  return(expand)
}

FoldRow <- function(image, layer = 4){
  # fold the layers vertically
  m <- ncol(image)
  n <- nrow(image)
  fold <- array(dim = m*n)
  dim(fold) <- c(n/layer, m, layer) 
  for (i in 1:layer){
    fold[, , i] <- image[((i-1)*n/layer+1):(i*n/layer), ]
  }
  fold <- as.Image(fold)
  colorMode(fold) <- 2
  return(fold)
}

FoldCol <- function(image, layer = 4){
  # fold the layers horizontally
  m <- ncol(image)
  n <- nrow(image)
  fold <- array(dim = m*n)
  dim(fold) <- c(n, m/layer, layer) 
  for (i in 1:layer){
    fold[, , i] <- image[, ((i-1)*m/layer+1):(i*m/layer)]
  }
  fold <- as.Image(fold)
  colorMode(fold) <- 2
  return(fold)
}

X1 <- Denoise(imColorNoise, 40, p=0.2)

X2 <- DenoiseUnif(imColorNoise, 40, p=0.2)

X1_EFRow <- FoldRow( Denoise( ExpandRow( imColorNoise ), 40, p=0.2 ), 3 ) 

X1_EFCol <- FoldCol( Denoise( ExpandCol( imColorNoise ), 40, p=0.2 ), 3) 

X2_EFRow <- FoldRow( DenoiseUnif( ExpandRow( imColorNoise ), 40, p=0.2 ), 3)

X2_EFCol <- FoldCol( DenoiseUnif( ExpandCol( imColorNoise ), 40, p=0.2 ), 3)

display(X1)
display(imColorNoise)
display(X2_EFRow)
display(abs(im-X2_EFRow))
display(abs(im-imColorNoise))

Dist(im-imColorNoise)
Dist(im-X1)
Dist(im-X2)
Dist(im-X1_EFRow)
Dist(im-X1_EFCol)
Dist(im-X2_EFRow)
Dist(im-X2_EFCol)


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
  # X is an input m*n matrix
  # rank is the desired rank of the output matrix
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
  # X is an input m*n matrix
  # rank is the desired rank of the output matrix
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

