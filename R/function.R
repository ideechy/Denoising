positive <- function(x){ x * (x > 0) }

#' Robust PCA
#'
#' @param X handle to function that returns objective function values
#' @param lamnda tuning parameter
#' @param tol stopping tolerance, default 1e-7
#'
#' @export
RPCA <- function(X, lambda=1/sqrt(nrow(X)), tol=1e-7){
    dim <- dim(X)
    Akk <- Ak <- Ekk <- Ek <- matrix(0,dim[1], dim[2])
    tkk <- tk <- 1
    mu <- 0.99 * norm(X,"2")
    mumin <- 1e-5 * mu
    change <- 1
    while (change > 0){
        A <- Ak + (tkk-1)/tk * (Ak - Akk)
        E <- Ek + (tkk-1)/tk * (Ek - Ekk)
        YA <- A - (A + E - X)/2
        YE <- E - (A + E - X)/2
        svd <- svd(YA)
        Akk <- Ak
        Ak <- svd$u %*% diag(positive(svd$d - mu/2)) %*% t(svd$v)
        Ekk <- Ek
        Ek <- sign(YE) * positive(abs(YE) - lambda * mu/2)
        tkk <- tk
        tk <- (1 + sqrt(1 + 4 * tk^2))/2
        mu <- max(0.9 * mu, mumin)
        change <- norm(A - Ak + Ek - E, "F") -
            sqrt(2) * max(1, norm(cbind(Ak, Ek), "F")) * tol
    }
    return(list(A=Ak, E=Ek))
}

#' Partially Observed PCA
#'
#' @param X an input m*n matrix
#' @param rank the desired rank of the output matrix
#' @param p the corruption probability of every pixel, if p=0 then it is PCA.
#'
#' @export
POPCA <- function(X, rank, p=0){
    #m <- nrow(X)
    n <- ncol(X)
    # o is the probability of observe a uncorrupted pixel
    o <- 1-p
    center <- rowMeans(X)
    X <- X-center
    # Calculate the covariance matrix
    XtX = X %*% t(X)
    C <- 1/n * (1/o^2 * XtX + (1/o-1/o^2) * diag(diag(XtX)))
    # Get eigenvectors correspond to the biggest #=rank eigenvalues
    V <- eigen(C)$vector[,1:rank]
    # Project X to the column space of V
    P <- V %*% t(V) %*% X + center
    return(P)
}

#' Partially Observed PCA with Uniform noise
#'
#' @param X an input m*n matrix
#' @param rank the desired rank of the output matrix
#' @param p the corruption probability of every pixel, if p=0 then it is PCA.
#'
#' @export
POPCAUnif <- function(X, rank, p=0){
    #m <- nrow(X)
    n <- ncol(X)
    # o is the probability of observe a uncorrupted pixel
    o <- 1-p
    center <- rowMeans(X)
    X <- X-center
    # Calculate the covariance matrix
    Cdiag = (diag(X %*% t(X)) - n*p/3)/o
    Coffdiag = 1/o^2 * (X-p/2) %*% t(X-p/2)
    C <- 1/n * (Coffdiag - diag(diag(Coffdiag)) + diag(Cdiag))
    # Get eigenvectors correspond to the biggest #=rank eigenvalues
    V <- eigen(C)$vector[,1:rank]
    # Project X to the column space of V
    P <- V %*% t(V) %*% X + center
    return(P)
}

#' Noise Generator
#'
#' @param X an input m*n matrix
#' @param p the corruption probability of every pixel
#' @param type the type of noise, 0=white and 1=color
#'
#' @export
GenerateNoise <- function(X, p, type=0){
    m <- dim(X)[1]
    n <- dim(X)[2]
    NoiseLocation <- which(rbinom(m*n, 1, p)==1)
    NoiseNumber <- length(NoiseLocation)
    if (is.na(dim(X))){
        X[NoiseLocation] <- runif(n=NoiseNumber, max=type)
    } else {
        for (d in 1:dim(X)[3]){
            X[NoiseLocation+(d-1)*m*n] <- runif(n=NoiseNumber, max=type)
        }
    }
    return(X)
}
