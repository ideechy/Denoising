library(EBImage)
im <-  readImage("NCState.png")

X=matrix(1,10,10)
X=GenerateNoise(X, 0.2, 1)

par(mfrow=c(3,1))
contour(X-1, main="origin")
contour(RPCA(X)$E, main="RPCA")
contour(X-POPCAUnif(X,1,0.2), main="POPCA")

