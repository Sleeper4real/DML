# Comparison between Cross-fitting and Naive DML in a Linear Model #
library(mvtnorm)
library("hdm")
library("sandwich")
# DML2.for.PLM function taken and DML2.Naive.PLM modified from:
# https://github.com/VC2015/DMLonGitHub/blob/master/DMLCodeandTutorial.zip
N = 250 # Number of observations #
m = 300 # Number of simulations #
folds = 5 # Number of folds in cross-fitting DML #
# Model #
# Y = theta * D + X %*% eta + N(0, 1) where (D, X) ~ N(0, Sig)
K = 100 # Number of nuisance parameters, i.e. variables in X #
theta = 2 # True value of effect of d on y #
eta = c(1/(1:100))
Sig = rbind(rep(0.2, K + 1), cbind(rep(0.2, K), matrix(rep(0.1, K*K), nrow = K, byrow = TRUE)))
diag(Sig) = 1
# Results #
thetaHat = matrix(NA, m, 2)
thetaHatSE = matrix(NA, m, 2)
colnames(thetaHat) = c("Cross-fitting", "No Cross-fitting")
colnames(thetaHatSE) = c("Cross-fitting", "No Cross-fitting")
# DML Functions #
DML2.for.PLM <- function(x, d, y, dreg, yreg, nfold=2, print=FALSE) {
# this implements DML2 algorithm, where there moments are estimated via DML, before constructing
# the pooled estimate of theta randomly split data into folds
nobs <- nrow(x)
foldid <- rep.int(1:nfold,times = ceiling(nobs/nfold))[sample.int(nobs)]
I <- split(1:nobs, foldid)
# create residualized objects to fill
ytil <- dtil <- rep(NA, nobs)
# obtain cross-fitted residuals
if(print){ cat("fold: ") }
for(b in 1:length(I)){
  dfit <- dreg(x[-I[[b]],], d[-I[[b]]])  #take a fold out
  yfit <- yreg(x[-I[[b]],], y[-I[[b]]])  # take a folot out
  dhat <- predict(dfit, x[I[[b]],], type="response")  #predict the fold out
  yhat <- predict(yfit, x[I[[b]],], type="response")  #predict the fold out
  dtil[I[[b]]] <- (d[I[[b]]] - dhat) #record residual
  ytil[I[[b]]] <- (y[I[[b]]] - yhat) #record residial
  if(print){  cat(b," ") }
  }
  rfit <- lm(ytil ~ dtil)               #estimate the main parameter by regressing one residual on the other
  coef.est <- coef(rfit)[2]             #extract coefficient 
  se <- sqrt(vcovHC(rfit)[2,2])         #record standard error
  if(print){ cat(sprintf("\ncoef (se) = %g (%g)\n", coef.est , se)) }
  return( list(coef.est =coef.est , se=se, dtil=dtil, ytil=ytil) )
}
DML2.Naive.PLM <- function(x, d, y, dreg, yreg, print=FALSE) {
  dfit <- dreg(x, d)
  yfit <- yreg(x, y)
  dhat <- predict(dfit, x, type="response")
  yhat <- predict(yfit, x, type="response")
  dtil <- d - dhat #record residual
  ytil <- y - yhat #record residial
  rfit <- lm(ytil ~ dtil)               #estimate the main parameter by regressing one residual on the other
  coef.est <- coef(rfit)[2]             #extract coefficient 
  se <- sqrt(vcovHC(rfit)[2,2])         #record standard error
  if(print){ cat(sprintf("coef (se) = %g (%g)\n", coef.est , se)) }
  return( list(coef.est =coef.est , se=se, dtil=dtil, ytil=ytil) )
}
Dreg <- function(x,d){ rlasso(x, d) }
Yreg <- function(x,y){ rlasso(x, y) }
# Estimator Distribution Simulation #
for(i in 1:m)
{
	temp = rmvnorm(N, sigma = Sig)
	D = temp[,1]
	X = temp[,2:(K + 1)]
	Y = theta*D + X %*% eta + rnorm(N)
	# Cross-fitting #
	DML_CF = DML2.for.PLM(X, D, Y, Dreg, Yreg, nfold = folds)
	thetaHat[i, 1] = DML_CF$coef.est
	thetaHatSE[i, 1] = DML_CF$se
	# No Cross-fitting #
	DML_Naive = DML2.Naive.PLM(X, D, Y, Dreg, Yreg)
	thetaHat[i, 2] = DML_Naive$coef.est
	thetaHatSE[i, 2] = DML_Naive$se
}
print("Mean of Theta Hat")
colMeans(thetaHat)
print("SE of Theta Hat")
colMeans(thetaHatSE)
plot(density(thetaHat[, 1]), xlim = c(1, 3), ylim = c(0, 5))
lines(density(thetaHat[, 2]), col = 2)
# Comparison Between Number of Folds #
start = 2
interval = 4
iterations = 5
print("Naive DML")
DML_Naive = DML2.Naive.PLM(X, D, Y, Dreg, Yreg, TRUE)
print("Cross-fitting DML")
while(iterations >= 0)
{
	sprintf("Folds = %i", start)
	DML_CF = DML2.for.PLM(X, D, Y, Dreg, Yreg, nfold = start, TRUE)
	start = start + interval
	iterations = iterations - 1
}


