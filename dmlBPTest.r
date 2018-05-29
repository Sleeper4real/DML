# Code modified from "Cross-Fitting Double Machine Learning estimator" by Gabriel Vasconcelos #
# https://www.r-bloggers.com/cross-fitting-double-machine-learning-estimator/ #

library(clusterGeneration)
library(mvtnorm)
library(randomForest)
set.seed(123)
N = 500 # Number of observations #
K = 10 # Number of nuisance parameters, i.e. variables in z #
# True value of effect of d on y #
theta1 = 0.2
theta2 = 0.8
BP = 0
b = 1/(1:K)
sigma = genPositiveDefMat(K,"unifcorrmat")$Sigma # genPositiveDefMat gives eigenvalues and the matrix, we only need the matrix here #
sigma = cov2cor(sigma)
s = 100 # Number of simulations #
thetaHat = matrix(NA, s, 3)
colnames(thetaHat) = c("OLS", "Naive DML", "Cross-fitting DML")
for(i in 1:s)
{
	z = rmvnorm(N, sigma = sigma) # z ~ N(0, sigma) #
	g = as.vector(cos(z%*%b)^2)
	m = as.vector(sin(z%*%b) + cos(z%*%b))
	d = m + rnorm(N)
	y = c()
	for(j in 1:N)
	{
		if(d[j] <= 0)
		{
			y[j] = theta1*d[j] + g[j] + rnorm(1)
		}
		else
		{
			y[j] = theta2*d[j] + g[j] + rnorm(1)
		}
	}
	# OLS #
	theta_ols = coef(lm(y~d))[2]
	thetaHat[i, 1] = theta_ols
	# Naive DML #
	model = randomForest(z, y, maxnodes = 20)
	gHat = predict(model, z)
	model_d = randomForest(z, d, maxnodes = 20)
	mHat = predict(model_d, z)
	yTilda = y - gHat
	vHat = d - mHat
	theta_nv = mean(vHat*yTilda)/mean(vHat*d)
	thetaHat[i, 2] = theta_nv
	# Cross-fitting DML #
	I = sort(sample(1:N, N/2))
	IC = setdiff(1:N, I)
	model1 = randomForest(z[IC,], y[IC], maxnodes = 10)
	model2 = randomForest(z[I,], y[I], maxnodes = 10)
	gHat1 = predict(model1, z[I,])
	gHat2 = predict(model2, z[IC,])
	model_d1 = randomForest(z[IC,], d[IC], maxnodes = 10)
	model_d2 = randomForest(z[I,], d[I], maxnodes = 10)
	mHat1 = predict(model_d1, z[I,])
	mHat2 = predict(model_d2, z[IC,])
	yTilda1 = y[I] - gHat1
	yTilda2 = y[IC] - gHat2
	vHat1 = d[I] - mHat1
	vHat2 = d[IC] - mHat2
	theta_cf1 = mean(vHat1*yTilda1)/mean(vHat1*d[I])
	theta_cf2 = mean(vHat2*yTilda2)/mean(vHat2*d[IC])
	theta_cf = mean(c(theta_cf1, theta_cf2))
	thetaHat[i,3] = theta_cf
}
colMeans(thetaHat)
plot(density(thetaHat[, 1]), xlim = c(0.2, 0.8), ylim = c(0, 10))
lines(density(thetaHat[, 2]), col = 2)
lines(density(thetaHat[, 3]), col = 4)
abline(v = 0.2, lty = 2, col = 3) # Dashed line is the true value #
abline(v = 0.8, lty = 2, col = 3) # Dashed line is the true value #
legend("topleft", legend = c("OLS","Naive DML","Cross-fitting DML"), col = c(1,2,4), lty = 1, cex = 0.7, seg.len = 0.7, bty = "n")