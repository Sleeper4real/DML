library(mvtnorm)
library("hdm")
library("sandwich")
library(randomForest)
# DML1.segmented function modified from:
# https://github.com/VC2015/DMLonGitHub/blob/master/DMLCodeandTutorial.zip
N = 250 # Number of observations #
m = 300 # Number of simulations #
folds = 2 # Number of folds in DML #

# Model #
# Y = D * theta_1 + Z %*% eta_y + N(0, 1) if D <= tau
# Y = D * theta_2 + Z %*% eta_y + N(0, 1) if D > tau
# D = Z %*% eta_d + N(0, 1)
# Z = N(0.5, Sig)
K = 100 # Number of nuisance parameters, i.e. variables in Z #
theta_1 = 0
theta_2 = 10
tau = 0
eta_y = c(1/(1:K))
eta_d = c(1/(1:K))
Sig = matrix(rep(0.1, K*K), nrow = K, byrow = TRUE)
diag(Sig) = 1

DML.segmented <- function(Z, D, Y, N, gHat, Dreg, nfold=2, A, B, step) {
	noBreaks <- trunc((B - A)/step) + 1 # Number of possible breakpoints w/ A, B
	breaks <- seq(from = A, to = B, length.out = noBreaks) # Generate all possible breakpoints w/ A, B
	breaks <- breaks[2:(noBreaks - 1)]
	foldid <- rep.int(1:nfold, times = ceiling(N/nfold))[sample.int(N)] #fold IDs
	I <- split(1:N, foldid)
	# Create residualized objects to fill
	What <- Vhat <- rep(NA, N)
	tau.est <- rep(NA, length(I)) # Estimated DML1, DML2 threshold locations for each fold
	coef <- coef.naive <- matrix(NA, nrow = 2, ncol = length(I)) # Estimated DML coefficients for each fold
	thetaHat.orc <- thetaHat <- thetaHat.naive <- list(rep(NA, 2)) # Estimated oracle, DML, naive ML coefficients (mean of all folds)
	# Obtain cross-fitted residuals & Estimate threshold
	## Estimate oracle threshold
	tauHat.orc <- NA # Set oracle threshold location to null
	minRSS.orc <- Inf 
	for(thresh in breaks) {
		model.orc <- lm(W ~ D*(D<thresh) + D*(D>=thresh))
		RSS.orc <- as.numeric(summary(model.orc)[6])
		if(RSS.orc < minRSS.orc) {
			minRSS.orc <- RSS.orc
			tauHat.orc <- thresh # Update oracle threshold location
		}
	}
	## Estimate DML threshold
	cat("fold: ")
	for(b in 1:length(I)){
		# Calculate V hat
		dfit <- Dreg(Z[-I[[b]],], D[-I[[b]]]) # m model
		mhat <- predict(dfit, Z[I[[b]],], type="response") # m hat
		Vhat[I[[b]]] <- (D[I[[b]]] - mhat) # V hat
		# Calculate W hat for DML
		for(i in I[[b]]) {
			gHat[i] = g[i] + 1/N^(0.5)
		}
		What[I[[b]]] <- (Y[I[[b]]] - gHat[I[[b]]]) # W hat for DML1
		# Estimate threshold
		minRSS <- Inf # Set minRSS of DML1, DML2 estimator to infinity
		threshIndx = floor(N/10)
		for(thresh in breaks){
			# Calculate RSS for DML1 corresponding to thresh
			model <- lm(What[I[[b]]] ~ D[I[[b]]]*(D[I[[b]]]<thresh) + D[I[[b]]]*(D[I[[b]]]>=thresh))
			RSS <- as.numeric(summary(model)[6])
			if(RSS < minRSS){
				minRSS <- RSS
				tau.est[b] <- thresh
			}
		}
	}
	tauHat <- mean(tau.est) # Estimated threshold location for DML1 (mean of all folds)
	# Estimate coefficients using the estimated thresholds
	## Estimate oracle coefficients
	tauHatIndx.orc = 1
	while(D[tauHatIndx.orc] < tauHat.orc && tauHatIndx.orc < (N - ceiling(N/10))){
		tauHatIndx.orc = tauHatIndx.orc + 1
	}
	thetaHat.orc[1] <- mean(V[1:(tauHatIndx.orc-1)]*W[1:(tauHatIndx.orc-1)])/mean(V[1:(tauHatIndx.orc-1)]*D[1:(tauHatIndx.orc-1)])
	thetaHat.orc[2] <- mean(V[tauHatIndx.orc:N]*W[tauHatIndx.orc:N])/mean(V[tauHatIndx.orc:N]*D[tauHatIndx.orc:N])
	## Estimate DML1, DML2 coefficients
	for(b in 1:length(I)){
		# Estimate DML1 coefficients
		tauHatIndx = 1
		while(D[I[[b]]][tauHatIndx] < tauHat && tauHatIndx < (N - ceiling(N/10))){
			tauHatIndx = tauHatIndx + 1
		}
		I1 = I[[b]][1:(tauHatIndx-1)]
		I2 = I[[b]][tauHatIndx:length(I[[b]])]
		coef[1,b] <- mean(Vhat[I1]*What[I1])/mean(Vhat[I1]*D[I1])
		coef[2,b] <- mean(Vhat[I2]*What[I2])/mean(Vhat[I2]*D[I2])
		coef.naive[1,b] <- mean(D[I1]*What[I1])/mean(D[I1]*D[I1])
		coef.naive[2,b] <- mean(D[I2]*What[I2])/mean(D[I2]*D[I2])
	}
	thetaHat[1] <- mean(coef[1,])
	thetaHat[2] <- mean(coef[2,])
	thetaHat.naive[1] <- mean(coef.naive[1,])
	thetaHat.naive[2] <- mean(coef.naive[2,])
	return( list(thetaHat.orc = thetaHat.orc, thetaHat = thetaHat, thetaHat.naive = thetaHat.naive, tauHat.orc = tauHat.orc, tauHat = tauHat) )
}


Dreg <- function(Z,D){ rlasso(Z, D) } 
Yreg <- function(Z,Y){ rlasso(Z, Y) }

# Estimator Distribution Simulation #
thetaHat1 = matrix(NA, nrow = m, ncol = 3)
thetaHat2 = matrix(NA, nrow = m, ncol = 3)
tauHat = matrix(NA, nrow = m, ncol = 2)
colnames(thetaHat1) = c("Oracle", "DML", "Naive")
colnames(thetaHat2) = c("Oracle", "DML", "Naive")
colnames(tauHat) = c("Oracle", "DML/Naive")
for(c in 1:m)
{
	# Generate and order (Y_i, D_i, Z_i)_{i = 1}^N such that D_i is ascending
	Z <- matrix(, nrow = N, ncol = K)
	D <- Y <- rep(NA, N)
	V <- W <- rep(NA, N)
	g <- gHat <- rep(NA, N)
	# D = Z %*% eta_d + rnorm(N)
	# Y = Z %*% eta_y + rnorm(N)
	for(i in 1:N){
		Z[i,] = 0.5 + rmvnorm(1, rep(0, K), Sig)
		V[i] = rnorm(1)
		D[i] = Z[i,] %*% eta_d + V[i]
		W[i] = rnorm(1)
		if(D[i] <= tau) {
			W[i] = W[i] + theta_1*D[i]
#			cat("<= tau: W = ", W[i], " D = ", D[i], "\n")
		} else {
			W[i] = W[i] + theta_2*D[i]
#			cat("> tau: W = ", W[i], " D = ", D[i], "\n")
		}
		g[i] = Z[i,] %*% eta_y
		Y[i] = g[i] + W[i]
	}
	orderD = order(D)
	Z <- Z[orderD,]
	Y <- Y[orderD]
	D <- D[orderD]
	V <- V[orderD]
	W <- W[orderD]
	g <- g[orderD]
	# Set boundaries of threshold range
	A = floor(D[ceiling(N/10)])
	B = ceiling(D[N - floor(N/10)])

#	plot(D, W, xlim = c(-7, 7), ylim = c(-20, 65))
#	plot(D, W, xlim = c(-7, 7), ylim = c(-20, 65), col = "darkblue")
#	points(D, Y, col= adjustcolor("forestgreen", alpha.f = 0.6) )

	step = 0.25
	DML_test = DML.segmented(Z, D, Y, N, gHat, Dreg, nfold=2, A, B, step)
	thetaHat1[c, 1] = DML_test$thetaHat.orc[[1]][1]
	thetaHat1[c, 2] = DML_test$thetaHat[[1]][1]
	thetaHat1[c, 3] = DML_test$thetaHat.naive[[1]][1]
	thetaHat2[c, 1] = DML_test$thetaHat.orc[[2]][1]
	thetaHat2[c, 2] = DML_test$thetaHat[[2]][1]
	thetaHat2[c, 3] = DML_test$thetaHat.naive[[2]][1]
	tauHat[c, 1] = DML_test$tauHat.orc
	tauHat[c, 2] = DML_test$tauHat
}


colMeans(thetaHat1)
colMeans(thetaHat2)
colMeans(tauHat)

plot(density(thetaHat1[, 1]), xlim = c(-1, 1), ylim = c(0, 4))
lines(density(thetaHat1[, 2]), col = 2)
lines(density(thetaHat1[, 3]), col = 4)

plot(density(thetaHat2[, 1]), xlim = c(9.6, 10.4), ylim = c(0, 25))
lines(density(thetaHat2[, 2]), col = 2)
lines(density(thetaHat2[, 3]), col = 4)

plot(density(tauHat[, 1]), xlim = c(-0.5, 0.5), ylim = c(0, 12))
lines(density(tauHat[, 2]), col = 2)
