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
# Z = N(0, Sig)
K = 100 # Number of nuisance parameters, i.e. variables in Z #
theta_1 = 0
theta_2 = 10
tau = 0
eta_y = c(1/(1:K))
eta_d = c(1/(1:K))
Sig = matrix(rep(0.1, K*K), nrow = K, byrow = TRUE)
diag(Sig) = 1

DML.segmented <- function(Z, D, Y, N, Dreg, Yreg1, Yreg2, nfold=2, A, B, step) {
	noBreaks <- trunc((B - A)/step) + 1 # Number of possible breakpoints w/ A, B
	breaks <- seq(from = A, to = B, length.out = noBreaks) # Generate all possible breakpoints w/ A, B
	breaks <- breaks[2:(noBreaks - 1)]
	foldid <- rep.int(1:nfold, times = ceiling(N/nfold))[sample.int(N)] #fold IDs
	I <- split(1:N, foldid)
	# Create residualized objects to fill
	What1 <- What2 <- Vhat <- rep(NA, N)
	tau.est1 <- tau.est2 <- rep(NA, length(I)) # Estimated DML1, DML2 threshold locations for each fold
	coef.est1 <- coef.est2 <- matrix(NA, nrow = 2, ncol = length(I)) # Estimated DML1, DML2, oracle coefficients for each fold
	thetaHat.orc <- thetaHat.est1 <- thetaHat.est2 <- list(rep(NA, 2)) # Estimated oracle, DML1, DML2 coefficients (mean of all folds)
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
	## Estimate DML1, DML2 thresholds
	cat("fold: ")
	for(b in 1:length(I)){
		# Calculate V hat
		dfit <- Dreg(Z[-I[[b]],], D[-I[[b]]]) # m model
		mhat <- predict(dfit, Z[I[[b]],], type="response") # m hat
		Vhat[I[[b]]] <- (D[I[[b]]] - mhat) # V hat
		# Calculate W hat for DML1
		yfit1 <- Yreg1(cbind(D[-I[[b]]], Z[-I[[b]],]), Y[-I[[b]]]) # g model for DML1
		ghat1 <- predict(yfit1, cbind(D[I[[b]]], Z[I[[b]],]), type="response") # g hat for DML1
		What1[I[[b]]] <- (Y[I[[b]]] - ghat1) # W hat for DML1
		# Estimate threshold
		minRSS.est1 <- minRSS.est2 <- Inf # Set minRSS of DML1, DML2 estimator to infinity
		threshIndx = floor(N/10)
		for(thresh in breaks){
			# Calculate RSS for DML1 corresponding to thresh
			model.est1 <- lm(What1[I[[b]]] ~ D[I[[b]]]*(D[I[[b]]]<thresh) + D[I[[b]]]*(D[I[[b]]]>=thresh))
			RSS.est1 <- as.numeric(summary(model.est1)[6])
			if(RSS.est1 < minRSS.est1){
				minRSS.est1 <- RSS.est1
				tau.est1[b] <- thresh
			}
			# Calculate W hat & RSS for DML2 corresponding to thresh
			while(D[I[[b]]][threshIndx] < thresh && threshIndx < (N - ceiling(N/10))){
				threshIndx = threshIndx + 1
			}
			I1.est2 = I[[b]][1:(threshIndx-1)]
			I2.est2 = I[[b]][threshIndx:length(I[[b]])]
			yfit2_1 <- Yreg2(cbind(D[-I1.est2], Z[-I1.est2,]), Y[-I1.est2]) # E[Y|Z, D<=tresh] model for DML1
			ghat2_1 <- predict(yfit2_1, cbind(D[I1.est2], Z[I1.est2,]), type="response") # g hat (<=thresh) for DML1
			yfit2_2 <- Yreg2(cbind(D[-I2.est2], Z[-I2.est2,]), Y[-I2.est2]) # E[Y|Z, D>tresh] model for DML1
			ghat2_2 <- predict(yfit2_2, cbind(D[I2.est2], Z[I2.est2,]), type="response") # g hat (>thresh) for DML1
			What2[I1.est2] <- (Y[I1.est2] - ghat2_1) # W hat for DML1
			What2[I2.est2] <- (Y[I2.est2] - ghat2_2)
			model.est2 <- lm(What2[I[[b]]] ~ D[I[[b]]]*(D[I[[b]]]<thresh) + D[I[[b]]]*(D[I[[b]]]>=thresh))
			RSS.est2 <- as.numeric(summary(model.est2)[6])
			if(RSS.est2 < minRSS.est2){
				minRSS.est2 <- RSS.est2
				tau.est2[b] <- thresh
			}
		}
	}
	tauHat.est1 <- mean(tau.est1) # Estimated threshold location for DML1 (mean of all folds)
	tauHat.est2 <- mean(tau.est2) # Estimated threshold location for DML2 (mean of all folds)
	# Estimate coefficients using the estimated thresholds
	## Estimate oracle coefficients
	tauHatIndx.orc = 1
	while(D[tauHatIndx.orc] < tauHat.orc && tauHatIndx.orc < (N - ceiling(N/10))){
		tauHatIndx.orc = tauHatIndx.orc + 1
	}
	thetaHat.orc[1] <- mean(V[1:(tauHatIndx.orc-1)]*W[1:(tauHatIndx.orc-1)])/mean(V[1:(tauHatIndx.orc-1)]*V[1:(tauHatIndx.orc-1)])
	thetaHat.orc[2] <- mean(V[tauHatIndx.orc:N]*W[tauHatIndx.orc:N])/mean(V[tauHatIndx.orc:N]*D[tauHatIndx.orc:N])
	## Estimate DML1, DML2 coefficients
	for(b in 1:length(I)){
		# Estimate DML1 coefficients
		tauHatIndx.est1 = 1
		while(D[I[[b]]][tauHatIndx.est1] < tauHat.est1 && tauHatIndx.est1 < (N - ceiling(N/10))){
			tauHatIndx.est1 = tauHatIndx.est1 + 1
		}
		yfit1 <- Yreg1(Z[-I[[b]],], Y[-I[[b]]]) # g model for DML1
		ghat1 <- predict(yfit1, Z[I[[b]],], type="response") # g hat for DML1
		What1[I[[b]]] <- (Y[I[[b]]] - ghat1) # W hat for DML1
		I1.est1 = I[[b]][1:(tauHatIndx.est1-1)]
		I2.est1 = I[[b]][tauHatIndx.est1:length(I[[b]])]
		coef.est1[1,b] <- mean(Vhat[I1.est1]*What1[I1.est1])/mean(Vhat[I1.est1]*D[I1.est1])
		coef.est1[2,b] <- mean(Vhat[I2.est1]*What1[I2.est1])/mean(Vhat[I2.est1]*D[I2.est1])
		# Estimate DML2 coefficients
		tauHatIndx.est2 = 1
		while(D[I[[b]]][tauHatIndx.est2] < tauHat.est2 && tauHatIndx.est2 < (N - ceiling(N/10))){
			tauHatIndx.est2 = tauHatIndx.est2 + 1
		}
		I1.est2 = I[[b]][1:(tauHatIndx.est2-1)]
		I2.est2 = I[[b]][tauHatIndx.est2:length(I[[b]])]
		yfit2_1 <- Yreg2(Z[-I1.est2,], Y[-I1.est2]) # E[Y|Z, D<=tresh] model for DML1
		ghat2_1 <- predict(yfit2_1, Z[I1.est2,], type="response") # g hat (<=thresh) for DML1
		yfit2_2 <- Yreg2(Z[-I2.est2,], Y[-I2.est2]) # E[Y|Z, D>tresh] model for DML1
		ghat2_2 <- predict(yfit2_2, Z[I2.est2,], type="response") # g hat (>thresh) for DML1
		What2[I1.est2] <- (Y[I1.est2] - ghat2_1) # W hat for DML1
		What2[I2.est2] <- (Y[I2.est2] - ghat2_2)
		coef.est2[1,b] <- mean(Vhat[I1.est2]*What2[I1.est2])/mean(Vhat[I1.est2]*D[I1.est2])
		coef.est2[2,b] <- mean(Vhat[I2.est2]*What2[I2.est2])/mean(Vhat[I2.est2]*D[I2.est2])
	}
	thetaHat.est1[1] <- mean(coef.est1[1,])
	thetaHat.est1[2] <- mean(coef.est1[2,])
	thetaHat.est2[1] <- mean(coef.est2[1,])
	thetaHat.est2[2] <- mean(coef.est2[2,])
	return( list(thetaHat.orc = thetaHat.orc, thetaHat.est1 = thetaHat.est1, thetaHat.est2 = thetaHat.est2, tauHat.orc = tauHat.orc, tauHat.est1 = tauHat.est1, tauHat.est2 = tauHat.est2) )
}


Dreg <- function(Z,D){ rlasso(Z, D) } 
Yreg1 <- function(Z,Y){ randomForest(Z, Y, maxnodes = 10) }
Yreg2 <- function(Z,Y){ rlasso(Z, Y) }

# Estimator Distribution Simulation #
thetaHat1 = matrix(NA, nrow = m, ncol = 3)
thetaHat2 = matrix(NA, nrow = m, ncol = 3)
tauHat = matrix(NA, nrow = m, ncol = 3)
colnames(thetaHat1) = c("Oracle", "DML1", "DML2")
colnames(thetaHat2) = c("Oracle", "DML1", "DML2")
colnames(tauHat) = c("Oracle", "DML1", "DML2")
for(c in 1:m)
{
	# Generate and order (Y_i, D_i, Z_i)_{i = 1}^N such that D_i is ascending
	Z <- matrix(, nrow = N, ncol = K)
	D <- Y <- rep(NA, N)
	V <- W <- rep(NA, N)
	# D = Z %*% eta_d + rnorm(N)
	# Y = Z %*% eta_y + rnorm(N)
	for(i in 1:N){
		Z[i,] = rmvnorm(1, rep(0, K), Sig)
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
		Y[i] = Z[i,] %*% eta_y + W[i]
	}
	orderD = order(D)
	Z <- Z[orderD,]
	Y <- Y[orderD]
	D <- D[orderD]
	V <- V[orderD]
	W <- W[orderD]
	# Set boundaries of threshold range
	A = floor(D[ceiling(N/10)])
	B = ceiling(D[N - floor(N/10)])

#	plot(D, W, xlim = c(-7, 7), ylim = c(-20, 65))
#	plot(D, W, xlim = c(-7, 7), ylim = c(-20, 65), col = "darkblue")
#	points(D, Y, col= adjustcolor("forestgreen", alpha.f = 0.6) )

	step = 0.25
	DML_test = DML.segmented(Z, D, Y, N, Dreg, Yreg1, Yreg2, nfold=2, A, B, step)
	thetaHat1[c, 1] = DML_test$thetaHat.orc[[1]][1]
	thetaHat1[c, 2] = DML_test$thetaHat.est1[[1]][1]
	thetaHat1[c, 3] = DML_test$thetaHat.est2[[1]][1]
	thetaHat2[c, 1] = DML_test$thetaHat.orc[[2]][1]
	thetaHat2[c, 2] = DML_test$thetaHat.est1[[2]][1]
	thetaHat2[c, 3] = DML_test$thetaHat.est2[[2]][1]
	tauHat[c, 1] = DML_test$tauHat.orc
	tauHat[c, 2] = DML_test$tauHat.est1
	tauHat[c, 3] = DML_test$tauHat.est2
}


colMeans(thetaHat1)
colMeans(thetaHat2)
colMeans(tauHat)

plot(density(thetaHat1[, 1]), xlim = c(-4, 4), ylim = c(0, 5))
lines(density(thetaHat1[, 2]), col = 2)
lines(density(thetaHat1[, 3]), col = 4)

plot(density(thetaHat2[, 1]), xlim = c(4, 18), ylim = c(0, 0.5))
lines(density(thetaHat2[, 2]), col = 2)
lines(density(thetaHat2[, 3]), col = 4)







for(b in 1:length(I)){
		# Estimate DML1 coefficients
		tauHatIndx.est1 = 1
		while(D[I[[b]]][tauHatIndx.est1] < tauHat.est1){
			tauHatIndx.est1 = tauHatIndx.est1 + 1
		}
		I1.est1 = I[[b]][1:(tauHatIndx.est1-1)]
		I2.est1 = I[[b]][tauHatIndx.est1:length(I[[b]])]
		coef.est1[1,b] <- mean(Vhat[I1.est1]*What1[I1.est1])/mean(Vhat[I1.est1]*D[I1.est1])
		coef.est1[2,b] <- mean(Vhat[I2.est1]*What1[I2.est1])/mean(Vhat[I2.est1]*D[I2.est1])
		# Estimate DML2 coefficients
		tauHatIndx.est2 = 1
		while(D[I[[b]]][tauHatIndx.est2] < tauHat.est2){
			tauHatIndx.est2 = tauHatIndx.est2 + 1
		}
		I1.est2 = I[[b]][1:(tauHatIndx.est2-1)]
		I2.est2 = I[[b]][tauHatIndx.est2:length(I[[b]])]
		yfit2_1 <- Yreg2(cbind(D[-I1.est2], Z[-I1.est2,]), Y[-I1.est2]) # E[Y|Z, D<=tresh] model for DML1
		ghat2_1 <- predict(yfit2_1, cbind(D[I1.est2], Z[I1.est2,]), type="response") # g hat (<=thresh) for DML1
		yfit2_2 <- Yreg2(cbind(D[-I2.est2], Z[-I2.est2,]), Y[-I2.est2]) # E[Y|Z, D>tresh] model for DML1
		ghat2_2 <- predict(yfit2_2, cbind(D[I2.est2], Z[I2.est2,]), type="response") # g hat (>thresh) for DML1
		What2[I1.est2] <- (Y[I1.est2] - ghat2_1) # W hat for DML1
		What2[I2.est2] <- (Y[I2.est2] - ghat2_2)
		coef.est2[1,b] <- mean(Vhat[I1.est2]*What2[I1.est2])/mean(Vhat[I1.est2]*D[I1.est2])
		coef.est2[2,b] <- mean(Vhat[I2.est2]*What2[I2.est2])/mean(Vhat[I2.est2]*D[I2.est2])
	}
	thetaHat.est1[1] <- mean(coef.est1[1,])
	thetaHat.est1[2] <- mean(coef.est1[2,])
	thetaHat.est2[1] <- mean(coef.est2[1,])
	thetaHat.est2[2] <- mean(coef.est2[2,])
