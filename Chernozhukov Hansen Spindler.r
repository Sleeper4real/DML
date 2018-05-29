library(clusterGeneration)
library(mvtnorm)
library(Matrix)
library(hdm)
set.seed(123)
n = 200 # Number of observations #
p_x = 200 # Number of controls # 
p_z = 150 # Number of instruments #
v = 4/9
for(i in 5:p_z)
{
	v <- (v + 1)
}
alpha = 0
beta = vector(, p_x)
beta[1:4] = 1/(9*v)
beta[5:p_x] = 1/((5:p_x)^2*v)
delta = 3/((1:p_z)^2)
Pi = cbind(diag(p_z), matrix(0, p_z, p_x - p_z))
sigma_eu = matrix(c(1, 0.6, 0.6, 1), nrow = 2, ncol = 2)
sigma_zeta = diag(p_z)
sigma_x = diag(p_x)
for(i in 1:(p_x - 1))
{
	for(j in (i + 1):p_x)
	{
		sigma_x[i, j] <- sigma_x[i, j - 1]*0.5
		sigma_x[j, i] <- sigma_x[i, j]
	}
}
sigma_zeta_x = rbind(cbind(sigma_zeta, matrix(0, p_z, p_x)), cbind(matrix(0, p_x, p_z), sigma_x))
Sigma = rbind(cbind(sigma_eu, matrix(0, 2, p_z + p_x)), cbind(matrix(0, p_z + p_x, 2), sigma_zeta_x))
rm(sigma_eu, sigma_x, sigma_zeta, sigma_zeta_x)
m = 200 # Number of simulations #
alpha_hat = matrix(, m, 2)
colnames(alpha_hat) = c("Oracle", "Non-orthogonal")
for(s in 1:m)
{
	e = vector(, n)
	u = vector(, n)
	zeta = matrix(, n, p_z)
	x = matrix(, n, p_x)
	z = matrix(, n, p_z)
	for(i in 1:n)
	{
		temp = rmvnorm(2 + p_z + p_x, sigma = Sigma)
		e[i] = temp[1]
		u[i] = temp[2]
		for(j in 1:p_z)
		{
			zeta[i, j] = temp[j + 2]
		}
		for(j in 1:p_x)
		{
			x[i, j] = temp[j + p_z + 2]
		}
		rm(temp)
		for(j in 1:p_z)
		{
			z[i, j] = x[i, j] + 0.125*zeta[i, j]
		}
	}
	d = c(x %*% beta + z %*% delta + u)
	y = c(alpha*d + x %*% beta + 2*e)
	# Oracle #
	y_po = y - x %*% (alpha*(beta + t(Pi) %*% delta)) # y_i - E[y_i|x_i] #
	d_po = d - x %*% (beta + t(Pi) %*% delta) # d_i - E[d_i|x_i] #
	IV = zeta %*% delta # n x 1 #
	alpha_hat[s, 1] = solve(t(IV) %*% d) %*% (t(IV) %*% y) # = (t(IV) %*% y)/(t(IV) %*% d) #
	# Non-orthogonal #
	lasso_d_xz = rlasso(d~cbind(x, z), post = FALSE, intercept = FALSE)$coefficients
	lasso_y_x = rlasso(y~x, post = FALSE, intercept = FALSE)$coefficients
	xSelect = 0
	for(i in 1:p_x)
	{
		if(lasso_d_xz[i] != 0 || lasso_y_x[i] != 0)
		{
			xSelect = xSelect + 1
		}
	}
	zSelect = 0
	for(i in (p_x + 1):(p_x + p_z))
	{
		if(lasso_d_xz[i] != 0)
		{
			zSelect = zSelect + 1
		}
	}
	alpha_hat[s, 2] = tsls(y = y, d = d, x = x[1:n, 1:xSelect], z = z[1:n, 1:zSelect], intercept = FALSE)$coefficients[1]
}
colMeans(alpha_hat)
plot(density(alpha_hat_large[, 1]), xlim = c(-1.5, 1.5), ylim = c(0, 2))
lines(density(alpha_hat_large[, 2]), col = 2)