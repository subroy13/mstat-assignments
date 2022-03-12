setwd('D:/Academic/Assignments/Statistical Inference 2_TS/Assignment 2/')

# Data Generation
theta <- 2 * c(1:10)   # different means
sigma.sq <- 1
n <- 25
k <- 10

set.seed(1911)
errors <- rnorm(n * k, mean = 0, sd = sqrt(sigma.sq))
y <- matrix(rep(theta, each = n) + errors, nrow = k, ncol = n, byrow = TRUE)


#################

# Perform Gibbs Sampling
do.Gibbs <- function(y, a1, a2, b1, b2, mu0, sigma0.sq, n_samples = 100, burnin = 1000, nthin = 10) {
    k <- nrow(y)
    n <- ncol(y)
    
    # initialize Gibbs sampler
    theta <- rep(1, k)
    mu_pi <- 0
    sigma_pi.sq <- 1
    sigma.sq <- 1
    
    nchain = burnin + (n_samples * nthin)
    theta_post <- matrix(NA, nrow = n_samples, ncol = k)
    pb <- txtProgressBar(max = nchain, style = 3)
    
    for (t in 1:nchain) {
        # generate theta
        mu1 <- (sigma_pi.sq / (sigma_pi.sq + sigma.sq / n)) * rowMeans(y) + 
            ((sigma.sq / n) / (sigma_pi.sq + sigma.sq / n)) * mu_pi
        Sigma1 <- diag(k) * (sigma_pi.sq * sigma.sq / n) / (sigma_pi.sq + sigma.sq / n)
        theta <- MASS::mvrnorm(mu = mu1, Sigma = Sigma1)
        
        # generate sigma^2
        ss <- sum((y - matrix(theta, nrow = k, ncol = n, byrow = FALSE))^2)
        sigma.sq <- 1 / rgamma(n = 1, shape = (a1 + n/2), rate = (b1 + ss/2))
        
        # generate mu_pi
        mu2 <- (sigma0.sq/(sigma0.sq + sigma_pi.sq/k)) * mean(theta) + 
            ((sigma0.sq/k)/(sigma0.sq + sigma_pi.sq/k)) * mu0
        sigma2 <- (sigma0.sq * sigma_pi.sq / k) / (sigma0.sq + sigma_pi.sq / k)
        mu_pi <- rnorm(n = 1, mean = mu2, sd = sqrt(sigma2))
        
        # generate sigma_pi.sq
        ss2 <- sum((theta - mu_pi)^2)
        sigma_pi.sq <- 1 / rgamma(n = 1, shape = (a2 + k/2), rate = (b2 + ss2/2))
        
        if (t > burnin) {
            if ((t - burnin) %% nthin == 0) {
                index <- (t - burnin) / nthin
                theta_post[index, ] <- theta
            }
        }
        
        setTxtProgressBar(pb, value = t)
    }
    
    close(pb)
    
    return(theta_post)
}


set.seed(1234)
params <- c(1000, 1000, 1000, 1000, 0, 1000)
res <- do.Gibbs(y, params[1], params[2], params[3], params[4], params[5], params[6], 
                n_samples = 1000, burnin = 1000, nthin = 10)
theta_est <- colMeans(res)

paste(paste("(", paste(params, collapse = ", "), ")", sep = ""),
      paste("(", paste(round(theta_est, 3), collapse = ", "), ")", sep = ""), sep = " & ")



# Convergence Plots
plot(1000 + 10 * c(1:1000), res[, 1], type = "l", col = "red", xlab = "iteration", 
     ylab = expression(theta[1]))

plot(1000 + 10 * c(1:1000), res[, 2], type = "l", col = "red", xlab = "iteration", 
     ylab = expression(theta[2]))

plot(1000 + 10 * c(1:1000), res[, 5], type = "l", col = "red", xlab = "iteration", 
     ylab = expression(theta[5]))














