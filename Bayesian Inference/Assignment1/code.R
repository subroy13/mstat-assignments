# Question 2
HPD.Beta <- function(a, b, alpha = 0.95) {
    # function to find HPD credible set with shape parameters a and b
    
    credint <- function(right_endpoint) {
        left_endpoint <- qbeta(pbeta(right_endpoint, a, b) - alpha, a, b) 
        diff <- dbeta(left_endpoint, a, b) - dbeta(right_endpoint, a, b)
        return(diff^2)
    }
    
    root <- optimise(credint, interval = c(qbeta(alpha + 1e-5, a, b), 1))$minimum
    int <- c(qbeta(pbeta(root, a, b) - alpha, a, b), root)  # the credible interval
    return(int)
}


simQ2 <- function(theta) {
    set.seed(1911)  # set a seed for reproducibility
    B <- 1000   # perform 1000 simulations
    n <- 20   # sample size is 20
    pb <- txtProgressBar(min = 0, max = B, style = 3)  # set a progress bar
    
    bayesCoverage <- logical(B)
    freqCoverage <- logical(B)   # holds the indicator for coverage of the true theta by 95% CI
    
    for (b in 1:B) {
        data <- rbinom(n, size = 1, prob = theta)   # simulate the data
        Xsum <- sum(data)
        
        credible_int <- HPD.Beta(Xsum + 0.5, n - Xsum + 0.5)
        if (credible_int[1] <= theta & credible_int[2] >= theta) {
            bayesCoverage[b] <- TRUE
        }
        
        theta_hat <- Xsum / n
        
        if ((theta_hat + qnorm(0.025) * sqrt(theta_hat * (1 - theta_hat)/n) <= theta )
            & (theta_hat - qnorm(0.025) * sqrt(theta_hat * (1 - theta_hat)/n) >= theta )) {
            freqCoverage[b] <- TRUE
        }
        
        setTxtProgressBar(pb, value = b)
    }
    close(pb)
    
    cat("Frequentist CI Coverage", sum(freqCoverage)/B, "\n")
    cat("Bayesian CR Coverage", sum(bayesCoverage)/B, "\n")
}

simQ2(theta = 7/8)  # pass theta = 1/8, 1/4, 1/2, 3/4, 7/8


#############################
# Question 6

HPD.CS <- function(mcmc_samples, alpha = 0.95) {
    order_theta <- sort(mcmc_samples)    # perform step 2
    B <- length(order_theta)   
    b <- floor(alpha * B)
    
    lengths <- numeric(B - b)   # array to store the lengths of the intervals
    
    for (j in 1:(B-b)) {
        lengths[j] <- order_theta[j+b] - order_theta[j]  # compute the lengths (step 3)
    }
    
    jstar <- which.min(lengths)  # find j* (step 4)
    return(c(order_theta[jstar], order_theta[jstar + b]))  # return the HPD interval
}

MH.samples <- function(X_samples, burnin = 1000, B = 1000, thin = 20) {
    mcmc_samples <- numeric(B)  # array to hold mcmc_samples
    nchain <- burnin + (B * thin) - 1
    
    # initialize theta_0
    theta <- mean(X_samples)
    pb <- txtProgressBar(min = 0, max = nchain, style = 3)  # set a progress bar
    
    for (i in 1:nchain) {
        accept <- FALSE
        while (!accept) {
            theta_prime <- rnorm(1, mean = theta, sd = 1)  # get theta'
            post_ratio <- sum(log(1 + (X_samples - theta)^2 ) 
                              - log(1 + (X_samples - theta_prime)^2 ))
            A <- min(1, exp(post_ratio))
            check <- rbinom(1, size = 1, prob = A) 
            
            if (check == 1) {
                # if accepted, update theta_t
                theta <- theta_prime
                accept <- TRUE
            }
        }
        
        # if we have thin-th sample after burnin, take that as iid posterior sample
        if ((i >= burnin) & ((i-burnin) %% thin == 0)) {
            index <- (i - burnin) %/% thin + 1
            mcmc_samples[index] <- theta
        }
        
        setTxtProgressBar(pb, value = i)  # update progress bar
    }
    
    close(pb)
    return(mcmc_samples)
}


# now, for n = 5
set.seed(1911)
n <- 5
theta <- 10
X_samples <- rcauchy(n, location = theta)
mcmc_samples <- MH.samples(X_samples)
HPD.CS(mcmc_samples)
HPD.CS(mcmc_samples, alpha = 0.99)


# now, for n = 25
n <- 25
X_samples <- rcauchy(n, location = theta)
mcmc_samples <- MH.samples(X_samples)
HPD.CS(mcmc_samples)
HPD.CS(mcmc_samples, alpha = 0.99)

# now, for n = 100
n <- 100
X_samples <- rcauchy(n, location = theta)
mcmc_samples <- MH.samples(X_samples)
HPD.CS(mcmc_samples)
HPD.CS(mcmc_samples, alpha = 0.99)





