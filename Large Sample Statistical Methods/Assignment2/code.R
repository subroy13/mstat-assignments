# Problem 6

tanhinv <- function(r) {
    return(1/2 * log((1+r)/(1-r)))
}

estimate.coverage <- function(n, rho) {
    # set a seed for reproducibility
    set.seed(1911)
    
    # create the parameters
    mu_vec <- c(0, 0)
    sigma_mat <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
    coverage <- logical(1000)  # a logical array to see if the coverage happens
    
    # loop through each simulation
    for (i in 1:1000) {
        data <- MASS::mvrnorm(n = n, mu = mu_vec, Sigma = sigma_mat)
        r <- cor(data[, 1], data[, 2])   # sample correlation coefficient
        z <- tanhinv(r)  # obtain tanh inverse transformation of sample correlation
        z.upper <- tanh(z + qnorm(0.975) / sqrt(n))
        z.lower <- tanh(z - qnorm(0.975) / sqrt(n))
        
        if (rho <= z.upper & rho >= z.lower) {
            coverage[i] <- TRUE
        }
    }
    
    # finally output the proportion of simulation where the coverage occurs
    return(sum(coverage) / 1000)
}

estimate.coverage(20, 0.2)
estimate.coverage(20, 0.4)
estimate.coverage(20, 0.6)
estimate.coverage(20, 0.8)

estimate.coverage(30, 0.2)
estimate.coverage(30, 0.4)
estimate.coverage(30, 0.6)
estimate.coverage(30, 0.8)


# r_n

gen.r <- function(n, rho) {
    # set a seed for reproducibility
    set.seed(1911)
    
    # create the parameters
    mu_vec <- c(0, 0)
    sigma_mat <- matrix(c(1, rho, rho, 1), nrow = 2, byrow = TRUE)
    
    r <- sapply(1:1000, FUN = function(i) {
        data <- MASS::mvrnorm(n = n, mu = mu_vec, Sigma = sigma_mat)
        return(cor(data[, 1], data[, 2]))
    })
    
    return(r)
}

r <- gen.r(30, 0.4)
z <- tanhinv(r)

par(mfrow = c(1, 2))
hist(r, breaks = 20, col = "salmon", freq = FALSE, xlab = expression(r[n]), main = expression(paste("Histogram of ", r[n])))
hist(z, breaks = 20, col = "slateblue", freq = FALSE, xlab = expression(paste(tanh^-1,"(", r[n], ")" )), 
     main = expression(paste("Histogram of ", tanh^-1,"(", r[n], ")" )))



# Problem 13

set.seed(1911)
samp <- rcauchy(25, location = 1, scale = 1)

median(samp)

update.theta <- function(init.theta) {
    
    theta_t <- init.theta  # initialize
    error <- Inf   # initialize error to some large value
    n_iter <- 0    # store number of iterations
    
    while(error > 1e-10) {
        # update theta by iteration
        new_theta <- theta_t + (2/25) * sum( 2 * (samp - theta_t)/(1 + (samp - theta_t)^2 ) )
        
        error <- (new_theta - theta_t) / theta_t    # calculate relative error
        n_iter <- n_iter + 1    # increase counter of iteration
        
        theta_t <- new_theta   # looping
    }
    
    return(list(MLE = theta_t, iteration = n_iter))
}

update.theta(init.theta = median(samp))
update.theta(init.theta = 0)






