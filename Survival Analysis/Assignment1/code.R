setwd("D:/Academic/Assignments/Survival Analysis_AG/Assignment 1/")

library(dplyr)
library(tidyr)
library(ggplot2)

# Comparing different Confidence interval performances
# Simulation for Type 1 and Type 2 censoring
sim.type1 <- function(n, T0, theta = 1, seed = 1234) {
    set.seed(seed)  # set a seed
    t_i <- rexp(n, rate = 1/theta)
    x_i <- pmin(t_i, T0)
    del_i <- ifelse(t_i < T0, 1, 0)
    
    return(data.frame("Observation" = x_i, "Indicator" = del_i))
}


sim.type2 <- function(n, r, theta = 1, seed = 1234) {
    set.seed(seed)  # set a seed
    t_i <- rexp(n, rate = 1/theta)
    t_r <- sort(t_i)[r]
    x_i <- pmin(t_i, t_r)
    del_i <- ifelse(t_i < t_r, 1, 0)
    
    return(data.frame("Observation" = x_i, "Indicator" = del_i))
}


# Function to create the CI's
ci.type2 <- function(data, alpha = 0.05) {
    x_i <- data$Observation
    del_i <- data$Indicator
    r <- sum(del_i)
    theta_mle <- sum(x_i) / r
    ci1 <- c(2*r*theta_mle/ qchisq(1-alpha/2, df = 2*r), 2*r*theta_mle/ qchisq(alpha/2, df = 2*r))
    
    return(ci1)
}

ci.type1 <- function(data, T0, alpha = 0.05) {
    x_i <- data$Observation
    del_i <- data$Indicator
    r <- sum(del_i)
    n <- length(x_i)
    V <- sum(x_i)
    theta_mle <- V / r
    
    # asymptotic CI based on I(theta)
    ci1 <- c(theta_mle + qnorm(alpha/2) * theta_mle / sqrt(n * (1 - exp(-T0/theta_mle) )),
             theta_mle - qnorm(alpha/2) * theta_mle / sqrt(n * (1 - exp(-T0/theta_mle) )) )
    
    # asymptotic CI based on observed Fisher's information
    ci2 <- c(theta_mle + qnorm(alpha/2) * theta_mle / sqrt(r),
             theta_mle - qnorm(alpha/2) * theta_mle / sqrt(r) )
    
    # Cox's CI
    ci3 <- c(2*r*theta_mle/ qchisq(1-alpha/2, df = 2*r + 1), 
             2*r*theta_mle/ qchisq(alpha/2, df = 2*r + 1))
    
    # LRT based CI
    ci4.upper <- tryCatch({
        uniroot(f = function(theta) { 2*r*log(theta) + 2*V/theta - (qchisq(1-alpha, df = 1) + 2*r*log(theta_mle) + 2*r)}, interval = c(theta_mle, V))$root
    }, error = function(cond){
        return(theta_mle)
    })
    
    ci4.lower <- tryCatch({
        uniroot(f = function(theta) { 2*r*log(theta) + 2*V/theta - (qchisq(1-alpha, df = 1) + 2*r*log(theta_mle) + 2*r)}, interval = c(1e-5, theta_mle))$root
    }, error = function(cond){
        return(theta_mle)
    })
    
    ci4 <- c(ci4.lower, ci4.upper)
    
    return(list("CI1" = ci1, "CI2" = ci2, "CI3" = ci3, "CI4" = ci4))
    
}


do.sim_ci <- function(B, n, T0, theta, seed = 1234) {
    # B = number of resamples
    set.seed(seed)
    data_seeds <- sample(1:1e6, size = B)   # seeds to use to generate data
    
    r <- round(n * (1 - exp(-T0/theta)))  # E(r) under type1 is used as parameter for type 2
    pb <- txtProgressBar(min = 0, max = B, style = 3)  # initialize a progress bar
    
    # matrix to contain output of length of CI and coverage by the CI
    lenMat <- matrix(NA, nrow = B, ncol = 5)
    coverageMat <- matrix(FALSE, nrow = B, ncol = 5)  
    
    for (b in 1:B) {
        df.type1 <- sim.type1(n, T0, theta, seed = data_seeds[b])
        df.type2 <- sim.type2(n, r, theta, seed = data_seeds[b])
        
        out <- ci.type1(df.type1, T0)
        if ((out$CI1[1] <= theta) & (out$CI1[2] >= theta)) { coverageMat[b, 1] <- TRUE }
        if ((out$CI2[1] <= theta) & (out$CI2[2] >= theta)) { coverageMat[b, 2] <- TRUE }
        if ((out$CI3[1] <= theta) & (out$CI3[2] >= theta)) { coverageMat[b, 3] <- TRUE }
        if ((out$CI4[1] <= theta) & (out$CI4[2] >= theta)) { coverageMat[b, 4] <- TRUE }
        lenMat[b, 1] <- out$CI1[2] - out$CI1[1]
        lenMat[b, 2] <- out$CI2[2] - out$CI2[1]
        lenMat[b, 3] <- out$CI3[2] - out$CI3[1]
        lenMat[b, 4] <- out$CI4[2] - out$CI4[1]
        
        out <- ci.type2(df.type2)
        if ((out[1] <= theta) & (out[2] >= theta)) { coverageMat[b, 5] <- TRUE }
        lenMat[b, 5] <- out[2] - out[1]
        
        setTxtProgressBar(pb, value = b)  # update progress bar
    }
    
    close(pb) # close the connection to progress bar
    
    result <- data.frame(Method = c("CI1", "CI2", "CI3", "CI4", "CI5"),
                         `Average Length` = round(colMeans(lenMat, na.rm = TRUE),3),
                         `Coverage Probability` = round(colMeans(coverageMat, na.rm = TRUE), 3))
    return(result)
}


# Proportion of Censoring 5%
theta <- 1
T0 <- theta * log(1/0.05)
out <- do.sim_ci(B = 1000, n = 100, T0 = T0, theta = theta, seed = 1911)  # change value of n as required
out
paste(out$Average.Length, collapse = " & ")
paste(out$Coverage.Probability, collapse = " & ")

# Proportion of Censoring 10%
theta <- 1
T0 <- theta * log(1/0.1)
out <- do.sim_ci(B = 1000, n = 100, T0 = T0, theta = theta, seed = 1911)  # change value of n as required
out
paste(out$Average.Length, collapse = " & ")
paste(out$Coverage.Probability, collapse = " & ")

# Proportion of Censoring 20%
theta <- 1
T0 <- theta * log(1/0.2)
out <- do.sim_ci(B = 1000, n = 5, T0 = T0, theta = theta, seed = 1911)  # change value of n as required
out
paste(out$Average.Length, collapse = " & ")
paste(out$Coverage.Probability, collapse = " & ")


############################################
library(dplyr)
library(tidyr)
library(ggplot2)

# Comparing Performances of different tests
# Random Censoring under Exponential Model
sim.RandomCensor <- function(n, lambda, seed = 1234, censoring_dist = "exp", ...) {
    set.seed(seed)
    
    Ti <- rexp(n, rate = lambda)
    if (censoring_dist == "exp") {
        Ci <- rexp(n, ...)
    } else if (censoring_dist == "weibull") {
        Ci <- rweibull(n, ...)
    } else if (censoring_dist == "unif") {
        Ci <- runif(n, ...)
    } else {
        stop("Unknown Censoring distribution")
    }
    
    Xi <- pmin(Ti, Ci)
    del_i <- ifelse(Ti < Ci, 1, 0)
    
    return(data.frame("Observation" = Xi, "Indicator" = del_i))
}

# Obtain power of 4 asymp. tests for H0: lambda = lambda0 vs H1: lambda != lambda0
do.sim_test <- function(B, n, lambda0, lambda_alt, seed = 1234, censoring_dist = "exp", censor_param = NULL) {
    # B = number of resamples
    set.seed(seed)
    data_seeds <- sample(1:1e6, size = B)   # seeds to use to generate data
    
    pb <- txtProgressBar(min = 0, max = B, style = 3)  # initialize a progress bar
    
    lambda_test <- c(lambda0, lambda_alt)
    power_mat <- array(0, dim = c(B, length(lambda_test), 3))  # a 3d array to hold the indicator of test result
    crit_val <- qchisq(p = 0.95, df = 1)
    
    for (b in 1:B) {
        for (i in 1:length(lambda_test)) {
            if (censoring_dist == "exp") {
                tau <- censor_param$rate[i]
                df <- sim.RandomCensor(n, lambda = lambda_test[i], seed = data_seeds[b], censoring_dist = "exp", rate = tau)
            } else if (censoring_dist == "unif") {
                theta <- censor_param$theta[i]
                df <- sim.RandomCensor(n, lambda = lambda_test[i], seed = data_seeds[b], censoring_dist = "unif", max = theta)
            } else if (censoring_dist == "weibull") {
                rate <- censor_param$rate[i]
                df <- sim.RandomCensor(n, lambda = lambda_test[i], seed = data_seeds[b], censoring_dist = "weibull", shape = 2, scale = rate)
            }
            d <- sum(df$Indicator)
            if (d == 0) { d <- 1e-8 }
            V <- sum(df$Observation)
            
            Ts <- (d - V*lambda0)^2 * ((d/V)^2 / (d * lambda0^2))
            if (Ts > crit_val) { power_mat[b, i, 1] <- 1 }
            Tw <- (d - V*lambda0)^2/d
            if (Tw > crit_val) { power_mat[b, i, 2] <- 1 }
            Tlr <- 2 * (lambda0 - (d/V))*V - 2*d*log(lambda0/(d/V))
            if (Tlr > crit_val) { power_mat[b, i, 3] <- 1 }
            
        }
        
        setTxtProgressBar(pb, value = b)  # update progress bar
    }
    
    close(pb)
    
    # compute the power curve
    out <- sapply(1:3, FUN = function(i) { colMeans(power_mat[, , i]) })
    out <- data.frame("lambda" = lambda_test, "Score Test" = out[, 1], "Wald Test" = out[, 2], "LRT" = out[, 3])
    out <- out %>% pivot_longer(2:4, names_to = "Test", values_to = "Power")
    
    return(out)
}

# Obtain area and size metric
metric <- function(out) {
    approx.size <- out %>% filter(lambda == 1) %>% select(Test, Size = Power)
    approx.area <- out %>% filter(lambda > 1) %>% group_by(Test) %>% 
        summarise(Area = (2 * sum(Power) - first(Power) - last(Power))/(2 * max(lambda)), .groups = 'drop')
    
    res <- approx.size %>% left_join(approx.area, by = "Test")
    print(paste(round(c(res$Size, res$Area),3), collapse = " & "))
    
    return(res)
}


# censoring under G ~ Exp(tau)
lambda0 <- 1
lambda_alt <- 2:15
censor_prop <- 0.4  # change as required
sample_size <- 5   # change as required
censor_param <- list("rate" = c(lambda0, lambda_alt) * censor_prop / (1 - censor_prop))

out <- do.sim_test(1000, sample_size, lambda0, lambda_alt, seed = 1911, censoring_dist = "exp", censor_param = censor_param)
metric(out)
ggplot(out, aes(x = lambda, y = Power, color = Test)) + geom_line() + geom_point() + xlab(expression(lambda)) + ylab("Power of the Test") + theme_bw() + 
    ggsave(paste0(paste('ExpTest',censor_prop, sample_size, sep = "-"), '.jpg'))


# censoring under G ~ Unif(theta)
lambda0 <- 1
lambda_alt <- 2:15
censor_prop <- 0.4  # change as required
sample_size <- 5  # change as required
censor_param <- list("theta" = 2 * (c(lambda0, lambda_alt)-censor_prop)/c(lambda0, lambda_alt)^2 )

out <- do.sim_test(1000, sample_size, lambda0, lambda_alt, seed = 1911, censoring_dist = "unif", censor_param = censor_param)
metric(out)
ggplot(out, aes(x = lambda, y = Power, color = Test)) + geom_line() + geom_point() + xlab(expression(lambda)) + ylab("Power of the Test") + theme_bw() 


# censoring under G ~ Weibull(tau, 2)
lambda0 <- 1
lambda_alt <- 2:15
censor_prop <- 0.4  # change as required
sample_size <- 5  # change as required

censor_param <- list("rate" = sapply(c(lambda0, lambda_alt), FUN = function(lambda) {
    uniroot(f = function(tau) {
        (1 - censor_prop) - sqrt(pi)/ tau * exp(lambda^2/(2 * tau^2)) * pnorm(lambda/(2*tau),lower.tail = FALSE)
    }, interval = c(0.5, 100))$root
}))

out <- do.sim_test(1000, sample_size, lambda0, lambda_alt, seed = 1911, censoring_dist = "weibull", censor_param = censor_param)
metric(out)
ggplot(out, aes(x = lambda, y = Power, color = Test)) + geom_line() + geom_point() + xlab(expression(lambda)) + ylab("Power of the Test") + theme_bw() 








