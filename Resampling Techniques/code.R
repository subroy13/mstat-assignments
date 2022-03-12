#####
# Model 1:
# y_i = 5x_i^2 + 7x_i + e_i; x_i, e_i ~ Unif[-1, 1] and Unif[-1/2, 1/2]

# Model 2:
# y_i = 5x_i^2 + 7x_i + (x_i + 1)e_i; x_i, e_i ~ Unif[-1, 1] and Unif[-1/2, 1/2]

# K(x) = 3/4 (1 - x)^2 I{|x| < 1}
# n = 50
# lambda = [n * log(n)]^(-1/5) = 0.35
# B = 1000
# simulation_size = 100


sim.model <- function(n, seed = 1234, model = 1) {
    set.seed(seed)
    x <- runif(n, min = (-1), max = 1)
    error <- runif(n, min = (-1/2), max = (1/2))
    if (model == 1) {
        y <- 5 * x^2 + 7 * x + error
    } else {
        y <- 5 * x^2 + 7 * x + ((x + 1) * error)
    }
    return(data.frame(x = x, y = y))
}

kernel.estimate.norm <- function(x, data, lambda = 0.35, conf = 0.9) {
    f <- ecdf(data$x)
    norm.x <- (f(x) - f(data$x)) / lambda     # apply ecdf
    kern.weights <- (3/4) * (1 - norm.x^2) * ifelse(abs(norm.x) > 1, 0, 1)
    
    estimate <- weighted.mean(data$y, kern.weights)
    
    var_est <- (weighted.mean(data$y^2, kern.weights) - estimate^2) * (3/5)   # (3/5) = integral(K(u)du)
    len_est <- qnorm((1 - conf)/2, lower.tail = FALSE) * sqrt(var_est) / sqrt(nrow(data) * lambda)
    
    result <- c(estimate, len_est ) 
    names(result) <- c("Estimate", "Halflength")
    
    return(result)
}

get.normal.result <- function(x, n = 50, seed = 1234, model = 1, sim_size = 100, lambda = 0.35, conf = 0.9) {
    set.seed(seed)
    data_seeds <- sample(1:1e9, size = sim_size)   # generate some seeds to be used to generate data
    
    cp <- 0   # coverage indicator
    el <- numeric(sim_size)   # expected length
    
    for (i in 1:sim_size) {
        data <- sim.model(n = n, seed = data_seeds[i], model = model)  # simulate the data
        res <- kernel.estimate.norm(x, data, lambda = lambda, conf = conf)
        true_val <- 5 * x^2 + 7 * x
        
        if ( abs(true_val - res[1]) < res[2] ) {
            cp <- cp + 1   # since it is covered
        } 
        el[i] <- (2 * res[2])
    }
    
    return(list("CP" = cp/sim_size, "EL" = mean(el)))
}

kernel.bootstrap <- function(x, data, seed = 1234, B = 1e3, lambda = 0.35, conf = 0.9) {
    set.seed(seed)
    boot_samples <- sample(1:1e9, size = B)   # generate some seeds to be used to generate data
    
    boot_ests <- numeric(B)
    for (b in 1:B) {
        index <- sample(1:nrow(data), size = nrow(data), replace = TRUE)
        boot_data <- data[index, ]
        
        f <- ecdf(boot_data$x)
        norm.x <- (f(x) - f(boot_data$x)) / lambda     # apply ecdf
        kern.weights <- (3/4) * (1 - norm.x^2) * ifelse(abs(norm.x) > 1, 0, 1)
        
        boot_ests[b] <- weighted.mean(boot_data$y, kern.weights)
    }
    
    return(quantile(boot_ests, probs = c( (1-conf)/2, 1 - (1-conf)/2 ) ))
}

get.boot.result <- function(x, n = 50, seed = 1234, model = 1, sim_size = 100, lambda = 0.35, conf = 0.9, B = 1e3) {
    set.seed(seed)
    data_seeds <- sample(1:1e9, size = sim_size)   # generate some seeds to be used to generate data
    boot_seeds <- sample(1:1e9, size = sim_size)
    
    cp <- 0   # coverage indicator
    el <- numeric(sim_size)   # expected length
    
    pb <- txtProgressBar(max = sim_size, style = 3)
    for (i in 1:sim_size) {
        data <- sim.model(n = n, seed = data_seeds[i], model = model)  # simulate the data
        
        f <- ecdf(data$x)
        norm.x <- (f(x) - f(data$x)) / lambda     # apply ecdf
        kern.weights <- (3/4) * (1 - norm.x^2) * ifelse(abs(norm.x) > 1, 0, 1)
        m_n0 <- weighted.mean(data$y, kern.weights)
        
        res <- kernel.bootstrap(x, data, seed = boot_seeds[i], B = B, lambda = lambda, conf = conf)
        conf.int <- c(2*m_n0 - res[2], 2*m_n0 - res[1])
        
        true_val <- 5 * x^2 + 7 * x
        
        if ( (true_val < conf.int[2]) & (true_val > conf.int[1]) ) {
            cp <- cp + 1   # since it is covered
        } 
        el[i] <- conf.int[2] - conf.int[1]
        
        setTxtProgressBar(pb, value = i)
    }
    close(pb)
    return(list("CP" = cp/sim_size, "EL" = mean(el)))
}


###########################################################

get.normal.result(x = -0.7, model = 1, seed = 1911)
get.normal.result(x = -0.4, model = 1, seed = 1911)
get.normal.result(x = -0.2, model = 1, seed = 1911)
get.normal.result(x = 0, model = 1, seed = 1911)
get.normal.result(x = 0.2, model = 1, seed = 1911)
get.normal.result(x = 0.4, model = 1, seed = 1911)
get.normal.result(x = 0.7, model = 1, seed = 1911)


get.normal.result(x = -0.7, model = 2, seed = 1911)
get.normal.result(x = -0.4, model = 2, seed = 1911)
get.normal.result(x = -0.2, model = 2, seed = 1911)
get.normal.result(x = 0, model = 2, seed = 1911)
get.normal.result(x = 0.2, model = 2, seed = 1911)
get.normal.result(x = 0.4, model = 2, seed = 1911)
get.normal.result(x = 0.7, model = 2, seed = 1911)

###########################################################

get.boot.result(x = -0.7, model = 1, seed = 1911)
get.boot.result(x = -0.4, model = 1, seed = 1911)
get.boot.result(x = -0.2, model = 1, seed = 1911)
get.boot.result(x = 0, model = 1, seed = 1911)
get.boot.result(x = 0.2, model = 1, seed = 1911)
get.boot.result(x = 0.4, model = 1, seed = 1911)
get.boot.result(x = 0.7, model = 1, seed = 1911)







