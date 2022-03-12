################
# QUESTION 4

FFT_put_price <- function(K, maturity, alpha, param_list) {
    eta <- param_list[["eta"]]
    N <- param_list[["N"]]
    lambda <- (2 * pi)/(N * eta)
    vj <- eta * (1:N - 1)
    
    # VG process char function of log(S(t))
    char_VG <- function(u, t, param_list) {
        theta <- param_list[["theta"]]
        nu <- param_list[["nu"]]
        sigma <- param_list[["sigma"]]
        S0 <- param_list[["S0"]]
        r <- param_list[["r"]]
        q <- param_list[["q"]]
        
        phi_u <- (1 - (0+1i)*u*theta*nu + sigma^2 * u^2 * nu/2)
        phi_u <- (1/phi_u)^(t/nu)
        phi_i <- as.complex(1 - theta * nu - sigma^2 * nu/2)
        phi_i <- (1 / phi_i)^(t/nu)
        ch <- exp((0+1i)*u*(log(S0) + (r-q)*t )) * phi_u / phi_i^((0 + 1i)*u)
        
        return(ch)
    }
    
    xvec <- rep(eta, N)
    xvec[1] <- xvec[1]/2
    xvec <- xvec * exp(-param_list[["r"]] * maturity) / ((alpha + (0+1i)*vj) * (alpha + (0+1i)*vj + 1))
    xvec <- xvec * exp(-(0+1i)*(log(K) - lambda * N/2)*vj)
    xvec <- xvec * char_VG(vj - (alpha + 1)*(0+1i), maturity, param_list)
    
    yvec <- fft(xvec)
    p_seq <- log(K) - lambda*N/2 + (0:(N-1))*lambda
    price <- Re(yvec) * exp(-alpha*p_seq)/pi
    
    price <- price[N/2]
    
    return(price)
}

simulate_put_price <- function(K, maturity, param_list, seed = 1234) {
    B <- 1e6  # number of resamples to take the average
    
    S0 <- param_list[["S0"]]
    r <- param_list[["r"]]
    q <- param_list[["q"]]
    sigma <- param_list[["sigma"]]
    nu <- param_list[["nu"]]
    theta <- param_list[["theta"]]
    
    set.seed(seed)
    gamma_g <- rgamma(B, shape = maturity/nu, scale = nu)
    Wt <- rnorm(B, mean = 0, sd = sqrt(gamma_g))
    Xt <- theta * gamma_g + sigma * Wt
    w <- (1/nu) * log(1 - sigma^2 * nu/2 - theta * nu)
    St <- S0 * exp((r - q) * maturity + Xt + w*maturity)
    
    values <- ifelse(K > St, (K - St), 0)
    price <- mean(values) * exp(-r*maturity)
    return(price)
}

param_list <- list("S0" = 100, "r" = 0.0475, "q" = 0.0125, "sigma" = 0.25, "nu" = 0.5, "theta" = -0.3,
                   "eta" = 0.05, "N" = 2^12)

FFT_put_price(K = 105, maturity = 1, alpha = -1.25, param_list = param_list)
FFT_put_price(K = 105, maturity = 1, alpha = -1.5, param_list = param_list)
FFT_put_price(K = 105, maturity = 1, alpha = -2, param_list = param_list)

simulate_put_price(K = 105, maturity = 1, param_list = param_list, seed = 1911)


#################
# Question 3

# uniform mesh grid on S
upoutPriceType1 <- function(K, B, sigma, param_list) {
    delS <- param_list[["delS"]]
    delT <- param_list[["delT"]]
    S0 <- param_list[["S0"]]
    r <- param_list[["r"]]
    q <- param_list[["q"]]
    
    Smin <- 0
    Smax <- B
    N <- (Smax - Smin)/delS
    S <- Smin + 1:(N-1) * delS
    U <- pmax((S - K), 0)
    
    maturity <- 1
    nT <- maturity / delT
    for (i in 1:nT) {
        A_mat <- matrix(0, nrow = (N-1), ncol = (N-1))
        
        # fill up the rows
        for (j in 1:(N-1)) {
            if (j > 1) {
                A_mat[j,(j-1)] <- (-sigma^2*S[j]^2/(2 * delS^2) + (r-q)*S[j]/(2 * delS))*delT
            }
            if (j < (N-1)) {
                A_mat[j,(j+1)] <- (-sigma^2*S[j]^2/(2 * delS^2) - (r-q)*S[j]/(2 * delS))*delT
            }
            A_mat[j,j] <- (1/delT + sigma^2*S[j]^2/delS^2 + r)*delT
        }
        
        # update U vector by implicit scheme
        U <- solve(A_mat, U)
    }
    
    # the price corresponds to U(S0, 0)
    S0index <- (S0 - Smin)/delS
    return(U[S0index])
}

# uniform mesh grid on xi
upoutPriceType2 <- function(K, B, sigma, param_list) {
    delS <- param_list[["delS"]]
    delxi <- param_list[["delxi"]]
    delT <- param_list[["delT"]]
    S0 <- param_list[["S0"]]
    r <- param_list[["r"]]
    q <- param_list[["q"]]
    
    Smin <- 0
    Smax <- B
    alpha <- (Smax - Smin)/20
    c1 <- asinh((Smax - B)/alpha)
    c2 <- asinh((Smin - B)/alpha)
    N <- 1/delxi
    xi <- 1:(N-1) * delxi  # uniform grid on xi
    S <- B + alpha * sinh(c1 * xi + c2 * (1-xi))
    Sprime <- alpha * (c1-c2) * cosh(c1 * xi + c2 * (1-xi))  # dS/dxi
    Sdprime <- alpha * (c1 - c2)^2 * sinh(c1 * xi + c2 * (1-xi))  # d2S/dxi2
    U <- pmax((S - K), 0)
    
    maturity <- 1
    nT <- maturity / delT
    for (i in 1:nT) {
        A_mat <- matrix(0, nrow = (N-1), ncol = (N-1))
        # fill up the rows
        for (j in 1:(N-1)) {
            alphajk <- (r-q)*S[j]/Sprime[j] - (sigma^2*S[j]^2*Sdprime[j]/(2 * Sprime[j]^3))
            betajk <- S[j]^2 * sigma^2 / (2 * Sprime[j]^2)

            if (j > 1) {
                A_mat[j,(j-1)] <- delT * (alphajk/(2 * delxi) - betajk / delxi^2 )
            }
            if (j < (N-1)) {
                A_mat[j,(j+1)] <- delT * (-alphajk/(2*delxi) - betajk / delxi^2)
            }
            A_mat[j,j] <- delT * (r + 1/delT + 2 * betajk/delxi^2)
        }
        
        # update U vector by implicit scheme
        U <- solve(A_mat, U)
    }
    
    # the price corresponds to U(S0, 0)
    S0index <- round((asinh((S0-B)/alpha) - c2)/((c1 - c2) * delxi))
    return(U[S0index])
}


# closed form price
upoutPriceClosed <- function(K, B, sigma, param_list) {
    S0 <- param_list[["S0"]]
    r <- param_list[["r"]]
    q <- param_list[["q"]]
    maturity <- 1
    sigmat <- sigma * sqrt(maturity)
    
    b <- (r - q)
    mu <- (b - sigma^2/2)/sigma^2
    lambda <- sqrt(mu^2 + 2*r/sigma^2)
    x1 <- log(S0/K)/sigmat + (1 + 2*mu)*sigmat
    x2 <- log(S0/B)/sigmat + (1+2*mu)*sigmat
    y1 <- log(B^2/(S0 * K))/sigmat + (1 + 2*mu)*sigmat
    y2 <- log(B/S0)/sigmat + (1 + 2*mu)*sigmat
    
    A <- S0*exp((b-r)*maturity)*pnorm(x1) - K*exp(-r*maturity)*pnorm(x1 - sigmat)
    B_new <- S0*exp((b-r)*maturity)*pnorm(x2) - K*exp(-r*maturity)*pnorm(x2 - sigmat)
    C <- S0*exp((b-r)*maturity)*((B/S0)^(2*(mu+1)))*pnorm(-y1) - K*exp(-r*maturity)*((B/S0)^(2*mu))*pnorm(-y1 + sigmat)
    D <- S0*exp((b-r)*maturity)*((B/S0)^(2*(mu+1)))*pnorm(-y2) - K*exp(-r*maturity)*((B/S0)^(2*mu))*pnorm(-y2 + sigmat)
    
    price <- abs(A - B_new + C - D)
    
    return(price)
}

param_list <- list("delS" = 5, "delxi" = 1/20, "delT" = 1/365, "S0" = 100, "r" = 0.0475, "q" = 0.0175)
upoutPriceClosed(K = 115, B = 125, sigma = 0.15, param_list = param_list)
upoutPriceType1(K = 115, B = 125, sigma = 0.15, param_list = param_list)
upoutPriceType2(K = 115, B = 125, sigma = 0.15, param_list = param_list)

upoutPriceClosed(K = 115, B = 125, sigma = 0.30, param_list = param_list)
upoutPriceType1(K = 115, B = 125, sigma = 0.30, param_list = param_list)
upoutPriceType2(K = 115, B = 125, sigma = 0.30, param_list = param_list)

upoutPriceClosed(K = 115, B = 125, sigma = 0.50, param_list = param_list)
upoutPriceType1(K = 115, B = 125, sigma = 0.50, param_list = param_list)
upoutPriceType2(K = 115, B = 125, sigma = 0.50, param_list = param_list)







