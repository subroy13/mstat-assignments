library(ggplot2)

######################################################
# Question 2
set.seed(1911)   # my roll number
B <- 1000
X <- 1
Z <- rnorm(B)

# do this for B samples to converge to stationary distribution
for (t in 2:B) {
    X <- 0.1 * X + Z[t] + 0.4 * Z[t-1]
}

# generate and store ARMA(1, 1) samples
arma_samples <- numeric(40)
arma_samples[1] <- X
Z <- rnorm(40)
for (t in 2:40) {
    arma_samples[t] <- 0.1 * arma_samples[t-1] + Z[t] + 0.4 * Z[t-1] 
}

arma_samples

# plot the sample values
ggplot() + 
    geom_line(aes(x = 1:40, y = arma_samples), size = 1) +
    geom_point(aes(x = 1:40, y = arma_samples), size = 4, color = "red", pch = 12) + 
    theme_bw() + xlab("Time") + ylab("Sample")

# Auto correlation values computation
sample_mean <- sum(arma_samples) / 40
sample_mean

autocovar <- numeric(8)   # obtain auto-covariances
for (h in 0:7) {
    autocovar[h+1] <- sum((arma_samples[1:(40-h)] - sample_mean) * (arma_samples[(h+1):40] - sample_mean))/40
}
autocovar

autocor <- autocovar[2:8] / autocovar[1]
autocor   # sample autocorrelation

# calculate w_11 using Bartlett's formula
w11 <- sum( (autocor[2:7] * c(1, autocor[1:5]) - 2 * autocor[1] * autocor[1:6])^2 )
w11

autocor[1] - 1.96 * sqrt(w11/40)  # lower bound
autocor[1] + 1.96 * sqrt(w11/40)  # upper bound


# True value of autocorrelation
(0.4 + 0.1) * (1 + 0.4 * 0.1) / (1 + 2 * 0.1 * 0.4 + 0.4^2)




###################################
# QUESTION 3

set.seed(1911)
iid_seq <- sample(0:9, size = 40, replace = TRUE)   # generate 40 samples and treat at iid seq
iid_seq

# plot the iid sample values
ggplot() + 
    geom_line(aes(x = 1:40, y = iid_seq), size = 1) +
    geom_point(aes(x = 1:40, y = iid_seq), size = 4, color = "red", pch = 12) + 
    theme_bw() + xlab("Time") + ylab("(IID) Samples")


sample_mean <- sum(iid_seq) / 40
sample_mean

autocovar <- numeric(8)   # obtain auto-covariances
for (h in 0:7) {
    autocovar[h+1] <- sum((iid_seq[1:(40-h)] - sample_mean) * (iid_seq[(h+1):40] - sample_mean))/40
}
autocovar

autocor <- autocovar[2:8] / autocovar[1]
autocor   # sample autocorrelation

# Portmanteau's test
Q <- 40 * sum(autocor^2)
qchisq(0.95, df = 7)    

# Ljung and Box's test
QLB <- 40 * (40 + 2) * sum(autocor^2 / c(39:33) )
    
# Mcleod and Li's test
sq.autocovar <- numeric(8)   # obtain auto-covariances
for (h in 0:7) {
    sq.autocovar[h+1] <- sum((iid_seq[1:(40-h)]^2 - mean(iid_seq^2) ) * 
                              (iid_seq[(h+1):40]^2 - mean(iid_seq^2) ))/40
}
sq.autocor <- sq.autocovar[2:8] / sq.autocovar[1]
sq.autocor
QML <- 40 * (40 + 2) * sum(sq.autocor^2 / c(39:33) )


# Rank test
P <- sum(sapply(1:39, FUN = function(i) { 
        sum(iid_seq[i] < iid_seq[(i+1):40])  # count the number of j's from which i-th number is smaller
    } ))
muP <- 40 * (40 - 1) / 4
sigmaP.sq <- 40 * (40 - 1) * (2 * 40 + 5) / 72

abs((P - muP)/ sqrt(sigmaP.sq))










