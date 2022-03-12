data <- read.csv('./writers_data.csv')

ctab <- as.matrix(data[, -1])
rownames(ctab) <- data[, 1]
ctab

library(FactoMineR)
fit <- CA(ctab, ncp = 14, graph = TRUE)

head(fit$eig)

plot(1:14, fit$eig[,3], type = "b", lwd = 2, xlab = "Number of Eigenvalues to keep", 
     ylab = "Percentage of Chi-square explained")
abline(h = 80, col = "blue", lty = 2)
abline(h = 90, col = "purple", lty = 2)
abline(h = 95, col = "red", lty = 2)
