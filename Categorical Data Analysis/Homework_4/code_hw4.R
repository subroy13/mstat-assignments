# Question 1

data <- matrix(c(14, 3, 3, 13, 1, 2, 6, 1, 1), nrow = 3, byrow = T, 
               dimnames = list("Gender" = c("Male","Female","Other"),
                               "Handedness" = c("Right","Left","Both")))
chisq.test(data)

DescTools::GTest(data)


meltdata <- reshape2::melt(data)
meltdata <- as.data.frame(lapply(meltdata, rep, meltdata$value))[, -3]
meltdata$Gender <- factor(meltdata$Gender)
meltdata$Handedness <- factor(meltdata$Handedness)
meltdata

chisq.stat <- function(tab) {
    n <- sum(tab)
    rowmars <- apply(tab, 1, sum)
    colmars <- apply(tab, 2, sum)
    exp.tab <- rowmars %*% t(colmars)/n
    
    return(sum((tab - exp.tab)^2 / exp.tab))
}

lrt.stat <- function(tab) {
    n <- sum(tab)
    rowmars <- apply(tab, 1, sum)
    colmars <- apply(tab, 2, sum)
    exp.tab <- rowmars %*% t(colmars)/n
    
    return(2*sum(tab * log(tab / exp.tab), na.rm = T))
    
}

perm.stat <- function(data, stat, nperms = 10e4, seed = 1234) {
    
    if (! (stat %in% c("chisq","lrt"))) {
        stop("Not Implemented yet! Only 'chisq' and 'lrt' is available.")
    }
    else {
        set.seed(seed)
        x <- sapply((1:nperms), function(a){ sample(data[, 1], size = nrow(data)) })
        x <- t(x)
        x <- unique(x)
        
        if (stat == "chisq") {
            y <- apply(x, 1, function(a) {
                tab <- table(a, data[, 2])
                return(chisq.stat(tab))
            })
        }
        else if (stat == "lrt") {
            y <- apply(x, 1, function(a) {
                tab <- table(a, data[, 2])
                return(lrt.stat(tab))
            })
        }
        
        return(y)
    }
}


y <- perm.stat(meltdata, stat = "chisq", nperms = 5000, seed = 1911)
head(y)
paste("Approximated p-Value is", sum(y > chisq.stat(data))/5000)

library(ggplot2)

ggplot(as.data.frame(y), aes(x = y)) + 
    geom_histogram(aes(y=..density..), color = "black", fill="white", binwidth = 0.5)+
    geom_density(alpha=.2, fill="#FF1111") +
    ylab("Density") + xlab("Chi Square statistic under Null hypothesis of independence") +
    geom_vline(xintercept = chisq.stat(data), color = "red", size = 1, linetype = "dashed")
    

y <- perm.stat(meltdata, stat = "lrt", nperms = 5000, seed = 1911)
head(y)
paste("Approximated p-Value is", sum(y > lrt.stat(data))/5000)


nperms <- c(5, 10, 50, 100, 500, 1000, 5000, 10000)

for (i in nperms) {
    print(i)
    y <- perm.stat(meltdata, stat = "lrt", nperms = i, seed = 1234)
    print(paste("Approximated p-Value for Chi-square is", sum(y > lrt.stat(data))/i))
    
    y <- perm.stat(meltdata, stat = "lrt", nperms = i, seed = 1234)
    print(paste("Approximated p-Value for Chi-square is", sum(y > lrt.stat(data))/i))
    
}




# Question 2

library(elrm)


df <- data.frame(Gender = c("M","F","M","F"), CollegeEd = c("No", "No","Yes","Yes"),
                 WhiteJob = c(1, 1, 7, 6), Total = c(8, 6, 10, 6))

fit <- elrm(WhiteJob / Total ~ Gender + CollegeEd, interest = ~ Gender, dataset = df, 
            iter = 10e5, burnIn = 100)

summary(fit)

fit <- elrm(WhiteJob / Total ~ Gender + CollegeEd, interest = ~ CollegeEd, dataset = df, 
            iter = 10e5, burnIn = 100)

summary(fit)


######

library(logistiX)
Gender <- factor(c(rep("M", 8), rep("F", 6), rep("M", 10), rep("F", 6)))
CollegeEd <- factor(c(rep("No", 14), rep("Yes", 16)))
WhiteJob <- c(1, rep(0, 7), 1, rep(0, 5), rep(1, 7), rep(0, 3), rep(1, 6))
longdf <- data.frame(Gender = Gender, CollegeEd = CollegeEd, WhiteJob = WhiteJob)

x <- model.matrix(WhiteJob ~ Gender + CollegeEd, longdf)[, -1]
fit <- logistiX(x, y = longdf$WhiteJob)
summary(fit)


## 

fit <- glm(cbind(WhiteJob, Total - WhiteJob) ~ Gender + CollegeEd, data = df, family = "binomial")
summary(fit)







