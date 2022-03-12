library(readxl)
library(dplyr)
library(readr)
library(tidyr)
library(reshape2)
library(survival)
library(ggplot2)
library(survminer)

setwd('./Assignment 2/')

##############################
# READING THE RAW FILES
rawFilePaths2018 <- list.files('./datasets/J Thomas raw/2018/')
rawFilePaths2019 <- list.files('./datasets/J Thomas raw/2019/')

rawFiles2018 <- lapply(rawFilePaths2018, FUN = function(x){
    read_xls(paste0('./datasets/J Thomas raw/2018/', x), sheet = 1)
})
names(rawFiles2018) <- substr(rawFilePaths2018, 5, 6)

rawFiles2019 <- lapply(rawFilePaths2019, FUN = function(x){
    read_xls(paste0('./datasets/J Thomas raw/2019/', x), sheet = 1)
})
names(rawFiles2019) <- substr(rawFilePaths2019, 5, 6)

df2018 <- bind_rows(rawFiles2018, .id = "Week")
df2019 <- bind_rows(rawFiles2019, .id = "Week")
df_raw <- df2018 %>% 
    mutate(Year = 2018) %>%
    bind_rows(df2019 %>% mutate(Year = 2019))


#############################
# BASIC PREPROCESSING

sourcedf <- read_csv('./datasets/processed/TeaSource.csv')
df <- df_raw %>% 
    mutate(Volume = Net * Pkgs) %>%         # Volume = Net weight of packages x Number of packages
    select(Year, Week, Garden, Grade, Volume, Val, Price) %>%
    mutate(Week = as.numeric(Week), Val = as.numeric(Val), Price = as.numeric(Price)) %>%
    filter(!is.na(Price)) %>%    # remove unsold packages
    filter(Val > 10, Price > 10) %>%   # remove data where valuation or price wrongly recorded as < 10 Rs.
    left_join(sourcedf, by = "Garden")

# Use the grade based clustering mentioned in https://arxiv.org/pdf/2005.02814.pdf
GradeClusterList <- c(
    "OD" = 1, "OD(S)" = 1, "OPD1" = 1,
    "OCD" = 2, "OCD1" = 2, "OD1" = 2,
    "D(F)" = 3, "CD1" = 3, "CHD1" = 3, "RD1" = 3,
    "D" = 4, "D(SPL)" = 4, "CD" = 4, "CHD" = 4, "CHU" = 4, "CHUR" = 4, "PD" = 4, "PD(SPL)" = 4,
    "PD(FINE)" = 5, "GTDUST" = 5,
    "OPD" = 6, "OPD(CLONAL)" = 6, "ORD" = 6, "D1" = 6, "D1(SPL)" = 6, "PD1" = 6, "PD1(SPL)" = 6
)
df$GradeCluster <- GradeClusterList[df$Grade]

# Use the source based clustering mentioned in https://arxiv.org/pdf/2005.02814.pdf
SourceClusterList <- c(
    "Darjeeling" = 1, "Cooch Bihar" = 1, "Uttar Dinajpur" = 1, "Jalpaiguri" = 1, "Alipurduar" = 1,
    "Karimganj" = 2, "Hailakandi" = 2,
    "Bongaigaon" = 3, "Cachar" = 3, "Udalguri" = 3, "Darrang" = 3, "Dima Hasao" = 3,
    "Lakhimpur" = 4, "Nagaon" = 5,
    "Sivasagar" = 6, "Tinsukia" = 6, "Golaghat" = 6, "Jorhat" = 6,
    "Baksa" = 7, "Dibrugarh" = 7, "Sonitpur" = 7
)
df$SourceCluster <- SourceClusterList[df$District]
df <- df[!((df$GradeCluster == 3)&(df$SourceCluster==4)), ]
df <- df[!((df$GradeCluster == 5)&(df$SourceCluster==1)), ]

write_csv(df, './datasets/processed/J Thomas 2018-2019.csv')



####################################
# STATISTICAL ANALYSIS STARTS HERE

df <- read_csv('./datasets/processed/J Thomas 2018-2019.csv')
df <- df %>% mutate(GradeCluster = as.character(GradeCluster), 
                    SourceCluster = as.character(SourceCluster))
str(df)  # 9 columns, approximately 26000 samples

# Uncensored are those where Price <= Valuation
df$Uncensored <- ifelse(df$Val >= df$Price, 1, 0)

##################################
# Correlation analysis
cor(df$Val, df$Price)  # Very high correlation
cor(df$Val, df$Price, method = "spearman")  # Very high rank correlation

cordf <- df %>% 
    group_by(GradeCluster, SourceCluster) %>% 
    summarise(Observation = n(), Correlation = cor(Val, Price))

ggplot(cordf, aes(x = GradeCluster, y = SourceCluster, fill=Correlation)) + 
    geom_raster() + geom_point(aes(size = Observation)) + theme_bw() +
    scale_fill_viridis_c()

########################################
# Survival Analysis

f1 <- survfit(Surv(Price, Uncensored) ~ SourceCluster, df)
ggsurvplot(f1, df, facet.by = "GradeCluster", palette = "jco", conf.int = TRUE, 
           ylab = "Bid raising probability", xlab = "Price (or Bid)", 
           legend.title = "SourceCluster", surv.median.line = "hv", censor = FALSE)

f2 <- survfit(Surv(Price, Uncensored) ~ SourceCluster + GradeCluster, df)
f2


##########
# Cox PH (Model 1)

coxmod <- coxph(Surv(Price, Uncensored) ~ log(Val) + log(Volume) + GradeCluster + SourceCluster, df)
coxmod

base_cumhaz <- basehaz(coxmod)
risk <- predict(coxmod, type = "risk")  # gives exp(zi * beta)
Lambda0 <- base_cumhaz$hazard
names(Lambda0) <- base_cumhaz$time

coxsnell_resid <- Lambda0[as.character(df$Price)] * risk  # Cox-snell residual

# Obtain its Nelson-Aalen estimate and plot it
f1 <- survfit(Surv(coxsnell_resid, Uncensored) ~ 1, df)
ggsurvplot(f1, df, fun = "cumhaz", xlab = "Cox snell residuals", xlim = c(0, 4), break.x.by = 1)$plot +
    geom_abline(aes(intercept = 0, slope = 1), color = "blue", linetype = "dashed")


##########
# Cox PH (Model 2)
coxmod <- coxph(Surv(Price, Uncensored) ~ log(Val) + log(Volume) + strata(GradeCluster, SourceCluster), df)
coxmod

base_cumhaz <- basehaz(coxmod)
risk <- predict(coxmod, type = "risk")  # gives exp(zi * beta)
coxsnell_resid <- numeric(nrow(df))   # computing Cox-Snell residuals
df$Strata <- paste(df$GradeCluster, df$SourceCluster, sep = ", ")
for (i in 1:nrow(df)) {
    coxsnell_resid[i] <- risk[i] * base_cumhaz$hazard[(base_cumhaz$time == df$Price[i]) & (base_cumhaz$strata == df$Strata[i])]
}

# Obtain its Nelson-Aalen estimate and plot it
f1 <- survfit(Surv(coxsnell_resid, Uncensored) ~ 1, df)
ggsurvplot(f1, df, fun = "cumhaz", xlab = "Cox snell residuals", xlim = c(0, 4), break.x.by = 1)$plot +
    geom_abline(aes(intercept = 0, slope = 1), color = "blue", linetype = "dashed")


############################
# Semi-parametric AFT model, Buckley James estimator (Model 3)
library(rms)

bjmod <- bj(Surv(Price, Uncensored) ~ log(Val) + log(Volume) + GradeCluster + SourceCluster, df, control = list("iter.max" = 100), x = TRUE, y = TRUE)
bjmod
latex(bjmod)

# Check model assumptions
X <- model.matrix(formula(~log(Val) + log(Volume) + GradeCluster + SourceCluster), df)
fitted_vals <- X %*% bjmod$coefficients
resid_vals <- log(df$Price) - fitted_vals

ggplot() + geom_point(aes(x = fitted_vals, y = resid_vals + 0.1), alpha = 0.1) +
    xlim(4, 6) + ylim(-0.5, 0.5) + theme_bw() + 
    xlab(expression(Fitted~Values~alpha+Z[i]~beta)) + 
    ylab(expression(Residuals~log~(T)~-~alpha~-~Z[i]~beta)) +
    geom_abline(intercept = 0, slope = 0, color = "red", linetype = "dashed")







