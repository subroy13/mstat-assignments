library(gRbase)
library(gRim)

data("reinis")
str(reinis)

dm.sat <- dmod( ~ . ^ ., reinis) # fit saturated model
dm.null <- dmod( ~ . ^ 1, reinis)   # fit independence model

# forward model selection using AIC
fit.forward <- stepwise(dm.null, criterion = "aic", direction = "forward", type = "unrestricted", k = 2)
summary(fit.forward)
iplot(fit.forward)



# using bic
fit.forward <- stepwise(dm.null, criterion = "aic", direction = "forward", 
                        type = "unrestricted", k = log(sum(reinis)))


iplot(fit.forward)


# backward model selection
fit.backward <- stepwise(dm.sat, criterion = "aic", direction = "backward", type = "unrestricted", k = log(sum(reinis)))
summary(fit.backward)

iplot(fit.backward)


####################################################


deathpenalty <- array(c(19, 132, 11, 52, 0, 9, 6, 97), dim = c(2,2,2),
                      dimnames = list("Death.penalty" = c("Yes", "No"),
                                      "Defendant.race"= c("White","Black"),
                                      "Victim.race" = c("White","Black")))

dm.sat <- dmod( ~ . ^ ., as.table(deathpenalty)) # fit saturated model
dm.null <- dmod( ~ . ^ 1, as.table(deathpenalty))   # fit independence model

# forward model selection using AIC
fit.forward <- stepwise(dm.null, criterion = "aic", direction = "forward", type = "unrestricted", k = 2)
iplot(fit.forward)

# backward model selection using AIC
fit.backward <- stepwise(dm.sat, criterion = "aic", direction = "backward", type = "unrestricted", k = 2)
iplot(fit.backward)

# forward model selection using BIC
fit.forward <- stepwise(dm.null, criterion = "aic", direction = "forward", type = "unrestricted", k = log(sum(deathpenalty)))
iplot(fit.forward)

# backward model selection using BIC
fit.backward <- stepwise(dm.sat, criterion = "aic", direction = "backward", type = "unrestricted", k = log(sum(deathpenalty)))
iplot(fit.backward)



IPS <- function(form, data, maxiter = 1000, tol = 1e-05) {
    # initialize an array of 1's only
    tempdata <- array(1, dim = dim(data), dimnames = dimnames(data))
    form.vars <- labels(terms(formula(form)))
    
    current_error <- Inf
    current_iter <- 0
    while (current_error > tol & current_iter < maxiter) {
        # store the current table for computing the error later
        current_tab <- tempdata
        
        # for each given margin, perform the updation
        for (var in rev(form.vars))   {
            true_margins <- ar_marg(data, formula(paste0("~", var)))   # compute the true marginals
            current_margins <- ar_marg(tempdata, formula(paste0("~", var)))   # compute the current margins
            
            # expand them to higher dimension as original so that multiplication can be performed
            true_margins <- ar_expand(true_margins, dimnames(data))    
            current_margins <- ar_expand(current_margins, dimnames(data))
            
            tempdata <- (tempdata * true_margins)/current_margins
        }
        
        current_iter <- current_iter+ 1        # increase number of iteration
        current_error <- max(abs(tempdata - current_tab))   # compute the error
    }
    
    return(list("MLE Table" = tempdata, "Iteration" = current_iter))
}

IPS( ~ Victim.race + Death.penalty + Defendant.race + Victim.race:Death.penalty + 
        Death.penalty:Defendant.race + Victim.race:Defendant.race, deathpenalty)


IPS( ~ Victim.race + Death.penalty + Defendant.race + Victim.race:Defendant.race +
         Victim.race:Death.penalty, deathpenalty)


a = dmod(~ Victim.race:Death.penalty + Death.penalty:Defendant.race + Victim.race:Defendant.race, 
         as.table(deathpenalty))
logLik.iModel(a)

b = dmod(~ Victim.race:Death.penalty + Victim.race:Defendant.race, 
         as.table(deathpenalty))
logLik.iModel(b)











