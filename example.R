#######################################################################
##### Example for the simulation study and the real data analysis #####
#######################################################################
library(dplyr)
library(mnormt)
library(Rlab)
library(nlme)
library(CompQuadForm)
library(glmnet)
library(expm)
library(huge)
library(optparse)
library(doParallel)
library(foreach)
library(rags2ridges)
library(corpcor)
library(abind)
library(PDSCE)

# set the working directory
source("pairedtest.R")
source("generate_sim_data.R")

# ------------------------------- setup parameters -------------------------------- #
##### data generation (for simulation, skipped for the real data analysis)
n <- 15
p <- 200
q <- 50


type_t <- 3
type_s <- 1
corr_type_s <- 1
corr_type_t <- 0
mod_s <- 7
mod_t <- 15
df <- "infty"
rho_input <- 0.2


##### algorithm
m <- 0.5
alpha <- 0.01
z <- seq(0 : ((floor(sqrt(4 * log(p))) + 2) * 100))/100
len_z <- length(z)
lambda_sets <- seq(3, 20)

indpt <- 1 # 0: independent case; 1: paired case
method_rm <- 2


# ------------------------------- running algorithm -------------------------------- #
###################
## generate data ## 
###################
## simulation
list_info1 <- InfoGenerate(p, q, type = c(type_t, type_s))
list_info2 <- InfoGenerate(p, q, type = c(type_t + 1, 0))

list_Omega_s <- RandOnSpaceInfo(list_info1$Omega_s, m = m)
Omega1 <- list_Omega_s$Omega1
Omega2 <- list_Omega_s$Omega2

Sigma1 <- list_info1$Sigma_t
Sigma2 <- list_info2$Sigma_t

Z <- DataGenerate(n, list(Sigma1, Sigma2), list(Omega1, Omega2), corr_type_s = corr_type_s, corr_type_t = corr_type_t, rho_input = rho_input, mod_t = mod_t, mod_s = mod_s, df = df)
X1 <- Z[[1]]
X2 <- Z[[2]]

# truth
Wtrue <- Omega1/outer(diag(Omega1), diag(Omega1), "/") - Omega2/outer(diag(Omega2), diag(Omega2), "/")
Wtrue[abs(Wtrue) <= 0.001] <- 0

## real data
# load("img_conversion.RData")
# n <- length(convert_img)
# q <- dim(convert_img[[1]][[1]])[1]
# p <- dim(convert_img[[1]][[1]])[2] - 1
# X1 <- lapply(convert_img, function(x) c(t(do.call(x[[1]], what = cbind)[, -1]))) %>% do.call(what = cbind)
# X1 <- (X1 - rowMeans(X1)) %>% array(dim = c(p, q, n))
# X2 <- lapply(convert_img, function(x) c(t(do.call(x[[2]], what = cbind)[, -1]))) %>% do.call(what = cbind)
# X2 <- (X2 - rowMeans(X2)) %>% array(dim = c(p, q, n))     


#######################
## testing procedure ## 
#######################
# remove the time info
Y1 <- TimeRemoval(X1, method_rm)
Y2 <- TimeRemoval(X2, method_rm)


A_time <- TimeSignChange(InfoGenerate_Time(q, type_t = corr_type_t), mod_t)
result <- sapply(lambda_sets, function(lambda){
    print(paste0("Lambda: ", lambda))
    # Test Statistics
    Wstat <- TestStat(Y1, Y2, lambda, indpt, X1, X2)$W
    
    # rejection threshold
    rr <- sapply(seq(len_z), function(k){
        p * (p - 1) * (1 - pnorm(z[k]))/max(sum(abs(Wstat) >= z[k])/2, 1)
    })
    
    if(rr[len_z] <= alpha){
        t_threshold <- min(which((rr <= alpha)))
    }else{
        t_threshold <- len_z
    }
    
    
    # empirical size (only for the simulation data since the truth is known)
    ffp <- sum(abs(Wstat[which(Wtrue == 0, arr.ind = TRUE)]) >= z[t_threshold])/max(sum(abs(Wstat) >= z[t_threshold]), 1)
    # select lambda
    sssp <- sapply(seq(10), function(k){
        (sum(abs(Wstat) >= qnorm(1 - k/floor(10/(1 - pnorm(sqrt(log(p)))))))/2/(k/floor(10/(1 - pnorm(sqrt(log(p))))) * (p * (p-1))) - 1)^2
    }) %>% sum
    # empirical power (only for the simulation data since the truth is known)
    fpp <- (max(sum(abs(Wstat) >= z[t_threshold]), 1) - sum(abs(Wstat[which(Wtrue == 0, arr.ind = TRUE)]) >= z[t_threshold]))/((sum(Wtrue != 0)))
    
    c(ffp, sssp, fpp)
})

best_index <- which.min(result[2, ])
print(paste0("Empirical size: ", result[1, best_index], "; empirical power: ", result[3, best_index]))
