###############################################################################################
################################# Data Generation for simulation ##############################
###############################################################################################

InfoGenerate <- function(p, q, type){
    Sigma_t <- InfoGenerate_Time(q, type[1])
    Omega_s <- InfoGenerate_Space(p, type[2])
    list(Sigma_t = Sigma_t, Omega_s = Omega_s)
}

InfoGenerate_Time <- function(q, type_t){
    ## Time Sigma
    if(type_t == 0){
        Sigma_t <- diag(rep(1, q))
    }else if(type_t == 1){
        Sigma_t <- 0.4^(outer(seq(q), seq(q), FUN = function(x, y) abs(x - y)))
    }else if(type_t == 2){
        Sigma_t <- 0.5^(outer(seq(q), seq(q), FUN = function(x, y) abs(x - y)))
    }else if(type_t == 3){
        Sigma_t <- sapply(seq(q), function(j){
            sapply(seq(q), function(i){
                if(i <= q - 5 & j %in% seq(i, i + 3)){
                    1/(abs(i - j) + 1)
                }else{
                    0
                }
            })
        })
        Sigma_t <- Sigma_t + t(Sigma_t)
        diag(Sigma_t) <- 1
    }else if(type_t == 4){
        Sigma_t <- sapply(seq(q), function(j){
            sapply(seq(q), function(i){
                if(i <= q - 6 & j %in% seq(i, i + 4)){
                    1/(abs(i - j) + 1)
                }else{
                    0
                }
            })
        })
        Sigma_t <- Sigma_t + t(Sigma_t)
        diag(Sigma_t) <- 1
    }else if(type_t == 5){
        Sigma_t <- 0.1^(outer(seq(q), seq(q), FUN = function(x, y) abs(x - y)))
    }
    
    if(type_t %in% c(3, 4)){
        eigen_val <- eigen(Sigma_t)
        de <- abs(min(eigen_val$values)) + 0.05
        Sigma_t <- (Sigma_t + diag(rep(de, q)))/(1 + de)
    }
    
    Sigma_t
}

InfoGenerate_Space <- function(p, type_s){
    ## Space Sigma
    if(type_s == 1){
        Omega <- huge.generator(d = p, graph = "band", g = 3, verbose = FALSE)$omega
    }else if(type_s == 2){
        Omega <- huge.generator(d = p, verbose = FALSE)$omega
    }else if(type_s == 3){
        Omega <- huge.generator(d = p, prob = 0.01, verbose = FALSE)$omega
    }else if(type_s == 4){
        Omega <- huge.generator(d = p, graph = "hub", g = 20, verbose = FALSE)$omega
    }else if(type_s == 5){
        Omega <- createS(n = 1, p = p, topology="small-world", precision = TRUE, banded.n = 5)
        Omega[abs(Omega) <= 0.001] <- 0
        eigen_val <- eigen(Omega)
        de <- abs(min(eigen_val$values)) + 0.05
        Omega_s <- Omega + diag(rep(de, p))
    }
    
    if(type_s %in% seq(4)){
        Omega[abs(Omega) <= 0.001] <- 0
        Omega_s <- Omega
    }else if(type_s == 0){
        Omega_s <- NULL
    }
    
    if(!is.null(Omega_s)){Omega_s <- (Omega_s + t(Omega_s))/2}
    Omega_s
}

RandOnSpaceInfo <- function(Omega, m = 0.5){
    Omega_temp <- Omega
    Omega_random <- Omega
    
    if(m > 0){
        Omega_temp[row(Omega_temp) >= col(Omega_temp)] <- 0
        nonzero <- which(Omega_temp != 0)
        nn <- round(ifelse(m < 1, length(nonzero) * m, round(m)))
        index_random <- sample(nonzero, nn)
        
        Omega_random[index_random] <- 0
        Omega_random[row(Omega_random) > col(Omega_random)] <- t(Omega_random)[row(Omega_random) > col(Omega_random)]
    }
    eigen_val <- eigen(Omega_random)
    de <- abs(min(eigen_val$values)) + 0.05
    
    Omega1 <- Omega + diag(rep(de, p))
    Omega2 <- Omega_random + diag(rep(de, p))
    
    list(Omega1 = Omega1, Omega2 = Omega2)
}


DataGenerate <- function(n, Sigma_t_list, Omega_s_list, corr_type_s = 1, corr_type_t = 0, rho_input, mod_t = 1, mod_s = 1, df = "infty"){
    Sigma_t1 <- Sigma_t_list[[1]]
    Sigma_t2 <- Sigma_t_list[[2]]
    
    Sigma_s1 <- solve(Omega_s_list[[1]])
    Sigma_s2 <- solve(Omega_s_list[[2]])
    
    p <- dim(Sigma_s1)[1]
    q <- dim(Sigma_t1)[1]
    
    ## Spatial cross covariance
    Sigma_s12 <- SigmaCross_Space(Sigma_s1, Sigma_s2, corr_type_s, mod_s)
    
    ## Temporal cross covariance
    Sigma_t12 <- SigmaCross_Time(Sigma_t1, Sigma_t2, corr_type_t, mod_t) # this is inconsistent with the paper. It in fact produces P_{T_{1,2}} in the paper instead Sigma_{T_{1,2}}. But it affects little to the results.
    
    
    #### data generation
    ## spatial part
    Sigma_s1_sqrt <- mysqrtm(Sigma_s1)
    Sigma_s2_sqrt <- mysqrtm(Sigma_s2)
    tilde_S12 <- solve(Sigma_s1_sqrt, Sigma_s12) %*% solve(Sigma_s2_sqrt)
    tilde_S2112 <- t(tilde_S12) %*% tilde_S12
    eigen_S2112 <- eigen(tilde_S2112, symmetric = TRUE)
    Qs <- eigen_S2112$vectors
    lam_s <- eigen_S2112$values
    
    ## temporal part
    Sigma_t1_sqrt <- mysqrtm(Sigma_t1)
    Sigma_t2_sqrt <- mysqrtm(Sigma_t2)
    tilde_T12 <- solve(Sigma_t1_sqrt, Sigma_t12) %*% solve(Sigma_t2_sqrt)
    tilde_T2112 <- t(tilde_T12) %*% tilde_T12
    eigen_T2112 <- eigen(tilde_T2112, symmetric = TRUE)
    Qt <- eigen_T2112$vectors
    lam_t <- eigen_T2112$values
    
    rho <- sign(rho_input) * min(abs(rho_input), 0.95/sqrt((max(abs(lam_t)) * max(abs(lam_s)))))
    
    # step 1:
    if(df == "infty"){
        Z_temp <- rnorm(n = n * 2 * p * q) 
    }else{
        Z_temp <- rt(n = n * 2 * p * q, df = df)
    }
    Z_temp <- Z_temp * sqrt(c(rep(1, p * q), 1 - rho^2 * c(lam_s %o% lam_t))) %>% matrix(nrow = 2 * p * q, ncol = n) %>% t
    
    Z <- sapply(seq(n), function(k){
        Z1 <- matrix(Z_temp[k, seq(p * q)], p, q)
        Z2 <- matrix(Z_temp[k, -seq(p * q)], p, q)
        # step 2
        Z2 <- Qs %*% Z2 %*% t(Qt)
        # step 3
        Z2 <- rho * t(tilde_S12) %*% Z1 %*% tilde_T12 + Z2
        # step 4
        Z1 <- Sigma_s1_sqrt %*% Z1 %*% Sigma_t1_sqrt
        Z2 <- Sigma_s2_sqrt %*% Z2 %*% Sigma_t2_sqrt
        
        return(rbind(Z1, Z2))
    }, simplify = "array")
    
    list(X1 = Z[seq(p),,], X2 = Z[-seq(p),,], rho = rho)
}

DataGenerate2 <- function(n, Sigma_t_list, Omega_s_list, corr_type_s = 1, corr_type_t = 0, rho_input, mod_t = 1, mod_s = 1, distort_prop = 0., distort_level = 1.){
    ######
    # A slower but more general generator for arbitrary off-diagonal covariance matrix
    ######
    
    Sigma_t1 <- Sigma_t_list[[1]]
    Sigma_t2 <- Sigma_t_list[[2]]
    
    Sigma_s1 <- solve(Omega_s_list[[1]])
    Sigma_s2 <- solve(Omega_s_list[[2]])
    
    p <- dim(Sigma_s1)[1]
    q <- dim(Sigma_t1)[1]
    
    ## Spatial cross covariance
    Sigma_s12 <- SigmaCross_Space(Sigma_s1, Sigma_s2, corr_type_s, mod_s)
    
    ## Temporal cross covariance
    Sigma_t12 <- SigmaCross_Time(Sigma_t1, Sigma_t2, corr_type_t, mod_t) # this is inconsistent with the paper. It in fact produces P_{T_{1,2}} in the paper instead Sigma_{T_{1,2}}. But it affects little to the results.
    
    
    Sigma_s1_sqrt <- mysqrtm(Sigma_s1)
    Sigma_s2_sqrt <- mysqrtm(Sigma_s2)
    Sigma_t1_sqrt <- mysqrtm(Sigma_t1)
    Sigma_t2_sqrt <- mysqrtm(Sigma_t2)
    
    
    Sigma1_sqrt <- Sigma_t1_sqrt %x% Sigma_s1_sqrt
    Sigma2_sqrt <- Sigma_t2_sqrt %x% Sigma_s2_sqrt
    Sigma1 <- Sigma_t1 %x% Sigma_s1
    Sigma2 <- Sigma_t2 %x% Sigma_s2
    Sigma12 <- Sigma_t12 %x% Sigma_s12
    
    if(distort_prop > 0 && distort_level > 0){
        sample_idx <- sample(x = seq(length(Sigma12)), size = round(length(Sigma12) * distort_prop), replace = FALSE)
        sd12 <- sd(c(Sigma12))
        Sigma12[sample_idx] <- Sigma12[sample_idx] + rnorm(length(sample_idx), mean = 0, sd = sd12 * distort_level)
    }
    
    shur_comp12_half <- (solve(Sigma_t1_sqrt) %x% solve(Sigma_s1_sqrt)) %*% Sigma12 %*% (solve(Sigma_t2_sqrt) %x% solve(Sigma_s2_sqrt))
    max_eig_shur_comp12_half <- RSpectra::eigs_sym(t(shur_comp12_half) %*% shur_comp12_half, 1)
    rho <- sign(rho_input) * min(abs(rho_input), 0.95/max_eig_shur_comp12_half$values[1])
    
    
    shur_comp12 <- Sigma2 - rho^2 * t(Sigma12) %*% (solve(Sigma_t1) %x% solve(Sigma_s1)) %*% Sigma12
    shur_comp12_sqrt <- RSpectra::eigs_sym(shur_comp12, dim(shur_comp12)[1])
    
    
    list(Omega1 = Omega_s_list[[1]], 
         Omega2 = Omega_s_list[[2]],
         Sigma_s1_sqrt = Sigma_s1_sqrt,
         Sigma_t1_sqrt = Sigma_t1_sqrt,
         Sigma12 = Sigma12,
         shur_comp12_sqrt = shur_comp12_sqrt,
         rho = rho)
}


SigmaCross_Space<- function(Sigma1, Sigma2, corr_type = 1, mod = 1){
    d <- nrow(Sigma1)
    if(corr_type == 1){
        Sigma_tmp <- Sigma1
    }else if(corr_type > 1){
        Sigma_tmp <- Sigma1
        diag(Sigma_tmp)[(seq(d) %% mod) %in% seq(1,5,by=2)] <- -diag(Sigma_tmp)[(seq(d) %% mod) %in% seq(1,5,by=2)]
        if (corr_type == 2){
            Sigma_tmp <- diag(diag(Sigma_tmp))
        }else if(corr_type == 3){
            Sigma_tmp <- Sigma_tmp * (-1)^(row(Sigma_tmp) + col(Sigma_tmp))
        }
    }
    
    return(Sigma_tmp)
}

SigmaCross_Time <- function(Sigma1, Sigma2, corr_type = 0, mod = 1){
    d <- nrow(Sigma1)
    Sigma_tmp <- InfoGenerate_Time(d, type_t = corr_type)
    
    
    if(mod != 1){
        Sigma_tmp <- TimeSignChange(Sigma_tmp, mod)
    }
    Sigma_tmp <- mysqrtm(Sigma1) %*% Sigma_tmp %*% mysqrtm(Sigma2)
    
    return(Sigma_tmp)
}

TimeSignChange <- function(Sigma, mod){
    d <- nrow(Sigma)
    Sigma <- Sigma * (-1)^(row(Sigma) + col(Sigma))
    if(mod > 0){
        diag(Sigma)[(seq(d) %% mod) %in% seq(1,5,by=2)] <- -diag(Sigma)[(seq(d) %% mod)%in% seq(1,5,by=2)]
    }else if(mod == 0){
        diag(Sigma) <- diag(Sigma) * (-1)^(seq(d))
    }else{
        mod <- abs(mod)
        diag(Sigma)[(seq(d) %% mod) %in% c(seq(1,mod - 1,by=2), 0)] <- -diag(Sigma)[(seq(d) %% mod) %in% c(seq(1,mod - 1,by=2), 0)]
    }
    Sigma
}

