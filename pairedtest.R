TimeRemoval <- function(X, method, ...){
  # Returns an array with time information removed
  #  Args:
  #   X: a 3d-array with three dimension: p, the space dimension; 
  #      q, the time dimension; n, the number of subjects
  #   method: the method to remove time information. To be specifice,
  #           "data-driven", "banded-estimation", "oracle", "band_xia"
  #  Dependent Packages:
  #   expm 0.999-0
  #   dplyr
  #   PDSCE
  optional_par <- list(...)
  dimen <- dim(X)
  p <- dimen[1]
  q <- dimen[2]
  n <- dimen[3]
  ## root of the time covariance
  if(method == 1){
    # "data_driven"
    Sigma_t <- EmpSigma(X, NULL, 2)
  }else if(method == 2){
    #  "banded_estimation"
    X_stack_p <- Stack(X, 1)
    Sigma_t <- EmpSigma(X, NULL, 2)
    lambda <- seq(10)
    
    Sigma_t <- Band_Est(X_stack_p, Sigma_t, lambda)
  }else if(method == 3){
    # "oracle"
    Sigma_t <- optional_par$oSigma_t
  }else if(method == 4){
    # "band_Xia"
    X_stack_p <- Stack(X, 1)
    Sigma_t <- EmpSigma(X, NULL, 2)
    band_par <- q
    lambda <- seq(10)
    
    ## root of the time covariance and remove the time information inside
    Y <- lapply(seq(p), function(k){
      CSigma_t <- t(X[k,,]) %*% X[k,,]/n
      Sigma_t <- Band_Est(t(X[k,,]), Sigma_t, lambda)
      de <- abs(min(eigen(Sigma_t)$values))+0.05
      Sigma_t <- (Sigma_t + de * diag(rep(1, q)))/(1 + de)
      
      rSigma_t <- mysqrtm(Sigma_t)
      solve(t(rSigma_t), X[k,,]) %>% t %>% c
    }) %>% 
      do.call(what = c) %>% 
      array(dim = c(n, q, p)) %>% 
      aperm(perm = c(3,2,1))
  }
  
  ## remove the time information
  if(method != 4){
    rSigma_t <- mysqrtm(Sigma_t)
    Y <- lapply(seq(n), function(k){
      solve(t(rSigma_t), t(X[,,k])) %>% t %>% c
    }) %>% do.call(what = c) %>% array(dim = c(p, q, n))
  }
  
  return(Y)
}

mysqrtm <- function(X, type = 1){
  if(type == 1){
    rlt <- Schur(X)
    Tmat <- rlt$T
    Tmat[which(abs(Tmat) < 1e-10, arr.ind = TRUE)] <- 0
    rlt$Q %*% sqrt(Tmat) %*% t(rlt$Q)
  }else if(type == 2){
    sqrtm(X)
  }
}

EmpSigma <- function(X1, X2 = NULL, d){
  # Returns the empirical covariance given the dimension
  #  Args:
  #   X: a 3d-array with three dimension: p, the space dimension; 
  #      q, the time dimension; n, the number of subjects
  if(is.null(X2)){
    X2 = X1
  }
  
  if(!all(dim(X1) == dim(X2))){
    stop("The dimension of X1 does not match that of X2.")
  }else{
    dimen <- dim(X1)
    p <- dimen[1]
    q <- dimen[2]
    n <- dimen[3]
    
    meanX1 <- matrix(c(X1), nrow = p * q, ncol = n) %>% rowMeans() %>% matrix(nrow = p, ncol = q)
    meanX2 <- matrix(c(X2), nrow = p * q, ncol = n) %>% rowMeans() %>% matrix(nrow = p, ncol = q)
    
    Sigma <- lapply(seq(n), function(k){
      if(d == 1){
        c((X1[,,k] - meanX1) %*% t((X2[,,k] - meanX2)))
      }else if(d == 2){
        c(t((X1[,,k] - meanX1)) %*% (X2[,,k] - meanX2))
      }
    }) %>% 
      do.call(what = cbind) %>% 
      rowMeans()
    
    if(d == 1){
      matrix(Sigma, nrow = p, ncol = p)/q 
    }else if(d == 2){
      matrix(Sigma, nrow = q, ncol = q)/p
    }
  }
}

Band_Est <- function(X, CSigma, lambda){
  # Returns the covaraince banded-estimated
  #  Args:
  #   X: a stacking matrix with dimension np * q (collasping the first dimension)
  #   CSigma: a comparison matrix to determine the banding paramter
  #   lambda: candidates of banding parameters
  dif <- sapply(lambda, function(k){
    (band.chol(X, k, centered = TRUE, method = "safe") - CSigma) %>% 
      abs %>%
      sum
  })
  
  thr <- max(lambda[which.min(dif)])
  band.chol(X, thr, centered = TRUE, method = "safe")
}

Stack <- function(X, d){
  # Stacking an 3-d array into a matrix given the collapsing dimension
  #  Args:
  #   X: a 3d-array with three dimension: p, the space dimension; 
  #      q, the time dimension; n, the number of subjects
  dimen <- dim(X)
  p <- dimen[1]
  q <- dimen[2]
  n <- dimen[3]
  if(d == 1){
    lapply(seq(p), function(k) t(X[k,,])) %>% do.call(what = rbind)
  }else if(d == 2){
    lapply(seq(q), function(k) t(X[,k,])) %>% do.call(what = rbind)
  }else if(d == 3){
    lapply(seq(p), function(k) X[k,,]) %>% do.call(what = rbind)
  }
}



InvReg <- function(Y, lambda){
  stackY <- Stack(Y, 2)
  sigmaY <- diag(cov(stackY))
  
  p <- dim(stackY)[2]
  n <- dim(stackY)[1]
  ## fit lasso for regression
  reg_rlt <- lapply(seq(p), function(i){ 
    y <- stackY[, i, drop = FALSE]
    resy <- stackY[, -i, drop = FALSE]
    glm_rlt<- glmnet(resy, y, family = "gaussian", intercept = TRUE, # need to think twice about the glmnet part
                     lambda = lambda/20 * sqrt(sigmaY[i] * log(p)/n))
    coef_rlt <- coef(glm_rlt, mode = "lambda") 
    beta <- matrix(coef_rlt[-1], ncol = 1)
    # xi
    xi <- y - mean(y) - scale(resy, scale = FALSE) %*% beta # need to think more
    
    c(beta, xi) 
  }) %>% do.call(what = cbind)
  Beta <- reg_rlt[seq(p-1), ]
  Xi <- reg_rlt[-seq(p-1), ]
  
  list(Beta = Beta, Xi = Xi)  
}

TestStat <- function(Y1, Y2, lambda, indpt = NULL, X1, X2, print_set){
  # Y1
  reg_rlt1 <- InvReg(Y1, lambda)
  Beta1 <- reg_rlt1$Beta
  Xi1 <- reg_rlt1$Xi
  n1 <- dim(Xi1)[1]
  p <- dim(Xi1)[2]
  R1 <- t(Xi1) %*% Xi1/n1
  s1 <- colMeans(Xi1^2)
  Rcrct1 <- RCorrect(R1, s1, Beta1)
  T1 <- TStat(Rcrct1)
  theta1 <- ThetaStat(Beta1, Rcrct1, n1)
  
  # Y2
  reg_rlt2 <- InvReg(Y2, lambda)
  Beta2 <- reg_rlt2$Beta
  Xi2 <- reg_rlt2$Xi
  n2 <- dim(Xi2)[1]
  R2 <- t(Xi2) %*% Xi2/n2
  s2 <- colMeans(Xi2^2)
  Rcrct2 <- RCorrect(R2, s2, Beta2)
  T2 <- TStat(Rcrct2)
  theta2 <- ThetaStat(Beta2, Rcrct2, n2)
  
  ## variance
  if(is.null(indpt)){
    #indpt_choices <- c(0, 1, 2, 3)
    indpt_choices <- c(0, 1, 6, 9)
  }else{
    indpt_choices <- indpt
  }
  Wstat <- lapply(seq(length(indpt_choices)), function(indpt_index){
    indpt <- indpt_choices[indpt_index]
    variance <- theta1 + theta2
    
    
    if(indpt != 0){
      if(n1 == n2){
        if(indpt != 10){
          corr4 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = indpt, X1, X2, print_set)
        }else{
          corr4_1 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = 1, X1, X2, print_set)
          corr4_6 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = 6, X1, X2, print_set)
          corr4 <- matrix(pmin(c(corr4_1), c(corr4_6)), nrow = nrow(corr4_1))
        }
        variance <- variance - 2 * corr4 
      }else{
        stop("Not repeated measures!")
      }
    }
    
    ## Ultimate stat
    W <- abs(T1 - T2)/sqrt(variance)
    diag(W) <- 0
    
    # temporary method to deal NA variance: using the variance of indpt = 1 for substitution
    if(sum(is.na(W) > 0)){
      tmp_variance <- theta1 + theta2
      if(indpt != 0){
        if(n1 == n2){
          tmp_corr4 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = 1, X1, X2, print_set)
          tmp_variance <- tmp_variance - 2 * tmp_corr4 
        }else{
          stop("Not repeated measures!")
        }
      }
      tmp_W <- abs(T1 - T2)/sqrt(tmp_variance)
      diag(tmp_W) <- 0
      W[is.na(W)] <- 0
      tmp_W[!is.na(W)] <- 0
      W <- W + tmp_W
    }
    
    list(W = W, variance = variance)
  })
  
  
  if(length(indpt_choices) == 1){
    return(Wstat[[1]])
  }else{
    return(Wstat)
  }
}

RCorrect <- function(R, s, Beta){
  p <- dim(R)[1]
  Rcrct <- lapply(seq(p), function(i){
    sapply(seq(p), function(j){
      if(i <= p - 1 && j > i){
        -(R[i, j] + s[i] * Beta[i, j] + s[j] * Beta[j - 1, i])
      }else{
        0
      }
    })
  }) %>% do.call(what = rbind)
  Rcrct + t(Rcrct) + diag(diag(R))
}

TStat <- function(Rcrct){
  Rcrct/outer(sqrt(diag(Rcrct)), sqrt(diag(Rcrct)))
}

ThetaStat <- function(Beta, Rcrct, n){
  p <- dim(Rcrct)[1]
  theta <- lapply(seq(p), function(i){
    sapply(seq(p), function(j){
      if(i <= p - 1 && j > i){
        (1 + Beta[i, j]^2 * Rcrct[i, i]/Rcrct[j, j])
      }else{
        0
      }
    })
  }) %>% do.call(what = rbind)
  
  (theta + t(theta))/n
}

Corr4 <- function(Z, lambda){
  reg_rlt <- InvReg(Z, lambda)
  Beta <- reg_rlt$Beta
  Xi <- reg_rlt$Xi
  n <- dim(Xi)[1]
  p <- dim(Xi)[2]/2
  R <- t(Xi) %*% Xi/n
  s <- colMeans(Xi^2)
  Rcrct <- RCorrect(R, s, Beta)
  Tz <- TStat(Rcrct)
  T11 <- Tz[seq(p), seq(p)]
  T12 <- Tz[seq(p), -seq(p)]
  T21 <- Tz[-seq(p), seq(p)]
  T22 <- Tz[-seq(p), -seq(p)]
  
  rho12 <- T11
  rho13 <- outer(diag(T12), rep(1, p))
  rho14 <- T12
  rho23 <- T21
  rho24 <- outer(rep(1, p), diag(T21))
  rho34 <- T22
  
  (rho13 * rho24 + rho14 * rho23)/n # it should be rho12 * rho34 + rho13 * rho24 + rho14 * rho23, but
  # rho12 * rho34 is cancelled by T1 * T2
}



VarCorrect <- function(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt, X1, X2, print_set){
  para_index <- print_set$para_index 
  inner_index <- print_set$inner_index
  rho <- print_set$rho
  irep <- print_set$irep
  lambda <- print_set$lambda
  
  p <- dim(Y1)[1]
  q <- dim(Y1)[2]
  n <- dim(Y1)[3]
  n_stack <- n * q
  R1 <- t(Xi1) %*% Xi1/n_stack
  R2 <- t(Xi2) %*% Xi2/n_stack
  CovR1R2 <- outer(sqrt(diag(R1)), sqrt(diag(R2)))
  
  ## empirical covaraince
  # temporal
  if(indpt %in% c(4,5,6)){
    #tmp_S1 <- EmpSigma(Y1, Y1, 1)
    #tmp_S2 <- EmpSigma(Y2, Y2, 1)
    #B_space <- EmpSigma(Y1, Y2, 1)/(outer(sqrt(diag(tmp_S1)), sqrt(diag(tmp_S2))))
    
    tmp_S1 <- EmpSigma(X1, X1, 1)
    tmp_S2 <- EmpSigma(X2, X2, 1)
    B_space <- EmpSigma(X1, X2, 1)/(outer(sqrt(diag(tmp_S1)), sqrt(diag(tmp_S2))))
    #if(indpt %in% c(4, 5, 6) & mean(abs(diag(B_space))) < 0.1){
    #  indpt <- 1
    #}else{
    # A_time <- EmpSigma(sign(diag(B_space)) * Y1, Y2, 2)
    
    
    A_time <- AlterSign(Y1, Y2, 2)
    ## thresholding (be careful)
    A_time[intersect(which(abs(A_time) < sqrt(log(q)/n/p)), which(abs(row(A_time) - col(A_time)) > 3))] <- 0
    #}
    print(paste0("File: ", para_index, "; inner_index: ", inner_index, "; rho: ", round(rho, digits = 6), 
                 "; Process: ", irep, ": lambda = ", lambda, ". Indpt: ", indpt, 
                 "; trace_abs: ", round(mean(abs(diag(B_space))), digits = 6), 
                 "; trace: ", round(mean(diag(B_space)), digits = 6),
                 "; factor:", q * sum(A_time^2)/(sum(abs(diag(A_time)))^2 + 1e-8),
                 "; true_factor:", q * sum(print_set$A_time^2)/(sum(abs(diag(print_set$A_time)))^2 + 1e-8),
                 "."))
  }else if(indpt == 7){
    tmp_S1 <- sqrt(diag(EmpSigma(Y1, Y1, 1)))
    tmp_S2 <- sqrt(diag(EmpSigma(Y2, Y2, 1)))
    CovY1Y2 <- tmp_S1 %o% tmp_S2
    
    shuffle_temp <- ShuffleSigma(Y1, Y2, CovY1Y2 = CovY1Y2, 2)
    max_index_temp <- shuffle_temp$max_index
    signs_temp <- shuffle_temp$signs
    A_time <- EmpSigma(Y1[max_index_temp[,1],,] * signs_temp/tmp_S1[max_index_temp[,1]], Y2/tmp_S2, 2)
    
  }else if(indpt == 8){
    # extension of indpt = 4
    # There is still problem to debug for indpt == 8!
    tmp_S1 <- sqrt(diag(EmpSigma(Y1, Y1, 1)))
    tmp_S2 <- sqrt(diag(EmpSigma(Y2, Y2, 1)))
    
    sign_mat <- cbind(rep(1, p), 
                      sapply(seq(round(p/2)), function(i){
                        tmp_vec <- rep(1, p)
                        tmp_vec[seq(2, 2 * i, 2)] <- -1
                        tmp_vec
                      }),
                      sapply(seq(round(p/2)), function(i){
                        tmp_vec <- rep(1, p)
                        tmp_vec[seq(1, 2 * i - 1, 2)] <- -1
                        tmp_vec
                      })
    )
    
    cutoff <- sapply(seq(ncol(sign_mat)), function(b){
      temp_mat <- EmpSigma(Y1 * sign_mat[, b]/tmp_S1, Y2/tmp_S2, 2)
      mean(abs(diag(temp_mat)))
    })
    
    signs_temp <- sign_mat[,which.max(cutoff)]
    A_time <- EmpSigma(Y1 * signs_temp/tmp_S1, Y2/tmp_S2, 2)
    
  }else if(indpt == 9){
    # semi-oracle case: true A
    A_time <- print_set$A_time
    # if(mean(abs(diag(A_time))) < 0.1){
    #   indpt <- 1
    # }else{
    indpt <- 6
    # }
  }
  
  # stack Y
  if(indpt == 7){
    shuffle_spt <- ShuffleSigma(Y1, Y2, CovY1Y2 = CovY1Y2, 1)
    max_index_spt <- shuffle_spt$max_index
    signs_spt <- shuffle_spt$signs
    A_time <- A_time[max_index_spt, ] * signs_spt
    Y1 <- ((Y1[,max_index_spt[,1],] %>% aperm(perm = c(2, 1, 3))) * signs_spt) %>% aperm(perm = c(2,1,3))
    if(abs(mean(diag(A_time))) < 0.1){
      indpt <- 1
    }
    print(paste0("File: ", para_index, "; inner_index: ", inner_index, "; rho:", round(rho, digits = 6), "; Process: ", irep, ": lambda = ", lambda, ". Indpt: ", indpt, "; trace: ", round(abs(mean(diag(A_time))), digits = 6), "."))
    
  }else if(indpt == 8){
    signs_spt <- sign(diag(A_time))
    A_time <- A_time * signs_spt
    Y1 <- (aperm(Y1, perm = c(2, 1, 3)) * signs_spt) %>% aperm(perm = c(2,1,3))
    if(abs(mean(diag(A_time))) < 0.1){
      indpt <- 1
    }
    print(paste0("File: ", para_index, "; inner_index: ", inner_index, "; rho:", round(rho, digits = 6), "; Process: ", irep, ": lambda = ", lambda, ". Indpt: ", indpt, "; trace: ", round(abs(mean(diag(A_time))), digits = 6), "."))
  }
  stackY1 <- Stack(Y1, 2)
  stackY2 <- Stack(Y2, 2)
  
  # spatial 
  if(indpt %in% c(seq(4), 7, 8)){
    emp_Sigma <- cov(stackY1, stackY2)
  }else if(indpt == 5){
    stackY1_pq <- Stack(Y1, 3)
    stackY2_pq <- Stack(Y2, 3) 
    
    temp <- stackY1_pq %*% t(stackY2_pq)/n
    temp <- lapply(seq(p), function(j){
      temp[,seq((j - 1) * q + 1, j * q)]
    }) %>% do.call(what = rbind) - rep(c(colMeans(stackY1) %o% colMeans(stackY2)), each = q)
    
    emp_Sigma <- matrix(c(t(temp)) * c(sign(t(A_time))), q^2) %>% colMeans() %>% matrix(nrow = p)
  }else if(indpt == 6){
    emp_Sigma <- cov(rep(sign(diag(A_time)), each = n) * stackY1, stackY2)
  }
  
  
  # beta estimation
  if(indpt != 3){
    if(indpt == 2){
      Beta1 <- BetaCorrect(Rcrct1)
      Beta2 <- BetaCorrect(Rcrct2)
    }
    
    temp1 <- rbind(rep(-1, p-1), matrix(Beta1, p))
    temp1 <- matrix(c(temp1, -1), p)
    
    temp2 <- rbind(rep(-1, p-1), matrix(Beta2, p))
    temp2 <- matrix(c(temp2, -1), p)
  }
  
  # empirical correlation 
  if(indpt != 3){
    corr_mat <- temp1 %*% emp_Sigma %*% temp2/CovR1R2
  }else if(indpt == 3){
    R12 <- t(Xi1) %*% Xi2/n_stack
    corr_mat <- R12/CovR1R2
  }
  
  # corr4 <- lapply(seq(p), function(j){
  #   sapply(seq(p), function(i){
  #     corr_mat[i, i] * corr_mat[j, j] + corr_mat[i, j] * corr_mat[j, i]
  #   })
  # }) %>% do.call(what = cbind)
  corr4 <- (diag(corr_mat) %o% diag(corr_mat)) + corr_mat * t(corr_mat)
  
  # output
  if(indpt <= 3){
    corr4/n_stack
  }else if(indpt %in% c(4, 7, 8)){
    corr4 * sum(A_time^2)/(sum(diag(A_time))^2 + 1e-8)/n
  }else if(indpt == 5){
    corr4 * q^2 * sum(A_time^2)/(sum(abs(A_time))^2 + 1e-8)/n
  }else if(indpt == 6){
    corr4 * sum(A_time^2)/(sum(abs(diag(A_time)))^2 + 1e-8)/n
  }
}

AlterSign <- function(Y1, Y2, num_alter){
  p <- dim(Y1)[1]
  q <- dim(Y1)[2]
  
  diag_sign_A <- rep(1, q)
  for(i in seq(num_alter)){
    tmp_S1 <- EmpSigma(rep(diag_sign_A, each = p) * Y1, rep(diag_sign_A, each = p) * Y1, 1)
    tmp_S2 <- EmpSigma(Y2, Y2, 1)
    B_space <- EmpSigma(rep(diag_sign_A, each = p) * Y1, Y2, 1)/(outer(sqrt(diag(tmp_S1)), sqrt(diag(tmp_S2))))
    
    A_time <- EmpSigma(sign(diag(B_space)) * Y1, Y2, 2)
    ## thresholding (be careful)
    A_time[intersect(which(abs(A_time) < sqrt(log(q)/n/p)), which(abs(row(A_time) - col(A_time)) > 3))] <- 0
    
    diag_sign_A <- sign(diag(A_time))
  }
  A_time
}

BetaCorrect <- function(Rcrct){
  p <- ncol(Rcrct)
  temp <- -diag(1/diag(Rcrct)) %*% Rcrct 
  lapply(seq(p), function(j) temp[-j, j]) %>% do.call(what = cbind)
}


ShuffleSigma <- function(Y1, Y2, CovY1Y2, d = 1){
  is_zero <- TRUE
  n <- dim(Y1)[3]
  if(d == 2){
    Y1 <- aperm(Y1, perm = c(2,1,3))
    Y2 <- aperm(Y2, perm = c(2,1,3))
  }
  dim1 <- dim(Y1)[1]
  dim2 <- dim(Y2)[2]
  
  stackY1 <- Stack(Y1, 2)
  stackY2 <- Stack(Y2, 2)
  
  CovY <- cov(t(sapply(seq(n), function(k) c(Y1[,,k]))), t(sapply(seq(n), function(k) c(Y2[,,k]))))
  
  # Frobenius norm of each cell
  # tmp_mat <- lapply(seq(dim2), function(l){
  #   CovY[, seq((l - 1) * dim1 + 1, l * dim1)]
  # }) %>% do.call(what = rbind) 
  # if(d == 1){
  #   normYY_Frob <- c(t(tmp_mat))/c(t(CovR1R2)) 
  # }else if(d == 2){
  #   normYY_Frob <- tmp_mat/rep(c(CovR1R2), each = dim1) %>% t %>% c
  # }
  # normYY_Frob <- matrix(normYY_Frob^2, nrow = dim1^2) %>% colMeans() %>% matrix(nrow = dim2)
  normYY_Frob <- sapply(seq(dim2), function(l2){
    sapply(seq(dim2), function(l1){
      if(d == 1){
        norm_factor <- CovY1Y2
      }else if(d == 2){
        norm_factor <- CovY1Y2[l1, l2]
      }
      (CovY[seq((l1 - 1) * dim1 + 1, l1 * dim1), seq((l2 - 1) * dim1 + 1, l2 * dim1)]/norm_factor)^2 %>% mean
    })
  })
  
  if(max(sqrt(normYY_Frob)) >= 0.1){
    is_zero <- FALSE
  }
  # max index
  max_index <- MaxIndex(normYY_Frob)
  
  # corresponding signs to max indexes
  max_block <- sapply(seq(dim2), function(l){
    if(d == 1){
      norm_factor <- CovY1Y2
    }else if(d == 2){
      norm_factor <- CovY1Y2[max_index[l, 1], max_index[l, 2]]
    }
    CovY[seq((max_index[l, 1] - 1) * dim1 + 1, max_index[l, 1] * dim1), 
         seq((max_index[l, 2] - 1) * dim1 + 1, max_index[l, 2] * dim1)]/norm_factor %>% c
  }) 
  sign_index <- rowMeans(abs(max_block)) %>% which.max
  signs <- sign(max_block[sign_index, ])
  
  # order max index by its second coordinate
  signs <- signs[order(max_index[, 2])]
  max_index <- max_index[order(max_index[,2]), ]
  
  list(is_zero = is_zero, max_index = max_index, signs = signs)
}


MaxIndex <- function(mat){
  if(length(mat) == 1){
    matrix(c(1,1), nrow = 1)
  }else{
    curr_index <- which(mat == max(mat), arr.ind = TRUE)[1,]
    next_index <- MaxIndex(mat[-curr_index[1], -curr_index[2]])
    rbind(curr_index, lapply(seq(nrow(next_index)), function(tmp_index) 
      c(ifelse(next_index[tmp_index, 1] < curr_index[1], next_index[tmp_index, 1], next_index[tmp_index, 1] + 1),
        ifelse(next_index[tmp_index, 2] < curr_index[2], next_index[tmp_index, 2], next_index[tmp_index, 2] + 1))) %>% do.call(what = rbind))
  }
}


SeqTest <- function(Xlist, alpha, indpt, method_rm = 1, Omegalist = NULL){
  nlist <- length(Xlist)
  p <- dim(Xlist[[1]])[1]
  
  z <- seq(0 : ((floor(sqrt(4 * log(p))) + 2) * 100))/100
  len_z <- length(z)
  
  cri = -2*log (-(8 * pi)^(1/2)*log(1 - alpha));
  cri = 4*log (p)-log (log (p)) + cri;
  
  
  Ylist <- lapply(seq(nlist), function(l){
    TimeRemoval(Xlist[[l]], method_rm)
  })
  
  ## truth
  if(!is.null(Omegalist)){
    Wlist <- lapply(seq(nlist - 1), function(l){
      Wtrue <- Omegalist[[l]]/outer(diag(Omegalist[[l]]), diag(Omegalist[[l]]), "/") - 
        Omegalist[[l+1]]/outer(diag(Omegalist[[l+1]]), diag(Omegalist[[l+1]]), "/")
      Wtrue[abs(Wtrue) <= 0.001] <- 0
      Wtrue
    })
  }
  
  ## global testing
  print_set <- list(para_index = 0, inner_index = 0, rho = 0, irep = 0, lambda  = 40, A_time = NULL)
  MBstat <- sapply(seq(nlist - 1), function(l){
    Wstat <- TestStat(Ylist[[l]], Ylist[[l + 1]], 40, indpt = indpt, Xlist[[l]], Xlist[[l + 1]], print_set)$W
    MB <- max(max(Wstat^2))
    as.numeric(MB > cri)
  })
  
  ## Multiple testing
  mul_rlt <- lapply(seq(nlist - 1), function(l){
    if(MBstat == 1){
      temp_rlt <- lapply(seq(3, 20), function(lambda){
        print_set <- list(para_index = 0, inner_index = 0, rho = 0, irep = 0, lambda  = lambda, A_time = NULL)
        
        # Test Statistics
        Wstat <- TestStat(Ylist[[l]], Ylist[[l + 1]], lambda, indpt = indpt, Xlist[[l]], Xlist[[l + 1]], print_set)$W
        
        # rejection threshold
        rr <- sapply(seq(len_z), function(l){
          p * (p - 1) * (1 - pnorm(z[l]))/max(sum(abs(Wstat) >= z[l])/2, 1)
        })
        
        if(rr[len_z] <= alpha){
          t_threshold <- min(which((rr <= alpha)))
        }else{
          t_threshold <- len_z
        }
        
        # different index
        diff_index <- which(abs(Wstat) >= z[t_threshold], arr.ind = TRUE)
        diff_index <- cbind(diff_index, 2 * (1 - pnorm(abs(Wstat[diff_index]))))
        # select lambda for fdr
        sssp <- sapply(seq(10), function(l){
          temp_alpha <- l/floor(10/(1 - pnorm(sqrt(log(p)))))
          (sum(abs(Wstat) >= qnorm(1 - temp_alpha))/2/(temp_alpha * (p * (p-1))) - 1)^2
        }) %>% sum
        if(!is.null(Omegalist)){
          # empirical size
          ffp <- sum(abs(Wstat[which(Wtrue == 0, arr.ind = TRUE)]) >= z[t_threshold])/max(sum(abs(Wstat) >= z[t_threshold]), 1)
          # empirical power
          fpp <- (max(sum(abs(Wstat) >= z[t_threshold]), 1) - sum(abs(Wstat[which(Wtrue == 0, arr.ind = TRUE)]) >= z[t_threshold]))/((sum(Wtrue != 0)))
        }else{
          ffp <- NULL
          fpp <- NULL
        }
        list(Wstat = Wstat, threshold = z[t_threshold], diff_index = diff_index, sssp = sssp, ffp = ffp, fpp = fpp)
      })
      best_index <- which.min(sapply(seq(length(temp_rlt)), function(lambda){
        temp_rlt[[lambda]]$sssp
      }))[1]
      temp_rlt <- temp_rlt[[best_index]]
      temp_rlt$lambda <- best_index
      temp_rlt
    }else{
      list()
    }
  })
  mul_rlt$MBstat <- MBstat
  mul_rlt
}


