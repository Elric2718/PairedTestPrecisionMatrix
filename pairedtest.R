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

TestStat <- function(Y1, Y2, lambda, indpt = NULL, X1, X2){
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
    indpt_choices <- c(0, 1)
  }else{
    indpt_choices <- indpt
  }
  Wstat <- lapply(seq(length(indpt_choices)), function(indpt_index){
    indpt <- indpt_choices[indpt_index]
    variance <- theta1 + theta2
    
    
    if(indpt != 0){
      if(n1 == n2){
        corr4 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = indpt, X1, X2)
        variance <- variance - 2 * corr4 
      }else{
        stop("Not repeated measures!")
      }
    }
    
    ## Ultimate stat
    W <- try(abs(T1 - T2)/sqrt(variance), silent = TRUE)
    diag(W) <- 0
    
    # temporary method to deal NA variance: using the variance of indpt = 1 for substitution
    if(sum(is.na(W) > 0)){
      tmp_variance <- theta1 + theta2
      if(indpt != 0){
        if(n1 == n2){
          tmp_corr4 <- VarCorrect(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt = 1, X1, X2)
          tmp_variance <- tmp_variance - 2 * tmp_corr4 
        }else{
          stop("Not repeated measures!")
        }
      }
      tmp_W <- try({abs(T1 - T2)/sqrt(tmp_variance)}, silent = TRUE)
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
  
  (rho13 * rho24 + rho14 * rho23)/n 
}



VarCorrect <- function(Y1, Y2, Beta1, Beta2, Xi1, Xi2, Rcrct1, Rcrct2, indpt, X1, X2){
  p <- dim(Y1)[1]
  q <- dim(Y1)[2]
  n <- dim(Y1)[3]
  n_stack <- n * q
  R1 <- t(Xi1) %*% Xi1/n_stack
  R2 <- t(Xi2) %*% Xi2/n_stack
  CovR1R2 <- outer(sqrt(diag(R1)), sqrt(diag(R2)))
  
  if(indpt == 0){
    return(matrix(0, p, p))
  }else if(indpt == 1){
    ## empirical covaraince
    # temporal
    tmp_S1 <- EmpSigma(X1, X1, 1)
    tmp_S2 <- EmpSigma(X2, X2, 1)
    B_space <- EmpSigma(X1, X2, 1)/(outer(sqrt(diag(tmp_S1)), sqrt(diag(tmp_S2))))
    
    A_time <- AlterSign(Y1, Y2, 2)
    ## thresholding (be careful)
    A_time[intersect(which(abs(A_time) < sqrt(log(q)/n/p)), which(abs(row(A_time) - col(A_time)) > 3))] <- 0
    
    stackY1 <- Stack(Y1, 2)
    stackY2 <- Stack(Y2, 2)
    emp_Sigma <- cov(rep(sign(diag(A_time)), each = n) * stackY1, stackY2)
    
    
    # beta estimation
    temp1 <- rbind(rep(-1, p-1), matrix(Beta1, p))
    temp1 <- matrix(c(temp1, -1), p)
    
    temp2 <- rbind(rep(-1, p-1), matrix(Beta2, p))
    temp2 <- matrix(c(temp2, -1), p)
    
    
    # empirical correlation 
    corr_mat <- temp1 %*% emp_Sigma %*% temp2/CovR1R2
    corr4 <- (diag(corr_mat) %o% diag(corr_mat)) + corr_mat * t(corr_mat)
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

