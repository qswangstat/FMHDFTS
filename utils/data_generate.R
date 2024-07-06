## This file contains data-generating functions used in all simulations
## I should have organized the code in a more concise form, such as defining a generic generator class, etc.,
## But at that time, I was not aware of the superiority of doing so
## Hence I resorted to the current naive approach presented here. Please bear with this limitation. 

## Paramter intro
# n: sample size
# p: data dimension
# r: number of factors
# s: dimension of basis expansion for factor process
# A: loading matrix
# varMat: coefficient matrix of VAR(1) model for factor basis coefficients
# bmat: basis funcion values
# delta: control factor signal strength
# q: control variance decaying of innovation terms in VAR(1) model for factor basis coefficients
# sigma: control variance amplititude of innovation terms in VAR(1) model for factor basis coefficients
# errdim: dimension of basis expansion for idiosyncratic process
# burn: to gurantee stationarity of time series, I drop the first "burn" generated terms
# err_ar_coef: in scenario 4-5, I consider serial dependence for idiosyncratic terms controlled by this param
# cross_depend: in scenario 5, I consider cross-correlation among p-dim idiosyncratic terms controlled by this param
data_gen_1 = function(n, p, r, s, A, varMat, bmat, delta = 1, q = 0.75,
                    sigma = 1, errdim = 20, err_ar_coef = 0.5, burn = 1000){
  ### factor process x_t(u) is r-dim
  ### each component is expanded under s basis
  ### first generate expansion coefficient
  ### coefficients are stored in a (r*s)*(2*n+burn)-dim matrix
  ### ((i-1)*r+1)~(i*r) rows correspond to the i-th component
  ### notice the coefficients of certain basis function are independent of those of other functions  
  coefMat = matrix(0, nrow = r*s, ncol = n + burn)
  for (i in 1:s) {
    coefMat_i = matrix(0, nrow = r, ncol = n + burn)
    coefMat_i[, 1] = runif(r, -1, 1)
    for(j in 2:ncol(coefMat_i)){
      coefMat_i[, j] = varMat %*% coefMat_i[, j-1] + rnorm(r, 0, sigma*i^(-q))
    }
    coefMat[((i-1)*r+1):(i*r), ] = coefMat_i
  }
  coefMat = coefMat[, -(1:burn)]
  
  data = array(0, dim = c(p, dim(bmat)[2], n))
  facto = array(0, dim = c(r, dim(bmat)[2], n))
  for (i in 1:(n)) {
    coefi = matrix(coefMat[, i], nrow = r, ncol = s)
    xn = coefi %*% bmat[1:s, ]
    facto[ , , i] = xn
    ### construct error process
    epscoef = matrix(rnorm(p*errdim, 0, 1), nrow = p, ncol = errdim)
    epscoef = t(t(epscoef) * 1/2^((1:errdim)+1))
    epsn = epscoef %*% bmat[1:errdim, ]
    data[, , i] = A %*% xn * delta + epsn
  }
  return(data)
}
# independent factor
data_gen_2 = function(n, p, r, s, A, varMat, bmat, delta = 1, q = 0.75,
                    sigma = 1, errdim = 20, err_ar_coef = 0.5, burn = 1000){
  coefMat = matrix(0, nrow = r*s, ncol = 2*n + burn)
  for (i in 1:s) {
    coefMat_i = matrix(0, nrow = r, ncol = 2*n + burn)
    coefMat_i[, 1] = runif(r, -1, 1)
    for(j in 2:ncol(coefMat_i)){
      coefMat_i[, j] = varMat %*% coefMat_i[, j-1] + rnorm(r, 0, sigma*i^(-q))
    }
    coefMat[((i-1)*r+1):(i*r), ] = coefMat_i
  }
  coefMat = coefMat[, -(1:burn)]
  
  data = array(0, dim = c(p, dim(bmat)[2], n))
  facto = array(0, dim = c(r, dim(bmat)[2], n))
  for (i in 1:(n)) {
    coefi = matrix(coefMat[, i], nrow = r, ncol = s)
    xn = coefi %*% bmat[1:s, ]
    facto[ , , i] = xn
    epscoef = matrix(rnorm(p*errdim, 0, 1), nrow = p, ncol = errdim)
    epscoef = t(t(epscoef) * 1/2^((1:errdim)+1))
    epsn = epscoef %*% bmat[1:errdim, ]
    
    epsn2 = t(as.vector(rnorm(4, 0, 1))) %*% bmat[1:4, ] * delta
    
    data[, , i] = t(t(A %*% xn + epsn)+as.vector(epsn2))
  }
  return(data)
}
# heterogeneous idiosyncratic terms
data_gen_3 = function(n, p, r, s, A, varMat, bmat, delta = 1, q = 0.75,
                    sigma = 1, errdim = 20, err_ar_coef = 0.5, burn = 1000){
  coefMat = matrix(0, nrow = r*s, ncol = 2*n + burn)
  for (i in 1:s) {
    coefMat_i = matrix(0, nrow = r, ncol = 2*n + burn)
    coefMat_i[, 1] = runif(r, -1, 1)
    for(j in 2:ncol(coefMat_i)){
      coefMat_i[, j] = varMat %*% coefMat_i[, j-1] + rnorm(r, 0, sigma*i^(-q))
    }
    coefMat[((i-1)*r+1):(i*r), ] = coefMat_i
  }
  coefMat = coefMat[, -(1:burn)]
  
  data = array(0, dim = c(p, dim(bmat)[2], n))
  facto = array(0, dim = c(r, dim(bmat)[2], n))
  for (i in 1:(n)) {
    coefi = matrix(coefMat[, i], nrow = r, ncol = s)
    xn = coefi %*% bmat[1:s, ]
    facto[ , , i] = xn
    ### construct error process
    epscoef = matrix(rnorm(p*errdim, 0, 1), nrow = p, ncol = errdim)
    epscoef = t(t(epscoef) * 1/2^((1:errdim)+1))
    epsn = epscoef %*% bmat[1:errdim, ]
    epsn = diag(sample(1:10, p, replace = TRUE)) %*% epsn / 5
    data[, , i] = A %*% xn * delta + epsn
  }
  return(data)
}
# serial dependent idiosyccratic terms
data_gen_4 = function(n, p, r, s, A, varMat, bmat, delta = 1, q = 0.75,
                      sigma = 1, errdim = 20, err_ar_coef = 0.5,cross_depend = 0.25,
                      burn = 1000){
 
  coefMat = matrix(0, nrow = r*s, ncol = n + burn)
  for (i in 1:s) {
    coefMat_i = matrix(0, nrow = r, ncol = n + burn)
    coefMat_i[, 1] = runif(r, -1, 1)
    for(j in 2:ncol(coefMat_i)){

      coefMat_i[, j] = varMat %*% coefMat_i[, j-1] + rnorm(r, 0, sigma*i^(-q))
    }
    coefMat[((i-1)*r+1):(i*r), ] = coefMat_i
  }
  coefMat = coefMat[, -(1:burn)]
  
  err_coefMat = matrix(0, nrow = p*errdim, ncol = n + burn)
  for (i in 1:errdim) {
    coefMat_i = matrix(0, nrow = p, ncol = n + burn)
    coefMat_i[, 1] = runif(p, -1, 1)
    for(j in 2:ncol(coefMat_i)){
      coefMat_i[, j] = err_ar_coef * coefMat_i[, j-1] + rnorm(p, 0, 1)
    }
    err_coefMat[((i-1)*p+1):(i*p), ] = coefMat_i
  }
  err_coefMat = err_coefMat[, -(1:burn)]
  
  data = array(0, dim = c(p, dim(bmat)[2], n))
  facto = array(0, dim = c(r, dim(bmat)[2], n))
  for (i in 1:(n)) {
    coefi = matrix(coefMat[, i], nrow = r, ncol = s)
    xn = coefi %*% bmat[1:s, ]
    facto[ , , i] = xn
    ### construct error process, serial dependence is allowed
    epscoef = matrix(err_coefMat[, i], nrow = p, ncol = errdim) # matrix(rnorm(p*errdim, 0, 1), nrow = p, ncol = errdim)
    epscoef = t(t(epscoef) * 1/2^((1:errdim)+1))
    epsn = epscoef %*% bmat[1:errdim, ]

    data[, , i] = A %*% xn * delta + epsn
  }
  return(data)
}
# cross-correlated + serial-dependent idiosyccratic terms
data_gen_5 = function(n, p, r, s, A, varMat, bmat, delta = 1, q = 0.75,
                      sigma = 1, errdim = 20, err_ar_coef = 0.5, cross_depend = 0.25,
                      burn = 1000){
   coefMat = matrix(0, nrow = r*s, ncol = n + burn)
  for (i in 1:s) {
    coefMat_i = matrix(0, nrow = r, ncol = n + burn)
    coefMat_i[, 1] = runif(r, -1, 1)
    for(j in 2:ncol(coefMat_i)){
      coefMat_i[, j] = varMat %*% coefMat_i[, j-1] + rnorm(r, 0, sigma*i^(-q))
    }
    coefMat[((i-1)*r+1):(i*r), ] = coefMat_i
  }
  coefMat = coefMat[, -(1:burn)]
  
  err_coefMat = matrix(0, nrow = p*errdim, ncol = n + burn)
  for (i in 1:errdim) {
    coefMat_i = matrix(0, nrow = p, ncol = n + burn)
    coefMat_i[, 1] = runif(p, -1, 1)
    for(j in 2:ncol(coefMat_i)){
        coefMat_i[, j] = err_ar_coef * coefMat_i[, j-1] + rnorm(p, 0, 1)
    }
    err_coefMat[((i-1)*p+1):(i*p), ] = coefMat_i
  }
  err_coefMat = err_coefMat[, -(1:burn)]
  
  data = array(0, dim = c(p, dim(bmat)[2], n))
  facto = array(0, dim = c(r, dim(bmat)[2], n))
  for (i in 1:(n)) {
    coefi = matrix(coefMat[, i], nrow = r, ncol = s)
    xn = coefi %*% bmat[1:s, ]
    facto[ , , i] = xn
    epscoef = matrix(err_coefMat[, i], nrow = p, ncol = errdim) # matrix(rnorm(p*errdim, 0, 1), nrow = p, ncol = errdim)
    epscoef = t(t(epscoef) * 1/2^((1:errdim)+1))
    epsn = epscoef %*% bmat[1:errdim, ]

    data[, , i] = A %*% xn * delta + toeplitz(cross_depend^(0:(p-1))) %*% epsn 
  }
  return(data)
}
data_gen_list = list(data_gen_1, data_gen_2, data_gen_3, data_gen_4, data_gen_5)

