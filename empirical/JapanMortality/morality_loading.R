load('mortality.RData')
sex = c('male', 'female', 'total')
all.data = list(data.sp.m, data.sp.f, data.sp)

library(doRNG)
library(foreach)
library(parallel)
library(doSNOW)
library(FMfts)
source('lag_cov_new.R')
source('mortality_func.R')

for(id in 1:3){
  #filename = paste('loading_', sex[id], '_new.RData', sep = '')
  data.train = all.data[[id]]

  k0 = 2
  ncol.Q = 12
  rhat = 2
  
  p          = dim(data.train)[1]
  t          = dim(data.train)[2]
  n.train    = dim(data.train)[3]
  h = 1/t
  ### centering
  ave        = apply(data.train, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train))
  data.train = data.train - aveArr
  ####################################
  ###### Sample Auto-Covariance ######
  ####################################
  Qmat       = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q)
  mMatArr.st = lag_cov_prod_int_scale(data.train, h, k0, Qmat)
  mMat       = apply(mMatArr.st, c(1,2), sum)
  mEig       = eigen(mMat)
  hatA       = mEig$vectors[, 1:rhat, drop=F]
  ##########################
  ###### thresholding ######
  ##########################
  data.full = data.train
  #id.train = sample(1:n, size = floor(n/2))
  id.train = 1:floor(n/2)
  data.tr = data.full[,,id.train]
  data.va = data.full[,,-id.train]
  k0 = 2
  hs = FMfts::lag_cov_hs(data.tr, h, k0)
  hs2 = FMfts::lag_cov_hs(data.va, h, k0)
  
  th_seq = seq(min(hs), median(hs), length.out = 5)
  cv.error = matrix(NA, nrow = k0, ncol = length(th_seq))
  
  for(k in 1:k0){
    hs.k = hs[,,k]
    cv.k = rep(NA, length(th_seq))
    for(i in 1:length(th_seq)){
      th = th_seq[i]
      thMat = 1*(hs.k >= rep(th, each = p^2))
      cv.k[i] = CVerror(data.tr, data.va, h, thMat, k)
    }
    cv.error[k, ] = cv.k
  }
  best.id = apply(cv.error, 1, which.min)
  
  sp = 23
  hs_t = FMfts::lag_cov_hs(data.train, h, k0)
  cat(min(hs_t), max(hs_t), '\n')
  thMat = array(NA, dim = dim(hs_t))
  for(k in 1:k0){
    thMat[,,k] = 1*(hs_t[,,k] >= th_seq[best.id[k]])
  }
  Qmat        = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
  mMatArr.th  = lag_cov_prod_int_scale_th(data.train, thMat, h, k0, Qmat)
  mMat        = apply(mMatArr.th, c(1,2), sum)
  hatA.sp     = gdefla(mMat, rep(sp, rhat))
  ###################
  filename = paste('loading_', sex[id],'_sp_',sp, '.RData', sep = '')
  
  save(hatA, hatA.sp, file = filename)
  cat('End', id, '\n')
}

