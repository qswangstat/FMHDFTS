#!/usr/bin/env Rscript
args = as.numeric(commandArgs(trailingOnly=TRUE))
type = args[1]
rhat = 2

load('mortality.RData')
for(sex in 1:2){
  if(sex == 1) {
    data = data.sp.m
    ave        = apply(data, c(1,2), mean)
    aveArr     = array(ave, dim = dim(data))
    data = data - aveArr
    if(type==1) filename = paste('mortality_male_sparse_r_', rhat, '_test_ahead_10.RData', sep = '')
    if(type==2) filename = paste('mortality_male_final_new_demean_ahead_10.RData')
    if(type==3) filename = paste('mortality_male_tnh_k=1_demean_ahead_10.RData')
  }
  if(sex == 2) {
    data = data.sp.f
    ave        = apply(data, c(1,2), mean)
    aveArr     = array(ave, dim = dim(data))
    data = data - aveArr
    if(type==1) filename = paste('mortality_female_sparse_r_', rhat, '_test_ahead_10.RData', sep = '')
    if(type==2) filename = paste('mortality_female_final_new_demean_ahead_10.RData')
    if(type==3) filename = paste('mortality_female_tnh_k=1_demean_ahead_10.RData')
    
  }
  if(sex == 3) {
    data = data.sp
    filename = paste('mortality_total_sparse_r_', rhat, '_test.RData', sep = '')
    filename = paste('mortality_total_final_new.RData')
    
  }
  library(doRNG)
  library(foreach)
  library(parallel)
  library(doSNOW)
  source('lag_cov_new.R')
  source('mortality_func.R')
  
  k0 = 2
  ncol.Q = 12
  n.train = 26
  n0 = 1
  n = dim(data)[3]
  C.seq = seq(0,3,by=0.1)
  
  begin = Sys.time()
  numCores = min(parallel::detectCores()-1, (n-n.train+1))
  clus = makeCluster(numCores, outfile=paste('sex_', sex, '.txt',sep=''))
  registerDoSNOW(clus)
  
  res_final = foreach(delta = 1:(n-n.train)) %dorng% {
    library(vars)
    library(fdapace)
    library(splines)
    library(forecast)
    library(ftsa)
    library(FMfts)
    cat('Begin', delta,'\n')
    id.train   = n0:(n.train + delta - 1)
    data.train = data[ , , id.train]
    ahead.max  = min(10, n-n.train+1-delta)
    if(type==1) res.delta  = fm_sparse_pred(data.train, h, k0, ncol.Q, tgrid, ahead.max, rhat.seq = rhat)
    if(type==2) res.delta  = pred(data.train, h, k0, ncol.Q, tgrid, ahead.max)
    if(type==3) res.delta  = tnh(data.train, tgrid, k.seq = 1, ahead.max = 10)
    #res.delta = tnh_select_k(data.train, tgrid, C.seq, J=10, P=5, kmax=10)
    
    cat('End', delta,'\n')
    
    res.delta
  }
  stopCluster(clus)

  end = Sys.time()
  elapse = end - begin
  print(elapse)
  attr(res_final, 'rng') = NULL
  data.train = data[,,n0:n.train]
  data.test  = data[,,-(n0:n.train)]
  cat(filename, '\n')
  
  save(elapse, data.train,data.test,res_final, file = paste('./', filename, sep=''))
}
