#!/usr/bin/env Rscript
args = as.numeric(commandArgs(trailingOnly=TRUE))
p = args[1]
type = args[2]

for(ncol.Q in c(12)){
    if(ncol.Q==12) k0_seq = 4
    else k0_seq = 1:5
    for(k0 in k0_seq){
        set.seed(2020)
library(foreach)
library(doSNOW)
library(doRNG)
#library(Rcpp)
#library(RcppArmadillo)
source('data_gen.R')
source('lag_cov.R')
#######################################
######  data generating process  ######
#######################################
data_gen = data_gen_list[[type]]


analysis = function(mMatArr, Q, r = NULL) {
  if(length(dim(mMatArr)) == 3)
    mMat = apply(mMatArr, c(1,2), sum)
  else
    mMat = mMatArr %*% t(mMatArr)
  p = dim(mMatArr)[1]
  ## threshold on each slice of mMatArr
  mEig = eigen(mMat)
  val = mEig$values
  vec = mEig$vectors
  
  R = p/2
  ratio = val[2:p] / val[1:p-1]
  rr = which.min(ratio[1:R])
  if(!is.null(r)) rhat = r
  hatA = vec[, 1:rhat]
  
  temp = hatA%*%t(hatA)%*%Q%*%t(Q)
  err = sqrt(1 - sum(diag(temp)) / max(ncol(Q), rhat))
  
  c(rr, err)
}
#######################################
######    start simulation        #####
######  basic parameter setting  ######
#######################################

sim = 120
n = 100
r = 4                       ## number of factors
s = 51                          ## number of basis functions
max_nb = 51
#k0 = 4                        ## upper bound of lag used
len = 21
u = seq(0, 1, length.out = len) ## time points, fully observed
u_est_id = seq(3,19,2)
h = 1 / (len-1)
CVnum = 30
bmat = matrix(0, nrow = max_nb, ncol = len)
for (i in 1:max_nb) {
  if(i == 1) bmat[i, ] = 1
  else if(i%%2 == 0) bmat[i, ] = sqrt(2)*sin(i*pi*u)
  else bmat[i, ] = sqrt(2)*cos((i-1)*pi*u)
  
} 
rho = 0.45
varMat = toeplitz(rho^(1:r)) ## coefficient matrix for VAR(1) model

#varMat = varMat / max(Mod(eigen(varMat)$values)) * 0.9
#eigen(varMat)$values

A = matrix(runif(p*(r), -1, 1), nrow = p, ncol = r) * sqrt(3)
#A = A / p^(0.25) ### weak signal
Q = qr.Q(qr(A))
if(type == 2) Q = qr.Q(qr(cbind(A, 1)))
R = qr.R(qr(A))
#nohup Rscript simulation.R 200 3 > ./print3.log 2> ./error3.log &

Qmat = matrix(runif(p*ncol.Q, -1, 1), nrow = p, ncol = ncol.Q) 
delta_seq = c(seq(0.1, 0.2,0.01),seq(0.25,0.5,0.025))
if(type == 3) delta_seq = c(seq(0.5,1,0.025), seq(1.25,2.5,0.25)) / 5
#delta_seq = c(0.15, 0.3, 0.45)
#if(type == 3) delta_seq = c(0.75, 1.5, 2.25)
## weak
#delta_seq = seq(0.3, 0.8, 0.025)

allres = list()
begin = Sys.time()
numCores = min(parallel::detectCores()-1, 30)
        
combineNestedLists <- function(x, y) {
  # 如果 x 为空（在第一个元素时），则创建一个列表
  if (is.null(x)) {
    list(y)
  } else {
    # 在非空列表上追加
    c(x, list(y))
  }
}
#delta = 1
#if(type == 2) delta = 0.2
for (delta in delta_seq) {
  print(delta)
  print(Sys.time())
  clus = makeCluster(numCores, outfile='')
  registerDoSNOW(clus)
  result = foreach(j = 1:sim, .combine = combineNestedLists) %dorng% {
    #for(j in 1:10) {
    cat("Simualtion", j, '\n')

    y = data_gen(n, p, r, s, A, varMat, bmat, delta = delta, q = 0.75, sigma = 1, errdim = 20, burn = 1000)
    
    ### 
    rhat = rep(NA, length(u_est_id))
    hatA_list = array(NA, c(p, r, length(u_est_id)))
    for(id in 1:length(u_est_id)){
        uid = u_est_id[id]
        mMatArr = array(NA, c(p, p, k0))

        for(k in 1:k0){
              mMat = matrix(0, p, p)
              yuk = y[ , uid, (k+1):n]
              yuk2 = y[ , uid, 1:(n-k)]

              yumean = apply(y[ , uid, ], 1, mean)
              yuk = yuk - yumean
              yuk2 = yuk2 - yumean

              sMat = (yuk %*% t(yuk2)) / (n-k)
              sigMat = sMat %*% t(sMat)

              mMatArr[ , , k] = sigMat
        }
        mMat = apply(mMatArr, c(1,2), sum)
        mEig = eigen(mMat)
        val = mEig$values
        vec = mEig$vectors
        R = p/2
        ratio = val[2:p] / val[1:p-1]
        rr = which.min(ratio[1:R])
        hatA = vec[, 1:r]
        
        rhat[id] = rr
        hatA_list[,,id] = hatA
    }
    
    cat("End", j, '\n')
    
    list(r_est = rhat, hatA_list=hatA_list)
  }
  
  stopCluster(clus)
  attr(result, 'rng') = NULL
  allres = c(allres, list(result))
  print(Sys.time())
}



end = Sys.time()
elapse = end - begin

elapse

save.image(paste("signal_strength_type_", type, "_p_", p, "_fixed_u.RData", sep = ''))

    }
}

