#!/usr/bin/env Rscript
args = as.numeric(commandArgs(trailingOnly=TRUE))
p = args[1]
type = args[2]


set.seed(2021)
library(foreach)
library(doSNOW)
library(doRNG)
library(Rcpp)
library(RcppArmadillo)
library(FMfts)
source('lag_covariance_estimate.R')
source('data_generate.R')

data_gen = data_gen_list[[1]]


sim = 120
n = 100
r = 4                          ## number of factors
s = 51                          ## number of basis functions
max_nb = 51
k0 = 4                          ## upper bound of lag used
len = 21
u = seq(0, 1, length.out = len) ## time points, fully observed
h = 1 / (len-1)
CVnum = 30
bmat = matrix(0, nrow = max_nb, ncol = len)
for (i in 1:max_nb) {
  if(i == 1) bmat[i, ] = 1
  else if(i%%2 == 0) bmat[i, ] = sqrt(2)*sin(i*pi*u)
  else bmat[i, ] = sqrt(2)*cos((i-1)*pi*u)
  
} ## basis function values at u

rho = 0.45
varMat = toeplitz(rho^(1:r)) ## coefficient matrix for VAR(1) model


if(type == 1){
  A = matrix(runif(p*r, -1, 1), nrow = p, ncol = r) * sqrt(3)
  A[sample(p, 0.8*p), ] = 0
} ### row sparsity, strong signal
if(type == 2){
  A = matrix(runif(p*r, -1, 1), nrow = p, ncol = r) / p^(0.125)
  A[sample(p, 0.8*p), ] = 0
} ### row sparsity, weak signal
if(type == 3){
  A = matrix(runif(p*r, -1, 1), nrow = p, ncol = r) * sqrt(3)
  for (i in 1:r) {
    A[sample(p, 0.8*p), i] = 0
  }
} ### column sparsity

Q = qr.Q(qr(A))
time.begin = Sys.time()

#######################################################################################################
#### cross validation: two steps
###############################################
#### first step: choose proper thresholding

analysis = function(mMatArr, r = NULL) {
  p = dim(mMatArr)[1]
  ## threshold on each slice of mMatArr
  mMat = apply(mMatArr, c(1,2), sum)
  mEig = eigen(mMat)
  val = mEig$values
  vec = mEig$vectors
  
  R = min(p/2, max(which(val > 1e-2)))
  ratio = val[2:p] / val[1:p-1]
  rhat_true = which.min(ratio[1:R])
  if(!is.null(r)) rhat = r
  
  hatA = vec[, 1:rhat]
  
  return(list(r = rhat, hatA = hatA))
}

ncol.Q = 12
Qmat = matrix(runif(p*ncol.Q, -1, 1), nrow = p, ncol = ncol.Q) * sqrt(3)

t_seq = seq(0.09, 0.45, 0.03) #type = 1,3
if(type == 2)  t_seq = seq(0.01, 0.5, 0.02) #type=2

cat('CV: choose thresholding value', '\n')
numCores = min(parallel::detectCores()-1, CVnum + 1)
clus = makeCluster(numCores, outfile='')
registerDoSNOW(clus)
res_th = foreach(j = 1:CVnum) %dorng% {
  library(FMfts)
  cat('Cross validation', j, 'starts.', '\n')

  data = data_gen(2*n, p, r, s, A, varMat, bmat)
  y_train = data[,,1:n]
  y_valid = data[,,-(1:n)]

  hs_t = lag_cov_hs(y_train, h, k0)
  mMatArr_v = lag_cov_prod_int_scale_th(y_valid, array(1, c(p, p, k0)), h, k0, Qmat)
  res_v = analysis(mMatArr_v, r)
  A_v = res_v$hatA
  r_v = res_v$r
  temp = Q%*%t(Q)%*%A_v%*%t(A_v)
  err_v = sqrt(1 - sum(diag(temp)) / max(r, r_v))

  resCV = rep(0, length(t_seq))
  resTrue = rep(0, length(t_seq))
  R_t = rep(0, length(t_seq))
  for (i in 1:length(t_seq)) {
    t = t_seq[i]
    th_val = rep(t, k0)
    thMat = 1*(hs_t >= rep(th_val, each = p^2))
    mMatArr_t = lag_cov_prod_int_scale_th(y_train, thMat, h, k0, Qmat)
    res_t = analysis(mMatArr_t, r)
    A_t = res_t$hatA
    r_t = res_t$r
    temp = A_v%*%t(A_v)%*%A_t%*%t(A_t)
    resCV[i] = sqrt(1 - sum(diag(temp)) / max(r_v, r_t))
    temp = Q%*%t(Q)%*%A_t%*%t(A_t)
    resTrue[i] = sqrt(1 - sum(diag(temp)) / max(r, r_t))
    R_t[i] = r_t
  }
  cat(r_v, '\n')
  #th_val = t_seq[which.min(resCV)]
  #cat('threhold:', th_val, '\n')
  cat('Cross validation', j, 'ends.', '\n')

  list(cv_err = resCV,
       cv_r = R_t,
       th_err = resTrue,
       res_v = c(r_v, err_v))
}
stopCluster(clus)
attr(res_th, 'rng') = NULL

cv_err = sapply(res_th, function(x) x$cv_err)
opt_id = which.min(apply(cv_err, 1, mean))
opt_th = t_seq[opt_id]
# 
# #######################################################################################################
# ###############################################
# #### second step: choose sparsity level
sp_seq = seq(0.05*p, 0.5*p, 5)
cat('CV: choose sparsity level', '\n')

numCores = min(parallel::detectCores()-1, CVnum + 1)
clus = makeCluster(numCores, outfile='')
registerDoSNOW(clus)
res_sp = foreach(j = 1:CVnum) %dorng% {
  library(FMfts)
  cat('Cross validation', j, 'starts.','\n')
  
  data = data_gen(2*n, p, r, s, A, varMat, bmat)
  y_train = data[,,1:n]
  y_valid = data[,,-(1:n)]
  
  hs_v = lag_cov_hs(y_valid, h, k0)
  thMat_v = 1*(hs_v >= rep(0, each = p^2))
  mMatArr_v = lag_cov_prod_int_scale_th(y_valid, thMat_v, h, k0, Qmat)
  mMat_v = apply(mMatArr_v, c(1,2), sum)
  res_v = analysis(mMatArr_v)
  r_v = res_v$r
  A_v = res_v$hatA
  temp = Q%*%t(Q)%*%A_v%*%t(A_v)
  err_v = sqrt(1 - sum(diag(temp)) / max(r, r_v))
  
  hs_t = lag_cov_hs(y_train, h, k0)
  thMat = 1*(hs_t >= rep(0, each = p^2))
  mMatArr_t = lag_cov_prod_int_scale_th(y_train, thMat, h, k0, Qmat)
  mMat = apply(mMatArr_t, c(1,2), sum)
  res_t = analysis(mMatArr_t)
  r_t = res_t$r
  
  resCV2 = rep(0, length(sp_seq))
  resTrue = resCV2
  
  for (i in 1:length(sp_seq)) {
    sp = sp_seq[i]
    Adefl = gdefla(mMat, rep(sp, r_t))
    temp = Adefl %*% t(Adefl) %*% A_v %*% t(A_v)
    resCV2[i] = sqrt(1 - sum(diag(temp)) / max(r_v, r_t))
    temp = Q%*%t(Q)%*%Adefl%*%t(Adefl)
    resTrue[i] = sqrt(1 - sum(diag(temp)) / max(r, r_t))
  }
  
  #sp_val = sp_seq[which.min(resCV2)]
  cat('Cross validation', j, 'ends.', '\n')
  
  list(cv_err = resCV2,
       sp_err = resTrue,
       cv_r = r_t, 
       res_v = c(r_v, err_v))
}
stopCluster(clus)
attr(res_sp, 'rng') = NULL

cv_err = sapply(res_sp, function(x) x$cv_err)
opt_id2 = which.min(apply(cv_err, 1, mean))
opt_sp = sp_seq[opt_id2]
# 
# #######################################################################################################
# ##### begin simulation
opt_th
opt_sp = 0.2 * p

cat('Now simulation begins', '\n')

numCores = min(parallel::detectCores()-1, 120)
clus = makeCluster(numCores, outfile='')
registerDoSNOW(clus)

res_final = foreach(id = 1:sim) %dorng% {
  library(FMfts)
  cat('Simulation', id, 'starts.', '\n')
  ## data generation
  data = data_gen(2*n, p, r, s, A, varMat, bmat)
  y_arr = data[,,1:n]
  y_valid = data[,,-(1:n)]
  
  ## estimate lag-k covariance
  mMatArr_no = lag_cov_prod_int_scale_th(y_arr, array(1, c(p,p,k0)), h, k0, Qmat)
  mMat_no = apply(mMatArr_no, c(1,2), sum)
  no = analysis(mMatArr_no, Q, r)
  rhat_no = no$r
  Adefl_no = gdefla(mMat_no, rep(opt_sp, rhat_no))
  temp_no = Adefl_no %*% t(Adefl_no) %*% Q %*% t(Q)
  trc_no = sum(diag(temp_no))
  meas_no = sqrt(1 - trc_no / max(r, rhat_no))
  sp_noth = list(sp_noth_dist = meas_no, sp_noth_vec = Adefl_no)
  
  # cat(' r = ', no$rhat, ', dist = ', round(no$dist, 3), '.', sep = '', '\n')
  
  hs = lag_cov_hs(y_arr, h, k0)
  thMat = 1*(hs >= rep(opt_th, each = p^2))
  mMatArr_yes = lag_cov_prod_int_scale_th(y_arr, thMat, h, k0, Qmat)
  ## yes only thresholding
  yes = analysis(mMatArr_yes, Q, r)
  mMat = apply(mMatArr_yes, c(1,2), sum)
  rhat = yes$r
  # cat(' r = ', yes$rhat, ', dist = ', round(yes$dist, 3), '.', sep = '', '\n')
  cat('Simulation', id, 'ends.', '\n')
  ## 
  Adefl = gdefla(mMat, rep(opt_sp, rhat))
  temp = Adefl %*% t(Adefl) %*% Q %*% t(Q)
  trc = sum(diag(temp))
  meas = sqrt(1 - trc / max(r, rhat))
  spth = list(spth_dist = meas, spth_vec = Adefl)
  c(no, yes, sp_noth, spth)
}
stopCluster(clus)
attr(res_final, 'rng') = NULL


time.end = Sys.time()
time.end - time.begin

op_no_dist = sapply(res_final, function(x) x[[2]])
op_th_dist = sapply(res_final, function(x) x[[6]])
sp_no_dist = sapply(res_final, function(x) x[[9]])
sp_th_dist = sapply(res_final, function(x) x[[11]])

save.image(paste('new_simulation_sparse_type_', type, '_p_', p, '.RData', sep = ''))

sink(paste('Simulation_output_type_', type, '_p_', p, '.txt', sep = ''))
cat('A. Simulation setting:', '\n\n')
cat('\t\t', '0. basic model: y_t = A*x_t + eps_t', '\n\n')
cat('\t\t', '1. basis are Fourier basis {1, sin(2*k*pi*x), cos(2*k*pi*x)}', '\n\n')
cat('\t\t', '2. idiosyncratic part eps_t can be expanded under first 20 Fourier basis
    and the i-th coefficient are N(0,1)*(0.5)^(i-1)', '\n\n')
if(type == 1) cat('\t\t', '3. matrix A is p*r and row-sparse with only 20% of rows nonzero, and nonzero elements
    are generated from unif[-sqrt(3), sqrt(3)] (strong factors)', '\n\n')
if(type == 2) cat('\t\t', '3. matrix A is p*r and row-sparse with only 20% of rows nonzero, and nonzero elements
    are generated from unif[-1, 1]/p^0.125 (weak factors)', '\n\n')
if(type == 3) cat('\t\t', '3. matrix A is p*r and column-sparse with only 20% of rows nonzero, and nonzero elements
    are generated from unif[-sqrt(3), sqrt(3)] (strong factors)', '\n\n')
cat('B. Tuning parameter choice:', '\n\n')
cat('\t', '0. run 30 cross-validtions, each containing a training set and a validation set', '\n\n')

cat('\t', '1. firstly choose thresholding', '\n\n')
cat('\t\t', 'training with threshold while validation without', '\n\n')
cat('\t\t', 'criterion is the distance bewteen column spaces of the recovered A', '\n\n')
cat('\t\t', 'for different lags, threshold values may vary, but here all are set equal', '\n\n')

cat('\t', '2. secondly choose sparsity level', '\n\n')
cat('\t\t', 'both training and validation sets are thresholded with the optimal threshold value', '\n\n')
cat('\t\t', 'apply sparse PCA on training set, while ordinary PCA on validation set', '\n\n')
cat('\t\t', 'criterion is the distance bewteen column spaces of the recovered A', '\n\n')
cat('\t\t', 'generalized deflation method for sparse PCA are applied', '\n\n')
cat('\t\t', 'again for different columns, sparsity may vary, but here all are set equal', '\n\n')
cat('\t\t', 'After averaging, level may not be multiple of 10, for simplicity, we artificially approximate', '\n\n')

cat('\t', '3. choose r based on eigen ratio of M', '\n\n')


cat('C. Simulation results:', '\n\n')

cat('\t', '0. Optimal thresholding =', opt_th, '\n\n')
cat('\t', '1. Optimal sparsity level =', opt_sp, '\n\n')


cat('Ordinary PCA without thresholding: mean =', mean(op_no_dist), ', sd =', sd(op_no_dist), '\n\n')
cat('Ordinary PCA with thresholding   : mean =', mean(op_th_dist), ', sd =', sd(op_th_dist), '\n\n')
cat('Sparse PCA without thresholding  : mean =', mean(sp_no_dist), ', sd =', sd(sp_no_dist), '\n\n')
cat('Sparse PCA with thresholding     : mean =', mean(sp_th_dist), ', sd =', sd(sp_th_dist), '\n\n')
cat('\t', '2. Each simulation ---- ', '\n\n')
for (i in 1:sim) {
  cat('\t\t', 'Simulation --', i, '\n\n')
  temp = res_final[[i]]
  cat('\t\t\t', 'ordinary PCA', '\n\n')
  cat('\t\t\t\t', 'no threhsold: rhat=', temp[[1]], 'dist =', temp[[2]], '\n\n')
  cat('\t\t\t\t', 'threhsolding: rhat=', temp[[5]], 'dist =', temp[[6]], '\n\n')
  cat('\t\t\t', 'sparse PCA', '\n\n')
  cat('\t\t\t\t', 'no threhsold: rhat=', temp[[1]], 'dist =', temp[[9]], '\n\n')
  cat('\t\t\t\t', 'threhsolding: rhat=', temp[[5]], 'dist =', temp[[11]], '\n\n')

}
sink()



