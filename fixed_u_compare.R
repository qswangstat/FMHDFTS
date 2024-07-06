p = 1
type = 1

set.seed(2020)
library(foreach)
library(doSNOW)
library(doRNG)

source('../utils/data_generate.R')
source('../utils/lag_covariance_estimate.R')
source('../utils/space_estimation_accuracy.R')

data_gen = data_gen_list[[type]]

#######################################
######    start simulation        #####
######  basic parameter setting  ######
#######################################

sim = 120
n = 100
r = 4                       
s = 51                          
max_nb = 51
k0 = 4                       
len = 21
u = seq(0, 1, length.out = len) 
u_est_id = seq(3,19,2)
h = 1 / (len-1)
ncol.Q = 12
bmat = matrix(0, nrow = max_nb, ncol = len)
for (i in 1:max_nb) {
  if(i == 1) bmat[i, ] = 1
  else if(i%%2 == 0) bmat[i, ] = sqrt(2)*sin(i*pi*u)
  else bmat[i, ] = sqrt(2)*cos((i-1)*pi*u)
  
} 
rho = 0.45
varMat = toeplitz(rho^(1:r)) ## coefficient matrix for VAR(1) model
A = matrix(runif(p*(r), -1, 1), nrow = p, ncol = r) * sqrt(3)
Q = qr.Q(qr(A))
delta_seq = c(seq(0.1, 0.2,0.01),seq(0.25,0.5,0.025))


allres = list()
begin = Sys.time()
numCores = min(parallel::detectCores()-1, 30)
        
combineNestedLists <- function(x, y) {
  if (is.null(x)) {
    list(y)
  } else {
    c(x, list(y))
  }
}
for (delta in delta_seq) {
  print(delta)
  print(Sys.time())
  clus = makeCluster(numCores, outfile='')
  registerDoSNOW(clus)
  result = foreach(j = 1:sim, .combine = combineNestedLists) %dorng% {
    cat("Simualtion", j, '\n')

    y = data_gen(n, p, r, s, A, varMat, bmat, delta = delta, q = 0.75, sigma = 1, errdim = 20, burn = 1000)
    Qmat = matrix(runif(p*ncol.Q, -1, 1), nrow = p, ncol = ncol.Q) 

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
print(elapse)

save.image(paste("output/signal_strength_type_", type, "_p_", p, "_fixed_u.RData", sep = ''))

