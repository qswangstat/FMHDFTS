library(foreach)
library(doSNOW)
library(doRNG)
source('../utils/data_generate.R')

simulation = function(p=100, type=1, if_weak=0, ncol.Q=12, k0=4,
                      if_verify=False) {
    set.seed(2020)
    #######################################
    ######  data generating process  ######
    #######################################
    data_gen = data_gen_list[[type]]
    #######################################
    ######    start simulation        #####
    ######  basic parameter setting  ######
    #######################################
    sim = 120
    n = 100
    r = 4                       
    s = 51                      
    max_nb = 51                        ## max number of basis functions
    len = 21                           ## number of observed time points
    u = seq(0, 1, length.out = len)    ## sequence of observed time points
    h = 1 / (len-1)                    ## time points interval
    bmat = matrix(0, nrow = max_nb, ncol = len) ## basis function matrix
    for (i in 1:max_nb) {
      if(i == 1) bmat[i, ] = 1
      else if(i%%2 == 0) bmat[i, ] = sqrt(2)*sin(i*pi*u)
      else bmat[i, ] = sqrt(2)*cos((i-1)*pi*u)
      
    } 
    rho = 0.45
    varMat = toeplitz(rho^(1:r)) 
    A = matrix(runif(p*(r), -1, 1), nrow = p, ncol = r) * sqrt(3)
    if(if_weak == 1) A = A / p^(0.25) 
    Q = qr.Q(qr(A))
    if(type == 3) Q = qr.Q(qr(cbind(A, 1))) ## scenario 2

    delta_seq = c(seq(0.1,0.2,0.01), seq(0.25,0.5,0.025))#c(seq(0.1,0.5,0.025))
    ## signal strength param adjustment
    if(type == 2) delta_seq = c(seq(0.5,1,0.025), seq(1.25,2.5,0.25)) / 5 
    if(if_weak == 1){
      delta_seq = seq(0.3, 1.25, 0.05)
    }
    allres = list()
    begin = Sys.time()
    numCores = min(parallel::detectCores()-1, 30)

    for (delta in delta_seq) {
      print(delta)
      print(Sys.time())
      clus = makeCluster(numCores, outfile='')
      registerDoSNOW(clus)
      result = foreach(j = 1:sim, .combine = 'rbind') %dorng% {
        cat("Simualtion", j, '\n')
        source('../utils/lag_covariance_estimate.R')
        source('../utils/space_estimation_accuracy.R')
        Qmat = matrix(runif(p*ncol.Q, -1, 1), nrow = p, ncol = ncol.Q) 
        y = data_gen(n, p, r, s, A, varMat, bmat, 
                     delta = delta,  err_ar_coef = 0.5,
                     cross_depend = 0.5, q = 0.75, 
                     sigma = 1, errdim = 20, burn = 1000)
        
        mMatArr1 = lag_cov_prod_int(y, h, k0)
        mMatArr2 = lag_cov_prod_int_s(y, h, k0)
        mMatArr3 = lag_cov_prod_int_scale(y, h, k0, Qmat)
        
        res1 = analysis(mMatArr1, Q, r=r)
        res2 = analysis(mMatArr2, Q, r=r)
        res3 = analysis(mMatArr3, Q, r=r)
        
        cat("End", j, '\n')
        
        c(res1, res2, res3)
      }
      
      stopCluster(clus)
      attr(result, 'rng') = NULL
      allres = c(allres, list(result))
      print(Sys.time())
    }
    end = Sys.time()
    elapse = end - begin
    print(elapse)

    filepath = paste("output/signal_strength_type_", type, "_p_", p, ".RData", sep = "")
    if(if_weak == 1){
      filepath = paste("output/weak_signal_strength_type_", type, "_p_", p, ".RData", sep = "")
    }
    if(if_verify==TRUE){
      paste("output/Q_k0_robustness/signal_strength_type_1_p_100_q_",ncol.Q,"_k0_",k0, ".RData", sep='')
    }
    save.image(filepath)
}


