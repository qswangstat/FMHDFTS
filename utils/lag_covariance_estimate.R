#######################################
######  matrix integral function  #####
#######################################
lag_cov_prod_int = function(y, h, k0){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  mMatArr = array(NA, c(p, p, k0))
  for(k in 1:k0){
    mMat = matrix(0, p, p)
    for(uid in 1:len){
      yuk = y[ , uid, (k+1):n]
      yumean = apply(y[ , uid, ], 1, mean)
      yuk = yuk - yumean
      for(vid in 1:len){
        yvk = y[ , vid, 1:(n-k)]
        yvmean = apply(y[ , vid, ], 1, mean)
        yvk = yvk - yvmean
        
        idinfo = (uid == 1 | uid == len) + (vid == 1 | vid == len)
        trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
        
        sMat = (yuk %*% t(yvk)) / (n-k)
        sigMat = sMat %*% t(sMat)
        
        mMat = mMat + sigMat * h^2 * trapcoef
      }
    }
    mMatArr[ , , k] = mMat;
  }
  mMatArr
}
lag_cov_prod_int_th = function(y, thMat, h, k0, Qmat){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  mMatArr = array(NA, c(p, p, k0))
  for(k in 1:k0){
    indMat = thMat[,,k]
    mMat = matrix(0, p, p)
    for(uid in 1:len){
      yuk = y[ , uid, (k+1):n]
      yumean = apply(y[ , uid, ], 1, mean)
      yuk = yuk - yumean
      for(vid in 1:len){
        yvk = y[ , vid, 1:(n-k)]
        yvmean = apply(y[ , vid, ], 1, mean)
        yvk = yvk - yvmean
        
        idinfo = (uid == 1 | uid == len) + (vid == 1 | vid == len)
        trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
        
        sMat = (yuk %*% t(yvk)) * indMat / (n-k)
        sigMat = sMat %*% t(sMat)
        
        mMat = mMat + sigMat * h^2 * trapcoef
      }
    }
    mMatArr[ , , k] = mMat;
  }
  mMatArr
}
lag_cov_prod_int_s = function(y, h, k0){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  mMatArr = array(NA, c(p, p, k0))
  
  for(k in 1:k0){
    mMat = matrix(0, p, p)
    for(uid in 1:len){
      yuk = y[ , uid, (k+1):n]
      yuk2 = y[ , uid, 1:(n-k)]
      
      yumean = apply(y[ , uid, ], 1, mean)
      yuk = yuk - yumean
      yuk2 = yuk2 - yumean
      
      idinfo = (uid == 1 | uid == len)
      trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
      
      sMat = (yuk %*% t(yuk2)) / (n-k)
      sigMat = sMat %*% t(sMat)
      
      mMat = mMat + sigMat * h * trapcoef
    }
    mMatArr[ , , k] = mMat;
  }
  mMatArr
}

lag_cov_prod_int_scale = function(y, h, k0, Q, D = NULL){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  if(is.null(D)) D = diag(rep(1, p))
  D = diag(1/sqrt(diag(D)))
  mMatArr = array(NA, c(p, p, k0))
  for(k in 1:k0){
    mMat = matrix(0, p, p)
    for(uid in 1:len){
      yuk = y[ , uid, (k+1):n]
      yumean = apply(y[ , uid, ], 1, mean)
      yuk = yuk - yumean
      for(vid in 1:len){
        yvk = y[ , vid, 1:(n-k)]
        yv0 = y[ , vid, ]
        yvmean = apply(yv0, 1, mean)
        yvk = yvk - yvmean
        yv0 = yv0 - yvmean
        
        idinfo = (uid == 1 | uid == len) + (vid == 1 | vid == len)
        trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
        
        sMat = (yuk %*% t(yvk)) / (n-k)
        sigma0v = (yv0 %*% t(yv0)) / n
        sigMat = sMat %*% D %*% Q %*% solve(t(Q) %*% D %*% sigma0v %*% D %*% Q) %*% t(Q) %*% D %*% t(sMat)
        #sigMat = sMat %*% D %*% Q %*% solve(t(Q) %*% sigma0v  %*% Q) %*% t(Q) %*% D %*%  t(sMat)
        
        mMat = mMat + sigMat * h^2 * trapcoef
      }
    }
    mMatArr[ , , k] = mMat;
  }
  mMatArr
}

lag_cov_prod_int_scale_th = function(y, thMat, h, k0, Q, D = NULL){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  if(is.null(D)) D = diag(rep(1, p))
  D = diag(1/sqrt(diag(D)))
  mMatArr = array(NA, c(p, p, k0))
  for(k in 1:k0){
    indMat = thMat[,,k]
    mMat = matrix(0, p, p)
    for(uid in 1:len){
      yuk = y[ , uid, (k+1):n]
      yumean = apply(y[ , uid, ], 1, mean)
      yuk = yuk - yumean
      for(vid in 1:len){
        yvk = y[ , vid, 1:(n-k)]
        yv0 = y[ , vid, ]
        yvmean = apply(yv0, 1, mean)
        yvk = yvk - yvmean
        yv0 = yv0 - yvmean
        
        idinfo = (uid == 1 | uid == len) + (vid == 1 | vid == len)
        trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
        
        sMat = (yuk %*% t(yvk)) * indMat / (n-k)
        sigma0v = (yv0 %*% t(yv0)) / n
        sigMat = sMat %*% D %*% Q %*% solve(t(Q) %*% D %*% sigma0v %*% D %*% Q) %*% t(Q) %*% D %*% t(sMat)
        #sigMat = sMat %*% D %*% Q %*% solve(t(Q) %*% sigma0v  %*% Q) %*% t(Q) %*% D %*%  t(sMat)
        
        mMat = mMat + sigMat * h^2 * trapcoef
      }
    }
    mMatArr[ , , k] = mMat;
  }
  mMatArr
}

cov_prod_int = function(y, h){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  mMat = matrix(0, p, p)
  for(uid in 1:len){
    yuk = y[ , uid, ]
    yumean = apply(y[ , uid, ], 1, mean)
    yuk = yuk - yumean
    for(vid in 1:len){
      yvk = y[ , vid, ]
      yvmean = apply(y[ , vid, ], 1, mean)
      yvk = yvk - yvmean
      
      idinfo = (uid == 1 | uid == len) + (vid == 1 | vid == len)
      trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
      
      sMat = (yuk %*% t(yvk)) / (n)
      sigMat = sMat %*% t(sMat)
      
      mMat = mMat + sigMat * h^2 * trapcoef
    }
  }
  mMat
}
cov_prod_int_s = function(y, h){
  ydim = dim(y)
  p = ydim[1]
  len = ydim[2]
  n = ydim[3]
  
  mMat = matrix(0, p, p)
  for(uid in 1:len){
    yuk = y[ , uid, ]

    yumean = apply(y[ , uid, ], 1, mean)
    yuk = yuk - yumean

    idinfo = (uid == 1 | uid == len)
    trapcoef = 1.0 / (idinfo * 2 + (idinfo == 0))
    
    sMat = (yuk %*% t(yuk)) / (n)
    sigMat = sMat %*% t(sMat)
    
    mMat = mMat + sigMat * h * trapcoef
  }
  
  mMat
}
lag_cov_prod_scale = function(y, k0, Q, D = NULL){
  ydim = dim(y)
  p = ydim[1]
  n = ydim[2]
  #if(is.null(D)) D = diag(rep(1, p))
  #D = diag(1/sqrt(diag(D)))
  mMatArr = array(NA, c(p, p, k0))
  ymean =  apply(y, 1, mean)
  y_cen = y - ymean
  sigma0 = (y_cen %*% t(y_cen)) / n
  #D = diag(1/sqrt(diag(sigma0)))
  #sigma0 = D %*% sigma0 %*% D
  for(k in 1:k0){
    sigmak = (y_cen[, (k+1):n] %*% t(y_cen[, 1:(n-k)])) / (n-k)
    #mMatArr[ , , k] = sigmak %*% D %*% Q %*% solve(t(Q) %*% D %*% sigma0 %*% D %*% Q) %*% t(Q) %*% D %*% t(sigmak)
    mMatArr[ , , k] = sigmak %*% Q %*% solve(t(Q) %*% sigma0 %*% Q) %*% t(Q) %*% t(sigmak)
  }
  mMatArr
}
lag_cov_prod = function(y, k0){
  ydim = dim(y)
  p = ydim[1]
  n = ydim[2]
  mMatArr = array(NA, c(p, p, k0))
  ymean =  apply(y, 1, mean)
  y_cen = y - ymean
  for(k in 1:k0){
    sigmak = (y_cen[, (k+1):n] %*% t(y_cen[, 1:(n-k)])) / (n-k)
    mMatArr[ , , k] = sigmak %*% t(sigmak)
  }
  mMatArr
}