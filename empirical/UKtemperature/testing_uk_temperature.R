load('UKstation.RData')
load('temperature.RData')
load('matM_temp.RData')
library(FMfts)
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
          mMat = mMat + sigMat * h^2 * trapcoef
        }
      }
      mMatArr[ , , k] = mMat;
    }
    mMatArr
}
## white noise test ref: Zhang(2016)
WN.test <- function(X, bsize="auto", M=500, nbasis=20, P=50){
    quan <- rep(NA,3)
    R <- dim(X)[1]
    N <- dim(X)[2]

    if(bsize=="auto")  block <- 1:floor(N^(1/2)) else block <- bsize
    basis <- create.bspline.basis(rangeval=c(0,1), nbasis=nbasis, norder=4)
    fd.data <- smooth.basis(0:(R-1)/(R-1), X, basis)$fd

    # Compute the spectra-based test

    data <- eval.fd(fd.data,1:P/P)
    covm <- array(NA, c(P,P,N-1))
    ncov <- rep(NA, N-1)
    for(h in 1:(N-1))
    {
      Gamma <- data[,(1+h):N]%*%t(data[,1:(N-h)])/N
      covm[,,h] <- Gamma+t(Gamma)
      ncov[h] <- sum(covm[,,h]^2)/P^2
    }
    stat <- sum(ncov/(1:(N-1))^2)/pi^2*N/8

    # Block bootstrap

    qb <- matrix(NA,length(block),3)
    pvalue <- rep(NA,length(block))
    for(cc in 1:length(block))
    {
      bN <- block[cc]
      MN <- floor((N-1)/bN)
      stat.b <- rep(NA,M)
      for(bt in 1:M)
      {
        b <- sample(MN, MN, replace=TRUE)
        index <- as.vector(t(matrix((b-1)*bN, MN, bN))+(1:bN))
        if(length(index)<N-1)
        {
          ab <- (sample(MN, 1, replace=TRUE)-1)*bN+(1:bN)
          index <- c(index,ab)[1:(N-1)]
        }
        index <- index+1
        ncovb <- rep(NA,N-1)
        covb1 <- matrix(0,P,P)
        for(hb in 1:(N-1))
        {
          ind <- index[(index-hb)>=1]
          Gb <- data[,ind]%*%t(data[,ind-hb])/N
          covb <- Gb+t(Gb)-covm[,,hb]*length(ind)/N
          ncovb[hb] <- sum(covb^2)/P^2
          covb1 <- covb^2/hb^2+covb1
        }
        stat.b[bt] <- sum(ncovb/(1:(N-1))^2)/pi^2*N/8
      }
      qb[cc,] <- quantile(stat.b,c(0.90,0.95,0.99))
      pvalue[cc] <- mean(stat.b > stat)
    }

    # Minimum volatility method

    if(bsize=="auto")
    {
      sb <- matrix(NA,length(block),3)
      for(cc in 1:length(block))
      {
        id <- max((cc-3),1):min((cc+3),length(block))
        sb[cc,] <- apply(qb[id,],2,sd)
      }
      for(L in 1:3) quan[L] <- qb[which.min(sb[,L]),L]
    } else
    {
      quan <- qb
    }

    ret <- matrix(stat>quan,1,3)
    colnames(ret) <- c("10%","5%","1%")
    if(bsize!="auto") ret<-list(rejection=ret, pvalue=c(pvalue))
    return(ret)
}

data        = aperm(temp22.sm10, c(2,3,1))
ave         = apply(data, c(1,2), mean)
aveArr      = array(ave, dim = dim(data))
data.demean = data - aveArr

p           = dim(data.train)[1]
t           = dim(data.train)[2]
n.train     = dim(data.train)[3]
k0          = 2
ncol.Q      = 12
h           = 1/12
set.seed(1909)
Qmat        = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
mMatArr.st  = lag_cov_prod_int_scale(data.train, h, k0, Qmat)
mMat       = apply(mMatArr.st, c(1,2), sum)
mEig        = eigen(mMat)
val         = mEig$values

ratio       = val[2:p] / val[1:p-1]

pv.mat = matrix(NA, 5, p)
for(rhat in 1:5){
  hatA       = mEig$vectors[, 1:rhat, drop=F]
  common.st  = array(NA, dim = c(rhat, t, n.train))
  residue    = array(NA, dim = dim(data.train))
  idm        = diag(1, nrow = p, ncol = p)
  for (i in 1:(t)) {
    common.st[,i,] = t(hatA) %*% data.train[,i,]
    residue[,i,]   = (idm - hatA %*% t(hatA)) %*% data.train[,i,]
  }
  
  pv = rep(NA, p)
  for(i in 1:p){
    res.i = WN.test(residue[i,,], bsize=1, nbasis = 10, M=200)
    pv[i] = res.i$pvalue
  }
  pv.mat[rhat,] = pv
}

pv.mat
pv.mat.adj = pv.mat
for(i in 1:5){
  pv.mat.adj[i,] = p.adjust(pv.mat[i,], 'BH')
}
pdf(paste('wn_pv_adj_folder/seed_', seed, '.pdf', sep = ''),
    width = cm(4), height = cm(2.5))
matplot(t(pv.mat.adj)[,1:3], type = "b", pch=1:3, col = c('blue', 'red', 'black'), lty=1:3, lwd=2,
        xlab = 'Dimension', ylab = 'P-value', main='Adjusted P-values (BH) of White Noise Tesing',
        ylim=c(0,0.7))
legend("topleft", legend = paste('r =', 1:3), col = c('blue', 'red', 'black'), pch=1:3, lty=1:3,ncol=3)
abline(h = 0.05, lty=4)
abline(h = 0.1, lty=4)
dev.off()
print(paste("FDR=0.05, r=1:",sum(pv.mat.adj[1,]<0.05),'r=3:',sum(pv.mat.adj[3,]<0.05)))
print(paste("FDR=0.1, r=1:",sum(pv.mat.adj[1,]<0.1),'r=3:',sum(pv.mat.adj[3,]<0.1)))