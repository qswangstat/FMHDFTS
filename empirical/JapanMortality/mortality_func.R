##' @param fts:   functional time series, stored in matrix of t*n 
##' @param tgrid: time grids
##' @param cvp:   cumulative percentage of variance, default 0.90
##' @details first perform functional FPCA, d FPCA scores are retained, 
##'          d is determined by cvp
fts_pred = function(fts, tgrid, ahead.max, cvp = 0.9){
  t          = dim(fts)[1]
  n.train    = dim(fts)[2]
  input      = MakeFPCAInputs(IDs = rep(1:n.train, each=t), tVec=rep(tgrid, n.train), fts)
  fpca       = FPCA(input$Ly, input$Lt)
  fpc.funct  = fpca$phi
  d          = which(fpca$cumFVE > cvp)[1]
  fpc.score  = fpca$xiEst[, 1:d, drop=F]
  fpc.mu     = fpca$mu
  score.fit = matrix(NA, nrow = ahead.max, ncol = d)
  for (idx in 1:d) {
    score.d = fpc.score[,idx]
    model.d = auto.arima(score.d, ic='aic')
    forc.d = forecast(model.d, h = ahead.max)
    score.fit[, idx] = forc.d$mean
  }
  
  # var.fit         = vars::VAR(ts(fpc.score), lag.max = 2, ic = "AIC")
  # var.pred        = predict(var.fit, n.ahead = ahead.max)
  # score.fit       = matrix(sapply(var.pred$fcst, function(x) x[, "fcst"]), nrow = ahead.max)
  pred = fpc.funct[, 1:d] %*% t(score.fit) + fpc.mu
  
  return(pred)
}
# component-wise prediction
cw_pred = function(data.train, tgrid, ahead.max, cvp=0.9){
  p          = dim(data.train)[1]
  t          = dim(data.train)[2]
  data.pred = array(NA, dim = c(p, t, ahead.max))
  for (i in 1:p) {
    data.pred[i,,] = fts_pred(data.train[i,,], tgrid, ahead.max, cvp)
  }
  return(data.pred)
}

gsy_pred = function(data.train, ahead.max, order = 2, r = 2){
  p          = dim(data.train)[1]
  data.list  = lapply(1:p, function(x) data.train[x, ,])
  hd.model   = hdfpca(data.list, order = order, r = r)
  #### level & B no use!
  hd.fore    = forecast(hd.model, h = ahead.max, level = 1, B = 1)
  data.pred  = abind::abind(hd.fore$forecast, along = 3)
  data.pred  = aperm(data.pred, c(3,1,2))
  return(data.pred)
}
## predict r factors
pred_r_fac = function(data.train, ave, tgrid, mEig, p, t, n.train, rhat){
  if(is.list(mEig)) hatA       = mEig$vectors[, 1:rhat, drop=F]
  else hatA = mEig
  common.st  = array(NA, dim = c(rhat, t, n.train))
  residue    = array(NA, dim = dim(data.train))
  idm = diag(1, nrow = p, ncol = p)
  for (i in 1:(t)) {
    common.st[,i,] = t(hatA) %*% data.train[,i,]
    residue[,i,]   = (idm - hatA %*% t(hatA)) %*% data.train[,i,]
  }
  
  factor.est = array(NA, dim = c(rhat, t, ahead.max))
  residu.est = array(NA, dim = c(p, t, ahead.max))
  for (i in 1:rhat) factor.est[i,,] = fts_pred(common.st[i,,], tgrid, ahead.max)
  for (i in 1:p) residu.est[i,,] = fts_pred(residue[i,,], tgrid, ahead.max)
  
  pred.fac = pred.fac.noise = array(NA, dim = c(p, t, ahead.max))
  for (i in 1:(t)) {
    pred.fac[,i,]  = hatA %*% matrix(factor.est[,i,], nrow=rhat)# + residu.est[,i,]
    pred.fac.noise[,i,] = hatA %*% matrix(factor.est[,i,], nrow=rhat) + residu.est[,i,]
  }
  
  pred.f = pred.fac + array(ave, dim = dim(pred.fac))
  pred.fac.noise = pred.fac.noise + array(ave, dim = dim(pred.fac.noise))
  
  return(list(pred_f = pred.f, pred_fn = pred.fac.noise, loading = hatA))
}

fm_pred = function(data.train.org, h, k0, ncol.Q, tgrid, ahead.max = 3, rhat.seq = NULL){
  ave        = apply(data.train.org, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train.org))
  data.train = data.train.org - aveArr
  
  p          = dim(data.train)[1]
  t          = dim(data.train)[2]
  n.train    = dim(data.train)[3]
  
  Qmat       = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
  mMatArr.st = lag_cov_prod_int_scale(data.train, h, k0, Qmat)
  mMat       = apply(mMatArr.st, c(1,2), sum)
  mEig       = eigen(mMat)
  val        = mEig$values
  
  ratio      = val[2:p] / val[1:p-1]
  R = p/2
  #if(is.null(rhat)) rhat = which.min(ratio[1:R])
  #hatA       = mEig$vectors[, 1:rhat, drop=F]
  res = list()
  for (rhat in rhat.seq) {
    tmp = pred_r_fac(data.train, ave, tgrid, mEig, p, t, n.train, rhat)
    res = c(res, list(tmp))
  }
  
  return(res)
}
# fm_2step_pred = function(data.train, h, k0, ncol.Q, tgrid, ahead.max){
#   p          = dim(data.train)[1]
#   t          = dim(data.train)[2]
#   n.train    = dim(data.train)[3]
#   
#   Qmat       = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
#   mMatArr.st = lag_cov_prod_int_scale(data.train, h, k0, Qmat)
#   mMat       = apply(mMatArr.st, c(1,2), sum)
#   mEig       = eigen(mMat)
#   val        = mEig$values
#   
#   ratio      = val[2:p] / val[1:p-1]
#   R = p/2
#   rhat       = which.min(ratio[1:R])
#   hatA       = mEig$vectors[, 1:rhat, drop=F]
#   
#   
#   ### obtain factor/residue series
#   common.st  = array(NA, dim = c(rhat, t, n.train))
#   residue    = array(NA, dim = dim(data.train))
#   idm = diag(1, nrow = p, ncol = p)
#   for (i in 1:(t)) {
#     common.st[,i,] = t(hatA) %*% data.train[,i,]
#     residue[,i,]   = (idm - hatA %*% t(hatA)) %*% data.train[,i,]
#   }
#   ### forecast first factor
#   factor.est = array(NA, dim = c(rhat, t, ahead.max))
#   residu.est = array(NA, dim = c(p, t, ahead.max))
#   
#   for (i in 1:rhat) factor.est[i,,] = fts_pred(common.st[i,,], tgrid, ahead.max)
#   for (i in 1:p) residu.est[i,,] = fts_pred(residue[i,,], tgrid, ahead.max)
#   
#   #### second step
#   Qmat       = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
#   mMatArr.st.res = lag_cov_prod_int_scale(residue, h, k0, Qmat)
#   mMat.res       = apply(mMatArr.st.res, c(1,2), sum)
#   mEig.res       = eigen(mMat.res)
#   val.res        = mEig.res$values
#   
#   ratio.res      = val.res[2:p] / val.res[1:p-1]
#   rhat.res       = which.min(ratio.res[1:R])
#   hatA.res       = mEig.res$vectors[, 1:rhat.res, drop=F]
#   ### obtain factor/residue series
#   common.st.res  = array(NA, dim = c(rhat.res, t, n.train))
#   for (i in 1:(t)) {
#     common.st.res[,i,] = t(hatA.res) %*% residue[,i,]
#   }
#   
#   factor.est2 = array(NA, dim = c(rhat.res, t, ahead.max))
# 
#   for (i in 1:rhat.res) factor.est2[i,,] = fts_pred(common.st.res[i,,], tgrid, ahead.max)
#   ### combine factor and residue
#   data.pred.first = data.pred.second = data.pred.full = array(NA, dim = c(p, t, ahead.max))
#   for (i in 1:(t)){
#     data.pred.first[,i,]  = hatA %*% matrix(factor.est[,i,], nrow=rhat)
#     data.pred.second[,i,]  = hatA %*% matrix(factor.est[,i,], nrow=rhat) + 
#       hatA.res %*% matrix(factor.est2[,i,], nrow=rhat.res)
#     data.pred.full[,i,]  = hatA %*% matrix(factor.est[,i,], nrow=rhat) + residu.est[,i,]
#   }
#   return(list(num = c(rhat, rhat.res), fm = data.pred.full,
#               fm_first = data.pred.first, fm_second = data.pred.second,
#               loading = cbind(hatA, hatA.res)))
# }
pred = function(data.train, h, k0, ncol.Q, tgrid, ahead.max, order = 2, r = 2){
  cw  = cw_pred(data.train, tgrid, ahead.max)
  #cat('end', '\n')
  gsy = gsy_pred(data.train, ahead.max, order, r)
  #cat('end', '\n')
  fm  = fm_pred(data.train, h, k0, ncol.Q, tgrid, ahead.max, rhat.seq = 1:4)
  
  
  #fm2 = fm_pred(data.train, h, k0, ncol.Q, tgrid, ahead.max)
  #cat('end', '\n')
  print('success')
  return(list(cw = cw, gsy = gsy, fm = fm))
}

fm_sparse_pred = function(data.train.org, h, k0, ncol.Q, tgrid, ahead.max = 3, rhat.seq = NULL){
  ave        = apply(data.train.org, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train.org))
  data.train = data.train.org - aveArr
  
  p          = dim(data.train)[1]
  t          = dim(data.train)[2]
  n.train    = dim(data.train)[3]
  rhat = rhat.seq
  
  hs_t = FMfts::lag_cov_hs(data.train, h, k0)
  res = list()
  #seq(0.035, 0.06, 0.005)
  #th_seq = seq(ceiling(min(hs_t)/ 0.005)*0.005, floor(max(hs_t) / 0.005 - 1)*0.005, 0.005)
  #th_seq = seq(min(hs_t), max(hs_t), length.out = 10)[1:5]
  th_seq = 0
  for(th in th_seq){
    res.th = list()
    thMat = 1*(hs_t >= rep(th, each = p^2))
    
    Qmat        = matrix(runif(p*ncol.Q, -1, 1) * sqrt(3), nrow = p, ncol = ncol.Q) 
    mMatArr.th  = lag_cov_prod_int_scale_th(data.train, thMat, h, k0, Qmat)
    mMat        = apply(mMatArr.th, c(1,2), sum)
    if(sum(mMat==0) > p^2*0.75) {cat('skip\n');next}
    else{
      for(sp in 47){
        hatA = gdefla(mMat, rep(sp, rhat))
        tmp = pred_r_fac(data.train, ave, tgrid, hatA, p, t, n.train, rhat)
        res.th = c(res.th, list(tmp))
        cat(th, sp, '\n')
      }
      res = c(res, list(res.th))
    }
      
  }
  
  
  return(res)
}
trapzRcpp.new = function(y, t){
  return(trapzRcpp(t, y))
}
tnh_select_k_sub = function(data.train.org, tgrid, c, kmax){
  ave        = apply(data.train.org, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train.org))
  data.train = data.train.org - aveArr # y_{it}
  
  ## N = 47, T <= 42, obtain T*T matrix F
  N = dim(data.train)[1]
  tlen = dim(data.train)[2]
  Tt = dim(data.train)[3]
  
  Fmat = matrix(NA, nrow = Tt, ncol = Tt)
  for(i in 1:Tt){
    for(j in 1:i){
      prod.ij = data.train[,, i] * data.train[,, j]
      Fmat[i, j] = sum(apply(prod.ij, 1, trapzRcpp.new, t=tgrid)) / N
    }
  }
  Fmat[upper.tri(Fmat)] = t(Fmat)[upper.tri(Fmat)]
  eig.F = eigen(Fmat)
  
  IC.seq = g.seq = V.seq = rep(NA, kmax)
  for(k in 1:kmax){
    lambda.hat = eig.F$values[1:k]/(sqrt(Tt))
    f.hat = sqrt(Tt) * eig.F$vectors[,1:k,drop=F]
    U.tilde = t(f.hat)
    
    B.tilde = array(NA, dim = c(N, tlen, k))
    for(l in 1:k){
      f.hat.arr = array(rep(f.hat[,l], each = N*tlen), dim = c(N, tlen, Tt))
      B.tilde[,, l] = apply(f.hat.arr * data.train, c(1,2), sum)/(Tt)
    }
    
    ## common component
    common.comp = array(NA, dim = dim(data.train))
    for (tid in 1:Tt) {
      u.t = U.tilde[,tid]
      common.comp[,, tid] = apply(array(rep(u.t, each = N*tlen), 
                                        dim = dim(B.tilde)) * B.tilde, c(1,2), sum)
    }
    ## idiosyncratic component
    idio.comp = data.train - common.comp
    V.value = 0
    for(tid in 1:Tt){
      prod.i = idio.comp[,, tid]^2
      V.value = V.value + sum(apply(prod.i, 1, trapzRcpp.new, t=tgrid))
    }
    gNT = c*k*((N+Tt)/(N*Tt))*log(min(N,Tt))
    
    g.seq[k] = gNT
    V.seq[k] = V.value
    IC.seq[k] = gNT + (V.value / (N*Tt))
  }
  k.best = which.min(IC.seq)[1]
  print(g.seq)
  print(V.seq)
  print(IC.seq)
  cat(k.best, '\n')
  return(k.best)
}
tnh_select_k = function(data.train.org, tgrid, C.seq, J, P, kmax){
  ave        = apply(data.train.org, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train.org))
  data.train = data.train.org - aveArr # y_{it}
  
  ## N = 47, T <= 42, obtain T*T matrix F
  N = dim(data.train)[1]
  tlen = dim(data.train)[2]
  Tt = dim(data.train)[3]
  k.best = array(NA, dim=c(length(C.seq), J, P))
  for (i in 1:length(C.seq)) {
    c = C.seq[i]
    for(j in 1:J){
      Nj = floor(4*N/5+N*(j-1)/45)
      for(perm in 1:P){
        index = sample(1:N, N)
        data.train.shuffle = data.train[index,,]
        data.train.sub = data.train[1:Nj,,]
        k.best[i,j,perm] = tnh_select_k_sub(data.train.sub, tgrid, c, kmax)
        #cat('c =',c,'j =',j,'p =',perm,' finish\n')
      }
    }
    #cat('c=', c, 'end\n')
  }
  return(list(kmat = k.best))
}

tnh_pred = function(data.train.org, tgrid, k = 2, ahead.max = 3){
  ave        = apply(data.train.org, c(1,2), mean)
  aveArr     = array(ave, dim = dim(data.train.org))
  data.train = data.train.org - aveArr # y_{it}
  
  ## N = 47, T <= 42, obtain T*T matrix F
  N = dim(data.train)[1]
  tlen = dim(data.train)[2]
  Tt = dim(data.train)[3]
  
  Fmat = matrix(NA, nrow = Tt, ncol = Tt)
  for(i in 1:Tt){
    for(j in 1:i){
      prod.ij = data.train[,, i] * data.train[,, j]
      
      Fmat[i, j] = sum(apply(prod.ij, 1, trapzRcpp.new, t=tgrid)) / N
    }
  }
  Fmat[upper.tri(Fmat)] = t(Fmat)[upper.tri(Fmat)]
  eig.F = eigen(Fmat)
  lambda.hat = eig.F$values[1:k]/(sqrt(Tt))
  f.hat = sqrt(Tt) * eig.F$vectors[,1:k,drop=F]
  
  U.tilde = t(f.hat)
  B.tilde = array(NA, dim = c(N, tlen, k))
  for(l in 1:k){
    f.hat.arr = array(rep(f.hat[,l], each = N*tlen), dim = c(N, tlen, Tt))
    B.tilde[,, l] = apply(f.hat.arr * data.train, c(1,2), sum)/(Tt)
  }
  
  ## common component
  common.comp = array(NA, dim = dim(data.train))
  for (tid in 1:Tt) {
    u.t = U.tilde[,tid]
    common.comp[,, tid] = apply(array(rep(u.t, each = N*tlen), 
                                      dim = dim(B.tilde)) * B.tilde, c(1,2), sum)
  }
  
  U.pred = matrix(NA, nrow = nrow(U.tilde), ncol = ahead.max)
  for(i in 1:nrow(U.pred)){
    arima.model = auto.arima(U.tilde[i,], ic='bic')
    U.pred[i,] = forecast(arima.model, h=ahead.max)$mean
  }
  
  common.comp.pred = array(NA, dim = c(N, tlen, ahead.max))
  for (tid in 1:ahead.max) {
    u.t = U.pred[,tid]
    common.comp.pred[,, tid] = apply(array(rep(u.t, each = N*tlen), 
                                      dim = dim(B.tilde)) * B.tilde, c(1,2), sum)
  }
  
  ## idiosyncratic component
  idio.comp = data.train - common.comp
  idio.comp.pred = cw_pred(idio.comp, tgrid, ahead.max, cvp = 0.95)
  
  ## final result
  final.pred = array(ave, dim = dim(idio.comp.pred)) + common.comp.pred + idio.comp.pred
  final.pred_wo_idio = array(ave, dim = dim(idio.comp.pred)) + common.comp.pred
  return(list(pred = final.pred, common = common.comp.pred,
              idio = idio.comp.pred, pred_wo_idio = final.pred_wo_idio))
}
tnh = function(data.train.org, tgrid, k.seq = 1:5, ahead.max = 3){
  res = list()
  for (k in k.seq) {
    tmp = tnh_pred(data.train.org, tgrid, k, ahead.max = ahead.max)
    res = c(res, list(tmp))
  }
  return(list(tnh=res))
}
