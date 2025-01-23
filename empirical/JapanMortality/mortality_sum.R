result = c()
for(sex in 1:1){
  ahead.max = 3
  rhat = 2
  result_mafe = result_msfe = matrix(NA, nrow = ahead.max, ncol = 5,
                                     dimnames = list(1:ahead.max, 
                                                     c('FFM', 'SFFM', 'TNH', 'GSY', 'ANH')))
  rownames(result_mafe) = rownames(result_msfe) = paste("h=", 1:ahead.max, sep = '')
  if (sex == 1) {
    load('./mortality_male_final.RData')
  } else if (sex == 2) {
    load('./mortality_female_final.RData')
  }
  len = length(res_final)
  for(lag in 1:ahead.max){
    data.pred.cw = sapply(res_final[1:(len-lag+1)], function(x) x$cw[,,lag], simplify = "array")
    diff.cw      = abs(data.pred.cw - data.test[,,lag:len])
    data.pred.gsy = sapply(res_final[1:(len-lag+1)], function(x) x$gsy[,,lag], simplify = "array")
    diff.gsy      = abs(data.pred.gsy - data.test[,,lag:len])
    
    data.pred.i = sapply(res_final[1:(len-lag+1)], 
                         function(x) x$fm[[rhat]]$pred_f[,,lag], simplify = "array") 
    diff.i      = abs(data.pred.i - data.test[,,lag:len])
    
    result_mafe[lag,c('FFM', 'GSY', 'ANH')] = 10*c(mean(diff.i), mean(diff.gsy), mean(diff.cw))
    result_msfe[lag,c('FFM', 'GSY', 'ANH')] = 10*c(mean(diff.i^2), mean(diff.gsy^2), mean(diff.cw^2))
    
  }
  if (sex == 1) {
    load('./mortality_male_sparse.RData')
  } else if (sex == 2) {
    load('./mortality_female_sparse.RData')
  }
  len = length(res_final)
  p = dim(data.test)[1]
  t = dim(data.test)[2]
  msfe_vec= mafe_vec = rep(NA, ahead.max)
  for (ahead in 1:ahead.max) {
    data.pred = array(NA, dim = c(p, t, (len-ahead+1)))
    sum = 0
    for(i in 1:(len-ahead+1)){
      true.i = data.test[ , , i + ahead - 1]
      loss = c()
      res.i = res_final[[i]]
      for (j in 1:length(res.i)) {
        tmp.ij = res.i[[j]]
        loss.j = sapply(tmp.ij, function(x){ mean((x$pred_f[,,ahead] - true.i)^2)}, simplify = 'vector')
        loss = c(loss, loss.j)
      }
      idx = which.min(loss)
      ll = length(tmp.ij)
      best.i = ceiling(idx / ll - 1e-5)
      best.j = idx - ll*best.i + ll
      data.pred[,,i] = res.i[[best.i]][[best.j]]$pred_f[,,ahead]
    }
    diff = data.pred - data.test[,,ahead:len]
    mafe_vec[ahead] = mean(abs(diff))*10
    msfe_vec[ahead] = mean(diff^2)*10
  }
  result_mafe[,'SFFM'] = mafe_vec
  result_msfe[,'SFFM'] = msfe_vec
  if (sex == 1) {
    load('./mortality_male_tnh.RData')
  } else if (sex == 2) {
    load('./mortality_female_tnh.RData')
  }
  len = length(res_final)
  mafe_vec = msfe_vec = rep(NA, ahead.max)
  for(lag in 1:ahead.max){
    
    data.pred.i = sapply(res_final[1:(len-lag+1)], 
                         function(x) x$tnh[[1]]$pred[,,lag], simplify = "array") 
    diff.i      = abs(data.pred.i - data.test[,,lag:len])
    mafe_vec[lag] = mean(diff.i) * 10
    msfe_vec[lag] = mean(diff.i^2) * 10
  }
  result_mafe[, 'TNH'] = mafe_vec
  result_msfe[, 'TNH'] = msfe_vec
  result = rbind(result, cbind(result_mafe,result_msfe))
}
