ahead.max = 3

filename = c('mortality_male3.RData','mortality_female3.RData')
all.error = matrix(NA, nrow = ahead.max, ncol = length(filename)*6*2)

for (i in 1:2) {
  load(filename[i])
  mafe = msfe = matrix(NA, nrow = ahead.max, ncol = 6,
                       dimnames = list(1:ahead.max, 
                                       c('FM-1', 'FM-2', 'FM-fix2', 'FM', 'GSY', 'CW')))
  rownames(mafe) = paste("Lag", 1:ahead.max)
  rownames(msfe) = paste("Lag", 1:ahead.max)
  
  len = length(res_final)
  for(lag in 1:ahead.max){
    data.pred4 = sapply(res_final[1:(len-lag+1)], function(x) x$fm$fm_first[,,lag], simplify = "array")
    data.pred5 = sapply(res_final[1:(len-lag+1)], function(x) x$fm$fm_second[,,lag], simplify = "array")
    data.pred6 = sapply(res_final[1:(len-lag+1)], function(x) x$fm2$fm[,,lag], simplify = "array")
    
    data.pred1 = sapply(res_final[1:(len-lag+1)], function(x) x$fm$fm[,,lag], simplify = "array")
    data.pred2 = sapply(res_final[1:(len-lag+1)], function(x) x$gsy[,,lag], simplify = "array")
    data.pred3 = sapply(res_final[1:(len-lag+1)], function(x) x$cw[,,lag], simplify = "array")
    
    diff4      = abs(data.pred4 - data.test[,,lag:len])
    diff5      = abs(data.pred5 - data.test[,,lag:len])
    diff6      = abs(data.pred6 - data.test[,,lag:len])
    diff1      = abs(data.pred1 - data.test[,,lag:len])
    diff2      = abs(data.pred2 - data.test[,,lag:len])
    diff3      = abs(data.pred3 - data.test[,,lag:len])
    
    mafe[lag, ] = c(mean(diff4), mean(diff5), mean(diff6), mean(diff1), mean(diff2), mean(diff3))
    msfe[lag, ] = c(mean(diff4^2), mean(diff5^2), mean(diff6^2), mean(diff1^2), mean(diff2^2), mean(diff3^2))
  }
  all.error[,(12*i-11):(12*i-6)] = mafe*10
  all.error[,(12*i-5):(12*i)] = msfe*10
}
kableExtra::kable(all.error, format = 'markdown', digits = 3,
                  col.names = rep(c('FM-1', 'FM-2', 'FM-fix2', 'FM', 'GSY', 'CW'), 
                                  2*length(filename)))
##############################
rm(list = ls()) 
ahead.max = 10

filename = c('revision/mortality_male_final_new_demean_ahead_10.RData',
             'revision/mortality_female_final_new_demean_ahead_10.RData')#, 'mortality_total_final_new.RData')
all.error = matrix(NA, nrow = ahead.max, ncol = length(filename)*6*2)

for (j in 1:length(filename)) {
  load(filename[j])
  mafe = msfe = matrix(NA, nrow = ahead.max, ncol = 6,
                       dimnames = list(1:ahead.max, 
                                       c('GSY','CW', paste('FM-', 1:4))))
  rownames(mafe) = paste("Lag", 1:ahead.max)
  rownames(msfe) = paste("Lag", 1:ahead.max)
  len = length(res_final)
  for(lag in 1:ahead.max){
    data.pred.cw = sapply(res_final[1:(len-lag+1)], function(x) x$cw[,,lag], simplify = "array")
    diff.cw      = abs(data.pred.cw - data.test[,,lag:len])
    data.pred.gsy = sapply(res_final[1:(len-lag+1)], function(x) x$gsy[,,lag], simplify = "array")
    diff.gsy      = abs(data.pred.gsy - data.test[,,lag:len])
    mafe_vec = c(mean(diff.gsy), mean(diff.cw))
    msfe_vec = c(mean(diff.gsy^2), mean(diff.cw^2))
    for(i in 1:4){
      data.pred.i = sapply(res_final[1:(len-lag+1)], 
                           function(x) x$fm[[i]]$pred_f[,,lag], simplify = "array") 
      diff.i      = abs(data.pred.i - data.test[,,lag:len])
      mafe_vec = c(mafe_vec, mean(diff.i))
      msfe_vec = c(msfe_vec, mean(diff.i^2))
    }
    
    
    mafe[lag, ] = mafe_vec
    msfe[lag, ] = msfe_vec
  }
  all.error[,(12*j-11):(12*j-6)] = mafe*10
  all.error[,(12*j-5):(12*j)] = msfe*10
}
library(kableExtra)
kable(all.error, format = 'html', digits = 2,
      col.names = rep(c('GSY', 'CW', paste('FM-', 1:4)), 
                      2*length(filename))) %>% kable_styling(font_size = 8)
for(j in 1:1){
  print(kableExtra::kable(all.error[,(12*j-11):(12*j)], format = 'latex', digits = 2,
                          col.names = rep(c('GSY', 'CW', paste('FMF-', 1:4)), 2)))
}
for(lag in 1:ahead.max){
  data.pred.cw = sapply(res_final[1:(len-lag+1)], function(x) x$cw[,,lag], simplify = "array")
  diff.cw      = abs(data.pred.cw - data.test[,,lag:len])
  cat(mean(diff.cw), mean(diff.cw^2), '\n')
}

########################################
######### sparse case ################
##################################
options(digits = 4)

ahead.max = 10
len = length(res_final)
p = dim(data.test)[1]
t = dim(data.test)[2]
for (ahead in 1:ahead.max) {
  len = length(res_final)
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
    #cat(best.i, best.j, '\n')
    data.pred[,,i] = res.i[[best.i]][[best.j]]$pred_f[,,ahead]
    #sum = sum + min(loss) * p * t
  }
  diff = data.pred - data.test[,,ahead:len]
  cat(c(mean(abs(diff)), mean(diff^2))*10, '\n')
  #cat(sum / (p*t*(len-ahead+1)))
}

######### TNH method
rm(list = ls()) 

ahead.max = 10

filename = c('revision/mortality_male_tnh_k=1_demean_ahead_10.RData',
             'revision/mortality_female_tnh_k=1_demean_ahead_10.RData')
all.error = matrix(NA, nrow = ahead.max, ncol = length(filename)*1*2)

for (j in 1:length(filename)) {
  load(filename[j])
  mafe = msfe = matrix(NA, nrow = ahead.max, ncol = 1,
                       dimnames = list(1:ahead.max, 
                                       c(paste('TNH-', 10))))
  rownames(mafe) = paste("Lag", 1:ahead.max)
  rownames(msfe) = paste("Lag", 1:ahead.max)
  len = length(res_final)
 
  for(lag in 1:ahead.max){
    mafe_vec = c()
    msfe_vec = c()
    for(i in 1:1){
      data.pred.i = sapply(res_final[1:(len-lag+1)], 
                           function(x) x$tnh[[i]]$pred[,,lag], simplify = "array") 
      diff.i      = abs(data.pred.i - data.test[,,lag:len])
      mafe_vec = c(mafe_vec, mean(diff.i))
      msfe_vec = c(msfe_vec, mean(diff.i^2))
    }
    mafe[lag, ] = mafe_vec
    msfe[lag, ] = msfe_vec
  }
  #all.error[,(10*j-9):(10*j-5)] = mafe*10
  all.error[,2*j-1] = mafe*10
  all.error[,2*j] = msfe*10
  #all.error[,(10*j-4):(10*j)] = msfe*10
}
library(kableExtra)
kable(all.error, format = 'html', digits = 2,
      col.names = c('TNH male MAPE', 'TNH male MSPE',
                    'TNH female MAPE', 'TNH female MSPE')) %>% kable_styling(font_size = 8)


