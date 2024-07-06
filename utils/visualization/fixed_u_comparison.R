load('output/signal_strength_type_1_p_100_fixed_u.RData')
prob = matrix(NA, nrow=length(allres), ncol=sim)
for(i in seq_along(allres)){
  res_i = allres[[i]]
  A = res_i$r_est
  rhat = floor(mean(A))
  
  prob[i,1] = rhat
  
  for(l in 3:length(res_i)){
    A = res_i[[l]]$r_est
    rhat = floor(mean(A))
    prob[i,l-1] = rhat
  }
}
prob1 = apply(prob==4, 1, mean)
for(i in 1:length(allres)){
  
  freq_table <- table(allres[[i]]$r_est)
  max_freq <- max(freq_table) 
  most_frequents <- names(freq_table[freq_table == max_freq])
  print(as.integer(max(most_frequents)))
}

for(i in seq_along(allres)){
  res_i = allres[[i]]
  mat = matrix(0, nrow=p, ncol=p)
  A = res_i$hatA_list
  for (j in 1:dim(A)[3]) {
    B <- A[,,j]
    mat = mat + B %*% t(B)
    #print(sum(diag(B%*%t(B)%*%Q%*%t(Q))))
  }
  mat = mat / dim(A)[3]
  eig = eigen(mat)
  hatA = eig$vectors[,1:r]
  temp = hatA%*%t(hatA)%*%Q%*%t(Q)
  dist[i,1] = sqrt(1 - sum(diag(temp)) / r)

  for(l in 3:length(res_i)){
    mat = matrix(0, nrow=p, ncol=p)
    A = res_i[[l]]$hatA_list
    for (j in 1:dim(A)[3]) {
      B <- A[,,j]
      mat = mat + B %*% t(B)
    }
    mat = mat / dim(A)[3]
    eig = eigen(mat)
    hatA = eig$vectors[,1:r]
    temp = hatA%*%t(hatA)%*%Q%*%t(Q)
    dist[i,l-1] = sqrt(1 - sum(diag(temp)) / r)
  }
}

dist1 = apply(dist, 1, mean)

load("outputu/signal_strength_type_1_p_100.RData")
dist2 = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
  apply(allres[[i]][1:100,c(2,4,6)], 2, mean)
}
dist2 = dist2[,3]

par(mfrow = c(1,3))
pdf("figure_fixed_u.pdf", width = cm(2), height = cm(2))
par(mai = c(2, 2, 1, 1)*0.5)
plot(delta_seq, delta_seq, ylim = c(0, 0.18), 'n', xlab = expression(Signal ~ strength ~ kappa[0]),
     ylab = 'Distance', cex.lab=1.5, cex.axis=1.25)
lines(delta_seq, dist1, lty = 2, lwd = 3, col = 'red')
lines(delta_seq, dist2, lty = 1, lwd = 3, col = 'black')
dev.off()
