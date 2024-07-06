library(foreach)
prob_mat = matrix(NA, ncol = 9,
                  nrow = 4)
l = 1
for(ncol.Q in c(12)){
  for(k0 in 2:5){
    load(paste("Q_k0_robustness/signal_strength_type_1_p_100_q_",ncol.Q,"_k0_",k0, ".RData", sep=''))
    prob = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
      apply(allres[[i]][,c(1,3,5)]==r, 2, mean)
    }
    prob_mat[l,] = prob[,3]
    l = l + 1
  }
}

pdf(paste('robustness_prob_q_', ncol.Q, '.pdf', sep = ''), width = cm(2), height = cm(2))
matplot(delta_seq[1:7], t(prob_mat[,1:7]), type = 'l', lty = c(4,5,2,1), lwd=3,col=4:1,
        xlab = expression(Signal ~ strength ~ kappa[0]),
        ylab = 'Proportion', ylim=c(0,1), cex.lab=1.5)
legend("bottomright", legend = expression(k[0] == 2, k[0] == 3, k[0] == 4, k[0] == 5),
       col = 4:1, lty = c(4,5,2,1), ncol=1, cex=1.5,lwd=3)
dev.off()

dist_mat = matrix(NA, ncol = 9,
                  nrow = 4)
l = 1
for(ncol.Q in c(12)){
  for(k0 in 2:5){
    load(paste("Q_k0_robustness/signal_strength_type_1_p_100_q_",ncol.Q,"_k0_",k0, ".RData", sep=''))
    dist = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
      apply(allres[[i]][,c(2,4,6)], 2, mean)
    }
    dist_mat[l,] = dist[,3]
    l = l + 1
  }
}
pdf(paste('robustness_dist_q_', ncol.Q, '.pdf', sep = ''), width = cm(2), height = cm(2))
matplot(delta_seq[1:7], t(dist_mat[,1:7]), type = 'l', lty = c(4,5,2,1), lwd=3,col=4:1,
        xlab = expression(Signal ~ strength ~ kappa[0]),
        ylab = 'Distance', cex.lab=1.5, ylim=c(0, 0.2))
legend("topright", legend = expression(k[0] == 2, k[0] == 3, k[0] == 4, k[0] == 5),
       col = 4:1, lty = c(4,5,2,1), ncol=1, cex=1.5,lwd=3)
dev.off()

prob_mat = matrix(NA, ncol = 9,
                  nrow = 4)
l = 1
for(ncol.Q in seq(8,20,4)){
  for(k0 in 4){
    load(paste("Q_k0_robustness/signal_strength_type_1_p_100_q_",ncol.Q,"_k0_",k0, ".RData", sep=''))
    prob = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
      apply(allres[[i]][,c(1,3,5)]==r, 2, mean)
    }
    prob_mat[l,] = prob[,3]
    l = l + 1
  }
}
pdf(paste('robustness_prob_k0_', k0, '.pdf', sep = ''), width = cm(2), height = cm(2))
matplot(delta_seq[1:7], t(prob_mat[,1:7]), type = 'l', lty = c(4,5,2,1), lwd=3,col=4:1,
        xlab = expression(Signal ~ strength ~ kappa[0]),
        ylab = 'Proportion', ylim=c(0, 1), cex.lab=1.5)
legend("bottomright",legend = expression(q==8, q==12, q==16, q==20), 
       col = 4:1, lty = c(4,5,2,1), ncol=1, cex=1.5,lwd=3)
dev.off()

dist_mat = matrix(NA, ncol = 9,
                  nrow = 4)
l = 1
for(ncol.Q in seq(8,20,4)){
  for(k0 in 4){
    load(paste("Q_k0_robustness/signal_strength_type_1_p_100_q_",ncol.Q,"_k0_",k0, ".RData", sep=''))
    dist = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
      apply(allres[[i]][,c(2,4,6)], 2, mean)
    }
    dist_mat[l,] = dist[,3]
    l = l + 1
  }
}
pdf(paste('robustness_dist_k0_', k0, '.pdf', sep = ''), width = cm(2), height = cm(2))
matplot(delta_seq[1:7], t(dist_mat[,1:7]), type = 'l', lty = c(4,5,2,1), lwd=3,col=4:1,
        xlab = expression(Signal ~ strength ~ kappa[0]),
        ylab = 'Distance', cex.lab=1.5, ylim=c(0, 0.2))
legend("topright", legend = expression(q==8, q==12, q==16, q==20),
       col = 4:1, lty = c(4,5,2,1), ncol=1, cex=1.5,lwd=3)
dev.off()
