library(foreach)
#####################################################################
plot_prob = function(p, type, if_weak){
  linet = c(2,4,1)
  color = c('red', 'cyan', 'black')
  filepath = paste("output/signal_strength_type_", type, "_p_", p, ".RData", sep = "")
  if(if_weak == 1){
    filepath = paste("output/weak_signal_strength_type_", type, "_p_", p, ".RData", sep = "")
  } 
  load(filepath)
  n_fac = ifelse(type==3, r+1, r)
  prob = foreach(i = 1:length(allres), .combine = 'rbind') %do% {
    apply(allres[[i]][,c(1,3,5)]==n_fac, 2, mean)
  }
  id = ifelse(if_weak==1,  max(which(delta_seq <= 1.25)), max(which(delta_seq <= 0.5)))
  if(type != 3)
    plot(delta_seq[1:id], delta_seq[1:id], ylim = c(0,1), 'n', xlab = expression(Signal ~ strength ~ kappa[0]),
         ylab = 'Proportion', cex.lab=1.5, cex.axis=1.25)
  else
    plot(delta_seq[1:id], delta_seq[1:id], ylim = c(0,1), 'n', xlab = expression(Signal ~ strength ~ kappa[1]),
         ylab = 'Proportion', cex.lab=1.5, cex.axis=1.25)
  for (j in 1:3) {
    lines(delta_seq[1:id], prob[1:id, j], lty = linet[j], lwd = 3, col = color[j])
  }
}

for(type in 1:5){
    p_seq = c(50,100,200)
    w_seq = c(0,1)
    if(type == 3) {w_seq=c(0)}
    if(type == 5) {p_seq = c(100) w_seq=c(0)}
    for(p in p_seq){
        for(if_weak in w_seq){
            figpath = paste("plots/figure_type_", i, "_p_", p, '.pdf', sep = '')
            if(if_weak == 1) figpath = paste("plots/figure_weak_type_", i, "_p_", p, '.pdf', sep = '')
            pdf(figpath, width = cm(2), height = cm(2))
            par(mai = c(2, 2, 1, 1)*0.5)
            plot_prob(p, type, if_weak)
            dev.off()
        }
    }
}

plot_dist = function(p, type, if_weak){
  linet = c(2,4,1)
  color = c('red', 'cyan', 'black')
  filepath = paste("output/signal_strength_type_", type, "_p_", p, ".RData", sep = "")
  if(if_weak == 1){
    filepath = paste("output/weak_signal_strength_type_", type, "_p_", p, ".RData", sep = "")
  } 
  load(filepath)
  dist= foreach(i = 1:length(allres), .combine = 'rbind') %do% {
    apply(allres[[i]][1:100,c(2,4,6)], 2, mean)
  } 
  if(if_weak == 1){
    delta_seq = delta_seq[-(1:2)]
    dist = dist[-(1:2),]
  } 
  if(type == 5)
    plot(delta_seq, delta_seq, ylim = c(0, max(dist)), 'n', xlab = expression(Signal ~ strength ~ kappa[0]),
         ylab = 'Distance', cex.lab=1.5, cex.axis=1.25)
  if(type == 4)
    plot(delta_seq, delta_seq, ylim = c(0, 0.4), 'n', xlab = expression(Signal ~ strength ~ kappa[0]),
         ylab = 'Distance', cex.lab=1.5, cex.axis=1.25)
  if(type == 3)
    plot(delta_seq, delta_seq, ylim = c(0, 0.1), 'n', xlab = expression(Signal ~ strength ~ kappa[1]),
         ylab = 'Distance', cex.lab=1.5, cex.axis=1.25)
  for (j in 1:3) {
    lines(delta_seq, dist[, j], lty = linet[j], lwd = 3, col = color[j])
  }
}

for(type in 1:5){
    p_seq = c(50,100,200)
    w_seq = c(0,1)
    if(type == 3) {w_seq=c(0)}
    if(type == 5) {p_seq = c(100) w_seq=c(0)}
    for(p in p_seq){
        for(if_weak in w_seq){
            figpath = paste("plots/figure_true_r_type_", i, "_p_", p, '.pdf', sep = '')
            if(if_weak == 1) figpath = paste("plots/figure_weak_true_r_type_", i, "_p_", p, '.pdf', sep = '')
            pdf(figpath, width = cm(2), height = cm(2))
            par(mai = c(2, 2, 1, 1)*0.5)
            plot_dist(p, type, if_weak)
            dev.off()
        }
    }
}