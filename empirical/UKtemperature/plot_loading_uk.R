load("UKstation.RData")
load("matM_temp.RData")

colfunc = colorRampPalette(rev(brewer.pal(11,'Spectral')))
num_col = 20
color = colfunc(num_col)

hatA = eigen(mMat)$vectors[,1:3, drop = FALSE]
max_v = max(abs(max(hatA, varimax(hatA)$loadings)), 
            abs(min(hatA, varimax(hatA)$loadings))) 
min_v = -max_v
# 
unit = (max_v - min_v) / num_col
loadMat = hatA
par(mfrow = c(2,3))

method = c('original', 'varimax', 'sparse')
plot_load = function(loadMat, type = 1){
  for (i in 1:ncol(loadMat)) {
    load  = loadMat[,i]
    idx = ceiling((load - min_v) / unit + 1e-5)
    idx[idx > num_col] = num_col
    
    pdf(paste("loading_UK_", method[type], "_factor", i, '.pdf', sep = ''), width = cm(2), height = cm(2))
    
    map('worldHires',
        c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
        xlim=c(-11,3), ylim=c(49,60.9), mar = c(5, 4, 4, 2) + 0.1)
    for (l in 1:nrow(loc_mat)) {
      points(loc_mat[l,2], loc_mat[l,1], col = color[idx[l]], pch = 19,
             cex = 2.5)
    }
    gradientLegend(
      valRange = c(min_v, -min_v),
      color = color,
      pos = c(0.95,0.1,1.05,0.4),
      fit.margin = TRUE,
      n.seg = 7,
      dec = 2,
    )
    dev.off()
  }
}

for (type in 1:2) {
  if(type == 1) loadMat = hatA
  if(type == 3) loadMat = hatA.sp
  if(type == 2) loadMat = varimax(hatA)$loadings
  plot_load(loadMat, type)
}
