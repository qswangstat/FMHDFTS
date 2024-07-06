sp = 23
sex = c('male', 'female', 'total')
for(id in 1:3){
  
  #filename = paste('loading_', sex[id], '_new.RData', sep = '')
  filename = paste('loading_', sex[id],'_sp_',sp, '.RData', sep = '')
  
  load(filename)
  library(NipponMap)
  library(plotfunctions)
  
  #colfunc = colorRampPalette(c("royalblue", "springgreen", "", "yellow", "red"))
  colfunc = colorRampPalette(rev(brewer.pal(11,'Spectral')))
  num_col = 20
  color = colfunc(num_col)
  
  
  max_v = max(abs(max(hatA, hatA.sp, varimax(hatA)$loadings)), 
              abs(min(hatA, hatA.sp, varimax(hatA)$loadings))) 
  min_v = -max_v
  # 
  unit = (max_v - min_v) / num_col
  
  library(NipponMap)
  library(plotfunctions)
  
  
  par(mfrow = c(1,3))
  
  method = c('original', 'varimax', 'sparse')
  plot_load = function(loadMat, type = 1){
    #unit = (max(loadMat) - min(loadMat)) / num_col
    
    for (i in 1:ncol(loadMat)) {
      load  = loadMat[,i]
      idx = ceiling((load - min_v) / unit + 1e-5)
      idx[idx > num_col] = num_col
      # unit = (max(load) - min(load)) / length(load)
      # idx = ceiling((load - min(load)) / unit + 1e-5)
      # idx[idx > length(load)] = length(load)
      pdf(paste("loading_", sex[id], "_sp_", sp, '_', 
                method[type], '_factor_', i, '.pdf', sep = ''), 
          width = cm(2), height = cm(2))
      
      JapanPrefMap(col = color[idx], border = gray(.8), axes = FALSE)
      gradientLegend(
        valRange = c(min_v, -min_v),
        color = color,
        pos = c(0.7,0.1,0.735,0.5),
        fit.margin = TRUE,
        n.seg = 7,
        dec = 2,
      )
      dev.off()
    }
  }
  for (type in 3) {
    if(type == 1) loadMat = hatA
    if(type == 3) loadMat = hatA.sp
    if(type == 2) loadMat = varimax(hatA)$loadings
    plot_load(loadMat, type)
  }
  
}