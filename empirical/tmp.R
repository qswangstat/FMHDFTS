library(RCurl)
p             = 47
f.year        = 1975
l.year        = 2017
n             = l.year - f.year + 1
max.age       = 110
data = data.m = data.f = array(NA, dim = c(p, max.age+1, n))
pf.names = rep(NA, p)
for (i in 1:p) {
  url.i       = getURL(paste("http://www.ipss.go.jp/p-toukei/JMD/", 
                             ifelse(i<10, paste('0', i, sep=''), i), "/STATS/Mx_1x1.txt", sep = ''), 
                       ssl.verifypeer = FALSE)
  data.i      = read.table(textConnection(url.i), header = FALSE, nrows = 1)
  pf.names[i] = data.i$V1
  print(i)
}
pf = pf.names
for(i in 1:p)
  pf[i] = substr(pf.names[i], 4, nchar(pf.names[i])-1)
for (i in 1:p) {
  if(i == 1) 
    pf[i] = str_extract(pf.names[i], '(?<=\\.).*(?=,)')
  else 
    pf[i] = str_extract(pf.names[i], '(?<=\\.).*(?=\\.)')
}
###################################################################
p          = dim(data.train)[1]
t          = dim(data.train)[2]
n.train    = dim(data.train)[3]

#ave        = apply(data.train, c(1,2), mean)
#aveArr     = array(ave, dim = dim(data.train))
#data.train = data.train - aveArr

#mMat       = apply(mMatArr.st, c(1,2), sum)
mEig       = eigen(mMat)
val        = mEig$values

ratio      = val[2:p] / val[1:p-1]
R = p/2
rhat = 3
if(is.null(rhat)) rhat = which.min(ratio[1:R])
hatA       = mEig$vectors[, 1:rhat, drop=F]

common.st  = array(NA, dim = c(rhat, t, n.train))
residue    = array(NA, dim = dim(data.train))
idm = diag(1, nrow = p, ncol = p)
for (i in 1:(t)) {
  common.st[,i,] = t(hatA) %*% data.train[,i,]
  residue[,i,]   = (idm - hatA %*% t(hatA)) %*% data.train[,i,]
}
h       = 1 / t
tgrid   = seq(0, 1, length.out = t)
par(mfrow = c(1,3))
#realt = seq(0,95,1)
#label.id = seq(1,96,10)
realt = seq(1,12,1)
label.id = seq(1,12,1)
for(i in 1:3){
  #pdf(paste("factor_male_", i, '.pdf', sep = ''), width = cm(2), height = cm(2))
  plot(tgrid, tgrid, ylim = c(min(common.st[i,,]), max(common.st[i,,])), 'n',
       xlab = 'Age', ylab = 'Factor', xaxt = "n")
  axis(1, at=tgrid[label.id], labels=realt[label.id])
  
  for(j in 1:n.train){
    lines(tgrid, common.st[i,,j], lwd = 0.5, col = j)
  }
  #dev.off()
}

par(mfrow = c(2,3))
rotA = varimax(hatA)
loadMat = rotA$loadings * (-10)
colfunc = colorRampPalette(c("red","yellow","springgreen","royalblue"))

library(NipponMap)
library(plotfunctions)
for (i in 1:ncol(loadMat)) {
  load  = loadMat[,i]
  grad = seq(min(load), max(load), length.out = length(load))
  unit = (max(load) - min(load)) / 46

  idx = ceiling((load - min(load)) / unit + 1e-5)
  pdf(paste("loading_male_", i, '.pdf', sep = ''), width = cm(2), height = cm(2))
  
  JapanPrefMap(col = colfunc(47)[idx], border = gray(.8), axes = FALSE)
  gradientLegend(
    valRange = load,
    color = colfunc(47),
    pos = c(0.7,0.1,0.75,0.4),
    fit.margin = TRUE,
    n.seg = 5,
    dec = 1
  )
  dev.off()
}

loadMat = matrix(loadMat, nrow = p)
rownames(loadMat) = pf
loadMat = data.frame(loadMat)
kable(loadMat, format = 'latex', digits = 3,
      col.names = paste('Factor', 1:6))


rownames(loadMat) = jpnprefs$prefecture
scaleMat = scale(loadMat)
distance = dist(scaleMat, method = 'euclidean')

hc = hclust(distance, method = "complete" )
clusDendro = as.dendrogram(hc)

plot(clusDendro, horiz = TRUE)

library(XML)

location = c('aberporthdata.txt',7,
             'armaghdata.txt',7,
             'bradforddata.txt',7,
             'cambridgedata.txt',7,
             'durhamdata.txt',7,
             'eastbournedata.txt',7,
             'eskdalemuirdata.txt',7,
             'heathrowdata.txt',7,
             'hurndata.txt',7,
             'lerwickdata.txt',7,
             'leucharsdata.txt',7,
             'newtonriggdata.txt',7,
             'oxforddata.txt',7,
             'paisleydata.txt',7,
             'rossonwyedata.txt',7,
             'shawburydata.txt',7,
             'sheffielddata.txt',7,
             'stornowaydata.txt',7,
             'suttonboningtondata.txt',7,
             'tireedata.txt',7,
             'valleydata.txt',7,
             'waddingtondata.txt',7
)

ind = matrix(location, ncol = 2,byrow = T)

loc_mat = matrix(NA, nrow = dim(ind)[1], ncol = 2)
for (i in 1:dim(ind)[1]){
  url.txt <- url(paste0("https://www.metoffice.gov.uk/pub/data/weather/uk/climate/stationdata/",ind[i,1]))
  tmp = read.table(url.txt, nrows = 1, skip = 1)         
  loc = c(as.numeric(tmp$V5), as.numeric(substr(tmp$V7, 1, nchar(tmp$V7)-1)))
  if(i == 2) loc = c(as.numeric(tmp$V7), as.numeric(substr(tmp$V9, 1, nchar(tmp$V9)-1)))
  loc_mat[i, ] = loc
}

map('worldHires',
    c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
    xlim=c(-11,3), ylim=c(49,60.9), mar = c(3,3,3,3))	

par(mfrow = c(1,3))
rotA = varimax(hatA)
loadMat = rotA$loadings * (10)
colfunc = colorRampPalette(c("red","yellow","springgreen","royalblue"))
min(loadMat)
max(loadMat)
library(plotfunctions)
for (i in 1:ncol(loadMat)) {
  load  = loadMat[,i]
  grad = seq(min(load), max(load), length.out = length(load))
  unit = (max(load) - min(load)) / 21
  
  idx = ceiling((load - min(load)) / unit + 1e-5)
  pdf(paste("loading_UK_", i, '.pdf', sep = ''), width = cm(2), height = cm(2))
  
  map('worldHires',
      c('UK', 'Ireland', 'Isle of Man','Isle of Wight'),
      xlim=c(-11,3), ylim=c(49,60.9), mar = c(3,3,3,3))
  for (l in 1:nrow(loc_mat)) {
    points(loc_mat[l,2], loc_mat[l,1], col = colfunc(22)[idx[l]], pch = 19,
           cex = 2.5)
  }
  gradientLegend(
    valRange = load,
    color = colfunc(47),
    pos = c(0.95,0.1,1.05,0.4),
    fit.margin = TRUE,
    n.seg = 5,
    dec = 1
  )
  dev.off()
}


plot_loading = function(loadMat, th, sp){
  colfunc = colorRampPalette(c("red","yellow","springgreen","royalblue"))
  for (i in 1:ncol(loadMat)) {
    load  = loadMat[,i]
    #unit = (max(load[load!=0]) - min(load[load!=0])) / sum(load != 0)
    #idx = ceiling((load - min(load[load!=0])) / unit + 1e-5)
    #idx[idx < 0] = 0
    #idx = idx + 1
    #pdf(paste("loading_male_", i, '.pdf', sep = ''), width = cm(2), height = cm(2))
    #color = c('black', colfunc(sum(load!=0)))
    unit = (max(load) - min(load)) / length(load)
    idx = ceiling((load - min(load)) / unit + 1e-5)
    idx[idx > length(load)] = length(load)
    color = colfunc(length(load))
    JapanPrefMap(col = color[idx], border = gray(.8), axes = FALSE)
    title(main = paste('Factor', i,'\n','Threshold:',th,'Sparsity:',sp))
    gradientLegend(
      valRange = load,
      color = color,
      pos = c(0.7,0.1,0.75,0.4),
      fit.margin = TRUE,
      n.seg = 5,
      dec = 1
    )
    #dev.off()
  }
  
}


for(th in seq(0.035, 0.06, 0.005)){
  pdf(paste('loading_male_th_', th, '.pdf', sep = ''))
  par(mfrow = c(2,2))
  for(sp in seq(5,45,5)){
    load(paste('loading/th_', th, '_sp_', sp, '.RData', sep = ''))
    plot_loading(hatA, th, sp)
  }
  dev.off()
}

