library(RCurl)
p             = 47
f.year        = 1975
l.year        = 2017
n             = l.year - f.year + 1
max.age       = 110
data = data.m = data.f = array(NA, dim = c(p, max.age+1, n))
for (i in 1:p) {
  url.i       = getURL(paste("http://www.ipss.go.jp/p-toukei/JMD/", 
                             ifelse(i<10, paste('0', i, sep=''), i), "/STATS/Mx_1x1.txt", sep = ''), 
                       ssl.verifypeer = FALSE)
  data.i      = read.table(textConnection(url.i), skip = 1, header = TRUE, fill = TRUE)
  data[i,,]   = matrix(data.i$Total, nrow = max.age+1)
  data.f[i,,] = matrix(data.i$Female, nrow = max.age+1)
  data.m[i,,] = matrix(data.i$Male, nrow = max.age+1)
}
mode(data)   = 'numeric'
mode(data.f) = 'numeric'
mode(data.m) = 'numeric'

t       = 96
h       = 1 / t
tgrid   = seq(0, 1, length.out = t)
data.sp = data.sp.m = data.sp.f = array(NA, dim = c(p, t, n))
bs.mat  = bs(tgrid, df = 9)

for (i in 1:p) {
  for(j in 1:n){
    resp = log(data.f[i,1:t,j])
    resp[resp == -Inf] = NA
    sp.fit = smooth.spline(tgrid[!is.na(resp)], resp[!is.na(resp)])
    sp.pred = predict(sp.fit, x = tgrid)
    sp.pred$y
    #bs.fit = lm(resp ~ bs.mat-1, na.action = na.omit)
    #bs.pred = bs.mat %*% bs.fit$coefficients
    
    data.sp.f[i,,j] = sp.pred$y
  }
}
for (i in 1:p) {
  for(j in 1:n){
    resp = log(data[i,1:t,j])
    resp[resp == -Inf] = NA
    sp.fit = smooth.spline(tgrid[!is.na(resp)], resp[!is.na(resp)])
    sp.pred = predict(sp.fit, x = tgrid)
    sp.pred$y
    #bs.fit = lm(resp ~ bs.mat-1, na.action = na.omit)
    #bs.pred = bs.mat %*% bs.fit$coefficients
    
    data.sp[i,,j] = sp.pred$y
  }
}
n.train = 26
n0 = 1
data.train = data.sp[,,n0:n.train]
data.test = data.sp[,,-(n0:n.train)]
save.image('mortality.RData')
