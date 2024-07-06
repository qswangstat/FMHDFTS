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
temp.list =list()
for (i in 1:dim(ind)[1]){
  url.txt <- url(paste0("https://www.metoffice.gov.uk/pub/data/weather/uk/climate/stationdata/",ind[i,1]))
  temp.list[[i]] = read.table(url.txt, header = F, skip = ind[i,2],na.strings = "---", col.names = paste0("V",seq_len(8)), fill = TRUE)[,1:4]         
                 
}

trunc = list()
for (i in 1:length(temp.list)){
  trunc[[i]]= temp.list[[i]][temp.list[[i]][,1]>=1959,]
  trunc[[i]][,3] = as.numeric(gsub("\\*", "", trunc[[i]][,3]))
  trunc[[i]][,4] = as.numeric(gsub("\\*", "", trunc[[i]][,4]))
}

temp22 = array(NA,dim = c(12,22,62))
for (i in 1:22){
  temp22[,i,] = apply(trunc[[i]][1:744,3:4],1,mean)
}

temp22 = aperm(temp22, c(3,2,1))

Y_temp = temp22
Tu = dim(Y_temp)[3]
n =dim(Y_temp)[1]
d = dim(Y_temp)[2]

if (Tu !=1){
  u = (0:(Tu-1)) / Tu # (1:Tu)/Tu #
}
phi <- cbind(
  rep(1, Tu) ,
  sqrt(2) * sin(2 * pi * u),
  sqrt(2) * cos(2 * pi * u),
  sqrt(2) * sin(4 * pi * u),
  sqrt(2) * cos(4 * pi * u),
  sqrt(2) * sin(6 * pi * u),
  sqrt(2) * cos(6 * pi * u),
  sqrt(2) * sin(8 * pi * u),
  sqrt(2) * cos(8 * pi * u),
  sqrt(2) * sin(10 * pi * u)
  # sqrt(2) * cos(10 * pi * u),
  # sqrt(2) * sin(12 * pi * u),
  # sqrt(2) * cos(12 * pi * u),
  # sqrt(2) * sin(14 * pi * u),
  # sqrt(2) * cos(14 * pi * u),
  # sqrt(2) * sin(16 * pi * u),
  # sqrt(2) * cos(16 * pi * u),
  # sqrt(2) * sin(18 * pi * u),
  # sqrt(2) * cos(18 * pi * u),
  # sqrt(2) * sin(20 * pi * u),
  # sqrt(2) * cos(20 * pi * u),
  # sqrt(2) * sin(22 * pi * u),
  # sqrt(2) * cos(22 * pi * u)
  
) 

trunc = 4
Ysmooth = array(dim  = c(n,d,Tu))
for (l in 1:n){
  for (k in 1:d){
    Y = Y_temp[l,k,]
    X = phi[,1:trunc]
    S = X%*%solve(t(X)%*%X)%*%t(X)
    Yest = S%*%Y
    Ysmooth[l,k,] = Yest
  }
  
}
temp22.sm10 = Ysmooth
save(temp22, temp22.sm10, file = "temperature.RData",version = 2)