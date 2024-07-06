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
