source("simulation.R")
for(ncol.Q in c(2:5)*4){
    for(k0 in 2:5){
        simulation(p=100, type=1, if_weak=0, ncol.Q=ncol.Q, k0=k0, if_verify=TRUE)
    }
}