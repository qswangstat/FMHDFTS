##########################################################################
######   pass the parameter via command line                        ######
###### nohup Rscript simulation.R 100 1 0 >./print.log 2>./error.log######
###### first param:   the dimension p of factor model               ######
###### second param:  the simulation scenario from 1-5              ######
###### third param:   strong/weak signal                            ######
##########################################################################
#!/usr/bin/env Rscript
args = as.numeric(commandArgs(trailingOnly=TRUE))
p = args[1]
type = args[2]
if_weak = args[3]

source("simulation.R")
simulation(p, type, if_weak)

## The following code collects results needed to generate figure 1/2/5/6/7
## But it is strongly recommended to run each combination of parameters individually as I do above due to computational resource limitation
# for(type in 1:5){
#     p_seq = c(50,100,200)
#     w_seq = c(0,1)
#     if(type == 3) {w_seq=c(0)}
#     if(type == 5) {p_seq = c(100) w_seq=c(0)}
#     for(p in p_seq){
#         for(if_weak in w_seq){
#             simulation(p, type, if_weak)
#         }
#     }
# }