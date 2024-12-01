# FMHDFTS
This is the repo for paper "Factor Modelling for High-dimensional Functional Time Series", and this README file provides procedure for reproducing numerical results presented in the paper.

## Simulation
Before you start, you install the package FMfts by running the following code in R.
```install.packages("youpath/FMHDFTS/FMfts_1.0.0.tar.gz", repos = NULL, type = "source")```
Besides, the dependencies also include the following packages.
```foreach, doRNG, doSNOW.```

To generate simualtion results, you should collect key outputs by running the following code in command line.
```nohup Rscript simulation.R 100 1 0 >./print.log 2>./error.log```

The first parameter (100) is the dimension p of factor model. The second parameter (1) denotes the simulation scenario from 1-5. The last parameter (0) indicates whether the simulation is for a strong signal or a weak signal. By using different combinations of the three parameters, you can gather all the necessary results to reproduce Figure 1/2/5/6/7 of the article.

If you are not familiar with command line, you can uncomment the following code and run it directly in R, but doing so will be more time-consuming and resource-intensive.
```
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
```
To generate Figure 8 or Table 4, you only need to run **additional/robustness_verify.R** and **utils/simulation_sparse.R** respectively.
Once you finish, you should turn to **utils/visualization** where three files contain key code to plot the figures. The code is self-explanatory, so we will not go into the details here.

## Empirical
For UK temperature data, **empirical/UKtemperature/data_collect.R** is code for data collection and **empirical/UKtemperature/plot_loading_uk.R** is for Figure 3.

For Japanese mortality data, **empirical/JapanMortality/plot_loading_japan.R** is for Figure 4/9, while **empirical/JapanMortality/mortality_predict_final.R** contains all code to compare various methods for forecasting, i.e. Table 2.

