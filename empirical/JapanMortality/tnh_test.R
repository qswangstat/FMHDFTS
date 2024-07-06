begin = Sys.time()
C.seq = c(0.01)
source('mortality_func.R')
res = tnh_select_k(data.train, tgrid, C.seq, J=10, P=1, kmax=10)
end = Sys.time()
elapse = end - begin
elapse

res_final[[1]]$fm[[i]]$pred_f[,,lag]

res_final[[1]]$tnh[[i]]$pred[,,lag]

res_final[[1]]$tnh[[i]]$idio[,,lag]

max(abs(res_final[[1]]$tnh[[i]]$pred_wo_idio[,,lag] - res_final[[1]]$tnh[[i]]$pred[,,lag]))
