rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 0
N_node = 100
spks_time_mlist = matrix(list(),N_node,1)

t_vec = res$data_param$t_vec

network_list = do.call(generate_data, args = res$data_param)
network_list = generate_data(SEED=38,
                             N_node=100,
                             N_spks = 1000,
                             conn_patt_sep = 1,
                             time_shift_rad = 0.05 )

spks_time_mlist = network_list$spks_time_mlist[1:100,1,drop=FALSE]
stim_onset_vec = network_list$stim_onset_vec[1:100]

center_density_array_true = network_list$center_density_array_true
mem_true_vec = network_list$mem_true_vec
clus_true_list = network_list$clus_true_list
v_true_list = network_list$v_vec_list

t_vec = network_list$t_vec

plot(t_vec,center_density_array_true[1,1,]+center_density_array_true[1,2,],type='l')
lines(density(spks_time_mlist[1,1][[1]],bw=0.015),col=2)

plot(t_vec,center_density_array_true[2,1,]+center_density_array_true[2,2,],type='l')




N_clus = 2
res = get_init(spks_time_mlist = spks_time_mlist, 
               stim_onset_vec = stim_onset_vec,
               N_clus = N_clus,
               N_component = 2,
               v0 = 0.5, v1 = 0.5,
               t_vec = t_vec, 
               fix_timeshift = FALSE)
clusters_list_init = res$clusters_list
v_vec_list_init = res$v_vec_list

print(clusters_list_init)
plot(v_true_list[[1]][1:100], v_vec_list_init[[1]][1:100]); abline(a=0,b=1,col=2)
plot(v_true_list[[2]][1:100], v_vec_list_init[[2]][1:100]); abline(a=0,b=1,col=2)


do_cluster_pdf(spks_time_mlist = spks_time_mlist,
               stim_onset_vec = stim_onset_vec,
               clusters_list_init = clusters_list_init,
               v_vec_list_init = v_vec_list_init,
               # v_vec_list_init = v_true_list,
               step_size = 0.0005,
               N_component = 2, 
               freq_trun = Inf,
               MaxIter = 10, 
               conv_thres = 0,
               gamma=0,
               t_vec=t_vec,
               fix_timeshift = FALSE)->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list
plot(v_true_list[[1]], tmp$v_vec_list[[1]]);abline(a=0, b=1, col='red')
plot(v_true_list[[2]], tmp$v_vec_list[[2]]);abline(a=0, b=1, col='red')

plot(v_vec_list_init[[1]], tmp$v_vec_list[[1]]);abline(a=0, b=1, col='red')
plot(v_vec_list_init[[2]], tmp$v_vec_list[[2]]);abline(a=0, b=1, col='red')


res_list = list(tmp)
id_res = 1
tmp3 = get_center_intensity_array(spks_time_mlist = spks_time_mlist,
                                 stim_onset_vec = stim_onset_vec,
                                 clusters_list = mem2clus(1:N_node),
                                 v_vec = tmp$v_vec_list[[1]],
                                 N_component = 1,
                                 freq_trun = 5,
                                 v0 = 0.5, 
                                 v1 = 0.5,
                                 t_vec = tmp$t_vec,
                                 bw=0.01)
tmp3$center_intensity_array-> node_intensity_array_shifted
tmp3$center_density_array-> node_density_array_shifted
tmp3$center_Nspks_mat[,1]-> node_Nspks_vec
tmp2 = plot_intensity_array(center_intensity_array = res_list[[id_res]]$center_intensity_array[,1,,drop=F]+
                             1*res_list[[id_res]]$center_intensity_array[,2,,drop=F],
                           center_Nspks_vec = rowSums(res_list[[id_res]]$center_Nspks_mat), 
                           clusters_list = res_list[[id_res]]$clusters_list, 
                           # node_intensity_array = node_intensity_array_shifted, 
                           # node_Nspks_vec = node_Nspks_vec, 
                           t_vec = res_list[[id_res]]$t_vec)
grid.arrange(tmp2$g)

