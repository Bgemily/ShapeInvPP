rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 0
reaction_time_vec = stim_onset_vec + 0.5
N_node = 100
spks_time_mlist_tmp = matrix(list(),N_node,1)

v_vec_true = rep(0,N_node)
for (i in 1:(N_node/2)) {
  for (j in 1:1) {
    v_vec_true[i] = runif(1,0,0.35)*0
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+rnorm(20,0.1,0.05)+v_vec_true[i])
  }
}

for (i in 1:(N_node/4)) {
  for (j in 1:1) {
    v_vec_true[i] = runif(1,0.05,0.2)*0
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+rnorm(40,0.1,0.05)+v_vec_true[i])
  }
}

for (i in (N_node/2+1):N_node) {
  for (j in 1:1) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(10,0,0.25))
  }
}




for (i in (1):N_node) {
  for (j in 1:1) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(10,0,0.5))
  }
}

t_vec = seq(0, 0.5, length.out=200)

res_list = list()

N_clus = 1
res = get_init(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec = stim_onset_vec,
               N_clus = N_clus,
               t_vec = t_vec,
               fix_timeshift = FALSE)
clusters_list_init = res$clusters_list
v_vec_init = res$v_vec
do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec, reaction_time_vec, 
               clusters_list_init = clusters_list_init, 
               v_vec_init = v_vec_init,
               freq_trun = Inf,
               MaxIter = 10,
               gamma=10,
               t_vec=t_vec,
               fix_timeshift = FALSE)->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list

res_list[[1]] = tmp



N_clus = 1
res = get_init(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec = stim_onset_vec,
               N_clus = N_clus,
               t_vec = t_vec,
               fix_timeshift = TRUE)
clusters_list_init = res$clusters_list
v_vec_init = res$v_vec
do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec, reaction_time_vec, 
               clusters_list_init = clusters_list_init, 
               v_vec_init = v_vec_init,
               freq_trun = Inf,
               MaxIter = 10,
               gamma=10,
               t_vec=t_vec,
               fix_timeshift = TRUE)->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list

res_list[[2]] = tmp



id_res = 1
tmp = plot_intensity_array(center_intensity_array = res_list[[id_res]]$center_intensity_array,
                     center_Nspks_vec = res_list[[id_res]]$center_Nspks_mat[,1], 
                     clusters_list = res_list[[id_res]]$clusters_list,
                     t_vec = res_list[[id_res]]$t_vec)
grid.arrange(tmp$g)
id_res = 2
tmp = plot_intensity_array(center_intensity_array = res_list[[id_res]]$center_intensity_array,
                           center_Nspks_vec = res_list[[id_res]]$center_Nspks_mat[,1], 
                           clusters_list = res_list[[id_res]]$clusters_list,
                           t_vec = res_list[[id_res]]$t_vec)
grid.arrange(tmp$g)


res_select_model = select_model(spks_time_mlist = spks_time_mlist_tmp, 
                                stim_onset_vec = stim_onset_vec, 
                                result_list = res_list)
res_select_model$log_lik_vec
res_select_model$log_lik_1_vec
res_select_model$log_lik_2_vec

ggplot()+
  geom_line(aes(x=N_clus_vec, 
                y=res_select_model$log_lik_vec
  ))

ggplot()+
  geom_line(aes(x=N_clus_vec, 
                y=res_select_model$log_lik_1_vec
  ))

ggplot()+
  geom_line(aes(x=N_clus_vec, 
                y=res_select_model$log_lik_2_vec
  ))

#######
library(ggplot2)
gamma_vec = factor(res_select_model$cand_gamma_vec)
ggplot()+
  geom_line(aes(x=res_select_model$cand_N_clus_vec, 
                y=res_select_model$log_lik_1_vec,
                group=gamma_vec,
                color=gamma_vec
                ))
ggplot()+
  geom_line(aes(x=res_select_model$cand_N_clus_vec, 
                y=res_select_model$log_lik_2_vec,
                group=gamma_vec,
                color=gamma_vec
  ))

ggplot()+
  geom_line(aes(x=res_select_model$cand_N_clus_vec, 
                y=res_select_model$log_lik_vec,
                group=gamma_vec,
                color=gamma_vec
  ))

ggplot()+
  geom_line(aes(x=res_select_model$cand_N_clus_vec, 
                y=res_select_model$ICL_vec,
                group=gamma_vec,
                color=gamma_vec
  ))

res_select_model$N_clus_best
res_select_model$gamma_best



