rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 1:1
reaction_time_vec = stim_onset_vec + 0.5
N_node = 100
spks_time_mlist_tmp = matrix(list(),N_node,1)

for (i in 1:(N_node/2)) {
  for (j in 1:1) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(sample(20+1:5,1),0,0.1))
  }
}

for (i in (N_node/2+1):N_node) {
  for (j in 1:1) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(sample(200+1:5,1),0,0.1))
  }
}


do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec, reaction_time_vec, 
               clusters_list_init = mem2clus(N_clus_min = 2, 
                                             membership = rep(1:2,times=c(N_node/2-N_node/5,N_node/2+N_node/5))), 
               v_vec_init = rep(0,100),
               freq_trun = Inf,
               MaxIter = 15,
               v0 = 0.15, v1=0, gamma=00000,
               t_vec=seq(0, 0.15, length.out=200),fix_timeshift = TRUE)->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list
summary(tmp$v_vec)

grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)


plot(res$loss_history,type='b')
center_intensity_array = res$center_intensity_array
grid.arrange(plot_intensity_array(center_intensity_array, res$clusters_list, res$t_vec)$g)

mem = clus2mem(res$clusters_list)
plot(jitter(mem), jitter(data_res$id_session_vec))

df = tibble(mem=mem, 
            resp_type=data_res$response_type_vec,
            brain_area = data_res$brain_area_vec,
            id_trial=data_res$id_trial_vec, 
            id_node=data_res$id_node_vec, 
            id_session=data_res$id_session_vec)
View(df)


