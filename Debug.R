rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 1:5
reaction_time_vec = stim_onset_vec + 0.5
N_node = 100
spks_time_mlist_tmp = matrix(list(),N_node,5)

for (i in 1:(N_node/2)) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(10,0,0.05))
  }
}

for (i in (N_node/2+1):N_node) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(10,0,0.1))
  }
}


do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec, reaction_time_vec, 
               clusters_list_init = mem2clus(N_clus_min = 2, 
                                             membership = rep(1:2,times=c(N_node/2-N_node/5,N_node/2+N_node/5))), 
               v_vec_init = rep(0,100),
               freq_trun = Inf,
               MaxIter = 15,
               v0 = 0.15, v1=0,
               t_vec=seq(0, 0.15, length.out=200))->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list
summary(tmp$v_vec)

grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)

res = get_center_intensity_array(spks_time_mlist = spks_time_mlist_tmp, 
                                                    stim_onset_vec = stim_onset_vec, 
                                                    reaction_time_vec = reaction_time_vec, 
                                                    clusters_list = tmp$clusters_list, 
                                                    v_vec = tmp$v_vec,
                                                    N_component = 1,
                                                    freq_trun = 8, 
                                                    t_vec = tmp$t_vec,
                                                    v0 = 0.15, v1 = 0,
                                                    rmv_conn_prob = TRUE)
center_intensity_array = res$center_intensity_array
grid.arrange(plot_intensity_array(center_intensity_array, clusters_list, t_vec)$g)
