rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 1:5
reaction_time_vec = stim_onset_vec + 0.5
spks_time_mlist_tmp = matrix(list(),200,5)

for (i in 1:50) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.2,10))
  }
}

for (i in 51:50) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.1,10))
  }
}

for (i in 101:150) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.3,10))
  }
}

for (i in 151:200) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.3,10))
  }
}

do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, stim_onset_vec, reaction_time_vec, 
               clusters_list_init = mem2clus(N_clus_min = 2, 
                                             membership = sample(1:2,200,replace=TRUE)) )->tmp
plot(tmp$loss_history,type='b')
grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)


cluster_kmeans_pdf(spks_time_mlist = spks_time_mlist_tmp, stim_onset_vec, reaction_time_vec, 
                   clusters_list = list(1:14,15:20)) -> tmp
tmp$clusters_list

tmp = res
tmp$t_vec = t_vec
grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)

