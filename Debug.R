rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 1:5
reaction_time_vec = stim_onset_vec + 0.5
spks_time_mlist_tmp = matrix(list(),20,5)

for (i in 1:10) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.2,1))
  }
}

for (i in 11:20) {
  for (j in 1:5) {
    spks_time_mlist_tmp[i,j] = list(rep(stim_onset_vec[j]+0.2,10))
  }
}


cluster_kmeans_pdf(spks_time_mlist = spks_time_mlist_tmp, stim_onset_vec, reaction_time_vec, 
                   clusters_list = list(1:11,12:20)) -> tmp
tmp$clusters_list

tmp = res
tmp$t_vec = t_vec
grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)

res = cluster::pam(x=poinproc_mat_2, k=N_clus, diss=FALSE, cluster.only=FALSE)
par(mfrow=c(4,1))
plot(res$medoids[1,],type='l')
plot(res$medoids[2,],type='l')
plot(res$medoids[3,],type='l')
plot(res$medoids[4,],type='l')

table(res$clustering)




