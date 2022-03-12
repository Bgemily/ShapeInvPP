rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

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
                   clusters_list = list(1:10,11:20)) -> tmp
tmp$clusters_list

plot(res$loss_history)
tmp=res
grid.arrange(plot_intensity_array(tmp$center_intensity_array, tmp$clusters_list, tmp$t_vec)$g)



