rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

stim_onset_vec = 1:1
reaction_time_vec = stim_onset_vec + 0.5
N_node = 100
spks_time_mlist_tmp = matrix(list(),N_node,1)

v_vec_true = rep(0,N_node)
for (i in 1:(N_node/2)) {
  for (j in 1:1) {
    v_vec_true[i] = runif(1,0,0.35)
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+rnorm(10,0.1,0.05)+v_vec_true[i])
  }
}

for (i in (N_node/2+1):N_node) {
  for (j in 1:1) {
    spks_time_mlist_tmp[i,j] = list(stim_onset_vec[j]+runif(7,0,0.25))
  }
}

t_vec = seq(0, 0.5, length.out=200)

tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist_tmp, 
                                 stim_onset_vec = stim_onset_vec, 
                                 clusters_list = mem2clus(1:N_node), 
                                 v_vec = rep(0,N_node),
                                 t_vec = t_vec , bw = 0.01)
node_density_array = tmp$center_density_array
node_Nspks_vec = tmp$center_Nspks_mat[,1]

N_sample = 5
grid.arrange(plot_intensity_array(center_intensity_array = node_density_array[sample(1:N_node,N_sample),,,drop=F], 
                                  center_Nspks_vec = node_Nspks_vec[sample(1:N_node,N_sample)],
                                  clusters_list = mem2clus(1:N_node)[sample(1:N_node,N_sample)], 
                                  t_vec = t_vec)$g)


res = get_init(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec = stim_onset_vec,
               N_clus = 2,
               t_vec = t_vec,
               v0 = v0, v1 = v1, 
               fix_timeshift = FALSE)

clusters_list_init = res$clusters_list
clusters_list_init
v_vec_init = res$v_vec
plot(v_vec_true, v_vec_init); abline(a=0,b=1,col=2)


do_cluster_pdf(spks_time_mlist = spks_time_mlist_tmp, 
               stim_onset_vec, reaction_time_vec, 
               clusters_list_init = clusters_list_init, 
               v_vec_init = v_vec_init,
               freq_trun = Inf,
               MaxIter = 5,
               v0 = 0.15, v1=0, gamma=0,
               t_vec=seq(0, 0.5, length.out=200),
               fix_timeshift = TRUE)->tmp
plot(tmp$loss_history,type='b')
tmp$clusters_list
summary(tmp$v_vec)
plot(v_vec_true, tmp$v_vec); abline(a=0,b=1,col=2)

grid.arrange(plot_intensity_array(center_intensity_array = tmp$center_intensity_array, 
                                  center_Nspks_vec = tmp$center_Nspks_mat[,1],
                                  clusters_list = tmp$clusters_list, t_vec = tmp$t_vec)$g)



i01 <- function(w,y=1, z=1) {
  list(y = y, z = z)
}

i02 <- function(x, w,...) {
  i01(w,...)
}

str(i02(x = 1, w=1,  z = 3,y=0))

