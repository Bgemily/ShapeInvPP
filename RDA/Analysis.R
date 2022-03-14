
# Default working directory: Sources/

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)
library(scales)
library(combinat)
library(cluster)
library(mclust)
library(reshape2)


# Load one session #####
id_session=1;new.path='../Data/Main/'
dat = readRDS(paste(new.path,"session",id_session,".rds",sep=''))

which(dat$scenario_num==1 & dat$reaction_type==1 & 
        dat$reaction_time-dat$stim_onset[,1]<0.5)->id_trials

spks_time_mlist = dat$spks_pp[,id_trials]
reaction_time_vec = dat$reaction_time[id_trials]
stim_onset_vec = dat$stim_onset[id_trials,1]


# Apply our algorithm ---------------------------------------------------------

N_clus_min = 1 # Number of clusters
N_clus_max = 10
N_clus_vec = c(4,8)
freq_trun = 7
MaxIter = 50
v0 = 0.2
v1 = 0.1
t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0)+0.01, by=0.01)
N_restart = 5

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for (ind_N_clus in 1:length(N_clus_vec)){
  N_clus = N_clus_vec[ind_N_clus]
  res = get_init(spks_time_mlist = spks_time_mlist, 
                 stim_onset_vec = stim_onset_vec,
                 reaction_time_vec = reaction_time_vec,
                 N_clus = N_clus,
                 freq_trun = freq_trun,
                 t_vec = t_vec,
                 v0 = v0, v1 = v1)
  
  clusters_list_init = res$clusters_list
  
  # Apply algorithm
  res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                       stim_onset_vec = stim_onset_vec,
                       reaction_time_vec = reaction_time_vec,
                       clusters_list_init = clusters_list_init,
                       freq_trun = freq_trun, 
                       MaxIter = MaxIter,
                       t_vec = t_vec,
                       v0 = v0, v1 = v1)
  # Restart
  loss_restart = c(tail(res$loss_history,1))
  if(N_restart>1){
    for (. in 1:(N_restart-1)) {
      mem_init = clus2mem(clusters_list_init)
      mem_init_2 = mem_init
      mem_init_2[sample(1:length(mem_init_2),length(mem_init_2)%/%5)] = sample(1:N_clus,length(mem_init_2)%/%5,replace=TRUE)
      clusters_list_init_2 = mem2clus(mem_init_2)
      res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                           stim_onset_vec = stim_onset_vec,
                           reaction_time_vec = reaction_time_vec,
                           clusters_list_init = clusters_list_init_2,
                           freq_trun = freq_trun, 
                           MaxIter = MaxIter,
                           t_vec = t_vec,
                           v0 = v0, v1 = v1)
      loss_restart = c(loss_restart, tail(res_2$loss_history,1))
      if(tail(res_2$loss_history,1) < tail(res$loss_history,1)){
        res = res_2
      }
    }
  }
  
  method = paste0("Model1")
  folder_path = paste0('../Results/Rdata/RDA/', method, '/', 'session',id_session,'/',now_trial)
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  save(res,loss_restart,
       spks_time_mlist,stim_onset_vec,reaction_time_vec,
       N_clus,freq_trun,v0,v1,t_vec,
       id_session,
       file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
}








