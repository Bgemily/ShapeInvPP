
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

# Parallel computing setup ------------------------------------------------

library(foreach)
library(doParallel)
N_cores = 20
registerDoParallel(cores=N_cores)


# Prepare data ----------------------------
id_session=1;
new.path='../Data/Main/'
spks_time_mlist = matrix(nrow=0,ncol=1)
stim_onset_vec = c()
id_session_vec = c()
id_trial_vec = c()
id_node_vec = c()
response_type_vec = c()
brain_area_vec = c()
for (id_session in 1:3) {
  ### Load data 
  dat = readRDS(paste(new.path,"session",id_session,".rds",sep=''))
  id_trial_vec_1 = which(dat$scenario_num==-1 
        # & dat$reaction_type==1 
        # & dat$contrast_left==1 & dat$contrast_right==0
        )
  id_node_vec_1 = which(dat$brain_region=='vis ctx')
  spks_time_mlist_1 = dat$spks_pp[id_node_vec_1, id_trial_vec_1]
  stim_onset_vec_1 = dat$stim_onset[id_trial_vec_1, 1]
  
  response_type_vec_1 = dat$response[id_trial_vec_1, 1]
  brain_area_1 = dat$brain_area[id_node_vec_1]
  

  ### Reshape: N_node x N_trial -> (N_node*N_trial) x 1
  N_node = length(id_node_vec_1)
  N_trial = length(id_trial_vec_1)
  spks_time_mlist_2 = matrix(t(spks_time_mlist_1),
                             nrow=length(spks_time_mlist_1),
                             ncol=1)
  stim_onset_vec_2 = rep(stim_onset_vec_1, N_node)
  id_trial_vec_2 = rep(id_trial_vec_1, N_node)
  response_type_vec_2 = rep(response_type_vec_1, N_node)
  id_node_vec_2 = rep(id_node_vec_1, each=N_trial)
  brain_area_2 = rep(brain_area_1, each=N_trial)
  
  ### Align (N_node*N_trial) spike trains by their stimuli onset time
  spks_time_mlist_3 = spks_time_mlist_2
  for (id_nodetrial in 1:nrow(spks_time_mlist_2)) {
    spks_time_mlist_3[id_nodetrial,1] = list(spks_time_mlist_2[id_nodetrial,1][[1]] - stim_onset_vec_2[id_nodetrial])
  }
  stim_onset_vec_3 = stim_onset_vec_2*0
  
  ### Stack processed data
  spks_time_mlist = rbind(spks_time_mlist, spks_time_mlist_3)
  stim_onset_vec = c(stim_onset_vec, stim_onset_vec_3)
  id_session_vec = c(id_session_vec, rep(id_session, length(stim_onset_vec_3)))
  id_node_vec = c(id_node_vec, id_node_vec_2)
  id_trial_vec = c(id_trial_vec, id_trial_vec_2)
  response_type_vec = c(response_type_vec, response_type_vec_2)
  brain_area_vec = c(brain_area_vec, brain_area_2)
}

### Remove neurons with small N_spks
spks_time_mlist_full = spks_time_mlist
stim_onset_vec_full = stim_onset_vec
id_session_vec_full = id_session_vec
id_node_vec_full = id_node_vec
id_trial_vec_full = id_trial_vec
response_type_vec_full = response_type_vec
brain_area_vec_full = brain_area_vec

N_node = nrow(spks_time_mlist_full)
N_trial = ncol(spks_time_mlist_full)
N_spks_vec = rep(0, N_node)
for (id_node in 1:N_node){
  for (id_trial in 1:N_trial){
    spks_time_tmp = unlist(spks_time_mlist_full[id_node,id_trial]) - stim_onset_vec_full[id_trial]
    spks_time_tmp = spks_time_tmp[which(spks_time_tmp<=0 & 
                                          spks_time_tmp>=-0.5)]
    N_spks_tmp = length(spks_time_tmp)
    N_spks_vec[id_node] = N_spks_vec[id_node] + N_spks_tmp
  }
}

subsample = which(N_spks_vec>=3)
# set.seed(831)
# subsample=sample(subsample,500)
spks_time_mlist = spks_time_mlist_full[subsample, ,drop=FALSE]
stim_onset_vec = stim_onset_vec_full[subsample]
id_session_vec = id_session_vec_full[subsample]
id_node_vec = id_node_vec_full[subsample]
id_trial_vec = id_trial_vec_full[subsample]
response_type_vec = response_type_vec_full[subsample]
brain_area_vec = brain_area_vec_full[subsample]


# Apply our algorithm ---------------------------------------------------------

method = paste0("Model4_est_timeshift")
signal_type = 'pre_stim'


N_clus_vec = c(2,3)
freq_trun_vec = c(Inf)
gamma_vec = c(0,1,10)
freq_trun = Inf
N_clus = 3
gamma = 1
MaxIter = 20
step_size = 5e-5
v0 = 0.0
v1 = 0.5
t_vec=seq(0-v1, v0, length.out=200)
fix_timeshift = FALSE
N_restart = 1

# now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
foreach(id_gamma = 1:length(gamma_vec)) %:% 
  foreach (ind_N_clus = 1:length(N_clus_vec)) %dopar% {
      gamma = gamma_vec[id_gamma]
      N_clus = N_clus_vec[ind_N_clus]
      set.seed(831)
      res = get_init(spks_time_mlist = spks_time_mlist, 
                     stim_onset_vec = stim_onset_vec,
                     N_clus = N_clus,
                     freq_trun = freq_trun,
                     t_vec = t_vec,
                     v0 = v0, v1 = v1, 
                     fix_timeshift = fix_timeshift)
      
      clusters_list_init = res$clusters_list
      v_vec_init = res$v_vec
      
      # Apply algorithm
      res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                           stim_onset_vec = stim_onset_vec,
                           clusters_list_init = clusters_list_init,
                           v_vec_init = v_vec_init,
                           freq_trun = freq_trun, 
                           MaxIter = MaxIter,
                           t_vec = t_vec,
                           v0 = v0, v1 = v1, 
                           gamma = gamma,
                           fix_timeshift = fix_timeshift,
                           step_size = step_size)
      # Restart
      loss_restart = c(tail(res$loss_history,1))
      if(N_restart>1){
        for (. in 1:(N_restart-1)) {
          mem_init = clus2mem(clusters_list_init)
          mem_init_2 = mem_init
          mem_init_2[sample(1:length(mem_init_2),length(mem_init_2)%/%5)] = sample(1:N_clus,length(mem_init_2)%/%5,replace=TRUE)
          clusters_list_init_2 = mem2clus(mem_init_2)
          v_vec_init_2 = v_vec_init
          res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                                 stim_onset_vec = stim_onset_vec,
                                 clusters_list_init = clusters_list_init_2,
                                 v_vec_init = v_vec_init_2,
                                 freq_trun = freq_trun, 
                                 MaxIter = MaxIter,
                                 t_vec = t_vec,
                                 v0 = v0, v1 = v1, 
                                 gamma = gamma,
                                 fix_timeshift = fix_timeshift)
          loss_restart = c(loss_restart, tail(res_2$loss_history,1))
          if(tail(res_2$loss_history,1) < tail(res$loss_history,1)){
            res = res_2
          }
        }
      }
      
      folder_path = paste0('../Results/Rdata/RDA', 
                           '/', method, "_", 
                           'gamma', gamma,
                           '/', 'session', 'ALL',
                           '/', signal_type)
      dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
      data_res = list(spks_time_mlist_analyzed = spks_time_mlist,
                      stim_onset_vec_analyzed = stim_onset_vec,
                      subsample = subsample,
                      spks_time_mlist = spks_time_mlist_full,
                      stim_onset_vec = stim_onset_vec_full,
                      id_trial_vec = id_trial_vec_full, 
                      id_node_vec = id_node_vec_full,
                      id_session_vec = id_session_vec_full,
                      response_type_vec = response_type_vec_full,
                      brain_area_vec = brain_area_vec_full)
      param_res = list(N_clus=N_clus,
                       freq_trun=freq_trun,
                       gamma=gamma,
                       v0=v0, v1=v1,
                       t_vec=t_vec,
                       N_restart=N_restart)
      save(res,
           loss_restart,
           data_res,
           param_res,
           file = paste0(folder_path, 
                         '/', "Nclus", N_clus, 
                         '_', "lr", step_size, 
                         '.Rdata'))
    }
  



