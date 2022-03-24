
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
        dat$contrast_left==1 & dat$contrast_right==0 &
        dat$reaction_time-dat$stim_onset[,1]<0.5)->id_trials

spks_time_mlist = dat$spks_pp[,id_trials]

N_node = nrow(spks_time_mlist)
N_trial = ncol(spks_time_mlist)
N_spks_vec = rep(0, N_node)
for (id_node in 1:N_node){
  for (id_trial in 1:N_trial){
    spks_time_tmp = unlist(spks_time_mlist[id_node,id_trial])
    N_spks_tmp = length(spks_time_tmp)
    N_spks_vec[id_node] = N_spks_vec[id_node] + N_spks_tmp
  }
}

id_nodes = which(N_spks_vec>=30)



# Apply our algorithm (vision) ---------------------------------------------------------
spks_time_mlist = dat$spks_pp[id_nodes,id_trials]
stim_onset_vec = dat$stim_onset[id_trials,1]
reaction_time_vec = dat$reaction_time[id_trials]

method = paste0("Model3_rmvsmlNspks_incstpsiz")
signal_type = 'vision'

N_clus_vec = c(6,5,4)
freq_trun_vec = c(Inf)
MaxIter = 20
v0 = 0.15
v1 = 0.2
t_vec=seq(0-v1, v0, length.out=200)
N_restart = 1

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for(freq_trun in freq_trun_vec){
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
    v_vec_init = res$v_vec
    
    # Apply algorithm
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         reaction_time_vec = reaction_time_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_init = v_vec_init,
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
        v_vec_init_2 = v_vec_init
        res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec,
                               reaction_time_vec = reaction_time_vec,
                               clusters_list_init = clusters_list_init_2,
                               v_vec_init = v_vec_init_2,
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
    
    folder_path = paste0('../Results/Rdata/RDA', 
                         '/', method, "_", 'freq_trun',freq_trun,
                         '/', 'session',id_session,
                         '/', signal_type)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    data_res = list(spks_time_mlist=spks_time_mlist,
                    stim_onset_vec=stim_onset_vec,
                    reaction_time_vec=reaction_time_vec,
                    id_trials=id_trials,
                    id_nodes=id_nodes,
                    id_session=id_session)
    param_res = list(N_clus=N_clus,
                     freq_trun=freq_trun,
                     v0=v0, v1=v1,
                     t_vec=t_vec,
                     N_restart=N_restart)
    save(res,loss_restart,
         data_res,
         param_res,
         file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
  }
  
}



# Apply our algorithm (reaction) ---------------------------------------------------------
spks_time_mlist = dat$spks_pp[id_nodes,id_trials]
stim_onset_vec = dat$reaction_time[id_trials]
reaction_time_vec = dat$reaction_time[id_trials]

method = paste0("Model3_rmvsmlNspks_incstpsiz")
signal_type = 'reaction'

N_clus_vec = c(6,5,4)
freq_trun_vec = c(Inf)
MaxIter = 20
v0 = 0.2
v1 = 0.15
t_vec=seq(0-v1, v0, length.out=200)
N_restart = 1

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for(freq_trun in freq_trun_vec){
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
    v_vec_init = res$v_vec
    
    # Apply algorithm
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         reaction_time_vec = reaction_time_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_init = v_vec_init,
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
        v_vec_init_2 = v_vec_init
        res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec,
                               reaction_time_vec = reaction_time_vec,
                               clusters_list_init = clusters_list_init_2,
                               v_vec_init = v_vec_init_2,
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
    
    folder_path = paste0('../Results/Rdata/RDA', 
                         '/', method, "_", 'freq_trun',freq_trun,
                         '/', 'session',id_session,
                         '/', signal_type)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    data_res = list(spks_time_mlist=spks_time_mlist,
                    stim_onset_vec=stim_onset_vec,
                    reaction_time_vec=reaction_time_vec,
                    id_trials=id_trials,
                    id_nodes=id_nodes,
                    id_session=id_session)
    param_res = list(N_clus=N_clus,
                     freq_trun=freq_trun,
                     v0=v0, v1=v1,
                     t_vec=t_vec,
                     N_restart=N_restart)
    save(res,loss_restart,
         data_res,
         param_res,
         file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
  }
  
}





# Apply our algorithm (gocue) ---------------------------------------------------------
spks_time_mlist = dat$spks_pp[id_nodes,id_trials]
stim_onset_vec = dat$gocue[id_trials,1] 
reaction_time_vec = dat$reaction_time[id_trials]

method = paste0("Model3_rmvsmlNspks_incstpsiz")
signal_type = 'gocue'

N_clus_vec = c(6,5,4)
freq_trun_vec = c(Inf)
MaxIter = 20
v0 = 0.1
v1 = 0.05
t_vec=seq(0-v1, v0, length.out=200)
N_restart = 1

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for(freq_trun in freq_trun_vec){
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
    v_vec_init = res$v_vec
    
    # Apply algorithm
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         reaction_time_vec = reaction_time_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_init = v_vec_init,
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
        v_vec_init_2 = v_vec_init
        res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec,
                               reaction_time_vec = reaction_time_vec,
                               clusters_list_init = clusters_list_init_2,
                               v_vec_init = v_vec_init_2,
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
    
    folder_path = paste0('../Results/Rdata/RDA', 
                         '/', method, "_", 'freq_trun',freq_trun,
                         '/', 'session',id_session,
                         '/',signal_type)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    data_res = list(spks_time_mlist=spks_time_mlist,
                    stim_onset_vec=stim_onset_vec,
                    reaction_time_vec=reaction_time_vec,
                    id_trials=id_trials,
                    id_nodes=id_nodes,
                    id_session=id_session)
    param_res = list(N_clus=N_clus,
                     freq_trun=freq_trun,
                     v0=v0, v1=v1,
                     t_vec=t_vec,
                     N_restart=N_restart)
    save(res,loss_restart,
         data_res,
         param_res,
         file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
  }
  
}




# Apply our algorithm (response) ---------------------------------------------------------
spks_time_mlist = dat$spks_pp[id_nodes,id_trials]
stim_onset_vec = dat$response_time[id_trials,1] 
reaction_time_vec = dat$reaction_time[id_trials]

method = paste0("Model3_rmvsmlNspks_incstpsiz")
signal_type = 'response'

N_clus_vec = c(6,5,4)
freq_trun_vec = c(Inf)
MaxIter = 20
v0 = 0.05
v1 = 0.1
t_vec=seq(0-v1, v0, length.out=200)
N_restart = 1

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for(freq_trun in freq_trun_vec){
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
    v_vec_init = res$v_vec
    
    # Apply algorithm
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         reaction_time_vec = reaction_time_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_init = v_vec_init,
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
        v_vec_init_2 = v_vec_init
        res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec,
                               reaction_time_vec = reaction_time_vec,
                               clusters_list_init = clusters_list_init_2,
                               v_vec_init = v_vec_init_2,
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
    
    folder_path = paste0('../Results/Rdata/RDA', 
                         '/', method, "_", 'freq_trun',freq_trun,
                         '/', 'session',id_session,
                         '/',signal_type)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    data_res = list(spks_time_mlist=spks_time_mlist,
                    stim_onset_vec=stim_onset_vec,
                    reaction_time_vec=reaction_time_vec,
                    id_trials=id_trials,
                    id_nodes=id_nodes,
                    id_session=id_session)
    param_res = list(N_clus=N_clus,
                     freq_trun=freq_trun,
                     v0=v0, v1=v1,
                     t_vec=t_vec,
                     N_restart=N_restart)
    save(res,loss_restart,
         data_res,
         param_res,
         file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
  }
  
}



# Apply our algorithm (feedback) ---------------------------------------------------------
spks_time_mlist = dat$spks_pp[id_nodes,id_trials]
stim_onset_vec = dat$feedback_time[id_trials,1] 
reaction_time_vec = dat$reaction_time[id_trials]

method = paste0("Model3_rmvsmlNspks_incstpsiz")
signal_type = 'feedback'

N_clus_vec = c(6,5,4)
freq_trun_vec = c(Inf)
MaxIter = 20
v0 = 0.1
v1 = 0.05
t_vec=seq(0-v1, v0, length.out=200)
N_restart = 1

now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
for(freq_trun in freq_trun_vec){
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
    v_vec_init = res$v_vec
    
    # Apply algorithm
    res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                         stim_onset_vec = stim_onset_vec,
                         reaction_time_vec = reaction_time_vec,
                         clusters_list_init = clusters_list_init,
                         v_vec_init = v_vec_init,
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
        v_vec_init_2 = v_vec_init
        res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                               stim_onset_vec = stim_onset_vec,
                               reaction_time_vec = reaction_time_vec,
                               clusters_list_init = clusters_list_init_2,
                               v_vec_init = v_vec_init_2,
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
    
    folder_path = paste0('../Results/Rdata/RDA', 
                         '/', method, "_", 'freq_trun',freq_trun,
                         '/', 'session',id_session,
                         '/',signal_type)
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    data_res = list(spks_time_mlist=spks_time_mlist,
                    stim_onset_vec=stim_onset_vec,
                    reaction_time_vec=reaction_time_vec,
                    id_trials=id_trials,
                    id_nodes=id_nodes,
                    id_session=id_session)
    param_res = list(N_clus=N_clus,
                     freq_trun=freq_trun,
                     v0=v0, v1=v1,
                     t_vec=t_vec,
                     N_restart=N_restart)
    save(res,loss_restart,
         data_res,
         param_res,
         file = paste0(folder_path, '/', "Nclus", N_clus, '.Rdata'))
  }
  
}




