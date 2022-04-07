
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
N_cores = 10
registerDoParallel(cores=N_cores)



# Loop over brain regions -------------------------------------------------

brain_region_vec = c("vis ctx", "hipp",
                     "frontal", "basal ganglia",
                      "midbrain", "thal"
                     )
brain_region_vec = c("vis ctx","hipp")
for (id_brain_region in 1:length(brain_region_vec)) {
  brain_region = brain_region_vec[id_brain_region]
  
  # Prepare data ----------------------------
  # brain_region='vis ctx'
  id_session=1;
  new.path='../Data/Main/'
  spks_time_mlist = matrix(nrow=0,ncol=1)
  stim_onset_vec = c()
  id_session_vec = c()
  id_trial_vec = c()
  id_node_vec = c()
  pre_feedback_type_vec = c()
  response_type_vec = c()
  is_passive_vec = c()
  brain_area_vec = c()
  for (id_session in 1:3) {
    ### Load data 
    dat = readRDS(paste(new.path,"session",id_session,".rds",sep=''))
    dat_P = readRDS(paste(new.path,"replay_session",id_session,".rds",sep=''))
    id_trial_vec_1 = which(dat$scenario_num==-1 
                           # & dat$reaction_type==1 
                           # & dat$contrast_left==1 & dat$contrast_right==0
    )
    id_trial_vec_1_P = which(dat_P$scenario=="passiveVisual")
    id_trial_vec_1_P = id_trial_vec_1_P[which(dat_P$contrast_left<dat_P$contrast_right)]
    id_node_vec_1 = which(dat$brain_region==brain_region)
    if(length(id_node_vec_1)>0){
      spks_time_mlist_1 = dat$spks_pp[id_node_vec_1, id_trial_vec_1]
      spks_time_mlist_1_P = dat_P$spks_pp[id_node_vec_1, id_trial_vec_1_P]
      spks_time_mlist_1 = cbind(spks_time_mlist_1, spks_time_mlist_1_P)
      
      stim_onset_vec_1 = dat$stim_onset[id_trial_vec_1, 1]
      stim_onset_vec_1_P = dat_P$stim_onset[id_trial_vec_1_P]
      stim_onset_vec_1 = c(stim_onset_vec_1, stim_onset_vec_1_P)
      
      feedback_type_vec_1 = dat$feedback_type[id_trial_vec_1, 1]
      feedback_type_vec_1_P = rep(NA,length(id_trial_vec_1_P))
      feedback_type_vec_1 = c(feedback_type_vec_1, feedback_type_vec_1_P)
      
      pre_feedback_type_vec_1 = c(NA, feedback_type_vec_1[-length(feedback_type_vec_1)])
      
      response_type_vec_1 = dat$response[id_trial_vec_1, 1]
      response_type_vec_1_P = rep(NA,length(id_trial_vec_1_P))
      response_type_vec_1 = c(response_type_vec_1, response_type_vec_1_P)
      
      is_passive_vec_1 = c(rep(0,length(id_trial_vec_1)), 
                           rep(1,length(id_trial_vec_1_P)) )
      
      id_trial_vec_1 = c(id_trial_vec_1, id_trial_vec_1_P+length(dat$scenario_num))
      
      brain_area_1 = dat$brain_area[id_node_vec_1]
      
      
      ### Reshape: N_node x N_trial -> (N_node*N_trial) x 1
      N_node = length(id_node_vec_1)
      N_trial = length(id_trial_vec_1)
      spks_time_mlist_2 = matrix(t(spks_time_mlist_1),
                                 nrow=length(spks_time_mlist_1),
                                 ncol=1)
      stim_onset_vec_2 = rep(stim_onset_vec_1, N_node)
      id_trial_vec_2 = rep(id_trial_vec_1, N_node)
      pre_feedback_type_vec_2 = rep(pre_feedback_type_vec_1, N_node)
      response_type_vec_2 = rep(response_type_vec_1, N_node)
      is_passive_vec_2 = rep(is_passive_vec_1, N_node)
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
      pre_feedback_type_vec = c(pre_feedback_type_vec, pre_feedback_type_vec_2)
      response_type_vec = c(response_type_vec, response_type_vec_2)
      is_passive_vec = c(is_passive_vec, is_passive_vec_2)
      brain_area_vec = c(brain_area_vec, brain_area_2)
    }
    
  }
  
  ### Remove neurons with small N_spks
  spks_time_mlist_full = spks_time_mlist
  stim_onset_vec_full = stim_onset_vec
  id_session_vec_full = id_session_vec
  id_node_vec_full = id_node_vec
  id_trial_vec_full = id_trial_vec
  pre_feedback_type_vec_full = pre_feedback_type_vec
  response_type_vec_full = response_type_vec
  is_passive_vec_full = is_passive_vec
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
  # subsample=sample(subsample,50)
  spks_time_mlist = spks_time_mlist_full[subsample, ,drop=FALSE]
  stim_onset_vec = stim_onset_vec_full[subsample]
  
  
  # Apply our algorithm ---------------------------------------------------------

  method = paste0("Model4_add_passive")
  signal_type = 'pre_stim'


  N_clus_vec = c(2,3,4)
  freq_trun_vec = c(Inf)
  gamma_vec = c(0)
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

  # foreach(id_gamma = 1:length(gamma_vec)) %:%
  #   foreach (id_N_clus = 1:length(N_clus_vec)) %dopar% {
  #     gamma = gamma_vec[id_gamma]
  #     N_clus = N_clus_vec[id_N_clus]
  #     set.seed(831)
  #     res = get_init(spks_time_mlist = spks_time_mlist,
  #                    stim_onset_vec = stim_onset_vec,
  #                    N_clus = N_clus,
  #                    freq_trun = freq_trun,
  #                    t_vec = t_vec,
  #                    v0 = v0, v1 = v1, rmv_conn_prob = TRUE,
  #                    fix_timeshift = fix_timeshift)
  # 
  #     clusters_list_init = res$clusters_list
  #     v_vec_init = res$v_vec
  # 
  #     # Apply algorithm
  #     res = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
  #                          stim_onset_vec = stim_onset_vec,
  #                          clusters_list_init = clusters_list_init,
  #                          v_vec_init = v_vec_init,
  #                          freq_trun = freq_trun,
  #                          MaxIter = MaxIter,
  #                          t_vec = t_vec,
  #                          v0 = v0, v1 = v1,
  #                          gamma = gamma,
  #                          fix_timeshift = fix_timeshift,
  #                          step_size = step_size)
  #     # Restart
  #     loss_restart = c(tail(res$loss_history,1))
  #     if(N_restart>1){
  #       for (. in 1:(N_restart-1)) {
  #         mem_init = clus2mem(clusters_list_init)
  #         mem_init_2 = mem_init
  #         mem_init_2[sample(1:length(mem_init_2),length(mem_init_2)%/%5)] = sample(1:N_clus,length(mem_init_2)%/%5,replace=TRUE)
  #         clusters_list_init_2 = mem2clus(mem_init_2)
  #         v_vec_init_2 = v_vec_init
  #         res_2 = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
  #                                stim_onset_vec = stim_onset_vec,
  #                                clusters_list_init = clusters_list_init_2,
  #                                v_vec_init = v_vec_init_2,
  #                                freq_trun = freq_trun,
  #                                MaxIter = MaxIter,
  #                                t_vec = t_vec,
  #                                v0 = v0, v1 = v1,
  #                                gamma = gamma,
  #                                fix_timeshift = fix_timeshift)
  #         loss_restart = c(loss_restart, tail(res_2$loss_history,1))
  #         if(tail(res_2$loss_history,1) < tail(res$loss_history,1)){
  #           res = res_2
  #         }
  #       }
  #     }
  # 
  #     folder_path = paste0('../Results/Rdata/RDA',
  #                          '/', method, "_",
  #                          'gamma', gamma,
  #                          '/', 'brain_region', brain_region,
  #                          '/', signal_type)
  #     dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  #     data_res = list(N_spks_vec = N_spks_vec,
  #                     brain_region = brain_region,
  #                     subsample = subsample,
  #                     spks_time_mlist = spks_time_mlist_full,
  #                     stim_onset_vec = stim_onset_vec_full,
  #                     id_trial_vec = id_trial_vec_full,
  #                     id_node_vec = id_node_vec_full,
  #                     id_session_vec = id_session_vec_full,
  #                     pre_feedback_type_vec = pre_feedback_type_vec_full,
  #                     response_type_vec = response_type_vec_full,
  #                     is_passive_vec = is_passive_vec_full,
  #                     brain_area_vec = brain_area_vec_full)
  #     param_res = list(N_clus=N_clus,
  #                      freq_trun=freq_trun,
  #                      gamma=gamma,
  #                      v0=v0, v1=v1,
  #                      t_vec=t_vec,
  #                      N_restart=N_restart,
  #                      MaxIter = MaxIter,
  #                      step_size = step_size,
  #                      fix_timeshift = fix_timeshift)
  #     save(res,
  #          loss_restart,
  #          data_res,
  #          param_res,
  #          file = paste0(folder_path,
  #                        '/', "Nclus", N_clus,
  #                        '.Rdata'))
  #   }


  # Select N_clus and gamma -------------------------------------------------

  ### Load analysis result 
  res_list = list()
  cand_N_clus_vec = c()
  cand_gamma_vec = c()
  cand_file_vec = c()
  for (id_gamma in 1:length(gamma_vec)) {
    gamma = gamma_vec[id_gamma]
    for (id_N_clus in 1:length(N_clus_vec)){
      N_clus = N_clus_vec[id_N_clus]
      folder_path = paste0('../Results/Rdata/RDA',
                           '/', method, "_",
                           'gamma', gamma,
                           '/', 'brain_region', brain_region,
                           '/', signal_type)
      file = paste0(folder_path,
                    '/', "Nclus", N_clus,
                    '.Rdata')
      load(file)
      
      res_list = c(res_list, list(res))
      cand_N_clus_vec = c(cand_N_clus_vec, N_clus)
      cand_gamma_vec = c(cand_gamma_vec, gamma)
      cand_file_vec = c(cand_file_vec, file)
    }
  }
  
  
  
  ### Get best model
  res_select_model = select_model(spks_time_mlist = spks_time_mlist, 
                                  stim_onset_vec = stim_onset_vec, 
                                  result_list = res_list)
  
  N_clus_best = cand_N_clus_vec[res_select_model$id_best_res]
  gamma_best = cand_gamma_vec[res_select_model$id_best_res]
  
  res_select_model = c(res_select_model,
                       list(N_clus_best = N_clus_best,
                            gamma_best = gamma_best,
                            cand_N_clus_vec = cand_N_clus_vec,
                            cand_gamma_vec = cand_gamma_vec,
                            cand_file_vec = cand_file_vec
                            )
                       )
  
  folder_path = paste0('../Results/Rdata/RDA',
                       '/', method, 
                       '/', 'brain_region', brain_region,
                       '/', signal_type)
  dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
  save(res_select_model,
       file = paste0(folder_path,
                     '/', "res_select_model",
                     '.Rdata'))
  
  file_best = cand_file_vec[res_select_model$id_best_res]
  rm(res, loss_restart, data_res, param_res)
  load(file_best)
  save(res,
       loss_restart,
       data_res,
       param_res,
       file = paste0(folder_path,
                     '/', "res_best_model",
                     '.Rdata'))
  
}
