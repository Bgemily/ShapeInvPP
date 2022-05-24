
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


# Get mouse info for all sessions -----------------------------------------
load("../Data/mouse_name.rdata")
target_session_vec = 1:39
session_vec = session_vec[target_session_vec]
mouse_name_vec = mouse_name_vec[target_session_vec]
brain_region_list = brain_region_list[target_session_vec]

mouse_name_set_vec = unique(mouse_name_vec)
id_session_mouse_list = list()
brain_region_mouse_list = list()
for (id_mouse in 1:length(mouse_name_set_vec)){
  mouse_name = mouse_name_set_vec[id_mouse]
  id_session_mouse_vec = session_vec[which(mouse_name_vec==mouse_name)]
  id_session_mouse_list[[id_mouse]] =  id_session_mouse_vec
  brain_region_mouse_vec = unique(unlist(brain_region_list[which(mouse_name_vec==mouse_name)]))
  brain_region_mouse_list[[id_mouse]] = brain_region_mouse_vec
}

# Loop over brain regions and mouse -------------------------------------------------
for (id_mouse in 1:length(mouse_name_set_vec)) {
  print(id_mouse)
  for (id_brain_region in 1:length(brain_region_mouse_list[[id_mouse]])) {
    
    mouse_name = mouse_name_set_vec[id_mouse]
    id_session_mouse_vec = id_session_mouse_list[[id_mouse]]
    brain_region = brain_region_mouse_list[[id_mouse]][id_brain_region]
    
    # Prepare data ----------------------------
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
    for (id_session in id_session_mouse_vec) {
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
        
        feedback_type_vec_active = dat$feedback_type[, 1]
        pre_feedback_type_vec_active = c(NA, feedback_type_vec_active[-length(feedback_type_vec_active)])
        pre_feedback_type_vec_1 = pre_feedback_type_vec_active[id_trial_vec_1]
        pre_feedback_type_vec_1_P = rep(NA,length(id_trial_vec_1_P))
        pre_feedback_type_vec_1 = c(pre_feedback_type_vec_1, pre_feedback_type_vec_1_P)
        
        id_trial_vec_1 = c(id_trial_vec_1, id_trial_vec_1_P+length(dat$scenario_num))
        
        brain_area_1 = dat$brain_area[id_node_vec_1]
        
        
        ### Reshape: N_node x N_trial -> (N_node*N_trial) x 1
        N_node = length(id_node_vec_1)
        N_trial = length(id_trial_vec_1)
        
        pre_feedback_type_vec_2 = rep(pre_feedback_type_vec_1, N_node)
        
        ### Stack processed data
        pre_feedback_type_vec = c(pre_feedback_type_vec, pre_feedback_type_vec_2)
        
      }
      
    }
    
    # Modify saved data in results ---------------------------------------------------------
    
    method = paste0("Model4_multi_mouse")
    signal_type = 'pre_stim'
    
    N_clus_vec = c(1,2,3)
    gamma_vec = c(0)
    
    for(id_gamma in 1:length(gamma_vec)) {
      for (id_N_clus in 1:length(N_clus_vec)) {
        gamma = gamma_vec[id_gamma]
        N_clus = N_clus_vec[id_N_clus]
        
        folder_path = paste0('../Results/Rdata/RDA',
                             '/', method, "_",
                             'gamma', gamma,
                             '/', 'mouse_', mouse_name,
                             '/', 'brain_region', brain_region,
                             '/', signal_type)
        
        load(file = paste0(folder_path,
                           '/', "Nclus", N_clus,
                           '.Rdata'))
        data_res$pre_feedback_type_vec = pre_feedback_type_vec
        # identical(data_res$is_passive_vec, is_passive_vec)
        
        save(res,
             loss_restart,
             data_res,
             param_res,
             file = paste0(folder_path,
                           '/', "Nclus", N_clus,
                           '.Rdata'))
      }
    }
    
    
    
  }
}



