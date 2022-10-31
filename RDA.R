#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)
library(fdapace)


# Prepare data ------------------------------------------------------------
new.path = '../Data/Main/'
id_session = 8
scenario_num = c(-1)
feedback_type = 1
dat = readRDS(paste(new.path, "session",id_session,".rds",sep=''))

id_trial_success_vec = which((dat$scenario_num %in% scenario_num) & (dat$feedback_type == feedback_type))
brain_region = 'vis ctx'
id_neuron_vis = which(dat$brain_region == brain_region)
N_neuron = length(id_neuron_vis)
N_trial = length(id_trial_success_vec)

t_vec = seq(-2, 2,length.out=200)

### Get the number of spikes for selected neurons and trials
N_spks_df = c()
for (i in 1:N_neuron){
  for (j in 1:N_trial){
    id_neuron = id_neuron_vis[i]
    id_trial = id_trial_success_vec[j]
    
    spks_vec = dat$spks_pp[id_neuron, id_trial][[1]]
    spks_shifted_vec = spks_vec - dat$gocue[id_trial]
    feedback_time_shifted = dat$feedback_time[id_trial] - dat$gocue[id_trial]
    stim_onset_time_shifted = dat$stim_onset[id_trial] - dat$gocue[id_trial]
    spks_shifted_vec = spks_shifted_vec[which((spks_shifted_vec <= min(max(t_vec), feedback_time_shifted)) &  
                                                (spks_shifted_vec >= max(min(t_vec), stim_onset_time_shifted)) )]
    
    N_spks_df = rbind(N_spks_df, c(id_neuron, id_trial, length(spks_shifted_vec), feedback_time_shifted))
  }
}
N_spks_df = as.data.frame(N_spks_df)
colnames(N_spks_df) = c("id_neuron", "id_trial", "N_spks", "feedback_time")

### Select neuron-trial pairs 
id_neuron_id_trial_selected = N_spks_df[(N_spks_df$N_spks >= 10) & (N_spks_df$feedback_time >= 0.2), c('id_neuron', 'id_trial')]
spks_time_mlist = matrix(list(), nrow = nrow(id_neuron_id_trial_selected), ncol = 1)
stim_onset_time_mat = matrix(nrow = nrow(id_neuron_id_trial_selected), ncol = 1)
for (id_subj in 1:nrow(id_neuron_id_trial_selected)) {
  id_neuron = id_neuron_id_trial_selected[id_subj, 'id_neuron']
  id_trial = id_neuron_id_trial_selected[id_subj, 'id_trial']
  
  spks_vec = dat$spks_pp[id_neuron, id_trial][[1]]
  spks_shifted_vec = spks_vec - dat$gocue[id_trial]
  feedback_time_shifted = dat$feedback_time[id_trial] - dat$gocue[id_trial]
  stim_onset_time_shifted = dat$stim_onset[id_trial] - dat$gocue[id_trial]
  spks_shifted_vec = spks_shifted_vec[which((spks_shifted_vec <= min(max(t_vec), feedback_time_shifted)) &  
                                              (spks_shifted_vec >= max(min(t_vec), stim_onset_time_shifted)) )]
  spks_time_mlist[id_subj, 1] = list(spks_shifted_vec)
  stim_onset_time_mat[id_subj, 1] = stim_onset_time_shifted
}
stim_onset_vec = 0


# Fit model for various cluster number ------------------------------------
N_clus_min = 3
N_clus_max = 3
N_component = 2
key_times_vec = c(min(t_vec), 0, max(t_vec))
N_start_kmean = 5
freq_trun = 10
fix_timeshift = FALSE
fix_comp1_timeshift_only = FALSE
use_true_timeshift = FALSE
v_trialwise_vec_list = rep(list(0), 2)
N_restart = 1

set.seed(1)
res_list = list()
for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
  res_list[[ind_N_clus]] = list()
  N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
  
  res_best = NA
  compl_log_lik_best = -Inf
  for (id_restart in 1:N_restart) {
    ### Get initialization -----------
    res = get_init(spks_time_mlist = spks_time_mlist, 
                   stim_onset_vec = stim_onset_vec,
                   N_clus = N_clus_tmp,
                   N_component = N_component,
                   t_vec = t_vec, 
                   key_times_vec = key_times_vec,
                   N_start_kmean = N_start_kmean,
                   freq_trun = freq_trun,
                   fix_timeshift = fix_timeshift, 
                   v_trialwise_vec_list = v_trialwise_vec_list,
                   rmv_conn_prob = TRUE)
    center_density_array_init = res$center_density_array
    center_Nspks_mat_init = res$center_Nspks_mat
    clusters_list_init = res$clusters_list
    v_mat_list_init = res$v_mat_list
    
    
    # Apply algorithm ---------
    time_start = Sys.time()
    ### Estimation z,v,f based on pdf
    res_new = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                             stim_onset_vec = stim_onset_vec,
                             v_trialwise_vec_list = v_trialwise_vec_list,
                             center_density_array_init = center_density_array_init,
                             center_Nspks_mat_init = center_Nspks_mat_init, 
                             clusters_list_init = clusters_list_init,
                             v_mat_list_init = v_mat_list_init,
                             N_component = N_component, 
                             freq_trun = freq_trun,
                             gamma=0,
                             t_vec=t_vec, 
                             key_times_vec = key_times_vec,
                             fix_timeshift = fix_timeshift, 
                             fix_comp1_timeshift_only = fix_comp1_timeshift_only )
    time_end = Sys.time()
    time_estimation = time_end - time_start
    time_estimation = as.numeric(time_estimation, units='secs')
    
    ### Calculate log likelihood
    res = select_model(spks_time_mlist, 
                       stim_onset_vec, 
                       N_component = N_component,
                       key_times_vec = key_times_vec,
                       result_list = list(res_new))
    compl_log_lik_new = res$compl_log_lik_vec[1]
    
    ### Update best estimation
    if(compl_log_lik_new > compl_log_lik_best){
      compl_log_lik_best = compl_log_lik_new
      res_best = res_new
    }
  }
  
  # Save results of N_clus_tmp ----------------------------------------------
  res_list[[ind_N_clus]] = res_best
  
}


# Select best cluster number using ICL ------------------------------------
res_select_model = select_model(spks_time_mlist = spks_time_mlist, 
                                stim_onset_vec = stim_onset_vec, 
                                N_component = N_component,
                                key_times_vec = key_times_vec,
                                result_list = res_list)
cand_N_clus_vec = N_clus_min:N_clus_max
N_clus_est = cand_N_clus_vec[res_select_model$id_best_res]
ICL_vec = res_select_model$ICL_vec 
compl_log_lik_vec = res_select_model$compl_log_lik_vec 
log_lik_vec = res_select_model$log_lik_vec
penalty_vec = res_select_model$penalty_vec


# Retrieve estimation results of the best cluster number ------------------
results = list(res_list = res_list, 
               cand_N_clus_vec = cand_N_clus_vec, 
               res_select_model = res_select_model,
               spks_time_mlist = spks_time_mlist,
               stim_onset_time_mat = stim_onset_time_mat,
               id_neuron_id_trial_selected = id_neuron_id_trial_selected)


# Save results ------------------------------------------------------------
top_level_folder = "../Results/Rdata"
setup = 'RDA_v2'
method = 'shape_inv_pp'
default_setting = paste0('Session ', id_session, 
                         ', ', brain_region, 
                         ', scenario_num = ', paste0(scenario_num, collapse = '_'),
                         ', feedback_type = ', feedback_type)
folder_path = paste0(top_level_folder,
                     '/', setup,
                     '/', method, 
                     '/', default_setting)
dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
save(results, file = paste0(folder_path, '/', 'res', '_', now_replicate, '.Rdata'))





