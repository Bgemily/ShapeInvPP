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
id_session = 13
scenario_num = c(1)
feedback_type = c(1)
brain_region = 'midbrain'

dat = readRDS(paste(new.path, "session",id_session,".rds",sep=''))
id_trial_selected = which((dat$scenario_num %in% scenario_num) & (dat$feedback_type %in% feedback_type))
id_neuron_selected = which(dat$brain_region == brain_region)
if (FALSE) {
  id_trial_selected = sample(id_trial_selected, 3)
  id_neuron_selected = sample(id_neuron_selected, 30)
}
N_neuron = length(id_neuron_selected)
N_trial = length(id_trial_selected)

trial_length = max(dat$trial_intervals[id_trial_selected,2] - dat$stim_onset[id_trial_selected]) + 0.2
if (identical(feedback_type, 1)) {
  t_vec = seq(0, trial_length, length.out=200)
} else {
  t_vec = seq(-2, 3,length.out=200)
}

### Select neuron-trial pairs 
trial_start_vec = dat$trial_intervals[,2] - trial_length
stim_onset_time_vec = (dat$stim_onset - trial_start_vec)[id_trial_selected]
gocue_time_vec = (dat$gocue - trial_start_vec)[id_trial_selected]
reaction_time_vec = (dat$reaction_time - trial_start_vec)[id_trial_selected]
feedback_time_vec = (dat$feedback_time - trial_start_vec)[id_trial_selected]
spks_time_mlist = matrix(list(), nrow = N_neuron, ncol = N_trial)
for (i in 1:N_neuron) {
  for (j in 1:N_trial) {
    id_neuron = id_neuron_selected[i]
    id_trial = id_trial_selected[j]
    
    spks_vec = dat$spks_pp[id_neuron, id_trial][[1]]
    spks_shifted_vec = spks_vec - trial_start_vec[id_trial]
    
    stim_onset_time_shifted = dat$stim_onset[id_trial] - trial_start_vec[id_trial]
    feedback_time_shifted = dat$feedback_time[id_trial] - trial_start_vec[id_trial]
    spks_shifted_vec = spks_shifted_vec[which( (spks_shifted_vec <= min(max(t_vec), Inf+feedback_time_shifted)) &  
                                                 (spks_shifted_vec >= max(min(t_vec), -Inf+stim_onset_time_shifted) ) )]
    
    spks_time_mlist[i, j] = list(spks_shifted_vec)
  }
}
N_spks_subjwise = rowSums(apply(spks_time_mlist, c(1,2), function(ls)length(unlist(ls))))
id_neuron_active = which(N_spks_subjwise/N_trial >= 1)
id_neuron_selected = id_neuron_selected[id_neuron_active]
N_neuron = length(id_neuron_selected)
spks_time_mlist = spks_time_mlist[id_neuron_active, ]


# Fit model for various cluster number ------------------------------------
N_clus_min = 4
N_clus_max = 4
N_component = 2
if (identical(feedback_type, 1)) {
  key_times_vec = c(min(stim_onset_time_vec), min(gocue_time_vec), trial_length)
} else {
  key_times_vec = c(-1.7, 0, 2.5)
}

N_start_kmean = 5
freq_trun = 10
fix_timeshift = FALSE
fix_comp1_timeshift_only = FALSE
v_true_mat_list = NULL
v_trialwise_vec_list = list(stim_onset_time_vec - min(stim_onset_time_vec), 
                            gocue_time_vec - min(gocue_time_vec))
N_restart = 5
MaxIter = 5 
conv_thres = 5e-6 
gamma = 0.01

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
                   N_clus = N_clus_tmp,
                   N_component = N_component,
                   t_vec = t_vec, 
                   key_times_vec = key_times_vec,
                   N_start_kmean = N_start_kmean,
                   freq_trun = freq_trun,
                   fix_timeshift = fix_timeshift, 
                   v_trialwise_vec_list = v_trialwise_vec_list,
                   fix_comp1_timeshift_only = fix_comp1_timeshift_only,
                   v_true_mat_list = v_true_mat_list,
                   rmv_conn_prob = TRUE)
    clusters_list_init = res$clusters_list
    v_mat_list_init = res$v_mat_list
    
    
    # Apply algorithm ---------
    time_start = Sys.time()
    res_new = do_cluster_pdf(spks_time_mlist = spks_time_mlist,
                             v_trialwise_vec_list = v_trialwise_vec_list,
                             clusters_list_init = clusters_list_init,
                             v_mat_list_init = v_mat_list_init,
                             N_component = N_component, 
                             freq_trun = freq_trun,
                             gamma = gamma,
                             t_vec=t_vec, 
                             key_times_vec = key_times_vec,
                             fix_timeshift = fix_timeshift, 
                             MaxIter = MaxIter, 
                             conv_thres = conv_thres, 
                             fix_comp1_timeshift_only = fix_comp1_timeshift_only )
    time_end = Sys.time()
    time_estimation = time_end - time_start
    time_estimation = as.numeric(time_estimation, units='secs')
    
    ### Calculate log likelihood
    res = select_model(spks_time_mlist = spks_time_mlist, 
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
               id_neuron_selected = id_neuron_selected,
               id_trial_selected = id_trial_selected)


# Save results ------------------------------------------------------------
top_level_folder = "../Results/Rdata"
setup = 'RDA_v2'
method = paste0('shape_inv_pp_v5.1_gamma',gamma)
default_setting = paste0('Session ', id_session, 
                         ', ', brain_region, 
                         ', scenario_num = ', paste0(scenario_num, collapse = '_'),
                         ', feedback_type = ', paste0(feedback_type, collapse = '_'))
folder_path = paste0(top_level_folder,
                     '/', setup,
                     '/', method, 
                     '/', default_setting)
dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
save(results, file = paste0(folder_path, '/', 'res', '_', now_replicate, '.Rdata'))





