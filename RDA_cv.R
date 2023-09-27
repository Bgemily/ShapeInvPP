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



for (id_session in c(13,28)) {
  # Prepare data ------------------------------------------------------------
  new.path = '../Data/Main/'
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
  N_clus_min = 2
  N_clus_max = 5
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
  N_restart = 1
  MaxIter = 10 
  conv_thres = 5e-6 
  # gamma = 0.007
  N_fold_cv = 5
  
  set.seed(1)
  for (gamma in c(0.01)) {
  res_list = list()
  compl_log_lik_vec = c()
  log_lik_vec = c()
  log_lik_tmp_1_vec = c()
  log_lik_tmp_2_vec = c()
  clus_entropy_vec = c()
  L2_loss_part_1_vec = c()
  L2_loss_part_2_vec = c()
  for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
    
    ### Split the trials into N_fold_cv folds -----------
    id_trial_shuffle = sample(1:N_trial, size = N_trial, replace = FALSE)
    fold_size = N_trial %/% N_fold_cv
    membership_trial = rep(1:N_fold_cv, each = fold_size)
    N_trial_remainder = N_trial - N_fold_cv*fold_size
    if (N_trial_remainder > 0) {
      membership_trial = c(membership_trial, 1:N_trial_remainder)
    }
    
    # Cross-validation --------
    compl_log_lik_cv = 0
    log_lik_cv = 0
    log_lik_tmp_1_cv = 0
    log_lik_tmp_2_cv = 0
    clus_entropy_cv = 0
    L2_loss_part_1_cv = 0
    L2_loss_part_2_cv = 0
    for (id_fold in 1:N_fold_cv) {
      # Get training set and testing set
      id_trial_training = id_trial_shuffle[membership_trial!=id_fold]
      spks_time_mlist_training = spks_time_mlist[ , id_trial_training, drop=FALSE]
      v_trialwise_vec_list_training = lapply(v_trialwise_vec_list, function(vec)vec[id_trial_training])
      id_trial_testing = id_trial_shuffle[membership_trial==id_fold]
      spks_time_mlist_testing = spks_time_mlist[ , id_trial_testing, drop=FALSE]
      v_trialwise_vec_list_testing = lapply(v_trialwise_vec_list, function(vec)vec[id_trial_testing])
      
      # Fit model on the rest (N_fold_cv-1) folds with multiple restarts -----------
      res_best = NA
      l2_loss_best = Inf
      l2_loss_history = c()
      for (id_restart in 1:N_restart) {
        ### Get initialization -----------
        res = get_init(spks_time_mlist = spks_time_mlist_training, 
                       N_clus = N_clus_tmp,
                       N_component = N_component,
                       t_vec = t_vec, 
                       key_times_vec = key_times_vec,
                       N_start_kmean = N_start_kmean,
                       freq_trun = freq_trun,
                       fix_timeshift = fix_timeshift, 
                       add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                       v_trialwise_vec_list = v_trialwise_vec_list_training,
                       fix_comp1_timeshift_only = fix_comp1_timeshift_only,
                       v_true_mat_list = v_true_mat_list,
                       rmv_conn_prob = TRUE)
        clusters_list_init = res$clusters_list
        v_mat_list_init = res$v_mat_list
        
        
        # Apply algorithm ---------
        time_start = Sys.time()
        res_new = do_cluster_pdf(spks_time_mlist = spks_time_mlist_training,
                                 v_trialwise_vec_list = v_trialwise_vec_list_training,
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
        
        ### Extract loss function value
        l2_loss_new = tail(res_new$loss_history, 2)[1]
        
        ### Update best estimation
        if(l2_loss_new < l2_loss_best){
          l2_loss_best = l2_loss_new
          res_best = res_new
        }
        l2_loss_history[id_restart] = l2_loss_new
      }
      res_best$l2_loss_history = l2_loss_history
      res_best$time_estimation = time_estimation
      
      # Save results of N_clus_tmp ----------------------------------------------
      res_list[[ind_N_clus]] = res_best
      
      # Apply fitted model on testing data ----------------------------------------------
      tmp = evaluate_model(spks_time_mlist = spks_time_mlist_testing, 
                           v_trialwise_vec_list = v_trialwise_vec_list_testing, 
                           N_component = N_component, 
                           key_times_vec = key_times_vec, 
                           model_fitted_list = res_best)
      compl_log_lik_cv = compl_log_lik_cv + tmp$compl_log_lik
      log_lik_cv = log_lik_cv + tmp$log_lik
      log_lik_tmp_1_cv = log_lik_tmp_1_cv + tmp$log_lik_tmp_1
      log_lik_tmp_2_cv = log_lik_tmp_2_cv + tmp$log_lik_tmp_2
      clus_entropy_cv = clus_entropy_cv + tmp$clus_entropy
      L2_loss_part_1_cv = L2_loss_part_1_cv + tmp$L2_loss_part_1
      L2_loss_part_2_cv = L2_loss_part_2_cv + tmp$L2_loss_part_2
    }
    compl_log_lik_vec[ind_N_clus] = compl_log_lik_cv
    log_lik_vec[ind_N_clus] = log_lik_cv
    log_lik_tmp_1_vec[ind_N_clus] = log_lik_tmp_1_cv
    log_lik_tmp_2_vec[ind_N_clus] = log_lik_tmp_2_cv
    clus_entropy_vec[ind_N_clus] = clus_entropy_cv
    L2_loss_part_1_vec[ind_N_clus] = L2_loss_part_1_cv
    L2_loss_part_2_vec[ind_N_clus] = L2_loss_part_2_cv
    
  }
  
  
  # Select best cluster number using compl_log_lik_vec ------------------------------------
  id_best_res = which.max(compl_log_lik_vec)
  cand_N_clus_vec = N_clus_min:N_clus_max
  N_clus_est = cand_N_clus_vec[id_best_res]

  
  # Retrieve estimation results of the best cluster number ------------------
  results = list(res_list = res_list, 
                 cand_N_clus_vec = cand_N_clus_vec, 
                 compl_log_lik_vec = compl_log_lik_vec,
                 log_lik_vec = log_lik_vec,
                 log_lik_tmp_1_vec = log_lik_tmp_1_vec,
                 log_lik_tmp_2_vec = log_lik_tmp_2_vec,
                 clus_entropy_vec = clus_entropy_vec,
                 L2_loss_part_1_vec = L2_loss_part_1_vec,
                 L2_loss_part_2_vec = L2_loss_part_2_vec,
                 spks_time_mlist = spks_time_mlist,
                 id_neuron_selected = id_neuron_selected,
                 id_trial_selected = id_trial_selected)
  
  
  # Save results ------------------------------------------------------------
  top_level_folder = "../Results/Rdata"
  setup = 'RDA_v2'
  method = paste0('shape_inv_pp_v6.1_gamma',gamma)
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
  
  }
}


