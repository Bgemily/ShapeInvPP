#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)

args <- commandArgs(trailingOnly = TRUE)
gamma <- as.numeric(args[1])
print(gamma)
print(class(gamma))

for (id_session in c(13)) {
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
  
  trial_length = 3.5
  t_vec = seq(0, trial_length, length.out=200)
  

  ### Select neuron-trial pairs
  trial_start_vec = ((dat$stim_onset-0.1)+(dat$feedback_time+2))/2 - trial_length/2
  stim_onset_time_vec = (dat$stim_onset - trial_start_vec)[id_trial_selected]
  gocue_time_vec = (dat$gocue - trial_start_vec)[id_trial_selected]
  reaction_time_vec = (dat$reaction_time - trial_start_vec)[id_trial_selected]
  feedback_time_vec = (dat$feedback_time - trial_start_vec)[id_trial_selected]
  
  spks_time_mlist = matrix(list(), nrow = length(id_neuron_selected), ncol = length(id_trial_selected))
  for (i in 1:length(id_neuron_selected)) {
    for (j in 1:length(id_trial_selected)) {
      id_neuron = id_neuron_selected[i]
      id_trial = id_trial_selected[j]

      spks_vec = dat$spks_pp[id_neuron, id_trial][[1]]
      if ((id_trial+1) <= ncol(dat$spks_pp)) {
        spks_vec_nexttrial = dat$spks_pp[id_neuron, id_trial+1][[1]]
        stim_onset_nexttrial = dat$stim_onset[id_trial+1] 
        spks_vec_nexttrial = spks_vec_nexttrial[spks_vec_nexttrial <= stim_onset_nexttrial]
        spks_vec = c(spks_vec, spks_vec_nexttrial)
      }
      spks_shifted_vec = spks_vec - trial_start_vec[id_trial]

      # Augment spikes if previous trial conclusion occurs after desired observation start time
      if ( dat$trial_intervals[id_trial,1] > trial_start_vec[id_trial] ) {
        spks_shifted_vec_augmented = spks_shifted_vec
        length_for_fill = 0.4
        spks_for_fill = spks_shifted_vec[spks_shifted_vec < (dat$trial_intervals[id_trial,1] - trial_start_vec[id_trial] +length_for_fill)]
        length_to_be_fill =  dat$trial_intervals[id_trial,1] - trial_start_vec[id_trial] 
        N_fill = ceiling(length_to_be_fill / length_for_fill )
        for (id_fill in 1:N_fill) {
          spks_shifted_vec_augmented = c(spks_for_fill - id_fill * length_for_fill, spks_shifted_vec_augmented)
        }
        spks_shifted_vec_augmented = spks_shifted_vec_augmented[spks_shifted_vec_augmented>0]
        spks_shifted_vec = spks_shifted_vec_augmented
      }
      # Augment spikes if next trial's stimuli onset occurs before desired observation end time
      if ( dat$stim_onset[id_trial+1] < (trial_start_vec[id_trial]+trial_length)  ) {
        spks_shifted_vec_augmented = spks_shifted_vec
        length_for_fill = 0.4
        spks_for_fill = spks_shifted_vec[spks_shifted_vec > (dat$stim_onset[id_trial+1]-trial_start_vec[id_trial]-length_for_fill)]
        length_to_be_fill = (trial_start_vec[id_trial]+trial_length) - (dat$stim_onset[id_trial+1])
        N_fill = ceiling(length_to_be_fill / length_for_fill )
        for (id_fill in 1:N_fill) {
          spks_shifted_vec_augmented = c(spks_shifted_vec_augmented, spks_for_fill + id_fill * length_for_fill)
        }
        spks_shifted_vec_augmented = spks_shifted_vec_augmented[spks_shifted_vec_augmented < trial_length]
        spks_shifted_vec = spks_shifted_vec_augmented
      }
      
      trial_end_time = max(t_vec)
      trial_start_time = min(t_vec)
      spks_shifted_vec = spks_shifted_vec[which( (spks_shifted_vec <= trial_end_time) &
                                                   (spks_shifted_vec >= trial_start_time ) )]
      spks_time_mlist[i, j] = list(spks_shifted_vec)
    }
  }
  N_spks_subjwise = rowSums(apply(spks_time_mlist, c(1,2), function(ls)length(unlist(ls))))
  id_neuron_active = which(N_spks_subjwise/length(id_trial_selected) >= 1)
  id_neuron_selected = id_neuron_selected[id_neuron_active]
  N_neuron = length(id_neuron_selected)
  spks_time_mlist = spks_time_mlist[id_neuron_active, ]

  N_neuron = length(id_neuron_selected)
  N_trial = length(id_trial_selected)
  
  # Fit model for various cluster number ------------------------------------
  N_clus_min = 2
  N_clus_max = 5
  cand_N_clus_vec = N_clus_min:N_clus_max
  N_component = 2
  key_times_vec = c(min(stim_onset_time_vec), min(gocue_time_vec), trial_length)
  
  N_start_kmean = 5
  freq_trun = 10
  fix_timeshift = FALSE
  fix_comp1_timeshift_only = FALSE
  v_true_mat_list = NULL
  v_trialwise_vec_list = list(stim_onset_time_vec - min(stim_onset_time_vec), 
                              gocue_time_vec - min(gocue_time_vec))
  N_restart = 1
  MaxIter = 10 
  conv_thres = 5e-3
  
  set.seed(1)
  res_list = list()
  compl_log_lik_vec = c()
  log_lik_vec = c()
  log_lik_tmp_1_vec = c()
  log_lik_tmp_2_vec = c()
  clus_entropy_vec = c()
  L2_loss_part_1_vec = c()
  L2_loss_part_2_vec = c()
  L2_loss_part_1_smoothdensity_vec = c()
  for (ind_N_clus in 1:length(cand_N_clus_vec)) {
    res_list[[ind_N_clus]] = list()
    N_clus_tmp = cand_N_clus_vec[ind_N_clus]
    
    res_best = NA
    l2_loss_best = Inf
    l2_loss_history = c()
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
                     add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                     v_trialwise_vec_list = v_trialwise_vec_list,
                     v_true_mat_list = v_true_mat_list,
                     rmv_conn_prob = TRUE)
      clusters_list_init = res$clusters_list
      v_mat_list_init = res$v_mat_list
      
      
      # Apply algorithm ---------
      time_start = Sys.time()
      res_new = apply_asimm(spks_time_mlist = spks_time_mlist,
                               v_trialwise_vec_list = v_trialwise_vec_list,
                               clusters_list_init = clusters_list_init,
                               v_mat_list_init = v_mat_list_init,
                               N_component = N_component, 
                               freq_trun = freq_trun,
                               gamma = gamma,
                               t_vec=t_vec, 
                               v_subjwise_max = 1.0,
                               key_times_vec = key_times_vec,
                               fix_timeshift = fix_timeshift, 
                               MaxIter = MaxIter, 
                               conv_thres = conv_thres )
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
    tmp = evaluate_model(spks_time_mlist = spks_time_mlist, 
                         v_trialwise_vec_list = v_trialwise_vec_list, 
                         N_component = N_component, 
                         key_times_vec = key_times_vec, 
                         model_fitted_list = res_best,
                         freq_trun = freq_trun)
    log_lik_vec[ind_N_clus] = tmp$log_lik
    log_lik_tmp_1_vec[ind_N_clus] = tmp$log_lik_tmp_1
    log_lik_tmp_2_vec[ind_N_clus] = tmp$log_lik_tmp_2
    clus_entropy_vec[ind_N_clus] = tmp$clus_entropy
    L2_loss_part_1_vec[ind_N_clus] = tmp$L2_loss_part_1
    L2_loss_part_2_vec[ind_N_clus] = tmp$L2_loss_part_2
    L2_loss_part_1_smoothdensity_vec[ind_N_clus] = tmp$L2_loss_part_1_smoothdensity
    compl_log_lik_vec[ind_N_clus] = tmp$compl_log_lik
  }
  
  
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
                 L2_loss_part_1_smoothdensity_vec = L2_loss_part_1_smoothdensity_vec,
                 spks_time_mlist = spks_time_mlist,
                 v_trialwise_vec_list = v_trialwise_vec_list, 
                 id_neuron_selected = id_neuron_selected,
                 id_trial_selected = id_trial_selected,
                 trial_length = trial_length,
                 trial_start_vec = trial_start_vec )
  
  
  # Save results ------------------------------------------------------------
  top_level_folder = "../Results/Rdata"
  setup = 'RDA_v3.2.8_Nspks_geq1_fill_spks_trials_Nrestart=1_v_subjwise_max=1'
  method = paste0('shape_inv_pp_v1_gamma',gamma)
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


