# Load functions -----
file_path = "../Code/Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(tidyverse)

# Load RDA results -----
id_session = 13
scenario_num = c(1)
feedback_type = c(1)
brain_region = 'midbrain'
gamma = 1e-4
method = paste0('shape_inv_pp_v1_gamma',gamma)

path_res = paste0("../Results/Rdata/RDA_v3.2.7_Nspks_geq1_fill_spks_trials_Nrestart=1_v_subjwise_max=1/", method, "/Session ",id_session,
                  ", ", brain_region, ", scenario_num = ",paste0(scenario_num, collapse = '_'),
                  ", feedback_type = ",paste0(feedback_type, collapse = '_'),
                  "/") 
file = tail(list.files(path = path_res, full.names = TRUE),1)
load(file)

# Extract result -----
res_list = results$res_list
cand_N_clus_vec = results$cand_N_clus_vec
res_select_model = results$res_select_model
spks_time_mlist = results$spks_time_mlist
stim_onset_time_mat = results$stim_onset_time_mat
id_neuron_selected = results$id_neuron_selected
id_trial_selected = results$id_trial_selected

id_res = 2
N_clus = cand_N_clus_vec[ id_res ]
res = res_list[[ id_res ]]
res$spks_time_mlist = spks_time_mlist
res$id_trial_selected = id_trial_selected
res$id_neuron_selected = id_neuron_selected

# Load original data ------
new.path = '../Data/Main/'
dat = readRDS(paste(new.path,"session",id_session,".rds",sep=''))

trial_length = 3.5
trial_start_vec = ((dat$stim_onset-0.1)+(dat$feedback_time+2))/2 - trial_length/2
stim_onset_time_vec = dat$stim_onset[id_trial_selected] - trial_start_vec[id_trial_selected]
reaction_time_vec = dat$reaction_time[id_trial_selected] - trial_start_vec[id_trial_selected]
gocue_time_vec = dat$gocue[id_trial_selected] - trial_start_vec[id_trial_selected]
response_time_vec = dat$response_time[id_trial_selected] - trial_start_vec[id_trial_selected]
feedback_time_vec = dat$feedback_time[id_trial_selected] - trial_start_vec[id_trial_selected]

res$stim_onset_time_vec = stim_onset_time_vec
res$gocue_time_vec = gocue_time_vec

key_times_vec = c(min(stim_onset_time_vec), min(gocue_time_vec), trial_length)

# Permute clusters -----
time_peak_density = c()
for (id_clus in 1:N_clus) {
  time_peak_density[id_clus] = which.max(res$center_density_array[id_clus, 1, ] + res$center_density_array[id_clus, 2, ])    
}
clus_permn_vec = order(rowSums(res$center_Nspks_mat)+res$center_intensity_baseline_vec*max(res$t_vec))

res$center_density_array = res$center_density_array[clus_permn_vec, , , drop = FALSE]
res$clusters_list = res$clusters_list[clus_permn_vec]
res$center_Nspks_mat = res$center_Nspks_mat[clus_permn_vec, , drop = FALSE]
res$center_intensity_array = res$center_intensity_array[clus_permn_vec, , , drop = FALSE]
res$center_intensity_baseline_vec = res$center_intensity_baseline_vec[clus_permn_vec]

# Extract subject-wise time shifts -----
v_trialwise_vec_list = list(res$stim_onset_time_vec - min(res$stim_onset_time_vec), 
                            res$gocue_time_vec - min(res$gocue_time_vec) )

timeshift_subjwise_list = list(res$v_mat_list[[1]][, 1] - v_trialwise_vec_list[[1]][1],
                               res$v_mat_list[[2]][, 1] - v_trialwise_vec_list[[2]][1])

# -----
center_density_array_update = res$center_density_array
center_Nspks_mat_update = res$center_Nspks_mat
center_intensity_array_update = res$center_intensity_array
v_mat_list_update = res$v_mat_list
v_subjwise_list_update = res$v_subjwise_vec_list
center_intensity_baseline_vec_update = res$center_intensity_baseline_vec

# Get next trial's stimuli onset time
res$stim_onset_time_nexttrial_vec = dat$stim_onset[res$id_trial_selected+1] - trial_start_vec[res$id_trial_selected]

# Re-estimate intensities for Cluster 1 & 2
set.seed(0)
train_trial = rbernoulli(n = ncol(res$spks_time_mlist), p = 1)
test_trial = !train_trial

loss_update_vec = c()
l2_loss_all_restarts_all_clus_list = list()
res_all_restarts_all_clus_list = list()
for (id_clus in (1:2)){
  spks_time_mlist_clus_tmp = res$spks_time_mlist[res$clusters_list[[id_clus]], train_trial]
  v_trialwise_vec_list_tmp = lapply(results$v_trialwise_vec_list, function(vec)vec[train_trial])
  
  N_component_tmp = 2
  N_restart = 20
  freq_trun = 10

  # Update intensities with restarts
  set.seed(0)
  res_update_best = NA
  metric_best = Inf
  l2_loss_all_restarts = c()
  timeshift_l1_loss_history = c()
  res_all_restarts = c()
  for (id_restart in 1:N_restart) {
    ### Get initialization -----------
    tmp = get_init(spks_time_mlist = spks_time_mlist_clus_tmp, 
                   N_clus = 1,
                   N_component = N_component_tmp,
                   t_vec = res$t_vec,
                   key_times_vec = key_times_vec,
                   N_start_kmean = 5,
                   freq_trun = freq_trun,
                   add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                   v_trialwise_vec_list = v_trialwise_vec_list_tmp,
                   rmv_conn_prob = TRUE)
    clusters_list_init = tmp$clusters_list
    v_mat_list_init = tmp$v_mat_list
    
    # Apply algorithm ---------
    tmp = apply_asimm(spks_time_mlist = spks_time_mlist_clus_tmp,
                         v_trialwise_vec_list = v_trialwise_vec_list_tmp,
                         clusters_list_init = clusters_list_init,
                         v_mat_list_init = v_mat_list_init,
                         N_component = N_component_tmp, 
                         freq_trun = freq_trun,
                         gamma = 0,
                         t_vec = res$t_vec,
                         key_times_vec = key_times_vec,
                         MaxIter = 10, 
                         conv_thres = 5e-3, 
                         alpha = 0 )
    
    ### Extract loss function value
    l2_loss_new = tail(tmp$loss_history, 2)[1]
    
    l2_loss_all_restarts[id_restart] = l2_loss_new
    res_all_restarts[[id_restart]] = tmp
    
    ### Update best estimation
    if(l2_loss_new < metric_best){
      metric_best = l2_loss_new
      res_update_best = tmp
    }
    
  }
  loss_update_vec[id_clus] = metric_best
  res_all_restarts_all_clus_list[[id_clus]] = res_all_restarts
  l2_loss_all_restarts_all_clus_list[[id_clus]] = l2_loss_all_restarts
  
  if( N_component_tmp > 2 ){
    center_density_array_update_new = array(NA, dim = c(dim(center_density_array_update)[1], 
                                                        N_component_tmp, 
                                                        dim(center_density_array_update)[3]) )
    center_Nspks_mat_update_new = matrix(nrow=nrow(center_Nspks_mat_update), ncol=N_component_tmp)
    center_intensity_array_update_new = array(NA, dim = c(dim(center_intensity_array_update)[1], 
                                                          N_component_tmp, 
                                                          dim(center_intensity_array_update)[3]) )
    center_density_array_update_new[ , 1:2, ] = center_density_array_update
    center_Nspks_mat_update_new[ , 1:2] = center_Nspks_mat_update
    center_intensity_array_update_new[ , 1:2, ] = center_intensity_array_update
    
    center_density_array_update = center_density_array_update_new
    center_Nspks_mat_update = center_Nspks_mat_update_new
    center_intensity_array_update = center_intensity_array_update_new
  }
  center_density_array_update[id_clus, 1:N_component_tmp, ] = res_update_best$center_density_array
  center_Nspks_mat_update[id_clus, 1:N_component_tmp] = res_update_best$center_Nspks_mat
  center_intensity_array_update[id_clus, 1:N_component_tmp, ] = res_update_best$center_intensity_array
  center_intensity_baseline_vec_update[id_clus] = res_update_best$center_intensity_baseline_vec
  
  
  for(id_component in 1:N_component_tmp){
    v_subjwise_list_update[[id_component]][res$clusters_list[[id_clus]]] = res_update_best$v_subjwise_vec_list[[id_component]]
  }
  
}

# Re-estimate intensities for Cluster 3

set.seed(0)
# train_trial = rbernoulli(n = ncol(res$spks_time_mlist), p = 1)
# test_trial = !train_trial
# 
# loss_update_vec = c()
# l2_loss_all_restarts_all_clus_list = list()
# res_all_restarts_all_clus_list = list()
for (id_clus in c(3)){
  spks_time_mlist_clus_tmp = res$spks_time_mlist[res$clusters_list[[id_clus]], train_trial]
  v_trialwise_vec_list_tmp = list(stim_onset_time_vec-min(stim_onset_time_vec))
  key_times_vec_tmp = c(min(stim_onset_time_vec), max(res$t_vec))
  
  N_component_tmp = 1
  N_restart = 1
  freq_trun = 10

  # Update intensities with restarts
  set.seed(0)
  res_update_best = NA
  metric_best = Inf
  l2_loss_all_restarts = c()
  timeshift_l1_loss_history = c()
  res_all_restarts = c()
  for (id_restart in 1:N_restart) {
    ### Get initialization -----------
    tmp = get_init(spks_time_mlist = spks_time_mlist_clus_tmp, 
                   N_clus = 1,
                   N_component = N_component_tmp,
                   t_vec = res$t_vec,
                   key_times_vec = key_times_vec_tmp,
                   N_start_kmean = 5,
                   freq_trun = freq_trun,
                   add_rand_to_init_timeshift = ifelse(id_restart>1, TRUE, FALSE),
                   v_trialwise_vec_list = v_trialwise_vec_list_tmp,
                   v_true_mat_list = NULL,
                   rmv_conn_prob = TRUE)
    clusters_list_init = tmp$clusters_list
    v_mat_list_init = tmp$v_mat_list
    
    # Apply algorithm ---------
    tmp = apply_asimm(spks_time_mlist = spks_time_mlist_clus_tmp,
                         v_trialwise_vec_list = v_trialwise_vec_list_tmp,
                         clusters_list_init = clusters_list_init,
                         v_mat_list_init = v_mat_list_init,
                         N_component = N_component_tmp, 
                         freq_trun = freq_trun,
                         gamma = 0,
                         t_vec = res$t_vec,
                         key_times_vec = key_times_vec_tmp,
                         MaxIter = 10, 
                         conv_thres = 5e-3 )
    
    ### Extract loss function value
    l2_loss_new = tail(tmp$loss_history, 2)[1]
    
    l2_loss_all_restarts[id_restart] = l2_loss_new
    res_all_restarts[[id_restart]] = tmp
    
    ### Update best estimation
    if(l2_loss_new < metric_best){
      metric_best = l2_loss_new
      res_update_best = tmp
    }
    
  }
  loss_update_vec[id_clus] = metric_best
  res_all_restarts_all_clus_list[[id_clus]] = res_all_restarts
  l2_loss_all_restarts_all_clus_list[[id_clus]] = l2_loss_all_restarts
  
  if( N_component_tmp > 2 ){
    center_density_array_update_new = array(NA, dim = c(dim(center_density_array_update)[1], 
                                                        N_component_tmp, 
                                                        dim(center_density_array_update)[3]) )
    center_Nspks_mat_update_new = matrix(nrow=nrow(center_Nspks_mat_update), ncol=N_component_tmp)
    center_intensity_array_update_new = array(NA, dim = c(dim(center_intensity_array_update)[1], 
                                                          N_component_tmp, 
                                                          dim(center_intensity_array_update)[3]) )
    center_density_array_update_new[ , 1:2, ] = center_density_array_update
    center_Nspks_mat_update_new[ , 1:2] = center_Nspks_mat_update
    center_intensity_array_update_new[ , 1:2, ] = center_intensity_array_update
    
    center_density_array_update = center_density_array_update_new
    center_Nspks_mat_update = center_Nspks_mat_update_new
    center_intensity_array_update = center_intensity_array_update_new
  }
  center_density_array_update[id_clus, 1:N_component_tmp, ] = res_update_best$center_density_array
  center_Nspks_mat_update[id_clus, 1:N_component_tmp] = res_update_best$center_Nspks_mat
  center_intensity_array_update[id_clus, 1:N_component_tmp, ] = res_update_best$center_intensity_array
  center_intensity_baseline_vec_update[id_clus] = res_update_best$center_intensity_baseline_vec
  if (N_component_tmp==1){
    center_density_array_update[id_clus, 2, ] = 0
    center_Nspks_mat_update[id_clus, 2] = 0
    center_intensity_array_update[id_clus, 2, ] = 0
  }
  
  for(id_component in 1:N_component_tmp){
    v_subjwise_list_update[[id_component]][res$clusters_list[[id_clus]]] = res_update_best$v_subjwise_vec_list[[id_component]]
  }
  
}


if (FALSE){
  res_reestimate = list(res_all_restarts_all_clus_list = res_all_restarts_all_clus_list,
                        l2_loss_all_restarts_all_clus_list = l2_loss_all_restarts_all_clus_list,
                        center_density_array_update = center_density_array_update,
                        center_Nspks_mat_update = center_Nspks_mat_update,
                        center_intensity_array_update = center_intensity_array_update,
                        center_intensity_baseline_vec_update = center_intensity_baseline_vec_update,
                        v_subjwise_list_update = v_subjwise_list_update,
                        clusters_list = res$clusters_list,
                        t_vec = res$t_vec,
                        spks_time_mlist = res$spks_time_mlist,
                        id_neuron_selected = res$id_neuron_selected,
                        id_trial_selected = res$id_trial_selected,
                        stim_onset_time_vec = res$stim_onset_time_vec,
                        gocue_time_vec = res$gocue_time_vec)
  saveRDS(res_reestimate, file = "../Results/Rdata/RDA_v3.2.7_Nspks_geq1_fill_spks_trials_Nrestart=1_v_subjwise_max=1/shape_inv_pp_v1_gamma1e-04/Session 13, midbrain, scenario_num = 1, feedback_type = 1/reestimate_res_20240131_233945_2.rds")
}
