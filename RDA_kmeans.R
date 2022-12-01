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
  key_times_vec = c(min(t_vec), min(gocue_time_vec), trial_length)
} else {
  key_times_vec = c(-1.7, 0, 2.5)
}

N_start_kmean = 10
freq_trun = Inf
fix_timeshift = TRUE
fix_comp1_timeshift_only = FALSE
use_true_timeshift = FALSE
rmv_conn_prob = FALSE
v_trialwise_vec_list = list(stim_onset_time_vec - min(stim_onset_time_vec), 
                            gocue_time_vec - min(gocue_time_vec))

set.seed(1)
res_list = list()
for (ind_N_clus in 1:length(N_clus_min:N_clus_max)) {
  res_list[[ind_N_clus]] = list()
  N_clus_tmp = c(N_clus_min:N_clus_max)[ind_N_clus]
  
  ### Apply k-means -----------
  res = get_init(spks_time_mlist = spks_time_mlist, 
                 N_clus = N_clus_tmp,
                 N_component = N_component,
                 t_vec = t_vec, 
                 key_times_vec = key_times_vec,
                 N_start_kmean = N_start_kmean,
                 freq_trun = freq_trun,
                 fix_timeshift = fix_timeshift, 
                 v_trialwise_vec_list = v_trialwise_vec_list,
                 rmv_conn_prob = rmv_conn_prob)
  clusters_list = res$clusters_list

  ### Segment densities to get components
  center_density_array = array(dim = c(N_clus_tmp, N_component, length(t_vec)))
  for (id_component in 1:N_component) {
    if (id_component == 1) {
      center_density_array[ , id_component, ] = res$center_density_mat * matrix(I(t_vec <= key_times_vec[2]), byrow = TRUE, 
                                                                                nrow = N_clus_tmp, ncol = length(t_vec))
    } else {
      center_density_array[ , id_component, ] = res$center_density_mat * matrix(I(t_vec > key_times_vec[2]), byrow = TRUE, 
                                                                                nrow = N_clus_tmp, ncol = length(t_vec))
    }
    
  }
  res$center_density_array = center_density_array
  
  center_Nspks_mat = matrix(0, nrow=N_clus_tmp, ncol=N_component)
  N_spks_mat = apply(spks_time_mlist, c(1,2), function(ls)length(unlist(ls)))
  for (id_clus in 1:N_clus_tmp) {
    N_spks_subjtrial_vec_q = c(N_spks_mat[clusters_list[[id_clus]], ])
    F_hat_q = mean(N_spks_subjtrial_vec_q) 
    density_q_mat = center_density_array[id_clus, , ]
    center_Nspks_mat[id_clus, ] = F_hat_q * rowSums(density_q_mat * (t_vec[2]-t_vec[1]))
  }
  res$center_Nspks_mat = center_Nspks_mat
  
  res$t_vec = t_vec
  
  # Save results of N_clus_tmp ----------------------------------------------
  res_list[[ind_N_clus]] = res
  
}


# Retrieve estimation results of the best cluster number ------------------
cand_N_clus_vec = N_clus_min:N_clus_max
results = list(res_list = res_list, 
               cand_N_clus_vec = cand_N_clus_vec, 
               spks_time_mlist = spks_time_mlist,
               id_neuron_selected = id_neuron_selected,
               id_trial_selected = id_trial_selected)

# Save results ------------------------------------------------------------
top_level_folder = "../Results/Rdata"
setup = 'RDA_v2'
method = ifelse(rmv_conn_prob, yes = 'kmeans_density', no = 'kmeans_intensity')
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





