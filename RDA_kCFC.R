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
  N_neuron = length(id_neuron_selected)
  N_trial = length(id_trial_selected)
  
  if (FALSE) {
    trial_length = max(dat$trial_intervals[id_trial_selected,2] - dat$stim_onset[id_trial_selected]) + 0.2
  } else {
    trial_length = 3.5
  }
  if (identical(feedback_type, 1)) {
    t_vec = seq(0, trial_length, length.out=200)
  } else {
    t_vec = seq(-2, 3,length.out=200)
  }
  
  ### Select neuron-trial pairs
  trial_start_vec = dat$stim_onset - 0.1
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
  id_neuron_active = which(N_spks_subjwise/N_trial >= 1)
  id_neuron_selected = id_neuron_selected[id_neuron_active]
  N_neuron = length(id_neuron_selected)
  spks_time_mlist = spks_time_mlist[id_neuron_active, ]
  
  # Prepare data for FPCA ######
  yList = list()
  tList = list()
  for (id_subj in 1:N_neuron){
    t_unit = t_vec[2]-t_vec[1]
    breaks = c(t_vec[1]-t_unit,t_vec)+t_unit/2
    emp_density_vec = hist(unlist(spks_time_mlist[id_subj, ]), breaks=breaks, plot=FALSE)$counts / t_unit / length(unlist(spks_time_mlist[id_subj, ]))
    yList[[id_subj]] = emp_density_vec
    tList[[id_subj]] = t_vec
    
    if (FALSE) {
      N_spks_curr_subj = mean(sapply(spks_time_mlist[id_subj, ], function(list_tmp)length(unlist(list_tmp))))
      yList[[id_subj]] = N_spks_curr_subj * yList[[id_subj]]
      
    }
    
  }
  
  # Apply kCFC -------
  if (TRUE) {
    N_clus = 3
    N_component = 2
    max_attempts <- 5  # You can adjust the number of maximum attempts as needed
    
    for (attempt in 1:max_attempts) {
      tryCatch({
        # Your original code with kSeed incremented by 1
        kcfcObj <- fdapace::kCFC(y = yList, t = tList, k = N_clus, 
                                 kSeed = 123 + attempt - 1, maxIter = 50,
                                 optnsSW = list(dataType='Dense', maxK=N_component), 
                                 optnsCS = list(dataType='Dense', maxK=N_component))
        
        # If no error, break out of the loop
        break
      }, error = function(e) {
        # Print the error (you can customize this part)
        cat(sprintf("Attempt %d failed with error: %s\n", attempt, conditionMessage(e)))
        
        if (attempt == max_attempts) {
          stop("Maximum number of attempts reached. Exiting.")
        }
        
        # Increment the kSeed for the next attempt
        cat("Retrying with a different kSeed...\n")
      })
    }
    
    
  }
  
  # Retrieve estimation results ------------------
  clusters_list_est = mem2clus(as.numeric(kcfcObj$cluster))
  kcfcObj$clusters_list = clusters_list_est
  
  results = list(kcfcObj = kcfcObj, 
                 spks_time_mlist = spks_time_mlist,
                 id_neuron_selected = id_neuron_selected,
                 id_trial_selected = id_trial_selected)
  
  
  # Save results ------------------------------------------------------------
  top_level_folder = "../Results/Rdata"
  setup = 'RDA_v3.2.7'
  method = 'kCFC_density'
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
  
  # }
}


