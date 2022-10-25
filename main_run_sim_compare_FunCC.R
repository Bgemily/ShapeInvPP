#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(mclust)
library(combinat)
library(FunCC)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)
library(parallel)

# User input setup --------------------------------------------------------

N_replicate_total = 10
N_split = 1

N_replicate = N_replicate_total/N_split


# Parallel computing setup ------------------------------------------------

N_cores = 10
doParallel::registerDoParallel(cores = N_cores)


# Experiment type ---------------------------------------------------------

select_tuning_parameter = FALSE
test_algorithm_performance = TRUE

# Select tuning parameter -------------------------------------------------
if (select_tuning_parameter) {
  N_component_true = 2
  clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5)
  N_replicate = 10
  
  # Define one replicate
  find_best_delta_tmp = function(SEED, clus_sep){
    # Generate synthetic data
    data_param = list(SEED = SEED,
                      N_subj = 100,
                      N_trial = 1,
                      N_clus = 4, 
                      u_1 = 1, u_0 = 1,
                      t_vec = seq(-1, 1, by = 0.01),
                      t_vec_extend = seq(-1, 1, by = 0.01),
                      N_spks_total = 100,
                      timeshift_subj_max_vec = c(1/4, 1/16),
                      clus_sep = clus_sep)
    data_generated = do.call(what = generate_data, args = data_param)
    
    ### Prepare data for FunCC 
    N_component_true = 2
    key_times_vec = c(-1, 0, 1)
    spks_time_mlist = data_generated$spks_time_mlist
    density_array = array(dim = c(data_param$N_subj, N_component_true, length(data_param$t_vec)))
    for (id_subj in 1:data_param$N_subj){
      res_smooth = density(spks_time_mlist[[id_subj]], bw = 'SJ', 
                           from = min(data_param$t_vec), to = max(data_param$t_vec),
                           n = length(data_param$t_vec))
      for (id_component in 1:N_component_true) {
        y_curr_comp = res_smooth$y * I((data_param$t_vec >= key_times_vec[id_component]) & (data_param$t_vec <= key_times_vec[id_component+1]))
        density_array[id_subj, id_component, ] = y_curr_comp
      }
    }
    
    # Get results for various delta
    res = FunCC_find_best_delta(fun_mat = density_array, 
                                delta_min = 0.001, delta_max = 0.1, num_delta = 20,
                                alpha = 0, beta = 0, theta = 1, 
                                shift.alignement = TRUE, shift.max = 0.25,
                                max.iter.align = 10, number = 10)
    return(list(SEED = SEED, res = res))
  }
  
  # Get results of multiple replicates
  for (clus_sep in clus_sep_list) {
    results <- foreach(j = 1:N_replicate) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(find_best_delta_tmp(SEED = SEED, clus_sep = clus_sep),
               error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    
    # Save results
    top_level_folder = "../Results/Rdata"
    setup = 'FunCC_tuning_param'
    method = 'FunCC_v4'
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method,
                         '/', 'clus_sep',
                         '/', as.character(clus_sep))
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
    now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'tuning_parameter_selection', 
                                '_', now_replicate, '.Rdata'))
    rm(results)
  }
}


# Run simulations ---------------------------------------------------------
test_N_component_1 = TRUE
test_N_component_2 = TRUE
test_N_clus_1 = TRUE
save_res_details = TRUE

top_level_folder = "../Results/Rdata"
setup = 'Compare_methods_v2.8.1'
method = 'funcc'

### Parameters' possible values:
timeshift_subj_max_vec_list = list(c(1/4, 1/16), c(1/4, 1/16)*1.5, c(1/4, 1/16)*2,
                              c(1/4, 1/16)*0.5, c(1/4, 1/16)*0.75, 
                              c(1/4, 1/16)*0.25, c(1/4, 1/16)*0.125,
                              c(1/4, 1/16)*1.25, c(1/4, 1/16)*1.75)
clus_sep_list = list(2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3)
N_subj_list = list(100, 150, 200, 250, 300)

if (test_algorithm_performance) {
  if (test_N_component_2) {
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=1.3,N_comp=2'
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_timeshift_subj_max_vec in 1:length(timeshift_subj_max_vec_list)) {
        timeshift_subj_max_vec = timeshift_subj_max_vec_list[[id_timeshift_subj_max_vec]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                               N_subj = 100,
                               N_clus = 4, 
                               N_component_true = 2,
                               t_vec = seq(-1,1,by=0.01),
                               N_spks_total = 100,
                               timeshift_subj_max_vec = timeshift_subj_max_vec,
                               ### params when N_clus==4:
                               clus_sep = 1.3,
                               ### Parameters for algorithms
                               delta = 0.03, 
                               theta = 1,
                               bw = 'SJ',
                               N_component = 2,
                               key_times_vec = c(-1, 0, 1),
                               save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "timeshift_subj_max_vec"
        param_value = paste0(timeshift_subj_max_vec, collapse = '_')
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
      for (id_clus_sep in 1:length(clus_sep_list)) {
        clus_sep = clus_sep_list[[id_clus_sep]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = 100,
                              N_clus = 4, 
                              N_component_true = 2,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = c(1/4, 1/16)*2,
                              ### params when N_clus==4:
                              clus_sep = clus_sep,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 2,
                              key_times_vec = c(-1, 0, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "clus_sep"
        param_value = clus_sep
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method,
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
      for (id_N_subj in 1:length(N_subj_list)) {
        N_subj = N_subj_list[[id_N_subj]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = N_subj,
                              N_clus = 4, 
                              N_component_true = 2,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = c(1/4, 1/16)*2,
                              ### params when N_clus==4:
                              clus_sep = 1.3,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 2,
                              key_times_vec = c(-1, 0, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_subj"
        param_value = N_subj
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
    }
  }
  
  if (test_N_component_1) {
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=4,clus_sep=1.3,N_comp=1'
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_timeshift_subj_max_vec in 1:length(timeshift_subj_max_vec_list)) {
        timeshift_subj_max_vec = timeshift_subj_max_vec_list[[id_timeshift_subj_max_vec]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = 100,
                              N_clus = 4, 
                              N_component_true = 1,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = timeshift_subj_max_vec,
                              ### params when N_clus==4:
                              clus_sep = 1.3,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 1,
                              key_times_vec = c(-1, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "timeshift_subj_max_vec"
        param_value = paste0(timeshift_subj_max_vec, collapse = '_')
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
      for (id_clus_sep in 1:length(clus_sep_list)) {
        clus_sep = clus_sep_list[[id_clus_sep]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = 100,
                              N_clus = 4, 
                              N_component_true = 1,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = c(1/4, 1/16)*2,
                              ### params when N_clus==4:
                              clus_sep = clus_sep,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 1,
                              key_times_vec = c(-1, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "clus_sep"
        param_value = clus_sep
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method,
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
      for (id_N_subj in 1:length(N_subj_list)) {
        N_subj = N_subj_list[[id_N_subj]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = N_subj,
                              N_clus = 4, 
                              N_component_true = 1,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = c(1/4)*2,
                              ### params when N_clus==4:
                              clus_sep = 1.3,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 1,
                              key_times_vec = c(-1, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_subj"
        param_value = N_subj
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
    }
  }
  
  if (test_N_clus_1) {
    default_setting = 'N_spks_total=100,N_subj=100,N_clus=1,N_comp=2'
    for (id_N_split in 1:N_split) {
      if (save_res_details & (id_N_split == 1)) {
        save_center_pdf_array = TRUE
      } else {
        save_center_pdf_array = FALSE
      }
      for (id_timeshift_subj_max_vec in 1:length(timeshift_subj_max_vec_list)) {
        timeshift_subj_max_vec = timeshift_subj_max_vec_list[[id_timeshift_subj_max_vec]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = 100,
                              N_clus = 1, 
                              N_component_true = 2,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = timeshift_subj_max_vec,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 2,
                              key_times_vec = c(-1, 0, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "timeshift_subj_max_vec"
        param_value = paste0(timeshift_subj_max_vec, collapse = '_')
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
      for (id_N_subj in 1:length(N_subj_list)) {
        N_subj = N_subj_list[[id_N_subj]]
        results <- foreach(j = 1:N_replicate) %dopar% {
          SEED = sample(1:1e7,1)
          tryCatch(main_funcc(SEED = SEED, 
                              N_subj = N_subj,
                              N_clus = 1, 
                              N_component_true = 2,
                              t_vec = seq(-1,1,by=0.01),
                              N_spks_total = 100,
                              timeshift_subj_max_vec = c(1/4, 1/16)*2,
                              ### Parameters for algorithms
                              delta = 0.03, 
                              theta = 1,
                              bw = 'SJ',
                              N_component = 2,
                              key_times_vec = c(-1, 0, 1),
                              save_center_pdf_array = save_center_pdf_array),
                   error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
        }
        param_name = "N_subj"
        param_value = N_subj
        folder_path = paste0(top_level_folder,
                             '/', setup,
                             '/', method, 
                             '/', default_setting,
                             '/', param_name, '/', param_value)
        dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)
        
        now_replicate = format(Sys.time(), "%Y%m%d_%H%M%S")
        save(results, file = paste0(folder_path, '/', 'N_replicate', N_replicate, '_', now_replicate, '.Rdata'))
        rm(results)
      }
    }
  }
}
