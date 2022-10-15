#!/usr/bin/env Rscript

# Import all functions ----------------------------------------------------

rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)
library(Matrix)
library(matrixcalc)
library(fdapace)

# Load libraries ----------------------------------------------------------

library(foreach)
library(doParallel)


# User input setup --------------------------------------------------------

N_trial_total = 20
split = 2

N_trial = N_trial_total/split


# Parallel computing setup ------------------------------------------------

N_cores = 10
registerDoParallel(cores=N_cores)


# Run simulations ---------------------------------------------------------

### Initialization ###########
### Parameters' possible values:
N_restart_list = list(1,5,10,20,40)

top_level_folder = "../Results/Rdata"
setup = 'init_Nclus4_v1.3.1'
default_setting = 'N_spks_total=50,N_node=100,clus_sep=2'
### Save estimated densities
for (. in 1:1) {
  method = 'rand_init'
  for (freq_trun in c(10)){
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(lapply(N_restart_list, function(N_restart){
        main_v5_pdf(SEED = SEED,
                    N_node = 100,
                    N_clus = 4,
                    u_1 = 1, u_0 = 1,
                    t_vec = seq(-1, 1, by=0.01),
                    t_vec_extend = seq(-3/2, 1, by=0.01),
                    ### params when N_clus==4:
                    N_spks_total = 50,
                    clus_sep = 2,
                    ### Parameters for algorithms
                    rand_init = TRUE,
                    N_restart = N_restart,
                    freq_trun = freq_trun,
                    fix_timeshift=FALSE,
                    save_center_pdf_array=TRUE )
      }),
      error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "N_restart_various"
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method, '_freqtrun', freq_trun,
                         '/', default_setting,
                         '/', param_name )
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)


  }
}

### NOT save estimated densities
for (. in 1:split) {
  method = 'rand_init'
  for (freq_trun in c(10)){
    results <- foreach(j = 1:N_trial) %dopar% {
      SEED = sample(1:1e7,1)
      tryCatch(lapply(N_restart_list, function(N_restart){
        main_v5_pdf(SEED = SEED,
                    N_node = 100,
                    N_clus = 4,
                    u_1 = 1, u_0 = 1,
                    t_vec = seq(-1, 1, by=0.01),
                    t_vec_extend = seq(-3/2, 1, by=0.01),
                    ### params when N_clus==4:
                    N_spks_total = 50,
                    clus_sep = 2,
                    ### Parameters for algorithms
                    rand_init = TRUE,
                    N_restart = N_restart,
                    freq_trun = freq_trun,
                    fix_timeshift=FALSE,
                    save_center_pdf_array=FALSE )
      }),
      error = function(e) print(paste0("SEED = ", SEED, " : ", e)) )
    }
    param_name = "N_restart_various"
    folder_path = paste0(top_level_folder,
                         '/', setup,
                         '/', method, '_freqtrun', freq_trun,
                         '/', default_setting,
                         '/', param_name )
    dir.create(path = folder_path, recursive = TRUE, showWarnings = FALSE)

    now_trial = format(Sys.time(), "%Y%m%d_%H%M%S")
    save(results, file = paste0(folder_path, '/', 'N_trial', N_trial, '_', now_trial, '.Rdata'))
    rm(results)


  }
}




