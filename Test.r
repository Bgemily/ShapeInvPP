# rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

N_restart_list = list(1,3)
tmp = lapply(N_restart_list, function(N_restart){
  main_v5_pdf(SEED = 831,
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
              freq_trun = 10,
              fix_timeshift=FALSE,
              save_center_pdf_array=TRUE )
})


tmp = plot_intensity_array(center_intensity_array = res$center_density_array_est_permn,
                           center_Nspks_mat = res$center_Nspks_mat_est_permn,
                           center_intensity_array_true = res$network_list$center_density_array_true,
                           clusters_list = res$clusters_list_est_permn, 
                           t_vec = res$data_param$t_vec_extend,
                           u_1 = res$data_param$u_1, u_0 = res$data_param$u_0,
                           N_component = 2)
grid.arrange(tmp$g)


tmp2 = plot_intensity_array(center_intensity_array = res_kernel$center_density_array_est_permn,
                           center_Nspks_mat = res_kernel$center_Nspks_mat_est_permn,
                           center_intensity_array_true = res_kernel$network_list$center_density_array_true,
                           clusters_list = res_kernel$clusters_list_est_permn, 
                           t_vec = res_kernel$data_param$t_vec_extend,
                           u_1 = res_kernel$data_param$u_1, u_0 = res_kernel$data_param$u_0,
                           N_component = 2)
grid.arrange(tmp2$g)



