# rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)


load("/Users/ztzhang/Documents/Academic/SC/ShapeInvPP/Results/Rdata/Nclus1/timeshifts_est_v4.1.2/N_node=10,N_spks_total=100/N_spks_ratio/0.2/N_trial10_20220626_113035.Rdata")

id_res = 8

results[[id_res]]$F_mse_squarel2_ratio_vec
results[[id_res]]$v_mean_sq_err_vec

results[[3]]$data_param
results[[3]]$data_param$SEED -> SEED

res = main_v5_pdf(SEED = SEED,
                  N_node = 100, 
                  N_clus = 1, 
                  u_1 = 1, u_0 = 1, 
                  t_vec = seq(-1, 1, by=0.01),
                  t_vec_extend = seq(-3/2, 1, by=0.01),
                  ### params when N_clus==1:
                  N_spks_total = 1000,
                  N_spks_ratio = 1,
                  sd_shrinkage = 2.5,
                  c_1 = 0.2, delta_1 = 0.4,
                  # c_2 = 0.1, delta_2 = 0.2,
                  # c_3 = 0.2, delta_3 = 0.4,
                  # identical_components = TRUE,
                  ### params when N_clus==2:
                  clus_mixture = 0,
                  ### Parameters for algorithms
                  fix_timeshift=FALSE,
                  use_true_timeshift = TRUE,
                  jitter_prop_true_timeshift = 0,
                  save_center_pdf_array=TRUE )

res$F_mse_squarel2_ratio_vec
res$F_mean_sq_err_vec
res$v_mean_sq_err_vec


tmp = plot_intensity_array(center_intensity_array = res$center_density_array_est_permn,
                           center_Nspks_mat = res$center_Nspks_mat_est_permn,
                           center_intensity_array_true = res$network_list$center_density_array_true,
                           clusters_list = res$clusters_list_est_permn, 
                           t_vec = res$data_param$t_vec_extend,
                           u_1 = res$data_param$u_1, u_0 = res$data_param$u_0,
                           N_component = 2)
grid.arrange(tmp$g)

res$center_density_array_est_permn[1,1,]



