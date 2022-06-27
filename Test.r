# rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)


load("/Users/ztzhang/Documents/Academic/SC/ShapeInvPP/Results/Rdata/Nclus1/timeshifts_jitter_true_v4.1.1/N_node=10,N_spks_total=100/N_spks_ratio/0.2/N_trial10_20220625_223748.Rdata")

id_res = 3

results[[3]]$data_param
results[[3]]$data_param$SEED -> SEED

res = main_v5_pdf(SEED = SEED,
                  N_node = 10, 
                  N_clus = 1, 
                  u_1 = 1, u_0 = 1,
                  ### params when N_clus==1:
                  N_spks_total = 100,
                  N_spks_ratio = 0.2,
                  sd_shrinkage = 1,
                  ### params when N_clus==2:
                  # clus_mixture = 0,
                  ### Parameters for algorithms
                  fix_timeshift=TRUE,
                  use_true_timeshift = TRUE,
                  jitter_prop_true_timeshift = 0.1,
                  save_center_pdf_array=TRUE )
res$F_mse_squarel2_ratio_vec
res$v_mean_sq_err_vec
N_clus = length(res$clusters_list)
N_node = 100

### Match clusters 
center_density_array_est = res$center_density_array
permn =1
center_density_array_est_permn = center_density_array_est[permn, , ,drop=FALSE]

center_density_array_true = res$network_list$center_density_array_true

tmp = plot_intensity_array(center_intensity_array = center_density_array_est_permn,
                           center_Nspks_mat = res$center_Nspks_mat,
                           center_intensity_array_true = center_density_array_true,
                           clusters_list = res$clusters_list[permn], 
                           t_vec = res$data_param$t_vec,
                           u_1 = res$data_param$u_1, u_0 = res$data_param$u_0,
                           N_component = 2)
grid.arrange(tmp$g)


res$N_iteration

