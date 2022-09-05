# rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)

SEED = 5789316
tmp = main_v5_pdf(SEED = SEED,
                  N_node = 100,
                  N_clus = 4,
                  u_1 = 1, u_0 = 1,
                  t_vec = seq(-1, 1, by=0.01),
                  t_vec_extend = seq(-3/2, 1, by=0.01),
                  ### params when N_clus==4:
                  N_spks_total = 50,
                  clus_sep = 1.5,
                  ### Parameters for algorithms
                  freq_trun = 10,
                  fix_timeshift=FALSE,
                  fix_membership = TRUE,
                  save_center_pdf_array=TRUE )

tmp$F_mse_squarel2_ratio
print(tmp$v_mean_sq_err)
