# rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

library(tidyverse)


load("/Users/ztzhang/Documents/Academic/SC/ShapeInvPP/Results/Rdata/Nclus2/timeshifts_est_v2.1/N_node=100/clus_mixture/0.3/N_trial10_20220622_183239.Rdata")

results[[1]]$data_param$SEED -> SEED

tmp = main_v5_pdf(SEED = SEED,
                  N_node = 100, 
                  N_clus = 2, 
                  u_1 = 1, u_0 = 1,
                  ### params when N_clus==2:
                  clus_mixture = 0,
                  ### Parameters for algorithms
                  fix_timeshift=FALSE,
                  save_center_pdf_array=FALSE )

