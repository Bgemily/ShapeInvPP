
### Generate synthetic networks
generate_data = function(SEED=NULL,
                         N_node=100, N_clus=2, 
                         clus_size_vec = rep(N_node/N_clus, N_clus),
                         N_spks = 20,
                         conn_patt_sep = 1,
                         v0 = 0.5, v1 = 0.5,
                         t_vec=seq(-v1, v0, by=0.01),
                         time_shift_rad=0.1 )
{
  if(!is.null(SEED)) set.seed(SEED)

# Generate cluster memberships ---------------------------------------
  mem_true_vec = rep(1:N_clus, clus_size_vec)
  clus_true_list = mem2clus(mem_true_vec, N_clus_min = N_clus)
  
  
# Generate time shifts ----------------------------------------------------
  v_vec_list = list()
  v_vec_list[[1]] = runif(n=N_node, 
                          min = -time_shift_rad,
                          max = time_shift_rad)
  v_vec_list[[2]] = runif(n = N_node,
                          min = 0,
                          max = 2*time_shift_rad )
  

# Generate spike density functions --------------------------------------------
  mean_mat = matrix(nrow = N_clus, ncol = 2)
  mean_mat[1,1] = -0.05
  mean_mat[1,2] = 0.15
  mean_mat[2,1] = -0.05
  mean_mat[2,2] = 0.15
  sd_mat = matrix(nrow = N_clus, ncol = 2)
  sd_mat[1,1] = 0.1
  sd_mat[1,2] = 0.05
  sd_mat[2,1] = 0.1*10
  sd_mat[2,2] = 0.05
  
  
  center_density_array_true = array(dim = c(N_clus, 2, length(t_vec)))
  if (N_clus==2) {
    center_density_array_true[1,1, ] = dnorm(x = t_vec, 
                                             mean = mean_mat[1,1], 
                                             sd = sd_mat[1,1]) * 1/(1+conn_patt_sep)
    center_density_array_true[1,2, ] = dnorm(x = t_vec, 
                                             mean = mean_mat[1,2], 
                                             sd = sd_mat[1,2]) * conn_patt_sep/(1+conn_patt_sep)
    
    center_density_array_true[2,1, ] = dnorm(x = t_vec, 
                                             mean = mean_mat[2,1], 
                                             sd = sd_mat[2,1]) * 1/(1+conn_patt_sep^(-1))
    center_density_array_true[2,2, ] = dnorm(x = t_vec, 
                                             mean = mean_mat[2,2], 
                                             sd = sd_mat[2,2]) * conn_patt_sep^(-1)/(1+conn_patt_sep^(-1))
    
  } else {
    stop("TODO: specify center_density_array_true")
  }
  
  
# Generate spike times ---------------------------------------------
  stim_onset_vec = 0
  spks_time_mlist = matrix(list(),N_node,1)
  for (id_node in clus_true_list[[1]]) {
    for (j in 1:1) {
      v_tmp_1 = v_vec_list[[1]][id_node]
      v_tmp_2 = v_vec_list[[2]][id_node]
      spks_time_mlist[id_node,j] = list(c( stim_onset_vec[j]+rnorm(# n = max(3,N_spks*1/(1+conn_patt_sep)),
                                                                   n = max(3,rpois(n=1,lambda=N_spks*1/(1+conn_patt_sep))),
                                                                  mean = mean_mat[1,1], sd = sd_mat[1,1])+v_tmp_1,
                                        stim_onset_vec[j]+rnorm(# n = max(3,N_spks*conn_patt_sep/(1+conn_patt_sep)),
                                                                n = max(3,rpois(n=1,lambda=N_spks*conn_patt_sep/(1+conn_patt_sep))),
                                                                mean = mean_mat[1,2], sd = sd_mat[1,2])+v_tmp_2 ))
    }
  }
  for (id_node in clus_true_list[[2]]) {
    for (j in 1:1) {
      v_tmp_1 = v_vec_list[[1]][id_node]
      v_tmp_2 = v_vec_list[[2]][id_node]
      spks_time_mlist[id_node,j] = list(c( stim_onset_vec[j]+rnorm(# n = max(3,N_spks*1/(1+conn_patt_sep^(-1))),
                                                                   n = max(3,rpois(n=1,lambda=N_spks*1/(1+conn_patt_sep^(-1)))),
                                                                   mean = mean_mat[2,1], sd = sd_mat[2,1])+v_tmp_1,
                                        stim_onset_vec[j]+rnorm(# n = max(3,N_spks*conn_patt_sep^(-1)/(1+conn_patt_sep^(-1))),
                                                                n = max(3,rpois(n=1,lambda=N_spks*conn_patt_sep^(-1)/(1+conn_patt_sep^(-1)))),
                                                                mean = mean_mat[2,2], sd = sd_mat[2,2])+v_tmp_2 ))
      # browser()
    }
  }
  
  
  

# Output ------------------------------------------------------------------

  
  return(list(spks_time_mlist=spks_time_mlist, 
              stim_onset_vec=stim_onset_vec,
              mem_true_vec=mem_true_vec, 
              clus_true_list=clus_true_list,
              v_vec_list=v_vec_list, 
              center_density_array_true=center_density_array_true,
              t_vec=t_vec
              ))
  
}

