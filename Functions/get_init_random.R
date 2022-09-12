### Initialize of cluster memberships and time shifts randomly
get_init_random = function(spks_time_mlist, 
                           stim_onset_vec,
                           N_clus,
                           N_component = 2,
                           v0 = 1, v1 = 1,
                           t_vec = seq(-v1, v0, by=0.01), 
                           freq_trun = 10,
                           bw = 0,
                           ### unused arguments:
                           N_restart = 1,
                           fix_timeshift = FALSE, 
                           fix_comp1_timeshift_only = FALSE,
                           use_true_timeshift = TRUE, 
                           jitter_prop_true_timeshift = 0, 
                           fix_membership = FALSE,
                           v_true_mat_list = NULL)
{
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  u_0 = v1
  u_1 = v0
  
  # Initialize clusters -----------------------------------------------------
  membership_vec = sample(1:N_clus, size = N_node, replace = TRUE)
  clusters_list = mem2clus(membership = membership_vec, N_clus_min = N_clus)
  
  
  # Initialize time shifts -------------------------------------------------------
  v_mat_list = list()
  v_mat_list[[1]] = runif(n=N_node*N_replicate, 
                          min = 0,
                          max = u_0/2)  
  v_mat_list[[1]] = matrix(v_mat_list[[1]], nrow = N_node, ncol = N_replicate)
  v_mat_list[[2]] = runif(n = N_node*N_replicate,
                          min = 0,
                          max = u_1-u_0/2 )
  v_mat_list[[2]] = matrix(v_mat_list[[2]], nrow = N_node, ncol = N_replicate)
  for(id_clus in 1:N_clus){
    v_mat_list[[1]][clusters_list[[id_clus]], ] = v_mat_list[[1]][clusters_list[[id_clus]], ] - 
      min(v_mat_list[[1]][clusters_list[[id_clus]], ])
  }
  
  
  

  return(list(v_mat_list=v_mat_list,
              membership_vec=membership_vec, 
              clusters_list=clusters_list))
}
