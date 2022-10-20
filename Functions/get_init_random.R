### Initialize of cluster memberships and time shifts
### Initialize time shifts by earliest edge time.
get_init_random = function(spks_time_mlist, stim_onset_vec, 
                    reaction_time_vec=NULL, 
                    N_clus,
                    N_component=1,
                    freq_trun=5, 
                    bw=0,
                    v0 = 0.15, v1 = 0.1,
                    t_vec=seq(0, v0, by=0.01),
                    key_times_vec = c(min(t_vec),0,max(t_vec)),
                    N_start_kmean = 5,
                    fix_timeshift=FALSE,
                    fix_comp1_timeshift_only=FALSE,
                    use_true_timeshift=FALSE, 
                    v_true_mat_list = NULL,
                    jitter_prop_true_timeshift=0,
                    fix_membership=FALSE,
                    rmv_conn_prob=FALSE,
                    default_timeshift=0
)
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  
  # Initialize time shifts -------------
  v_mat_list = rep(list(matrix(0, nrow = N_node, ncol = N_replicate)), N_component)
  for (id_node in 1:N_node) {
    for (id_replicate in 1:N_replicate) {
      spks_time_tmp = spks_time_mlist[id_node, id_replicate][[1]]-stim_onset_vec[id_replicate]
      for (id_component in 1:N_component){
        time_start_curr_comp = key_times_vec[id_component]
        time_end_curr_comp = key_times_vec[id_component + 1]
        spks_time_curr_comp_vec = spks_time_tmp[which(spks_time_tmp >= time_start_curr_comp &
                                                        spks_time_tmp <= time_end_curr_comp)]
        if (length(spks_time_curr_comp_vec) > 0) {
          v_mat_list[[id_component]][id_node, id_replicate] = quantile(spks_time_curr_comp_vec, 0.05) 
        }
      }
    }
  }
  
  ### Force minimum time shifts in each component to be zero
  for (id_component in 1:N_component){
    v_mat_list[[id_component]] = v_mat_list[[id_component]] - min(v_mat_list[[id_component]])
    v_mat_list[[id_component]] = round(v_mat_list[[id_component]]/t_unit)*t_unit
  }
  
  
  
  # Initialize densities -------------
  center_density_array = array(dim = c(N_clus, N_component, length(t_vec)) )
  center_Nspks_mat = matrix(nrow = N_clus, ncol = N_component)
  ind_random = sample(1:(N_node*N_replicate), size = N_clus, replace = FALSE)
  for (id_clus in 1:N_clus) {
    ind = ind_random[id_clus]
    spks_time_vec = spks_time_mlist[ind][[1]]
    tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                          freq_trun = freq_trun, 
                          t_vec = t_vec, 
                          bw=bw)
    intensity = tmp$intens_vec
    density = intensity/length(spks_time_vec)
    for (id_component in 1:N_component) {
      density_curr_component = I((t_vec >= key_times_vec[id_component]) &
                                   (t_vec <= key_times_vec[id_component+1])) * density
      center_density_array[id_clus, id_component, ] = density_curr_component
      center_Nspks_mat[id_clus, id_component] = length(spks_time_vec) * sum(density_curr_component * t_unit)
    }
    
  }
  
  
  return(list(v_vec=NA,
              v_mat_list=v_mat_list,
              membership_vec=NA, 
              clusters_list=NA,
              center_density_array = center_density_array,
              center_Nspks_mat = center_Nspks_mat,
              center_intensity_array = NA))
}

