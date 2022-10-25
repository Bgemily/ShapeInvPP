
# Estimate time shifts and cluster memberships, conditioning on densities
get_timeshift_and_clusters = function(spks_time_mlist,
                                      stim_onset_vec,
                                      v_trialwise_vec_list = NULL,
                                      center_density_array,
                                      center_Nspks_mat,
                                      v_mat_list,
                                      freq_trun,
                                      bw,
                                      v0, v1,
                                      t_vec,
                                      fix_timeshift,
                                      rand_init,
                                      fix_comp1_timeshift_only,
                                      gamma)
{
  t_unit = t_vec[2]-t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  N_clus = dim(center_density_array)[1]
  N_component = dim(center_density_array)[2]
  
  ### Get time shift between each node and each cluster -----
  if (is.null(v_mat_list)) {
    v_mat_list = rep(list(matrix(0, nrow = N_subj, ncol = N_replicate)), N_component)
    if (!is.null(v_trialwise_vec_list)) {
      for (id_component in 1:N_component) {
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_mat_list[[id_component]] = v_mat_list[[id_component]]  + matrix(v_trialwise_vec, byrow = TRUE, nrow = N_subj, ncol = N_replicate)
      }
    }
  }
  tmp = est_timeshift(spks_time_mlist = spks_time_mlist, 
                      stim_onset_vec = stim_onset_vec, 
                      v_trialwise_vec_list = v_trialwise_vec_list,
                      center_density_array = center_density_array,
                      v_mat_list = v_mat_list,
                      N_component = N_component,
                      freq_trun = freq_trun, 
                      bw = bw,
                      v0 = v0, v1 = v1,
                      t_vec = t_vec,
                      fix_timeshift = fix_timeshift,
                      rand_init = rand_init,
                      fix_comp1_timeshift_only = fix_comp1_timeshift_only )
  v_array_list_tmp = tmp$v_array_list
  dist_mat_tmp = tmp$dist_mat
  
  
  ### Get distance between each node and each cluster -----
  dist_mat = matrix(0, nrow=N_subj, ncol=N_clus)
  for (id_clus in 1:N_clus) {
    N_spks_mat = matrix(nrow=N_subj, ncol=N_replicate)
    for (id_node in 1:N_subj) {
      for (id_replicate in 1:N_replicate) {
        spks_time_mi_vec = spks_time_mlist[id_node, id_replicate][[1]] - stim_onset_vec[id_replicate] 
        spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                    spks_time_mi_vec>=min(t_vec))]
        N_spks_mi = length(spks_time_mi_vec)
        N_spks_mat[id_node, id_replicate] = N_spks_mi
      }
    }
    dist_1_vec = dist_mat_tmp[, id_clus]
    center_Nspks_q_scalar = sum(center_Nspks_mat[id_clus,1:N_component])
    dist_2_vec = gamma * rowSums( (center_Nspks_q_scalar+.Machine$double.eps)^(-1) * 
                                    (N_spks_mat - center_Nspks_q_scalar)^2 )
    dist_vec = dist_1_vec + dist_2_vec
    dist_mat[,id_clus] = dist_vec
  }
  
  
  ### Select memberships to minimize total distance -----
  membership = numeric(N_subj)
  dist_to_centr_vec = numeric(N_subj)
  for (i in 1:N_subj) {
    dist_vec_tmp = dist_mat[i, ]
    mem_tmp = which.min(dist_vec_tmp)
    if(length(mem_tmp)>1){
      mem_tmp = sample(mem_tmp, 1)
    }
    membership[i] = mem_tmp
    dist_to_centr_vec[i] = min(dist_vec_tmp)
  }
  clusters_list = mem2clus(membership = membership, N_clus_min = N_clus)
  l2_loss = sum(dist_to_centr_vec)
  
  
  ### Extract time shifts based on selected memberships -----
  for (id_clus in 1:N_clus) {
    for (id_component in 1:N_component) {
      v_mat_list[[id_component]][clusters_list[[id_clus]], 1:N_replicate] = v_array_list_tmp[[id_component]][clusters_list[[id_clus]], 1:N_replicate, id_clus]
    }
  }
  ### For each cluster, force the minimum time shifts to be trial-wise time shift
  if ( (!fix_timeshift) & (!fix_comp1_timeshift_only) ) {
    id_component = 1
    for (id_clus in 1:N_clus) {
      if (length(clusters_list[[id_clus]])>0) {
        id_replicate = 1
        v_subjwise_vec = v_mat_list[[id_component]][clusters_list[[id_clus]], id_replicate] - v_trialwise_vec_list[[id_component]][id_replicate]
        v_subjwise_vec = v_subjwise_vec - min(v_subjwise_vec)
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_mat_list[[id_component]][clusters_list[[id_clus]], ] = matrix(v_subjwise_vec, nrow = length(clusters_list[[id_clus]]), ncol = N_replicate) + matrix(v_trialwise_vec, byrow = TRUE, nrow = length(clusters_list[[id_clus]]), ncol = N_replicate)
      }
    }
  }
  
  
  
  
  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list = clusters_list, 
              v_mat_list = v_mat_list,
              l2_loss = l2_loss  ))
}


