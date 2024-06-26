
# Estimate time shifts and cluster memberships, conditioning on densities
get_timeshift_and_clusters = function(subjtrial_density_smooth_array,
                                      fft_subjtrial_density_unsmooth_array,
                                      N_spks_mat,
                                      v_trialwise_vec_list,
                                      center_density_array,
                                      center_Nspks_mat,
                                      v_mat_list,
                                      freq_trun,
                                      bw,
                                      t_vec,
                                      key_times_vec,
                                      fix_timeshift,
                                      rand_init,
                                      gamma,
                                      v_subjwise_max=NULL)
{
  t_unit = t_vec[2]-t_vec[1]
  N_subj = dim(subjtrial_density_smooth_array)[1]
  N_trial = dim(subjtrial_density_smooth_array)[2]
  N_clus = dim(center_density_array)[1]
  N_component = dim(center_density_array)[2]
  
  ### Get time shift between each subj and each cluster -----
  if (is.null(v_mat_list)) {
    v_mat_list = rep(list(matrix(0, nrow = N_subj, ncol = N_trial)), N_component)
    if (!is.null(v_trialwise_vec_list)) {
      for (id_component in 1:N_component) {
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_mat_list[[id_component]] = v_mat_list[[id_component]]  + matrix(v_trialwise_vec, byrow = TRUE, nrow = N_subj, ncol = N_trial)
      }
    }
  }
  tmp = est_timeshift(subjtrial_density_smooth_array = subjtrial_density_smooth_array,
                      fft_subjtrial_density_unsmooth_array = fft_subjtrial_density_unsmooth_array,
                      N_spks_mat = N_spks_mat,
                      v_trialwise_vec_list = v_trialwise_vec_list,
                      center_density_array = center_density_array,
                      v_mat_list = v_mat_list,
                      N_component = N_component,
                      freq_trun = freq_trun, 
                      bw = bw,
                      t_vec = t_vec,
                      v_subjwise_max = v_subjwise_max,
                      key_times_vec = key_times_vec,
                      fix_timeshift = fix_timeshift,
                      rand_init = rand_init )
  v_array_list_tmp = tmp$v_array_list
  dist_mat_tmp = tmp$dist_mat
  
  
  ### Get distance between each subj and each cluster -----
  dist_mat = matrix(0, nrow=N_subj, ncol=N_clus)
  l2_loss_mat_1 = matrix(0, nrow=N_subj, ncol=N_clus)
  l2_loss_mat_2 = matrix(0, nrow=N_subj, ncol=N_clus)
  for (id_clus in 1:N_clus) {
    dist_1_vec = dist_mat_tmp[, id_clus]
    center_Nspks_q_scalar = sum(center_Nspks_mat[id_clus,1:N_component])
    if (FALSE) {
      dist_2_vec = gamma * rowSums( (center_Nspks_q_scalar+.Machine$double.eps)^(-1) * 
                                      (N_spks_mat - center_Nspks_q_scalar)^2 )
    } else {
      dist_2_vec = gamma * rowSums( (N_spks_mat - center_Nspks_q_scalar)^2 )
    }
    dist_vec = dist_1_vec + dist_2_vec
    dist_mat[,id_clus] = dist_vec
    l2_loss_mat_1[,id_clus] = dist_1_vec
    l2_loss_mat_2[,id_clus] = dist_2_vec / gamma
  }
  
  
  ### Select memberships to minimize total distance -----
  membership = numeric(N_subj)
  dist_to_centr_vec = numeric(N_subj)
  l2_loss_vec_1 = numeric(N_subj)
  l2_loss_vec_2 = numeric(N_subj)
  for (i in 1:N_subj) {
    dist_vec_tmp = dist_mat[i, ]
    mem_tmp = which.min(dist_vec_tmp)
    if(length(mem_tmp)>1){
      mem_tmp = sample(mem_tmp, 1)
    }
    membership[i] = mem_tmp
    dist_to_centr_vec[i] = min(dist_vec_tmp)
    l2_loss_vec_1[i] = l2_loss_mat_1[i, mem_tmp]
    l2_loss_vec_2[i] = l2_loss_mat_2[i, mem_tmp]
  }
  clusters_list = mem2clus(membership = membership, N_clus_min = N_clus)
  l2_loss = sum(dist_to_centr_vec)
  l2_loss_part_1 = sum(l2_loss_vec_1)
  l2_loss_part_2 = sum(l2_loss_vec_2)
  
  ### Extract time shifts based on selected memberships -----
  for (id_clus in 1:N_clus) {
    for (id_component in 1:N_component) {
      v_mat_list[[id_component]][clusters_list[[id_clus]], 1:N_trial] = v_array_list_tmp[[id_component]][clusters_list[[id_clus]], 1:N_trial, id_clus]
    }
  }
  ### For Component 1, force the minimum subj-wise time shifts to be zero
  if ( (!fix_timeshift)  ) {
    id_component = 1
    for (id_clus in 1:N_clus) {
      if (length(clusters_list[[id_clus]])>0) {
        id_trial = 1
        v_subjwise_vec = v_mat_list[[id_component]][clusters_list[[id_clus]], id_trial] - v_trialwise_vec_list[[id_component]][id_trial]
        v_subjwise_vec = v_subjwise_vec - min(v_subjwise_vec)
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_mat_list[[id_component]][clusters_list[[id_clus]], ] = matrix(v_subjwise_vec, nrow = length(clusters_list[[id_clus]]), ncol = N_trial) + matrix(v_trialwise_vec, byrow = TRUE, nrow = length(clusters_list[[id_clus]]), ncol = N_trial)
      }
    }
  }
  
  
  
  
  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list = clusters_list, 
              v_mat_list = v_mat_list,
              l2_loss = l2_loss,
              l2_loss_part_1 = l2_loss_part_1,
              l2_loss_part_2 = l2_loss_part_2,
              dist_to_centr_vec = dist_to_centr_vec))
}


