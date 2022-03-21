
# Perform centering step once and clustering step once
cluster_kmeans_pdf = function(spks_time_mlist, stim_onset_vec, reaction_time_vec, 
                              clusters_list, 
                              v_vec,
                              N_component=1,
                              freq_trun=5, 
                              v0 = 0.15, v1 = 0.1,
                              t_vec=seq(0, v0, by=0.01),
                              # Unused arguments
                              order_list=NULL, 
                              opt_radius=max(t_vec)/2,
                              fix_timeshift=FALSE,
                              prob_err_mtplr=0.005,
                              gamma=0.001,
                              ...)
{
  
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  N_clus = length(clusters_list)
  
  # Update time shifts and intensities ------------------------------------------------------
  
  res = est_timeshift_density(spks_time_mlist = spks_time_mlist,
                                stim_onset_vec = stim_onset_vec,
                                reaction_time_vec = reaction_time_vec,
                                clusters_list = clusters_list,
                                v_vec = v_vec,
                                N_component = N_component,
                                freq_trun = freq_trun,
                                v0 = v0, v1 = v1,
                                t_vec = t_vec)
  v_vec = res$v_vec
  center_intensity_array = res$center_intensity_array
  
  # Update clusters--------------------------------------------------------------------------
  clusters = clusters_list

  ### Compute distance between each node and each cluster
  dist_mat = matrix(nrow=N_node, ncol=N_clus)
  for (id_node in 1:N_node) {
    for (id_clus in 1:N_clus) {
      dist_tmp = 0
      f_zi = center_intensity_array[id_clus, 1, ]
      f_zi = f_zi/(sum(f_zi)+.Machine$double.eps)
      for (id_trial in 1:N_trial) {
        spks_time_mi_vec = spks_time_mlist[id_node, id_trial][[1]] - stim_onset_vec[id_trial] - v_vec[id_node]
        spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                    spks_time_mi_vec>=min(t_vec))]
        N_spks_mi = length(spks_time_mi_vec)
        N_mi = hist(spks_time_mi_vec, breaks=c(t_vec[1]-t_unit,t_vec)+(t_unit/2),
                    plot=FALSE)$counts
        N_mi = N_mi/(sum(N_mi)+.Machine$double.eps)
        dist_tmp = dist_tmp + (sum(f_zi^2) - 2*sum(N_mi*f_zi))*N_spks_mi
      }
      dist_mat[id_node, id_clus] = dist_tmp
    }
  }

  ### Update memberships and clusters
  membership = numeric(N_node)
  dist_to_centr_vec = numeric(N_node)
  for (i in 1:N_node) {
    dist_vec = dist_mat[i, ]
    membership[i] = which.min(dist_vec)[1]
    dist_to_centr_vec[i] = min(dist_vec)
  }
  clusters = mem2clus(membership = membership, N_clus_min = N_clus)
  clusters_list = clusters
  l2_loss = sum(dist_to_centr_vec)
  
  
  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list=clusters_list, 
              v_vec=v_vec,
              l2_loss=l2_loss,
              t_vec=t_vec,
              center_intensity_array=center_intensity_array))
}


