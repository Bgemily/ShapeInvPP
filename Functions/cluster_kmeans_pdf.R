

# Perform centering step once and clustering step once
cluster_kmeans_pdf = function(spks_time_mlist, stim_onset_vec, reaction_time_vec, 
                              clusters_list, 
                              N_component=2,
                              freq_trun=5, 
                              v0 = 0.2, v1 = 0.1,
                              t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0)+0.01, by=0.01),
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
  
  # Update intensities ------------------------------------------------------------------------------
  center_intensity_array = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                                      stim_onset_vec = stim_onset_vec, 
                                                      reaction_time_vec = reaction_time_vec, 
                                                      clusters_list = clusters_list, 
                                                      N_component = N_component,
                                                      freq_trun = freq_trun, 
                                                      t_vec = t_vec,
                                                      v0 = v0, v1 = v1)
  
  
  # Update clusters--------------------------------------------------------------------------
  clusters = clusters_list

  node_intensity_array = get_node_intensity_array(spks_time_mlist = spks_time_mlist, 
                                                  stim_onset_vec = stim_onset_vec, 
                                                  reaction_time_vec = reaction_time_vec, 
                                                  N_component = N_component,
                                                  freq_trun = freq_trun,
                                                  v0 = v0, v1 = v1)
  
  ### Update clusters. Compare each node with each cluster
  membership = numeric(N_node)
  dist_to_centr_vec = numeric(N_node)
  
  ### Compute distance between each node and each cluster
  dist_mat = matrix(nrow=N_node, ncol=N_clus)
  t_vec_2 = t_vec+max(stim_onset_vec)-max(reaction_time_vec)
  for (id_node in 1:N_node) {
    for (id_clus in 1:N_clus) {
      dist_tmp = 0
      f_zi_vis = center_intensity_array[id_clus, 1, ]
      f_zi_act = center_intensity_array[id_clus, 2, ]
      for (id_trial in 1:N_trial) {
        spks_time_mi_vec = spks_time_mlist[id_node, id_trial][[1]] - stim_onset_vec[id_trial]
        spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec)&spks_time_mi_vec>=min(t_vec))]
        N_mi = hist(spks_time_mi_vec, breaks=c(t_vec[1]-t_unit,t_vec)+(t_unit/2),
                    plot=FALSE)$counts
        id_align = which.min(abs(t_vec_2+reaction_time_vec[id_trial]-stim_onset_vec[id_trial]))[1]
        if(id_align>1){
          f_zi_act_mi = c(f_zi_act[id_align:length(f_zi_act)],
                       rep(0,id_align-1))
        } else{
          f_zi_act_mi = f_zi_act
        }
        f_zi = f_zi_vis + f_zi_act_mi
        f_zi = f_zi/(sum(f_zi)+.Machine$double.eps)
        N_spks = sum(N_mi)
        N_mi = N_mi/(sum(N_mi)+.Machine$double.eps)
        dist_tmp = dist_tmp + (sum(f_zi^2) - 2*sum(N_mi*f_zi))*N_spks
      }
      dist_mat[id_node, id_clus] = dist_tmp
    }
  }
  
  
  
  ### Update memberships and clusters
  for (i in 1:N_node) {
    dist_vec = dist_mat[i, ]
    membership[i] = which.min(dist_vec)
    dist_to_centr_vec[i] = min(dist_vec)
  }
  clusters = mem2clus(membership = membership, N_clus_min = N_clus)
  clusters_list = clusters
  l2_loss = sum(dist_to_centr_vec)
  
  
  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list=clusters_list, 
              l2_loss=l2_loss,
              t_vec=t_vec,
              center_intensity_array=center_intensity_array))
}


