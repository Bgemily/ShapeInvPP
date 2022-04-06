
# Perform centering step once and clustering step once
cluster_kmeans_pdf = function(spks_time_mlist, stim_onset_vec, reaction_time_vec, 
                              clusters_list, 
                              v_vec,
                              N_component=1,
                              freq_trun=5, 
                              v0 = 0.15, v1 = 0.1,
                              t_vec=seq(0, v0, by=0.01),
                              fix_timeshift=FALSE,
                              gamma=0.06,
                              # Unused arguments
                              order_list=NULL, 
                              opt_radius=max(t_vec)/2,
                              prob_err_mtplr=0.005,
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
                              t_vec = t_vec,
                              fix_timeshift=fix_timeshift,
                              ...)
  # browser()
  v_vec = res$v_vec
  center_density_array = res$center_density_array
  center_Nspks_mat = res$center_Nspks_mat
  center_intensity_array = res$center_intensity_array
  
  # Update clusters--------------------------------------------------------------------------
  clusters = clusters_list
  
  N_spks_mat = matrix(nrow=N_node, ncol=N_trial)
  dNN_array = array(dim = c(N_node,N_trial,length(t_vec)))
  for (id_node in 1:N_node) {
    for (id_trial in 1:N_trial) {
      spks_time_mi_vec = spks_time_mlist[id_node, id_trial][[1]] - stim_onset_vec[id_trial] 
      spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                  spks_time_mi_vec>=min(t_vec))]
      spks_time_mi_vec = spks_time_mi_vec - v_vec[id_node]
      spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                  spks_time_mi_vec>=min(t_vec))]
      N_spks_mi = length(spks_time_mi_vec)
      dN_mi = hist(spks_time_mi_vec, breaks=c(t_vec[1]-t_unit,t_vec)+(t_unit/2),
                   plot=FALSE)$counts
      dNN_mi = dN_mi/(N_spks_mi+.Machine$double.eps)
      
      N_spks_mat[id_node, id_trial] = N_spks_mi
      dNN_array[id_node, id_trial, ] = dNN_mi
    }
  }

  ### Compute distance between each node and each cluster
  dist_mat = matrix(0,nrow=N_node, ncol=N_clus)
  for (id_trial in 1:N_trial) {
    ### Def of dist_1: N_spks_mi*(sum(density_zi^2)*t_unit - 2*sum(density_zi*dNN_mi))
    dist_1_mat = matrix(nrow=N_node, ncol=N_clus)
    N_spks_vec = N_spks_mat[ , id_trial]
    dNN_mat = dNN_array[ , id_trial, ]
    center_density_mat = center_density_array[, 1, ]
    if(N_clus==1){
      center_density_mat = matrix(center_density_mat, nrow=1)
    }
    inner_prod_mat = dNN_mat %*% t(center_density_mat)
    center_density_L2norm_vec = rowSums(center_density_mat^2)*t_unit
    dist_1_mat = N_spks_vec %*% t(center_density_L2norm_vec) - N_spks_vec * 2*inner_prod_mat
    
    N_spks_mat_tmp = matrix(N_spks_mat[ , id_trial], 
                            nrow=N_node, 
                            ncol=N_clus, 
                            byrow = FALSE)
    center_Nspks_mat_tmp = matrix(center_Nspks_mat[ ,1], 
                                  nrow=N_node, 
                                  ncol=N_clus, 
                                  byrow = TRUE)
    center_Nspks_invrs_mat_tmp = center_Nspks_mat_tmp^(-1)
    center_Nspks_invrs_mat_tmp[which(center_Nspks_mat_tmp==0)] = 0
    dist_2_mat = gamma * center_Nspks_invrs_mat_tmp * (N_spks_mat_tmp - center_Nspks_mat_tmp)^2
    
    dist_mat = dist_mat + dist_1_mat + dist_2_mat
  }
  
  ### (Debug) Evaluate l2 loss
  # dist_vec = unlist(sapply(1:length(clusters), function(id_clus)dist_mat[clusters[[id_clus]], id_clus]))
  # l2_loss = sum(dist_vec)
  # browser()
  

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
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              center_intensity_array=center_intensity_array,
              dist_1_mat=dist_1_mat,
              dist_2_mat=dist_2_mat))
}


