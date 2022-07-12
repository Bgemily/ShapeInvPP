
# Perform centering step once and clustering step once
cluster_kmeans_pdf = function(spks_time_mlist, 
                              stim_onset_vec, 
                              reaction_time_vec=NULL, 
                              clusters_list, 
                              v_vec=NULL,
                              v_mat_list=NULL,
                              N_component=1,
                              freq_trun=5, 
                              bw = 0,
                              v0 = 0.15, v1 = 0.1,
                              t_vec=seq(0, v0, by=0.01),
                              fix_timeshift=FALSE,
                              fix_comp1_timeshift_only=FALSE,
                              fix_membership=FALSE,
                              gamma=0.06,
                              # Unused arguments
                              order_list=NULL, 
                              opt_radius=max(t_vec)/2,
                              prob_err_mtplr=0.005,
                              ...)
{
  
  t_unit = t_vec[2]-t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  N_clus = length(clusters_list)
  
  if(N_component==1){
    # Update intensities ------------------------------------------------------
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                     stim_onset_vec = stim_onset_vec, 
                                     reaction_time_vec = reaction_time_vec, 
                                     clusters_list = clusters_list, 
                                     v_vec = v_vec,
                                     N_component = N_component,
                                     freq_trun = freq_trun, 
                                     t_vec = t_vec,
                                     v0 = v0, v1 = v1,
                                     fix_timeshift = fix_timeshift,
                                     align_density = TRUE)
    center_density_array = tmp$center_density_array
    center_Nspks_mat = tmp$center_Nspks_mat
    center_intensity_array = tmp$center_intensity_array
    v_vec = tmp$v_vec
    
    
    #Update time shifts and clusters--------------------------------------------------------------------------
    
    ### Get time shift between each node and each cluster
    tmp = est_timeshift(spks_time_mlist = spks_time_mlist, 
                        stim_onset_vec = stim_onset_vec, 
                        center_density_array = center_density_array,
                        v_vec = v_vec, 
                        N_component = N_component,
                        freq_trun = freq_trun, 
                        v0 = v0, v1 = v1,
                        t_vec = t_vec,
                        fix_timeshift = fix_timeshift,
                        ...)
    v_mat_tmp = tmp$v_mat
    
    ### Get distance between each node and each cluster
    dist_mat = matrix(0, nrow=N_node, ncol=N_clus)
    for (id_clus in 1:N_clus) {
      v_vec_tmp = v_mat_tmp[,id_clus]
      
      N_spks_mat = matrix(nrow=N_node, ncol=N_replicate)
      dNN_array = array(dim = c(N_node,N_replicate,length(t_vec)))
      for (id_node in 1:N_node) {
        for (id_replicate in 1:N_replicate) {
          spks_time_mi_vec = spks_time_mlist[id_node, id_replicate][[1]] - stim_onset_vec[id_replicate] 
          spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                      spks_time_mi_vec>=min(t_vec))]
          spks_time_mi_vec = spks_time_mi_vec - v_vec_tmp[id_node]
          spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                      spks_time_mi_vec>=min(t_vec))]
          N_spks_mi = length(spks_time_mi_vec)
          dN_mi = hist(spks_time_mi_vec, breaks=c(t_vec[1]-t_unit,t_vec)+(t_unit/2),
                       plot=FALSE)$counts
          dNN_mi = dN_mi/(N_spks_mi+.Machine$double.eps)
          
          N_spks_mat[id_node, id_replicate] = N_spks_mi
          dNN_array[id_node, id_replicate, ] = dNN_mi
        }
      }
      
      ### Compute distance between each node and current cluster
      dist_vec = rep(0, N_node)
      for (id_replicate in 1:N_replicate) {
        ### Def of dist_1: N_spks_mi*(sum(density_zi^2)*t_unit - 2*sum(density_zi*dNN_mi))
        dist_1_vec = rep(0, N_node)
        dist_2_vec = rep(0, N_node)
        
        N_spks_vec = N_spks_mat[ , id_replicate]
        dNN_mat = dNN_array[ , id_replicate, ]
        center_density_vec = center_density_array[id_clus, 1, ]
        inner_prod_vec = dNN_mat %*% center_density_vec
        center_density_L2norm = sum(center_density_vec^2)*t_unit
        dist_1_vec = N_spks_vec * center_density_L2norm - N_spks_vec * 2*inner_prod_vec
        
        N_spks_vec = N_spks_mat[ , id_replicate]
        center_Nspks_vec = rep(center_Nspks_mat[id_clus,1], N_node)
        center_Nspks_invrs_vec = center_Nspks_vec^(-1)
        center_Nspks_invrs_vec[which(center_Nspks_vec==0)] = 0
        dist_2_vec = gamma * center_Nspks_invrs_vec * (N_spks_vec - center_Nspks_vec)^2
        
        dist_vec = dist_vec + dist_1_vec + dist_2_vec
      }
      
      dist_mat[,id_clus] = dist_vec
    }
    
    ### Debug
    # dist_to_centr_vec = numeric(N_node)
    # mem = clus2mem(clusters_list)
    # for (i in 1:N_node) {
    #   mem_tmp = mem[i]
    #   dist_to_centr_vec[i] = dist_mat[i, mem_tmp]
    # }
    # l2_loss_2 = sum(dist_to_centr_vec)
    # print(l2_loss_2)
    
    
    ### Choose memberships 
    membership = numeric(N_node)
    dist_to_centr_vec = numeric(N_node)
    for (i in 1:N_node) {
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
    ### Debug:
    # print(l2_loss)
    
    for (id_clus in 1:N_clus) {
      v_vec[clusters_list[[id_clus]]] = v_mat_tmp[clusters_list[[id_clus]], id_clus]
    }
    
    
    
  } else if (N_component==2) {
    # Update intensities ------------------------------------------------------
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist, 
                                     stim_onset_vec = stim_onset_vec, 
                                     reaction_time_vec = reaction_time_vec, 
                                     clusters_list = clusters_list, 
                                     v_mat_list = v_mat_list,
                                     N_component = N_component,
                                     freq_trun = Inf, 
                                     bw = bw,
                                     t_vec = t_vec,
                                     v0 = v0, v1 = v1,
                                     fix_timeshift = fix_timeshift,
                                     align_density = FALSE)
    center_density_array = tmp$center_density_array
    center_Nspks_mat = tmp$center_Nspks_mat
    center_intensity_array = tmp$center_intensity_array
    v_mat_list = tmp$v_mat_list
    
    
    #Update time shifts and clusters--------------------------------------------------------------------------
    
    ### Get time shift between each node and each cluster
    tmp = est_timeshift(spks_time_mlist = spks_time_mlist, 
                        stim_onset_vec = stim_onset_vec, 
                        center_density_array = center_density_array,
                        v_mat_list = v_mat_list,
                        N_component = N_component,
                        freq_trun = freq_trun, 
                        bw = bw,
                        v0 = v0, v1 = v1,
                        t_vec = t_vec,
                        fix_timeshift = fix_timeshift,
                        fix_comp1_timeshift_only = fix_comp1_timeshift_only,
                        ...)
    v_array_list_tmp = tmp$v_array_list
    dist_mat_tmp = tmp$dist_mat
    
    
    ### Get distance between each node and each cluster
    dist_mat = matrix(0, nrow=N_node, ncol=N_clus)
    for (id_clus in 1:N_clus) {
      
      N_spks_mat = matrix(nrow=N_node, ncol=N_replicate)
      for (id_node in 1:N_node) {
        for (id_replicate in 1:N_replicate) {
          spks_time_mi_vec = spks_time_mlist[id_node, id_replicate][[1]] - stim_onset_vec[id_replicate] 
          spks_time_mi_vec = spks_time_mi_vec[which(spks_time_mi_vec<=max(t_vec) &
                                                      spks_time_mi_vec>=min(t_vec))]
          N_spks_mi = length(spks_time_mi_vec)
          N_spks_mat[id_node, id_replicate] = N_spks_mi
        }
      }
      
      ### Compute distance between each node and current cluster
      dist_1_vec = dist_mat_tmp[, id_clus]
      center_Nspks_q_scalar = sum(center_Nspks_mat[id_clus,1:N_component])
      dist_2_vec = gamma * rowSums( (center_Nspks_q_scalar+.Machine$double.eps)^(-1) * 
                                      (N_spks_mat - center_Nspks_q_scalar)^2 )
      dist_vec = dist_1_vec + dist_2_vec
      
      dist_mat[,id_clus] = dist_vec
    }
    

    ### Choose memberships 
    if(fix_membership==TRUE){
      dist_to_centr_vec = numeric(N_node)
      for(id_clus in 1:N_clus){
        dist_to_centr_vec[clusters_list[[id_clus]]] = dist_mat[clusters_list[[id_clus]], id_clus]
      }
      clusters_list = clusters_list
      l2_loss = sum(dist_to_centr_vec)
    } else {
      membership = numeric(N_node)
      dist_to_centr_vec = numeric(N_node)
      for (i in 1:N_node) {
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
    }
    
    ### Debug:
    # print(l2_loss)
    
    for (id_clus in 1:N_clus) {
      v_mat_list[[1]][clusters_list[[id_clus]], 1:N_replicate] = v_array_list_tmp[[1]][clusters_list[[id_clus]], 1:N_replicate, id_clus]
      v_mat_list[[2]][clusters_list[[id_clus]], 1:N_replicate] = v_array_list_tmp[[2]][clusters_list[[id_clus]], 1:N_replicate, id_clus]
    }
    
    ### For each cluster, force the minimum time shifts to be zero
    if ( (!fix_timeshift) & (!fix_comp1_timeshift_only) ) {
      for (id_clus in 1:N_clus) {
        v_mat_list[[1]][clusters_list[[id_clus]], 1:N_replicate] = v_mat_list[[1]][clusters_list[[id_clus]], 1:N_replicate] - 
                                                                    min(v_mat_list[[1]][clusters_list[[id_clus]], 1:N_replicate])
        v_mat_list[[1]] = round(v_mat_list[[1]]/t_unit)*t_unit
      }
    }
    
    
  }
  
  # Output -----------------------------------------------------------------------
  
  return(list(clusters_list=clusters_list, 
              v_vec=v_vec,
              v_mat_list=v_mat_list,
              l2_loss=l2_loss,
              t_vec=t_vec,
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              center_intensity_array=center_intensity_array
  ))
}


