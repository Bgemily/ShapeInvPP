### Initialize of cluster memberships and time shifts
### Initialize time shifts by earliest edge time.
get_init = function(spks_time_mlist, stim_onset_vec, 
                    reaction_time_vec=NULL, 
                    N_clus,
                    N_component=1,
                    freq_trun=5, 
                    v0 = 0.15, v1 = 0.1,
                    t_vec=seq(0, v0, by=0.01),
                    fix_timeshift=FALSE,
                    rmv_conn_prob=FALSE,
                    default_timeshift=0
                    )
{

  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)

  v_vec = v_vec_list = NULL
  if (N_component==1) {
    # Initialize time shifts --------------------------------------------------
    if (fix_timeshift) {
      v_vec = rep(default_timeshift, N_node)
    } else{
      v_vec = rep(0,N_node)
      spks_time_vec_list = list()
      for (id_node in 1:N_node) {
        spks_time_vec = c()
        for (id_trial in 1:N_trial) {
          spks_time_tmp = spks_time_mlist[id_node, id_trial][[1]]-stim_onset_vec[id_trial]
          spks_time_tmp = spks_time_tmp[which(spks_time_tmp<=max(t_vec) & spks_time_tmp>=min(t_vec))]
          spks_time_vec = c(spks_time_vec, spks_time_tmp)
        }
        if(length(spks_time_vec)>0){
          v_vec[id_node] = median(spks_time_vec) - median(t_vec)
        }
        spks_time_vec_list[id_node] = list(spks_time_vec)
      }
      # v_vec = v_vec - min(v_vec) - length(t_vec)%/%10*t_unit
    }
    
    
    # Initialize clusters -----------------------------------------------------
    tmp = get_center_intensity_array(spks_time_mlist = spks_time_mlist,
                                     stim_onset_vec = stim_onset_vec,
                                     reaction_time_vec = reaction_time_vec,
                                     clusters_list = mem2clus(1:N_node),
                                     v_vec = v_vec,
                                     N_component = N_component,
                                     freq_trun = freq_trun,
                                     v0 = v0, v1 = v1,
                                     t_vec = t_vec,
                                     rmv_conn_prob = TRUE)
    
    if (rmv_conn_prob){
      node_intensity_array = tmp$center_density_array
    } else{
      node_intensity_array = tmp$center_intensity_array
    }
    membership = kmeans(x=node_intensity_array[,1,], centers = N_clus, nstart = 5)$cluster
    
    clusters = mem2clus(membership = membership, N_clus_min = N_clus)
    clusters_list = clusters
    membership_vec = membership
    
  } else if (N_component==2) {
    # Initialize time shifts -------------
    if (fix_timeshift) {
      v_vec = rep(default_timeshift, N_node)
      v_vec_list = list(v_vec, v_vec)
    } else{
      v_vec_list = list(rep(0,N_node), rep(0,N_node))
      spks_time_vec_list = list()
      for (id_node in 1:N_node) {
        spks_time_vec_1 = c()
        spks_time_vec_2 = c()
        for (id_trial in 1:N_trial) {
          spks_time_tmp = spks_time_mlist[id_node, id_trial][[1]]-stim_onset_vec[id_trial]
          spks_time_tmp_1 = spks_time_tmp[which(spks_time_tmp<=0 & spks_time_tmp>=min(t_vec))]
          spks_time_vec_1 = c(spks_time_vec_1, spks_time_tmp_1)
          spks_time_tmp_2 = spks_time_tmp[which(spks_time_tmp<=max(t_vec) & spks_time_tmp>=0)]
          spks_time_vec_2 = c(spks_time_vec_2, spks_time_tmp_2)
        }
        if(length(spks_time_vec_1)>0){
          v_vec_list[[1]][id_node] = median(spks_time_vec_1) - min(t_vec)
        }
        if(length(spks_time_vec_2)>0){
          v_vec_list[[2]][id_node] = median(spks_time_vec_2) - 0
        }
      }
      v_vec_list[[1]] = v_vec_list[[1]] - median(v_vec_list[[1]])
      v_vec_list[[2]] = v_vec_list[[2]] - min(v_vec_list[[2]])
      
      v_vec_list[[1]] = round(v_vec_list[[1]]/t_unit)*t_unit
      v_vec_list[[2]] = round(v_vec_list[[2]]/t_unit)*t_unit
      
    }
    
    
    # Initialize clusters -------------
    node_intensity_array = array(dim=c(N_node, 1, length(t_vec)))
    node_density_array = array(dim=c(N_node, 1, length(t_vec)))
    node_Nspks_mat = matrix(nrow=N_node, ncol=1)
    for (id_node in 1:N_node) {
      intensity_tmp = rep(0, length(t_vec))
      density_tmp = rep(0, length(t_vec))
      F_hat_tmp = 0
      
      
      spks_time_vec_tmp = c()
      N_spks_nodetrial_vec_tmp = c()
      for (id_trial in 1:N_trial) {
        spks_time_tmp = unlist(spks_time_mlist[id_node,id_trial]) - stim_onset_vec[id_trial]
        spks_time_tmp_1 = spks_time_tmp[which(spks_time_tmp>=min(t_vec) & 
                                                spks_time_tmp<=0)]
        spks_time_tmp_1_shifted = spks_time_tmp_1 - v_vec_list[[1]][id_node]
        spks_time_tmp_1_shifted = spks_time_tmp_1_shifted[which(spks_time_tmp_1_shifted>=min(t_vec) &
                                                                  spks_time_tmp_1_shifted<=max(t_vec))]
        spks_time_tmp_2 = spks_time_tmp[which(spks_time_tmp>0 & 
                                                spks_time_tmp<=max(t_vec))]
        spks_time_tmp_2_shifted = spks_time_tmp_2 - v_vec_list[[2]][id_node]
        spks_time_tmp_2_shifted = spks_time_tmp_2_shifted[which(spks_time_tmp_2_shifted>=0 &
                                                                  spks_time_tmp_2_shifted<=max(t_vec))]
        
        spks_time_tmp = c(spks_time_tmp_1_shifted, spks_time_tmp_2_shifted)
        
        spks_time_vec_tmp = c(spks_time_vec_tmp, spks_time_tmp)
        N_spks_nodetrial_vec_tmp = c(N_spks_nodetrial_vec_tmp, length(spks_time_tmp))
      }
      
      
      if (length(spks_time_vec_tmp)>0) {
        fft_res = get_adaptive_fft(event_time_vec = spks_time_vec_tmp, 
                                     freq_trun_max = Inf, 
                                     t_vec = t_vec)
        fft_tmp = fft_res$fft_vec_best
      } else{
        fft_tmp = rep(0,length(t_vec))
      }
      intensity_tmp = Re(fft(fft_tmp, inverse = TRUE))
      
      N_spks_q = length(spks_time_vec_tmp)
      if (N_spks_q>0) {
        density_tmp = intensity_tmp / N_spks_q
      } else{
        density_tmp = intensity_tmp*0
      }
      
      F_hat_tmp = sqrt( mean(N_spks_nodetrial_vec_tmp^2) )
      
      intensity_tmp = density_tmp*F_hat_tmp
      
      
      
      node_intensity_array[id_node,1,] = intensity_tmp
      node_density_array[id_node,1,] = density_tmp
      node_Nspks_mat[id_node,1] = F_hat_tmp
      
    }
    
    
    if (rmv_conn_prob){
      membership = kmeans(x=node_density_array[,1,], centers = N_clus, nstart = 5)$cluster
    } else{
      membership = kmeans(x=node_intensity_array[,1,], centers = N_clus, nstart = 5)$cluster
    }
    
    
    clusters = mem2clus(membership = membership, N_clus_min = N_clus)
    clusters_list = clusters
    membership_vec = membership
    
    if (!fix_timeshift) {
      for (id_clus in 1:N_clus) {
        v_vec_list[[1]][clusters_list[[id_clus]]] = v_vec_list[[1]][clusters_list[[id_clus]]] - (quantile(v_vec_list[[1]][clusters_list[[id_clus]]], 0.5)-0)
        v_vec_list[[2]][clusters_list[[id_clus]]] = v_vec_list[[2]][clusters_list[[id_clus]]] - (quantile(v_vec_list[[2]][clusters_list[[id_clus]]], 0.0)-0)
        
        v_vec_list[[1]] = round(v_vec_list[[1]]/t_unit)*t_unit
        v_vec_list[[2]] = round(v_vec_list[[2]]/t_unit)*t_unit
        
      }
    }
    
  }

  
  return(list(membership_vec=membership_vec, 
              v_vec=v_vec,
              v_vec_list=v_vec_list,
              clusters_list=clusters_list))
}

