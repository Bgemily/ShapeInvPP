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
                    default_timeshift=0
                    )
{

  t_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)


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
  node_intensity_array = get_center_intensity_array(spks_time_mlist = spks_time_mlist,
                                                    stim_onset_vec = stim_onset_vec,
                                                    reaction_time_vec = reaction_time_vec,
                                                    clusters_list = mem2clus(1:N_node),
                                                    v_vec = v_vec,
                                                    N_component = N_component,
                                                    freq_trun = freq_trun,
                                                    v0 = v0, v1 = v1,
                                                    t_vec = t_vec,
                                                    rmv_conn_prob = TRUE)$center_intensity_array
  membership = kmeans(x=node_intensity_array[,1,], centers = N_clus, nstart = 5)$cluster
  
  clusters = mem2clus(membership = membership, N_clus_min = N_clus)
  clusters_list = clusters
  membership_vec = membership
  
  
  
  return(list(membership_vec=membership_vec, 
              v_vec=v_vec,
              clusters_list=clusters_list))
}

