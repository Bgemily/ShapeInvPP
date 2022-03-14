### Initialize of cluster memberships and time shifts
### Initialize time shifts by earliest edge time.
get_init = function(spks_time_mlist, stim_onset_vec, reaction_time_vec, 
                    N_clus,
                    N_component=2,
                    freq_trun=5, 
                    v0 = 0.2, v1 = 0.1,
                    t_vec=seq(0, max(reaction_time_vec-stim_onset_vec+v0), by=0.01)
                    )
{

  time_unit = t_vec[2] - t_vec[1]
  N_node = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)

  

  # Initialize clusters -----------------------------------------------------
  poinproc_mat = matrix(nrow=N_node, ncol=length(t_vec))
  for(id_node in 1:N_node){
    poinproc_mean = 0
    for (id_trial in 1:N_trial) {
      spks_time_tmp = spks_time_mlist[id_node, id_trial][[1]]-stim_onset_vec[id_trial]
      spks_time_tmp = spks_time_tmp[which(spks_time_tmp<=max(t_vec) & spks_time_tmp>=min(t_vec))]
      poinproc = hist(spks_time_tmp, breaks=t_vec, plot=FALSE)$counts
      poinproc = c(poinproc,0)
      poinproc_mean = poinproc_mean + poinproc
    }
    poinproc_mean = poinproc_mean/N_trial
    poinproc_mat[id_node, ] = poinproc_mean
  }
   
  poinproc_mat_2 = poinproc_mat/(rowSums(poinproc_mat)+.Machine$double.eps)
  membership = cluster::pam(x=poinproc_mat_2, k=N_clus, diss=FALSE, cluster.only=TRUE)

  clusters = mem2clus(membership = membership, N_clus_min = N_clus)
  clusters_list = clusters
  membership_vec = membership
  
  
  
  
  
  return(list(membership_vec=membership_vec, 
              poinproc_mat = poinproc_mat,
              clusters_list=clusters_list))
}

