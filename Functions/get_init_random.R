### Initialize of cluster memberships and time shifts
### Initialize time shifts by earliest edge time.
get_init_random = function(spks_time_mlist, 
                    reaction_time_vec=NULL, 
                    N_clus,
                    N_component=1,
                    freq_trun=5, 
                    bw=0,
                    t_vec=seq(0, 1, by=0.01),
                    key_times_vec = c(min(t_vec),0,max(t_vec)),
                    N_start_kmean = 5,
                    fix_timeshift=FALSE,
                    fix_comp1_timeshift_only=FALSE,
                    use_true_timeshift=FALSE, 
                    add_rand_to_init_timeshift=FALSE,
                    v_true_mat_list = NULL,
                    v_trialwise_vec_list = NULL,
                    jitter_prop_true_timeshift=0,
                    rmv_conn_prob=FALSE,
                    default_timeshift=0
)
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  # Initialize time shifts -------------
  if (fix_timeshift) {
    if (use_true_timeshift) {
      v_mat_list = v_true_mat_list
      ### Jitter true time shift
      if(jitter_prop_true_timeshift>0){
        for (id_component in 1:N_component) {
          v_mat_list[[id_component]] = jitter(v_mat_list[[id_component]], amount = jitter_prop_true_timeshift*((-min(t_vec))/2-0))
        }
      }
    } else{
      v_vec = matrix(default_timeshift, nrow = N_subj, ncol = N_trial)
      v_mat_list = rep(list(v_vec), N_component)
    }
  } else{
    ### Use earliest spike times as subject-wise time shifts
    v_subjwise_vec_list = rep(list(rep(0, N_subj)), N_component)
    for (id_subj in 1:N_subj) {
      for (id_component in 1:N_component){
        spks_time_shifted_vec = c()
        for (id_trial in 1:N_trial) {
          spks_time_tmp = spks_time_mlist[id_subj, id_trial][[1]]
          time_start_curr_comp = key_times_vec[id_component] 
          time_end_curr_comp = key_times_vec[id_component + 1]
          spks_time_curr_comp_vec = spks_time_tmp[which(spks_time_tmp >= time_start_curr_comp &
                                                          spks_time_tmp <= time_end_curr_comp)]
          if (length(spks_time_curr_comp_vec) > 0) {
            spks_time_curr_comp_vec = spks_time_curr_comp_vec - v_trialwise_vec_list[[id_component]][id_trial]
          }
          spks_time_shifted_vec = c(spks_time_shifted_vec, spks_time_curr_comp_vec)
        }
        
        if (TRUE) {
          v_subjwise_vec_list[[id_component]][id_subj] = runif(n = 1, min = min(t_vec), 
                                                               max = quantile(spks_time_shifted_vec, 0.0) )
        } else {
          if (length(spks_time_shifted_vec) > 0) {
            v_subjwise_vec_list[[id_component]][id_subj] = quantile(spks_time_shifted_vec, 0.0) 
            if (add_rand_to_init_timeshift){
              v_subjwise_vec_list[[id_component]][id_subj] = runif(n = 1, min = -0.05, max = 0.05) + v_subjwise_vec_list[[id_component]][id_subj]
            }
            v_subjwise_vec_list[[id_component]][id_subj] = v_subjwise_vec_list[[id_component]][id_subj] - key_times_vec[id_component] 
          } else {
            v_subjwise_vec_list[[id_component]][id_subj] = 0
          }
          
        }
      }
    }
    
    ### Force minimum time shifts in each component to be trial-wise time shift
    v_mat_list = rep(list(matrix(0, nrow = N_subj, ncol = N_trial)), N_component)
    for (id_component in 1:N_component){
      v_subjwise_vec = v_subjwise_vec_list[[id_component]]
      v_trialwise_vec = v_trialwise_vec_list[[id_component]]
      v_mat_list[[id_component]] = matrix(v_subjwise_vec, nrow = N_subj, ncol = N_trial) + matrix(v_trialwise_vec, byrow = TRUE, nrow = N_subj, ncol = N_trial)
      v_mat_list[[id_component]] = round(v_mat_list[[id_component]]/t_unit)*t_unit
    }
    
    ### Force time shifts of first component to be truth
    if( (!fix_timeshift) & fix_comp1_timeshift_only ){
      v_mat_list[[1]] = v_true_mat_list[[1]]
    }
  }
  
  
  # Initialize clusters ---------------------------------------
  ### Smooth the point process 
  if (FALSE) {
  subj_intensity_array = array(dim=c(N_subj, 1, length(t_vec)))
  subj_density_array = array(dim=c(N_subj, 1, length(t_vec)))
  subj_Nspks_mat = matrix(nrow=N_subj, ncol=1)
  for (id_subj in 1:N_subj) {
    spks_time_vec_tmp = unlist(spks_time_mlist[id_subj, ])
    tmp = get_smoothed_pp(event_time_vec = spks_time_vec_tmp, 
                          freq_trun = freq_trun, 
                          t_vec = t_vec, 
                          bw=bw)
    intensity_tmp = tmp$intens_vec
    N_spks_q = length(spks_time_vec_tmp)
    if (N_spks_q>0) {
      density_tmp = intensity_tmp / N_spks_q
    } else{
      density_tmp = intensity_tmp*0
    }
    N_spks_subjtrial_vec_tmp = sapply(spks_time_mlist[id_subj, , drop=FALSE], length)
    F_hat_tmp = sqrt( mean(N_spks_subjtrial_vec_tmp^2) )
    intensity_tmp = density_tmp*F_hat_tmp
    
    subj_intensity_array[id_subj,1,] = intensity_tmp
    subj_density_array[id_subj,1,] = density_tmp
    subj_Nspks_mat[id_subj,1] = F_hat_tmp
  }
  
  ### Apply k-means algorithm
  if (rmv_conn_prob){
    membership = kmeans(x=subj_density_array[,1,], centers = N_clus, nstart = N_start_kmean)$cluster
  } else{
    membership = kmeans(x=subj_intensity_array[,1,], centers = N_clus, nstart = N_start_kmean)$cluster
  }
}
  membership = sample(1:N_clus, size = N_subj, replace = TRUE)
  clusters = mem2clus(membership = membership, N_clus_min = N_clus)
  clusters_list = clusters
  membership_vec = membership
  
  ### Force minimum time shifts in each (cluster, component) to be zero
  if (!fix_timeshift) {
    for (id_clus in 1:N_clus) {
      for (id_component in 1:N_component){
        id_trial = 1
        v_subjwise_vec = v_mat_list[[id_component]][clusters_list[[id_clus]], id_trial] - v_trialwise_vec_list[[id_component]][id_trial]
        v_subjwise_vec = v_subjwise_vec - min(v_subjwise_vec)
        v_trialwise_vec = v_trialwise_vec_list[[id_component]]
        v_mat_list[[id_component]][clusters_list[[id_clus]], ] = matrix(v_subjwise_vec, nrow = length(clusters_list[[id_clus]]), ncol = N_trial) + matrix(v_trialwise_vec, byrow = TRUE, nrow = length(clusters_list[[id_clus]]), ncol = N_trial)
        v_mat_list[[id_component]] = round(v_mat_list[[id_component]]/t_unit)*t_unit
      }
    }
    ### Force time shifts of first component to be truth
    if( (!fix_timeshift) & fix_comp1_timeshift_only ){
      v_mat_list[[1]] = v_true_mat_list[[1]]
    }
  }
  
  
  return(list(v_vec=NA,
              v_mat_list=v_mat_list,
              membership_vec=membership_vec, 
              clusters_list=clusters_list ))
}

