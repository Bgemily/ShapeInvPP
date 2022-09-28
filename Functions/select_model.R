
select_model = function(spks_time_mlist, 
                        stim_onset_vec, 
                        N_component,
                        key_times_vec,
                        result_list)
{
  ICL_vec = rep(0,length(result_list))
  compl_log_lik_vec = rep(0,length(result_list))
  log_lik_vec = rep(0,length(result_list))
  log_lik_1_vec = rep(0,length(result_list))
  log_lik_2_vec = rep(0,length(result_list))
  clus_entropy_vec = rep(0,length(result_list))
  penalty_vec = rep(0,length(result_list))
  
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  
  for (id_res in 1:length(result_list)) {
    ### Retrieve estimates 
    res_tmp = result_list[[id_res]]
    clusters_list_tmp = res_tmp$clusters_list
    N_clus_tmp = length(clusters_list_tmp)
    clus_size_vec = sapply(clusters_list_tmp, length)
    v_mat_list_tmp = res_tmp$v_mat_list
    center_intensity_array_tmp = res_tmp$center_intensity_array
    center_Nspks_mat_tmp = res_tmp$center_Nspks_mat
    pi_vec = clus_size_vec / sum(clus_size_vec)
    tau_mat = matrix(0, nrow = N_node*N_replicate, ncol = N_clus_tmp)
    ### Let tau_{(i,r),q} == 1  if z_{i}=q
    for (q in 1:N_clus_tmp) {
      if (length(clusters_list_tmp[[q]]) >= 1){
        id_node_vec_tmp = rep((clusters_list_tmp[[q]]-1)*N_replicate, each=N_replicate) + rep(1:N_replicate, times=length(clusters_list_tmp[[q]]))
        tau_mat[id_node_vec_tmp, q] = 1
      }
    } 
    t_vec = res_tmp$t_vec_extend
    t_unit = t_vec[2]-t_vec[1]
    
    # First term of log likelihood: \sum_{i,r} ( -\sum_{q} (F_{q}(T)+G_{q}(T))*tau_{i,r,q} )
    F_q_T = rowSums(center_Nspks_mat_tmp)
    tau_F = tau_mat %*% F_q_T 
    log_lik_tmp_1 = sum(-tau_F)
    
    # Second term of log likelihood: \sum_{q}{ \sum_{i,r}\sum_{t} \log{S^{w_{i,r}}f_{q}(t)+S^{v_{i,r}}g_{q}(t)} *tau_{i,r,q} }
    log_lik_tmp_2 = 0
    for (id_clus in 1:N_clus_tmp) {
      for (id_node in clusters_list_tmp[[id_clus]]){
        for (id_replicate in 1:N_replicate) {
          ### Calculate estimated intensity for current (node, replicate)
          intensity_est = rep(0, length(t_vec))
          for (id_component in 1:N_component) {
            time_shift_tmp = v_mat_list_tmp[[id_component]][id_node, id_replicate]
            n0_shift_tmp = round(time_shift_tmp / t_unit)
            intensity_tmp = center_intensity_array_tmp[id_clus, id_component, ]
            intensity_shifted_curr_comp = c(rep(0, max(0, n0_shift_tmp) ),
                                      head(intensity_tmp, length(t_vec) - max(0, n0_shift_tmp)) )
            intensity_est = intensity_est + intensity_shifted_curr_comp
          }
          log_intensity_est = rep(0, length(t_vec))
          log_intensity_est[which(intensity_est>0)] = log(intensity_est[which(intensity_est>0)])
          
          ### Calculate observed intensity for current (node, replicate)
          event_time_vec_tmp = unlist(spks_time_mlist[id_node, id_replicate])
          if(length(event_time_vec_tmp)==0){
            next
          }
          breaks = c(t_vec[1]-t_unit, t_vec)
          counts_tmp = hist(event_time_vec_tmp, breaks=breaks, plot=FALSE, right=FALSE)$counts
          
          ### Add log-likelihood: \sum_{s} \log[ S^w_{i,r}f(t_{i,r,s})+S^v_{/i,r}g(t_{i,r,s}) ] 
          log_lik_tmp_2 = log_lik_tmp_2 + sum(log_intensity_est*counts_tmp)
        }
      }
    }
    log_lik_tmp = log_lik_tmp_1 + log_lik_tmp_2
    
    ### Fuzzy clustering entropy: \sum_{i}\sum_{q} \tau^{i,q} * \log(\pi_q)
    log_pi_vec = log(pi_vec)
    log_pi_vec[log_pi_vec==-Inf] = 0
    clus_entropy = sum(tau_mat %*% log_pi_vec)
    
    ### Compute penalty
    u_0 = -min(res_tmp$t_vec); u_1 = max(res_tmp$t_vec)
    degr_free_vec = c()
    for (id_component in 1:N_component) {
      ### TODO: Check whether degr_free_curr_comp is reasonable
      degr_free_curr_comp = length(which( (key_times_vec[id_component]<=t_vec) & (t_vec<=u_1-max(v_mat_list_tmp[[id_component]])) ))
      degr_free_vec[id_component] = degr_free_curr_comp
    }
    penalty_tmp = 1/2 * log(N_node*N_replicate) * 
      ( (N_clus_tmp - 1) + N_clus_tmp * sum(degr_free_vec) ) 
      
    ### Compute ICL
    compl_log_lik_tmp = log_lik_tmp + clus_entropy
    ICL_tmp = compl_log_lik_tmp - penalty_tmp 
    
    ### Store ICL, loglik, penalties under current res_tmp
    ICL_vec[id_res] = ICL_tmp
    compl_log_lik_vec[id_res] = compl_log_lik_tmp
    log_lik_vec[id_res] = log_lik_tmp
    log_lik_1_vec[id_res] = log_lik_tmp_1
    log_lik_2_vec[id_res] = log_lik_tmp_2
    clus_entropy_vec[id_res] = clus_entropy
    penalty_vec[id_res] = penalty_tmp
  }
  
  
  ### Retrieve index of best freq_trun and best cluster number
  id_best_res = which(ICL_vec==max(ICL_vec))
  
  # Output ------------------------------------------------------------------
  
  return(list(id_best_res = id_best_res, 
              ICL_vec = ICL_vec, 
              compl_log_lik_vec = compl_log_lik_vec, 
              log_lik_vec = log_lik_vec, 
              log_lik_1_vec = log_lik_1_vec,
              log_lik_2_vec = log_lik_2_vec,
              clus_entropy_vec = clus_entropy_vec,
              penalty_vec = penalty_vec))
}