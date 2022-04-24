
select_model = function(spks_time_mlist, 
                        stim_onset_vec, 
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
  N_trial = ncol(spks_time_mlist)
  
  for (id_res in 1:length(result_list)) {
    ### Retrieve estimates 
    res_tmp = result_list[[id_res]]
    clusters_list_tmp = res_tmp$clusters_list
    N_clus_tmp = length(clusters_list_tmp)
    clus_size_vec = sapply(clusters_list_tmp, length)
    v_vec_tmp = res_tmp$v_vec
    center_intensity_array_tmp = res_tmp$center_intensity_array
    center_Nspks_mat_tmp = res_tmp$center_Nspks_mat
    pi_vec = clus_size_vec / sum(clus_size_vec)
    tau_mat = matrix(0, nrow = N_node, ncol = N_clus_tmp)
    for (q in 1:N_clus_tmp) {
      tau_mat[clusters_list_tmp[[q]], q] = 1
    } 
    t_vec = res_tmp$t_vec
    
    # First term of log likelihood: \sum_{i,m} ( -\sum_{q} F_{q,k}(T)*tau_{i,m,q} )
    F_q_T = center_Nspks_mat_tmp[,1]
    tau_F = tau_mat %*% F_q_T 
    log_lik_tmp_1 = sum(-tau_F)
    
    # Second term of log likelihood: \sum_{q}{ \sum_{i,m,r} \log f_{q}[ (t_{i,m,r}-v_{vis,m}) - v_{i,m} ]*tau_{i,m,q} }
    log_lik_tmp_2 = 0
    for (q in 1:N_clus_tmp) {
      log_lik_q_vec = log(center_intensity_array_tmp[q,1,])
      log_lik_q_vec[is.na(log_lik_q_vec)] = 0
      adjst_edge_time_q_vec = c()
      for (id_node in clusters_list_tmp[[q]]) {
        id_trial = 1
        tmp = unlist(spks_time_mlist[id_node, id_trial]) - v_vec_tmp[id_node]
        adjst_edge_time_q_vec = c(adjst_edge_time_q_vec, tmp)
      }
      adjst_edge_time_q_vec = adjst_edge_time_q_vec[which(adjst_edge_time_q_vec<=max(t_vec) &
                                                            adjst_edge_time_q_vec>=min(t_vec))]
      
      if(length(adjst_edge_time_q_vec)==0){
        next
      }

      ### counts: number of spikes whose (adjusted) edge time is close to each breakpoint in t_vec
      breaks = c(2*t_vec[1]-t_vec[2], t_vec)
      counts = hist(adjst_edge_time_q_vec, breaks=breaks, plot=FALSE, right=FALSE)$counts
      
      ### Add compl_log_lik_tmp by \sum_{(i,m):z_{i,m}=q}{\log f_{q}[ (t_{i,m,r}-v_{vis,m}) - v_{i,m} ]}
      ind_tmp = which(counts > 0 & log_lik_q_vec>-Inf)
      log_lik_tmp_2 = log_lik_tmp_2 + sum(log_lik_q_vec[ind_tmp]*counts[ind_tmp])
    }
    
    
    log_lik_tmp = log_lik_tmp_1 + log_lik_tmp_2
    
    ### Fuzzy clustering entropy: \sum_{i}\sum_{q} \tau^{i,q} * \log(\pi_q)
    log_pi_vec = log(pi_vec)
    log_pi_vec[log_pi_vec==-Inf] = 0
    clus_entropy = sum(tau_mat %*% log_pi_vec)
    
    ### Compute penalty
    penalty_tmp = 1/2*(N_clus_tmp-1+2*N_clus_tmp)*log(N_node) 
    
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