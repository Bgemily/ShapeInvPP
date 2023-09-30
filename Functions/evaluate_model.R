# Apply fitted model on a given data with spike times and trialwise time shifts 
# and evaluate metrics including: log likelihood, l2 loss of distribution, l2 loss of N_spks, etc.
evaluate_model = function(spks_time_mlist, 
                          v_trialwise_vec_list,
                        N_component,
                        key_times_vec,
                        model_fitted_list,
                        freq_trun = 10)
{
  N_subj = nrow(spks_time_mlist)
  N_trial = ncol(spks_time_mlist)
  
  
  ### Retrieve estimates 
  clusters_list = model_fitted_list$clusters_list
  N_clus = length(clusters_list)
  clus_size_vec = sapply(clusters_list, length)
  v_subjwise_vec_list = model_fitted_list$v_subjwise_vec_list
  v_mat_list = rep(list(matrix(0, nrow = N_subj, ncol = N_trial)), N_component)
  for (id_component in 1:N_component){
    v_subjwise_vec = v_subjwise_vec_list[[id_component]] 
    v_trialwise_vec = v_trialwise_vec_list[[id_component]]
    v_mat_list[[id_component]] = matrix(v_subjwise_vec, nrow = N_subj, ncol = N_trial) + matrix(v_trialwise_vec, byrow = TRUE, nrow = N_subj, ncol = N_trial)
  }
  center_intensity_array = model_fitted_list$center_intensity_array
  center_Nspks_mat = model_fitted_list$center_Nspks_mat
  center_density_array = model_fitted_list$center_density_array
  pi_vec = clus_size_vec / sum(clus_size_vec)
  tau_mat = matrix(0, nrow = N_subj*N_trial, ncol = N_clus)
  ### Let tau_{(i,r),q} == 1  if z_{i}=q
  for (q in 1:N_clus) {
    if (length(clusters_list[[q]]) >= 1){
      id_subj_vec_tmp = rep((clusters_list[[q]]-1)*N_trial, each=N_trial) + rep(1:N_trial, times=length(clusters_list[[q]]))
      tau_mat[id_subj_vec_tmp, q] = 1
    }
  } 
  t_vec = model_fitted_list$t_vec_extend
  t_unit = t_vec[2]-t_vec[1]
  
  # First term of log likelihood: \sum_{i,r} ( -\sum_{q} (Lambda_{q}(T))*tau_{i,r,q} )
  F_q_T = rowSums(center_Nspks_mat)
  tau_F = tau_mat %*% F_q_T 
  log_lik_tmp_1 = sum(-tau_F)
  
  # Second term of log likelihood: \sum_{q}{ \sum_{i,r}\sum_{t} \log{lambda_{i,r}(t)} *tau_{i,r,q} }
  # L2_loss_part_1: sum_{i,r} N_{i,r}(T) * T^{-1} * \| y_{i,r}(t)/N_{i,r}(T) - lambda_{i,r}(t)/Lambda_{i,r}(T) \|^2
  # L2_loss_part_2: sum_{i,r} |Lambda_{i,r}(T)|^{-1} * |N_{i,r}(T)-Lambda_{i,r}(T)|^2
  log_lik_tmp_2 = 0
  L2_loss_part_1 = 0
  L2_loss_part_1_smoothdensity = 0
  L2_loss_part_2 = 0
  for (id_clus in 1:N_clus) {
    for (id_subj in clusters_list[[id_clus]]){
      for (id_trial in 1:N_trial) {
        ### Calculate estimated intensity for current (subj, trial), i.e., lambda_{i,r}
        intensity_est = rep(0, length(t_vec))
        for (id_component in 1:N_component) {
          time_shift_tmp = v_mat_list[[id_component]][id_subj, id_trial]
          n0_shift_tmp = round(time_shift_tmp / t_unit)
          if (n0_shift_tmp > length(t_vec)) {
            n0_shift_tmp = length(t_vec)
          }
          intensity_tmp = center_intensity_array[id_clus, id_component, ]
          intensity_shifted_curr_comp = c(tail(intensity_tmp, max(0, n0_shift_tmp)),
                                          head(intensity_tmp, length(t_vec) - max(0, n0_shift_tmp)) )
          intensity_est = intensity_est + intensity_shifted_curr_comp
        }
        log_intensity_est = rep(0, length(t_vec))
        log_intensity_est[which(intensity_est>0)] = log(intensity_est[which(intensity_est>0)])
        log_intensity_est[which(intensity_est<=0)] = log(min(intensity_est[which(intensity_est>0)])) 
        ### Calculate observed intensity for current (subj, trial)
        event_time_vec_tmp = unlist(spks_time_mlist[id_subj, id_trial])
        if(length(event_time_vec_tmp)==0){
          next
        }
        breaks = c(t_vec[1]-t_unit, t_vec)
        counts_tmp = hist(event_time_vec_tmp, breaks=breaks, plot=FALSE, right=FALSE)$counts
        
        ### Add log-likelihood: \sum_{s} \log{lambda_{i,r}(t_{i,r,s})}
        log_lik_tmp_2 = log_lik_tmp_2 + sum(log_intensity_est*counts_tmp)
        ### Add L2_loss_part_1_tmp: N_{i,r}(T) * T^{-1} * \| y_{i,r}(t)/N_{i,r}(T) - lambda_{i,r}(t)/Lambda_{i,r}(T) \|^2
        intensity_empirical = counts_tmp / t_unit
        L2_loss_part_1_tmp = length(event_time_vec_tmp) * max(t_vec)^{-1} * sum((intensity_empirical/sum(intensity_empirical*t_unit) - intensity_est/sum(intensity_est*t_unit))^2 * t_unit)
        L2_loss_part_1 = L2_loss_part_1 + L2_loss_part_1_tmp
        ### Add L2_loss_part_1_tmp_smoothdensity: N_{i,r}(T) * T^{-1} * \| \tilde{y}_{i,r}(t)/N_{i,r}(T) - lambda_{i,r}(t)/Lambda_{i,r}(T) \|^2
        intensity_empirical_fft = fft(intensity_empirical) / length(intensity_empirical)
        intensity_empirical_fft_trun = c(head(intensity_empirical_fft, freq_trun+1), 
                                     rep(0, length(t_vec)-2*freq_trun-1),
                                     tail(intensity_empirical_fft, freq_trun))
        intensity_empirical_smoothdensity = Re(fft(intensity_empirical_fft_trun, inverse = TRUE))
        L2_loss_part_1_tmp_smoothdensity = length(event_time_vec_tmp) * max(t_vec)^{-1} * sum((intensity_empirical_smoothdensity/sum(intensity_empirical*t_unit) - intensity_est/sum(intensity_est*t_unit))^2 * t_unit)
        L2_loss_part_1_smoothdensity = L2_loss_part_1_smoothdensity + L2_loss_part_1_tmp_smoothdensity
        ### Add L2_loss_part_2_tmp: |Lambda_{i,r}(T)|^{-1} * |N_{i,r}(T)-Lambda_{i,r}(T)|^2
        L2_loss_part_2_tmp = (sum(intensity_est*t_unit)+.Machine$double.eps)^(-1) * (length(event_time_vec_tmp) - sum(intensity_est*t_unit))^2 
        L2_loss_part_2 = L2_loss_part_2 + L2_loss_part_2_tmp
        
      }
    }
  }
  log_lik = log_lik_tmp_1 + log_lik_tmp_2
  
  ### Fuzzy clustering entropy: \sum_{i}\sum_{q} \tau^{i,q} * \log(\pi_q)
  log_pi_vec = log(pi_vec)
  log_pi_vec[log_pi_vec==-Inf] = 0
  clus_entropy = sum(clus_size_vec * log_pi_vec)
  
  ### Complete log likelihood
  compl_log_lik = log_lik + clus_entropy
  
  
  # Output ------------------------------------------------------------------
  
  return(list(compl_log_lik = compl_log_lik,
              log_lik = log_lik,
              log_lik_tmp_1 = log_lik_tmp_1,
              log_lik_tmp_2 = log_lik_tmp_2,
              clus_entropy = clus_entropy,
              L2_loss_part_1 = L2_loss_part_1,
              L2_loss_part_2 = L2_loss_part_2,
              L2_loss_part_1_smoothdensity = L2_loss_part_1_smoothdensity))
}
