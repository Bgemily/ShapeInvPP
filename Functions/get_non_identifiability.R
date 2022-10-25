get_non_identifiability = function(spks_time_mlist, 
                                   v_mat_list,
                                   N_component,
                                   t_vec,
                                   freq_trun) 
{
  non_identifiability = 0
  if (freq_trun == Inf) {
    freq_trun = (length(t_vec)-1) %/% 2
  }
  for (freq in 1:freq_trun) {
    # Get X matrix
    N_subj = nrow(spks_time_mlist)
    N_trial = ncol(spks_time_mlist)
    X_mat = matrix(nrow = N_subj * N_trial, ncol = N_component)
    for (id_component in 1:N_component) {
      timeshift_mat = v_mat_list[[id_component]]
      X_vec = exp( -1i*2*pi*freq*c(timeshift_mat)*1/(max(t_vec)-min(t_vec)) )
      X_mat[ , id_component] = X_vec
    }
    
    # Get W matrix
    N_spks_mat = matrix(nrow = N_subj, ncol = N_trial)
    for (id_subj in 1:N_subj) {
      for (id_trial in 1:N_trial) {
        spks_time_vec = spks_time_mlist[id_subj, id_trial][[1]]
        spks_time_vec = spks_time_vec[which(spks_time_vec >= min(t_vec) & spks_time_vec <= max(t_vec))]
        N_spks_mat[id_subj, id_trial] = length(spks_time_vec)
      }
    }
    W_mat = diag(c(N_spks_mat))
    
    # Calculate measurement of non-identifiability: sum of variances of coefficients
    variance_coef_mat = solve(t(Conj(X_mat)) %*% W_mat %*% X_mat)
    variance_coef_mat = Re(variance_coef_mat)
    non_identifiability = non_identifiability + sum(diag(variance_coef_mat)) 
  }
  
  return(non_identifiability)
}
