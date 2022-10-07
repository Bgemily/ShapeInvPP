get_non_identifiability = function(spks_time_mlist, 
                                   v_mat_list,
                                   N_component,
                                   t_vec) 
{
  # Get X matrix
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  X_mat = matrix(nrow = N_node * N_replicate, ncol = N_component)
  for (id_component in 1:N_component) {
    timeshift_mat = v_mat_list[[id_component]]
    X_vec = exp( -1i*2*pi*1*c(timeshift_mat)*1/(max(t_vec)-min(t_vec)) )
    X_mat[ , id_component] = X_vec
  }
  
  # Get W matrix
  N_spks_mat = matrix(nrow = N_node, ncol = N_replicate)
  for (id_node in 1:N_node) {
    for (id_replicate in 1:N_replicate) {
      spks_time_vec = spks_time_mlist[id_node, id_replicate][[1]]
      spks_time_vec = spks_time_vec[which(spks_time_vec >= min(t_vec) & spks_time_vec <= max(t_vec))]
      N_spks_mat[id_node, id_replicate] = length(spks_time_vec)
    }
  }
  W_mat = diag(c(N_spks_mat))
  
  # Calculate measurement of non-identifiability: sum of variances of coefficients
  variance_coef_mat = solve(t(Conj(X_mat)) %*% W_mat %*% X_mat)
  variance_coef_mat = Re(variance_coef_mat)
  non_identifiability = sum(diag(variance_coef_mat))
  
  return(non_identifiability)
}