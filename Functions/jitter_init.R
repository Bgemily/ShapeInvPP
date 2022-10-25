jitter_init_timeshift = function(v_mat_list_init, 
                                 jitter_level, 
                                 timeshift_min, 
                                 timeshift_max) 
{
  if (jitter_level == 0) {
    return(v_mat_list_init)
  } else {
    N_component = length(v_mat_list_init)
    N_subj = nrow(v_mat_list_init[[1]])
    N_replicate = ncol(v_mat_list_init[[1]])
    for (id_component in 1:N_component) {
      range = max(v_mat_list_init[[id_component]]) - min(v_mat_list_init[[id_component]])
      v_mat_list_init[[id_component]] = runif(n = length(v_mat_list_init[[id_component]]),
                                              min = -jitter_level * range, 
                                              max = jitter_level * range) + v_mat_list_init[[id_component]]
      v_mat_list_init[[id_component]][v_mat_list_init[[id_component]] < timeshift_min] = timeshift_min
      v_mat_list_init[[id_component]][v_mat_list_init[[id_component]] > timeshift_max] = timeshift_max
    }
    return(v_mat_list_init)
  }
}

jitter_init_membership = function(clusters_list_init, jitter_level) {
  N_clus = length(clusters_list_init)
  if ((jitter_level == 0) | (N_clus == 1)) {
    return(clusters_list_init)
  } else {
    membership_init = clus2mem(clusters = clusters_list_init)
    N_subj = length(membership_init)
    ind_switch_mem = sample(x = 1:N_subj, size = round(N_subj * jitter_level), replace = FALSE)
    membership_jitter = membership_init
    for (ind in ind_switch_mem) {
      membership_jitter[ind] = sample(setdiff(1:N_clus, membership_init[ind]), size = 1)
    }
    clusters_list_jitter = mem2clus(membership = membership_jitter, N_clus_min = N_clus)
    return(clusters_list_jitter)
  }
}
