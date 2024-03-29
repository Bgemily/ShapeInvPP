### Generate synthetic point processes with N_component = 1
generate_data_Ncomp_1 = function(SEED=NULL,
                                 N_subj=100, 
                                 N_trial=1,
                                 N_spks_total = 1000,
                                 N_clus=2, 
                                 clus_size_vec = rep(N_subj/N_clus, N_clus),
                                 t_vec = seq(-1,1,by=0.01),
                                 t_vec_extend = t_vec,
                                 timeshift_subj_max_vec = c(1/8),
                                 timeshift_trial_max = 1/8,
                                 ### params when N_clus==4:
                                 clus_sep = 2,
                                 ### unused params:
                                 N_spks_ratio = 3/2,
                                 sd_shrinkage = 1,
                                 c_1 = 0, delta_1 = 0,
                                 c_2 = 0, delta_2 = 0,
                                 c_3 = 0, delta_3 = 0,
                                 identical_components = FALSE,
                                 clus_mixture = 0)
{
  if(!is.null(SEED)) set.seed(SEED)
  t_unit = t_vec[2]-t_vec[1]
  N_component = 1
  
  # Generate cluster memberships ---------------------------------------
  mem_true_vec = rep(1:N_clus, clus_size_vec)
  clus_true_list = mem2clus(mem_true_vec, N_clus_min = N_clus)
  
  # Generate trial-wise time shifts -----------------------------------
  v_trialwise_vec_list = list()
  for (id_component in 1:1) {
    v_tmp = runif(n = N_trial, min = 0, max = timeshift_trial_max)
    v_tmp = v_tmp - min(v_tmp)
    v_trialwise_vec_list[[id_component]] = v_tmp
  }
  
  # Generate time shifts ----------------------------------------------------
  v_mat_list = list()
  v_mat_list[[1]] = runif(n=N_subj, 
                          min = 0,
                          max = timeshift_subj_max_vec[1])  
  v_mat_list[[1]] = matrix(v_mat_list[[1]], nrow = N_subj, ncol = N_trial)
  for(id_clus in 1:N_clus){
    v_mat_list[[1]][clus_true_list[[id_clus]], ] = v_mat_list[[1]][clus_true_list[[id_clus]], ] - min(v_mat_list[[1]][clus_true_list[[id_clus]], ])
  }
  v_mat_list[[1]] = v_mat_list[[1]] + matrix(v_trialwise_vec_list[[1]], byrow = TRUE, 
                                             nrow = N_subj, ncol = N_trial)
  
  
  # Generate expected number of spikes --------------------------------------------
  center_N_spks_mat = matrix(nrow=N_clus,ncol=N_component)
  if (N_clus==1){
    center_N_spks_mat[1,1] = N_spks_total*1
  } else if (N_clus==4){
    center_N_spks_mat[1,1] = N_spks_total*1
    center_N_spks_mat[2,1] = N_spks_total*1
    center_N_spks_mat[3,1] = N_spks_total*1
    center_N_spks_mat[4,1] = N_spks_total*1
  } 
  
  
  # Generate spike density functions --------------------------------------------  
  center_density_array_true = array(dim = c(N_clus, N_component, length(t_vec_extend)))
  if(N_clus==1){
    s_tmp = 3/10; mu_tmp = -0.5; 
    center_density_array_true[1,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
  } else if (N_clus==4) {
    ## Clus 1
    s_tmp = 1*(1/8)*(1); mu_tmp = -1*(1/2); 
    center_density_array_true[1,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    ## Clus 2
    s_tmp = 1*(1/8)*(sqrt(clus_sep)); mu_tmp = -1*(1/2); 
    center_density_array_true[2,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    ## Clus 3
    s_tmp = 1*(1/8)*(clus_sep^1); mu_tmp = -1*(1/2) 
    center_density_array_true[3,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    ## Clus 4
    s_tmp = 1*(1/8)*(clus_sep^2); mu_tmp = -1*(1/2); 
    center_density_array_true[4,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
  } 
  
  
  # Generate spike intensity functions --------------------------------------------      
  center_intensity_array_true = array(dim = c(N_clus, N_component, length(t_vec_extend)))
  for (id_clus in 1:N_clus){
    center_intensity_array_true[id_clus,1, ] = center_density_array_true[id_clus,1,]*sum(center_N_spks_mat[id_clus,])
  }
  
  
  # Generate spike times ---------------------------------------------
  rejection_sampling = function(density_vec, t_vec, N_sample){
    sample_tmp = c()
    while(length(sample_tmp)<N_sample){
      x_tmp = sample(1:length(t_vec), size=N_sample, replace=TRUE)
      y_tmp = runif(n = N_sample, min=0, max=max(density_vec)+1)
      sample_tmp = c(sample_tmp, t_vec[x_tmp][which(y_tmp <= density_vec[x_tmp])])
    }
    return(sample_tmp[1:N_sample])
  }
  
  
  stim_onset_vec = rep(0, N_trial)
  spks_time_mlist = matrix(list(), nrow = N_subj, ncol = N_trial)
  for (id_clus in 1:N_clus) {
    for (id_subj in clus_true_list[[id_clus]]) {
      for (id_trial in 1:N_trial) {
        v_tmp_1 = v_mat_list[[1]][id_subj, id_trial]
        spks_time_mlist[id_subj, id_trial] = list(c( rejection_sampling(density_vec = center_density_array_true[id_clus,1,], 
                                                                            t_vec = t_vec_extend, 
                                                                            N_sample = rpois(n=1, lambda=center_N_spks_mat[id_clus,1]) )+
                                                           stim_onset_vec[id_trial]+v_tmp_1 ))
        ### Only keep spike times during [min(t_vec), max(t_vec)] 
        spks_time_vec = spks_time_mlist[id_subj,id_trial][[1]]
        spks_time_mlist[id_subj,id_trial][[1]] = spks_time_vec[which(spks_time_vec >= min(t_vec) & 
                                                                           spks_time_vec <= max(t_vec))]
      }
    }
  }
  
  
  
  
  # Output ------------------------------------------------------------------
  
  
  return(list(spks_time_mlist=spks_time_mlist, 
              stim_onset_vec=stim_onset_vec,
              mem_true_vec=mem_true_vec, 
              clus_true_list=clus_true_list,
              v_mat_list=v_mat_list, 
              v_trialwise_vec_list = v_trialwise_vec_list,
              center_density_array_true=center_density_array_true,
              center_N_spks_mat=center_N_spks_mat,
              center_intensity_array_true=center_intensity_array_true,
              t_vec=t_vec,
              t_vec_extend=t_vec_extend
  ))
  
}
