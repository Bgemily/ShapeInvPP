### Generate synthetic point processes
generate_data = function(SEED=NULL,
                         N_subj=100, 
                         N_trial=1,
                         N_clus=2, 
                         clus_size_vec = rep(N_subj/N_clus, N_clus),
                         t_vec = seq(-1,1,by=0.01),
                         key_times_vec = c(-1,0,1),
                         N_spks_total = 1000,
                         timeshift_subj_max_vec = c(1/8, 1/32),
                         timeshift_trial_max = 1/8,
                         ### params when N_clus==4:
                         clus_sep = 1)
{
  if(!is.null(SEED)) set.seed(SEED)
  t_unit = t_vec[2]-t_vec[1]
  
  # Generate cluster memberships ---------------------------------------
  mem_true_vec = rep(1:N_clus, clus_size_vec)
  clus_true_list = mem2clus(mem_true_vec, N_clus_min = N_clus)
  
  # Generate trial-wise time shifts -----------------------------------
  v_trialwise_vec_list = list()
  for (id_component in 1:2) {
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
  v_mat_list[[2]] = runif(n = N_subj,
                          min = 0,
                          max = timeshift_subj_max_vec[2] )
  v_mat_list[[2]] = matrix(v_mat_list[[2]], nrow = N_subj, ncol = N_trial)
  for(id_clus in 1:N_clus){
    v_mat_list[[1]][clus_true_list[[id_clus]], ] = v_mat_list[[1]][clus_true_list[[id_clus]], ] - min(v_mat_list[[1]][clus_true_list[[id_clus]], ])
    v_mat_list[[2]][clus_true_list[[id_clus]], ] = v_mat_list[[2]][clus_true_list[[id_clus]], ] - min(v_mat_list[[2]][clus_true_list[[id_clus]], ])
  }
  
  v_mat_list[[1]] = v_mat_list[[1]] + matrix(v_trialwise_vec_list[[1]], byrow = TRUE, 
                                             nrow = N_subj, ncol = N_trial)
  v_mat_list[[2]] = v_mat_list[[2]] + matrix(v_trialwise_vec_list[[2]], byrow = TRUE, 
                                             nrow = N_subj, ncol = N_trial)
  
  
  # Generate expected number of spikes --------------------------------------------
  center_N_spks_mat = matrix(nrow=N_clus,ncol=2)
  if (N_clus==1){
    center_N_spks_mat[1,1] = N_spks_total*0.7*0.5
    center_N_spks_mat[1,2] = N_spks_total*0.7*0.5
  } else if (N_clus==4){
    center_N_spks_mat[1,1] = N_spks_total*0.7*(0.5)
    center_N_spks_mat[1,2] = N_spks_total*0.7*(0.5)
    center_N_spks_mat[2,1] = N_spks_total*0.8*(0.5-max(0,sqrt(clus_sep-0.5)/(2*sqrt(0.5))))
    center_N_spks_mat[2,2] = N_spks_total*0.8*(0.5+max(0,sqrt(clus_sep-0.5)/(2*sqrt(0.5))))
    if (TRUE){
      center_N_spks_mat[3,1] = N_spks_total*0.9*(0.5+(clus_sep)/4)
      center_N_spks_mat[3,2] = N_spks_total*0.9*(0.5-(clus_sep)/4)
    } else {
      center_N_spks_mat[3,1] = N_spks_total*0.9*0.5
      center_N_spks_mat[3,2] = N_spks_total*0.9*0.5
    }
    center_N_spks_mat[4,1] = N_spks_total*(0.5+clus_sep/2)
    center_N_spks_mat[4,2] = N_spks_total*(0.5-clus_sep/2)
  } 
  
  
  # Generate spike density functions --------------------------------------------  
  center_density_array_true = array(dim = c(N_clus, 2, length(t_vec)))
  if(N_clus==1){
    ## Clus 1
    if (TRUE) {
      center_density_array_true[1,1, ] = (2 - 2*cos(4*pi*(t_vec-0.4))) * I(0.4<=t_vec & t_vec<=0.9)
    } else {
      s_tmp = 1*(1/4)*(1); mu_tmp = -0.35; 
      center_density_array_true[1,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    }
    
    if (TRUE) {
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[1,2, ] = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
    } else {
      s_tmp = sqrt(1)*(1/2/sqrt(2))*(1); mu_tmp = s_tmp
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[1,2, ] = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
      
    }
    
    ### Add weights for two components
    center_density_array_true[1,1,] = center_density_array_true[1,1,]*center_N_spks_mat[1,1]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[1,2,] = center_density_array_true[1,2,]*center_N_spks_mat[1,2]/sum(center_N_spks_mat[1,1:2])
    
  } else if (N_clus==4) {
    new_version = TRUE
    ## Clus 1
    if (new_version) {
      center_density_array_true[1,1, ] = (2 - 2*cos(4*pi*(t_vec-0.4))) * I(0.4<=t_vec & t_vec<=0.9)
    } else {
      s_tmp = 1*(1/4)*(1); mu_tmp = -0.35; 
      center_density_array_true[1,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    }
    
    if (new_version) {
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[1,2, ] = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
    } else {
      s_tmp = sqrt(1)*(1/2/sqrt(2))*(1); mu_tmp = s_tmp
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[1,2, ] = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
    }
    
    ## Clus 2
    if (new_version) {
      center_density_array_true[2,1, ] = (2 - 2*cos(4*pi*(t_vec-0.4))) * I(0.4<=t_vec & t_vec<=0.9)
    } else {
      s_tmp = 1*(1/4)*(clus_sep*0+1); mu_tmp = -0.35; 
      center_density_array_true[2,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    }
    
    if (new_version) {
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[2,2, ] = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
    } else {
      s_tmp = sqrt(1)*(1/2/sqrt(2))*(sqrt(clus_sep*0+1)); mu_tmp = s_tmp
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[2,2, ] = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
    }
    
    ## Clus 3
    if (new_version) {
      center_density_array_true[3,1, ] = (2 - 2*cos(4*pi*(t_vec-0.4))) * I(0.4<=t_vec & t_vec<=0.9)
    } else {
      s_tmp = 1*(1/4)*(clus_sep^2*0+1); mu_tmp = -0.35
      center_density_array_true[3,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    }
    
    if (new_version) {
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[3,2, ] = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
    } else {
      s_tmp = sqrt(1)*(1/2/sqrt(2))*(clus_sep*0+1); mu_tmp = s_tmp
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[3,2, ] = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
    }
    
    
    ## Clus 4
    if (new_version) {
      center_density_array_true[4,1, ] = (2 - 2*cos(4*pi*(t_vec-0.4))) * I(0.4<=t_vec & t_vec<=0.9)
    } else {
      s_tmp = 1*(1/4); mu_tmp = -0.35; 
      center_density_array_true[4,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    }
    
    if (new_version) {
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[4,2, ] = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
    } else {
      s_tmp = sqrt(1)*(1/2/sqrt(2))*(clus_sep*0+1); mu_tmp = s_tmp
      t_vec_shift = t_vec - (key_times_vec[2]-0)
      center_density_array_true[4,2, ] = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
    }
    
    ### Add weights (prop to N_spks) for two components
    for (id_clus in 1:N_clus){
      center_density_array_true[id_clus,1,] = center_density_array_true[id_clus,1,]*center_N_spks_mat[id_clus,1]/sum(center_N_spks_mat[id_clus,1:2])
      center_density_array_true[id_clus,2,] = center_density_array_true[id_clus,2,]*center_N_spks_mat[id_clus,2]/sum(center_N_spks_mat[id_clus,1:2])
    }
    ### Adjust Cluster 2 intensity components                        
    if (TRUE){
      if (new_version) {
        t_vec_shift = 2*(t_vec - (0.8-0))
        tmp_density = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
        center_density_array_true[2,1, ] = center_density_array_true[2,1, ] + 0.8*tmp_density*min(clus_sep,0.5)
      } else {
        s_tmp = 1*(0.25); mu_tmp = s_tmp; 
        t_vec_shift = t_vec - (-0.2-0)
        tmp_density = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
        center_density_array_true[2,1, ] = center_density_array_true[2,1, ] + 0.4*tmp_density*min(clus_sep,0.5)
      }
      
      if (new_version) {
        t_vec_shift = 2*(t_vec - (key_times_vec[2]-0))
        tmp_density = (2 - 2*cos(2*pi*sqrt(abs(2*t_vec_shift)))) * I(0<=t_vec_shift & t_vec_shift<=0.5)
        center_density_array_true[2,2, ] = center_density_array_true[2,2, ] - 0.8*tmp_density*min(clus_sep,0.5)
      } else {
        s_tmp = 1*(0.25); mu_tmp = s_tmp; 
        t_vec_shift = t_vec - (key_times_vec[2]-0)
        tmp_density = 1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_shift)) + mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_shift & t_vec_shift<=(mu_tmp+s_tmp)^2) 
        center_density_array_true[2,2, ] = center_density_array_true[2,2, ] - 0.4*tmp_density*min(clus_sep,0.5)
      }
      
      N_spks_current_clus_total = sum(center_N_spks_mat[2,1:2])
      center_N_spks_mat[2,1] = sum(center_density_array_true[2,1, ]*t_unit) * N_spks_current_clus_total
      center_N_spks_mat[2,2] = sum(center_density_array_true[2,2, ]*t_unit) * N_spks_current_clus_total
    }
    
    
  } 
  
  # Generate spike intensity functions --------------------------------------------      
  center_intensity_array_true = array(dim = c(N_clus, 2, length(t_vec)))
  for (id_clus in 1:N_clus){
    center_intensity_array_true[id_clus,1, ] = center_density_array_true[id_clus,1,]*sum(center_N_spks_mat[id_clus,1:2])
    center_intensity_array_true[id_clus,2, ] = center_density_array_true[id_clus,2,]*sum(center_N_spks_mat[id_clus,1:2])    
  }
  
  ### Add baseline intensity for all components
  for (id_clus in 1:N_clus){
    if (N_clus==4) {
      intensity_baseline = 20
    } else if (N_clus==1){
      intensity_baseline = 20
    }
    
    center_intensity_array_true[id_clus,1, ] = intensity_baseline + center_intensity_array_true[id_clus,1, ] 
    center_N_spks_mat[id_clus,1] = sum(center_intensity_array_true[id_clus,1, ]*t_unit)
    
    center_density_array_true[id_clus,1,] = center_intensity_array_true[id_clus,1, ] / sum(center_N_spks_mat[id_clus,1:2])
    center_density_array_true[id_clus,2,] = center_intensity_array_true[id_clus,2, ] / sum(center_N_spks_mat[id_clus,1:2])
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
        v_tmp_2 = v_mat_list[[2]][id_subj, id_trial]
        spks_time_mlist[id_subj, id_trial] = list(c( rejection_sampling(density_vec = center_density_array_true[id_clus,1,], 
                                                                        t_vec = t_vec, 
                                                                        N_sample = 0*center_N_spks_mat[id_clus,1]+
                                                                          1*rpois(n=1, lambda=center_N_spks_mat[id_clus,1]) )+
                                                       stim_onset_vec[id_trial]+v_tmp_1,
                                                     rejection_sampling(density_vec = center_density_array_true[id_clus,2,], 
                                                                        t_vec = t_vec, 
                                                                        N_sample = 0*center_N_spks_mat[id_clus,2]+
                                                                          1*rpois(n=1, lambda=center_N_spks_mat[id_clus,2]) )+
                                                       stim_onset_vec[id_trial]+v_tmp_2 ))
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
              t_vec=t_vec
  ))
  
}


