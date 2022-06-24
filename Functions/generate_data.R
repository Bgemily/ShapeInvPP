
### Generate synthetic networks
generate_data = function(SEED=NULL,
                         N_node=100, 
                         N_clus=2, 
                         clus_size_vec = rep(N_node/N_clus, N_clus),
                         u_1 = 1, u_0 = 1,
                         t_vec=seq(-u_0*3/2, u_1, by=0.01),
                         ### params when N_clus==1:
                         N_spks_total = 60*10+40*10,
                         N_spks_ratio = 3/2,
                         sd_shrinkage = 1,
                         ### params when N_clus==2:
                         clus_mixture = 0)
{
  if(!is.null(SEED)) set.seed(SEED)
  t_unit = t_vec[2]-t_vec[1]
  
  # Generate cluster memberships ---------------------------------------
  mem_true_vec = rep(1:N_clus, clus_size_vec)
  clus_true_list = mem2clus(mem_true_vec, N_clus_min = N_clus)
  
  
  # Generate time shifts ----------------------------------------------------
  v_vec_list = list()
  v_vec_list[[1]] = runif(n=N_node, 
                          min = 0,
                          max = u_0/2)    
  v_vec_list[[2]] = runif(n = N_node,
                          min = 0,
                          max = u_1-u_0/2 )
  v_vec_list[[1]] = v_vec_list[[1]] - min(v_vec_list[[1]])
  
  
  
  
  # Generate expected number of spikes --------------------------------------------
  center_N_spks_mat = matrix(nrow=N_clus,ncol=2)
  if (N_clus==2) {
    center_N_spks_mat[1,1] = 60*10
    center_N_spks_mat[1,2] = 40*10
    center_N_spks_mat[2,1] = 30*10
    center_N_spks_mat[2,2] = 70*10
  } else if (N_clus==1){
    center_N_spks_mat[1,1] = N_spks_total*(N_spks_ratio/(N_spks_ratio+1))
    center_N_spks_mat[1,2] = N_spks_total*(1/(N_spks_ratio+1))
  } else{
    stop("TODO: specify center_density_array_true")
  }
  
  
  # Generate spike density functions --------------------------------------------  
  center_density_array_true = array(dim = c(N_clus, 2, length(t_vec)))
  if (N_clus==2) {
    mu_tmp = -u_0+u_0/10; s_tmp = u_0/10
    mu_tmp_prime = -u_0/2+u_0/3; s_tmp_prime = u_0/3   
    center_density_array_true[1,1, ] = (2/3)*1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) + 
      (1/3)*1/(2*s_tmp_prime)*( 1 + cos(((t_vec - mu_tmp_prime)/s_tmp_prime)*pi) ) * I(mu_tmp_prime-s_tmp_prime<=t_vec & t_vec<=mu_tmp_prime+s_tmp_prime)
    
    s_tmp = sqrt(u_1)/5
    center_density_array_true[1,2, ] = 1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec)) - s_tmp)/s_tmp) * pi) ) * I(0<=t_vec & t_vec<=(2*s_tmp)^2) 
    
    mu_tmp = -u_0*(1/2); s_tmp = u_0*(1/2)
    center_density_array_true[2,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) 
    
    mu_tmp = sqrt(u_1)/4+u_1/2; s_tmp = sqrt(u_1)/4
    center_density_array_true[2,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec & t_vec<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec & t_vec<=(mu_tmp+s_tmp)^2) 
    
    ### Add weights (prop to N_spks) for two components
    center_density_array_true[1,1,] = center_density_array_true[1,1,]*center_N_spks_mat[1,1]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[1,2,] = center_density_array_true[1,2,]*center_N_spks_mat[1,2]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[2,1,] = center_density_array_true[2,1,]*center_N_spks_mat[2,1]/sum(center_N_spks_mat[2,1:2])
    center_density_array_true[2,2,] = center_density_array_true[2,2,]*center_N_spks_mat[2,2]/sum(center_N_spks_mat[2,1:2])
    
  } else if(N_clus==1){
    s_tmp = u_0/10*sd_shrinkage; mu_tmp = -u_0+s_tmp; 
    s_tmp_prime = u_0/3; mu_tmp_prime = -u_0/2+s_tmp_prime; 
    center_density_array_true[1,1, ] = (2/3)*1/(2*s_tmp)*( 1 + cos(((t_vec - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec & t_vec<=mu_tmp+s_tmp) + 
      (1/3)*1/(2*s_tmp_prime)*( 1 + cos(((t_vec - mu_tmp_prime)/s_tmp_prime)*pi) ) * I(mu_tmp_prime-s_tmp_prime<=t_vec & t_vec<=mu_tmp_prime+s_tmp_prime)
    
    s_tmp = sqrt(u_1)/5*sqrt(sd_shrinkage)
    center_density_array_true[1,2, ] = 1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec)) - s_tmp)/s_tmp) * pi) ) * I(0<=t_vec & t_vec<=(2*s_tmp)^2) 
    
    ### Add weights for two components
    center_density_array_true[1,1,] = center_density_array_true[1,1,]*center_N_spks_mat[1,1]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[1,2,] = center_density_array_true[1,2,]*center_N_spks_mat[1,2]/sum(center_N_spks_mat[1,1:2])

  } else {
    stop("TODO: specify center_density_array_true")
  }
  
  
  # Generate spike intensity functions --------------------------------------------      
  center_intensity_array_true = array(dim = c(N_clus, 2, length(t_vec)))
  if (N_clus==2) {
    center_intensity_array_true[1,1, ] = center_density_array_true[1,1,]*sum(center_N_spks_mat[1,1:2])
    center_intensity_array_true[1,2, ] = center_density_array_true[1,2,]*sum(center_N_spks_mat[1,1:2])    
    center_intensity_array_true[2,1, ] = center_density_array_true[2,1,]*sum(center_N_spks_mat[2,1:2])   
    center_intensity_array_true[2,2, ] = center_density_array_true[2,2,]*sum(center_N_spks_mat[2,1:2])
  } else if (N_clus==1){
    center_intensity_array_true[1,1, ] = center_density_array_true[1,1,]*sum(center_N_spks_mat[1,1:2])
    center_intensity_array_true[1,2, ] = center_density_array_true[1,2,]*sum(center_N_spks_mat[1,1:2])    
  } else {
    stop("TODO: specify center_density_array_true")
  }
  
  
  # Mix two clusters for both components ---------------------------
  if(N_clus==2){
    center_intensity_array_true_mixed = center_intensity_array_true
    center_intensity_array_true_mixed[1,1,] = (1-clus_mixture)*center_intensity_array_true[1,1,] + clus_mixture*center_intensity_array_true[2,1,]
    center_intensity_array_true_mixed[2,1,] = clus_mixture*center_intensity_array_true[1,1,] + (1-clus_mixture)*center_intensity_array_true[2,1,]
    center_intensity_array_true_mixed[1,2,] = (1-clus_mixture)*center_intensity_array_true[1,2,] + clus_mixture*center_intensity_array_true[2,2,]
    center_intensity_array_true_mixed[2,2,] = clus_mixture*center_intensity_array_true[1,2,] + (1-clus_mixture)*center_intensity_array_true[2,2,]
    center_intensity_array_true = center_intensity_array_true_mixed
    for (id_clus in 1:N_clus) {
      for (id_component in 1:2) {
        center_N_spks_mat[id_clus, id_component] = sum(center_intensity_array_true[id_clus, id_component, ]*t_unit)
        N_spks_tmp = sum(center_intensity_array_true[id_clus, , ]*t_unit)
        center_density_array_true[id_clus, id_component, ] = center_intensity_array_true[id_clus, id_component, ] / N_spks_tmp
      }
    }
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
  
  
  stim_onset_vec = 0
  spks_time_mlist = matrix(list(),N_node,1)
  for (id_clus in 1:N_clus) {
    for (id_node in clus_true_list[[id_clus]]) {
      for (j in 1:1) {
        v_tmp_1 = v_vec_list[[1]][id_node]
        v_tmp_2 = v_vec_list[[2]][id_node]
        spks_time_mlist[id_node,j] = list(c( rejection_sampling(density_vec = center_density_array_true[id_clus,1,], 
                                                                t_vec = t_vec, 
                                                                N_sample = 0*center_N_spks_mat[id_clus,1]+
                                                                  1*rpois(n=1, lambda=center_N_spks_mat[id_clus,1]) )+
                                               stim_onset_vec[j]+v_tmp_1,
                                             rejection_sampling(density_vec = center_density_array_true[id_clus,2,], 
                                                                t_vec = t_vec, 
                                                                N_sample = 0*center_N_spks_mat[id_clus,2]+
                                                                  1*rpois(n=1, lambda=center_N_spks_mat[id_clus,2]) )+
                                               stim_onset_vec[j]+v_tmp_2 ))
      }
    }
  }
  
  
  
  
  # Output ------------------------------------------------------------------
  
  
  return(list(spks_time_mlist=spks_time_mlist, 
              stim_onset_vec=stim_onset_vec,
              mem_true_vec=mem_true_vec, 
              clus_true_list=clus_true_list,
              v_vec_list=v_vec_list, 
              center_density_array_true=center_density_array_true,
              center_N_spks_mat=center_N_spks_mat,
              center_intensity_array_true=center_intensity_array_true,
              t_vec=t_vec
  ))
  
}

