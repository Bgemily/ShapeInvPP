### Generate synthetic networks
generate_data = function(SEED=NULL,
                         N_node=100, 
                         N_replicate=1,
                         N_clus=2, 
                         clus_size_vec = rep(N_node/N_clus, N_clus),
                         u_1 = 1, u_0 = 1,
                         t_vec = seq(-u_0,u_1,by=0.01),
                         t_vec_extend = t_vec,
                         N_spks_total = 1000,
                         ### params when N_clus==4:
                         clus_sep = 2,
                         ### params when N_clus==1:
                         N_spks_ratio = 3/2,
                         sd_shrinkage = 1,
                         c_1 = 0, delta_1 = 0,
                         c_2 = 0, delta_2 = 0,
                         c_3 = 0, delta_3 = 0,
                         identical_components = FALSE,
                         ### params when N_clus==2:
                         clus_mixture = 0)
{
  if(!is.null(SEED)) set.seed(SEED)
  t_unit = t_vec[2]-t_vec[1]
  
  # Generate cluster memberships ---------------------------------------
  mem_true_vec = rep(1:N_clus, clus_size_vec)
  clus_true_list = mem2clus(mem_true_vec, N_clus_min = N_clus)
  
  
  # Generate time shifts ----------------------------------------------------
  v_mat_list = list()
  v_mat_list[[1]] = runif(n=N_node*N_replicate, 
                          min = 0,
                          max = u_0/2)  
  v_mat_list[[1]] = matrix(v_mat_list[[1]], nrow = N_node, ncol = N_replicate)
  v_mat_list[[2]] = runif(n = N_node*N_replicate,
                          min = 0,
                          max = u_1-u_0/2 )
  v_mat_list[[2]] = matrix(v_mat_list[[2]], nrow = N_node, ncol = N_replicate)
  for(id_clus in 1:N_clus){
    v_mat_list[[1]][clus_true_list[[id_clus]], ] = v_mat_list[[1]][clus_true_list[[id_clus]], ] - min(v_mat_list[[1]][clus_true_list[[id_clus]], ])
  }
  
  
  
  # Generate expected number of spikes --------------------------------------------
  center_N_spks_mat = matrix(nrow=N_clus,ncol=2)
  if (N_clus==2) {
    center_N_spks_mat[1,1] = N_spks_total*0.6
    center_N_spks_mat[1,2] = N_spks_total*0.4
    center_N_spks_mat[2,1] = N_spks_total*0.3
    center_N_spks_mat[2,2] = N_spks_total*0.7
  } else if (N_clus==1){
    center_N_spks_mat[1,1] = N_spks_total*(N_spks_ratio/(N_spks_ratio+1))
    center_N_spks_mat[1,2] = N_spks_total*(1/(N_spks_ratio+1))
  } else if (N_clus==4){
    center_N_spks_mat[1,1] = N_spks_total*0.5
    center_N_spks_mat[1,2] = N_spks_total*0.5
    center_N_spks_mat[2,1] = N_spks_total*0.5
    center_N_spks_mat[2,2] = N_spks_total*0.5
    center_N_spks_mat[3,1] = N_spks_total*0.5
    center_N_spks_mat[3,2] = N_spks_total*0.5
    center_N_spks_mat[4,1] = N_spks_total*0.5
    center_N_spks_mat[4,2] = N_spks_total*0.5
  } else{
    stop("TODO: specify center_density_array_true")
  }
  
  
  # Generate spike density functions --------------------------------------------  
  center_density_array_true = array(dim = c(N_clus, 2, length(t_vec_extend)))
  if (N_clus==2) {
    mu_tmp = -u_0+u_0/10; s_tmp = u_0/10
    mu_tmp_prime = -u_0/2+u_0/3; s_tmp_prime = u_0/3   
    center_density_array_true[1,1, ] = (2/3)*1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) + 
      (1/3)*1/(2*s_tmp_prime)*( 1 + cos(((t_vec_extend - mu_tmp_prime)/s_tmp_prime)*pi) ) * I(mu_tmp_prime-s_tmp_prime<=t_vec_extend & t_vec_extend<=mu_tmp_prime+s_tmp_prime)
    
    s_tmp = sqrt(u_1)/5
    center_density_array_true[1,2, ] = 1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp) * pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) 
    
    mu_tmp = -u_0*(1/2); s_tmp = u_0*(1/2)
    center_density_array_true[2,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    mu_tmp = sqrt(u_1)/4+u_1/2; s_tmp = sqrt(u_1)/4
    center_density_array_true[2,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_extend)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_extend & t_vec_extend<=(mu_tmp+s_tmp)^2) 
    
    ### Add weights (prop to N_spks) for two components
    center_density_array_true[1,1,] = center_density_array_true[1,1,]*center_N_spks_mat[1,1]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[1,2,] = center_density_array_true[1,2,]*center_N_spks_mat[1,2]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[2,1,] = center_density_array_true[2,1,]*center_N_spks_mat[2,1]/sum(center_N_spks_mat[2,1:2])
    center_density_array_true[2,2,] = center_density_array_true[2,2,]*center_N_spks_mat[2,2]/sum(center_N_spks_mat[2,1:2])
    
  } else if(N_clus==1){
    s_tmp = u_0/10*sd_shrinkage; mu_tmp = -u_0+s_tmp; 
    s_tmp_prime = u_0/3; mu_tmp_prime = -u_0/2+s_tmp_prime; 
    center_density_array_true[1,1, ] = (2/3)*1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) + 
      (1/3)*1/(2*s_tmp_prime)*( 1 + cos(((t_vec_extend - mu_tmp_prime)/s_tmp_prime)*pi) ) * I(mu_tmp_prime-s_tmp_prime<=t_vec_extend & t_vec_extend<=mu_tmp_prime+s_tmp_prime)
    
    s_tmp = sqrt(u_1)/5*sqrt(sd_shrinkage)
    center_density_array_true[1,2, ] = 1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp) * pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) 
    
    ### Add constant density to Component_1 before -u_0
    center_density_array_true[1,1, ] = (1-c_1)*center_density_array_true[1,1, ] + 
      c_1*c(rep(0,length(t_vec_extend[t_vec_extend < -u_0])-round(delta_1/t_unit)),
            rep(1/delta_1,round(delta_1/t_unit)), 
            rep(0,length(t_vec_extend[t_vec_extend >= -u_0])) )
    ### Add constant density to Component_1 before u_1
    center_density_array_true[1,1, ] = (1-c_2)*center_density_array_true[1,1, ] + 
      c_2*c(rep(0,length(t_vec_extend)-round(delta_2/t_unit)),
            rep(1/delta_2,round(delta_2/t_unit)) )
    ### Add constant density to Component_2 before u_1
    center_density_array_true[1,2, ] = (1-c_3)*center_density_array_true[1,2, ] + 
      c_3*c(rep(0,length(t_vec_extend)-round(delta_3/t_unit)),
            rep(1/delta_3,round(delta_3/t_unit)) )
    
    ### Force the two components to be identical up to a time shift
    if(identical_components == TRUE){
      center_density_array_true[1,1,] = c(rep(0, length(which(t_vec_extend < -u_0))),
                                          center_density_array_true[1,2,][which(t_vec_extend>=0)],
                                          rep(0, length(t_vec_extend) - 
                                                length(which(t_vec_extend < -u_0)) - 
                                                length(which(t_vec_extend >= 0)))
      )
    }
    
    ### Add weights for two components
    center_density_array_true[1,1,] = center_density_array_true[1,1,]*center_N_spks_mat[1,1]/sum(center_N_spks_mat[1,1:2])
    center_density_array_true[1,2,] = center_density_array_true[1,2,]*center_N_spks_mat[1,2]/sum(center_N_spks_mat[1,1:2])
    
    
    
  } else if (N_clus==4) {
    ## Clus 1
    s_tmp = u_0*(1/2)*(1/(clus_sep^2)); mu_tmp = -u_0*(1/2); 
    center_density_array_true[1,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    s_tmp = sqrt(u_1)*(1/2/sqrt(2))*(1/clus_sep); mu_tmp = u_1/2+s_tmp; 
    center_density_array_true[1,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_extend)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_extend & t_vec_extend<=(mu_tmp+s_tmp)^2) 
    ## Clus 2
    s_tmp = u_0*(1/2)*(1/clus_sep); mu_tmp = -u_0*(1/2); 
    center_density_array_true[2,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    s_tmp = sqrt(u_1)*(1/2/sqrt(2))*(1/sqrt(clus_sep)); mu_tmp = u_1/2+s_tmp; 
    center_density_array_true[2,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_extend)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_extend & t_vec_extend<=(mu_tmp+s_tmp)^2) 
    ## Clus 3
    s_tmp = u_0*(1/2)*1; mu_tmp = -u_0*(1/2) 
    center_density_array_true[3,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    s_tmp = sqrt(u_1)*(1/2/sqrt(2))*1; mu_tmp = u_1/2+s_tmp 
    center_density_array_true[3,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_extend)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_extend & t_vec_extend<=(mu_tmp+s_tmp)^2) 
    ## Clus 4
    s_tmp = u_0*(1/2)*(1/(clus_sep^2)); mu_tmp = -u_0*(1/2); 
    center_density_array_true[4,1, ] = 1/(2*s_tmp)*( 1 + cos(((t_vec_extend - mu_tmp)/s_tmp)*pi) ) * I(mu_tmp-s_tmp<=t_vec_extend & t_vec_extend<=mu_tmp+s_tmp) 
    
    s_tmp = sqrt(u_1)*(1/2/sqrt(2))*1; mu_tmp = u_1/2+s_tmp 
    center_density_array_true[4,2, ] = (1)*1/(2*s_tmp)^2*( 1 + cos(((sqrt(abs(t_vec_extend)) - s_tmp)/s_tmp)*pi) ) * I(0<=t_vec_extend & t_vec_extend<=(2*s_tmp)^2) +
      (0)*1/(2*s_tmp*2*mu_tmp)*( 1 + cos(((sqrt(abs(t_vec_extend)) - mu_tmp)/s_tmp)*pi) ) * I((mu_tmp-s_tmp)^2<=t_vec_extend & t_vec_extend<=(mu_tmp+s_tmp)^2) 
    
    ### Add weights (prop to N_spks) for two components
    for (id_clus in 1:N_clus){
      center_density_array_true[id_clus,1,] = center_density_array_true[id_clus,1,]*center_N_spks_mat[id_clus,1]/sum(center_N_spks_mat[id_clus,1:2])
      center_density_array_true[id_clus,2,] = center_density_array_true[id_clus,2,]*center_N_spks_mat[id_clus,2]/sum(center_N_spks_mat[id_clus,1:2])
    }
    
  } else {
    stop("TODO: specify center_density_array_true")
  }
  
  
  # Generate spike intensity functions --------------------------------------------      
  center_intensity_array_true = array(dim = c(N_clus, 2, length(t_vec_extend)))
  for (id_clus in 1:N_clus){
    center_intensity_array_true[id_clus,1, ] = center_density_array_true[id_clus,1,]*sum(center_N_spks_mat[id_clus,1:2])
    center_intensity_array_true[id_clus,2, ] = center_density_array_true[id_clus,2,]*sum(center_N_spks_mat[id_clus,1:2])    
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
  
  
  stim_onset_vec = rep(0, N_replicate)
  spks_time_mlist = matrix(list(), nrow = N_node, ncol = N_replicate)
  for (id_clus in 1:N_clus) {
    for (id_node in clus_true_list[[id_clus]]) {
      for (id_replicate in 1:N_replicate) {
        v_tmp_1 = v_mat_list[[1]][id_node, id_replicate]
        v_tmp_2 = v_mat_list[[2]][id_node, id_replicate]
        spks_time_mlist[id_node, id_replicate] = list(c( rejection_sampling(density_vec = center_density_array_true[id_clus,1,], 
                                                                            t_vec = t_vec_extend, 
                                                                            N_sample = 0*center_N_spks_mat[id_clus,1]+
                                                                              1*rpois(n=1, lambda=center_N_spks_mat[id_clus,1]) )+
                                                           stim_onset_vec[id_replicate]+v_tmp_1,
                                                         rejection_sampling(density_vec = center_density_array_true[id_clus,2,], 
                                                                            t_vec = t_vec_extend, 
                                                                            N_sample = 0*center_N_spks_mat[id_clus,2]+
                                                                              1*rpois(n=1, lambda=center_N_spks_mat[id_clus,2]) )+
                                                           stim_onset_vec[id_replicate]+v_tmp_2 ))
        ### Only keep spike times during [-u_0, u_1] 
        spks_time_vec = spks_time_mlist[id_node,id_replicate][[1]]
        spks_time_mlist[id_node,id_replicate][[1]] = spks_time_vec[which(spks_time_vec >= -u_0 & 
                                                                           spks_time_vec <= u_1)]
      }
    }
  }
  
  
  
  
  # Output ------------------------------------------------------------------
  
  
  return(list(spks_time_mlist=spks_time_mlist, 
              stim_onset_vec=stim_onset_vec,
              mem_true_vec=mem_true_vec, 
              clus_true_list=clus_true_list,
              v_mat_list=v_mat_list, 
              center_density_array_true=center_density_array_true,
              center_N_spks_mat=center_N_spks_mat,
              center_intensity_array_true=center_intensity_array_true,
              t_vec=t_vec,
              t_vec_extend=t_vec_extend
  ))
  
}