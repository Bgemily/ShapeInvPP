
### Obtain truncated fourier series (smoothed point process) for each cluster and intensity component
get_center_intensity_array = function(spks_time_mlist, 
                                      stim_onset_vec, 
                                      clusters_list, 
                                      v_vec=NULL,
                                      v_mat_list=NULL,
                                      N_component=1,
                                      freq_trun=5,
                                      v0 = 0.15, v1 = 0.1,
                                      t_vec=seq(0, v0, by=0.01),
                                      n0_mat_list=NULL,
                                      bw=0,
                                      align_density=FALSE,
                                      fix_timeshift=FALSE,
                                      # Unused arguments
                                      reaction_time_vec=NULL,
                                      rmv_conn_prob=FALSE)
{  
  t_unit = t_vec[2]-t_vec[1]
  N_clus = length(clusters_list)
  N_node = nrow(spks_time_mlist)
  N_replicate = ncol(spks_time_mlist)
  
  if(N_component==1){
    center_intensity_array = array(dim=c(N_clus, N_component, length(t_vec)))
    center_density_array = array(dim=c(N_clus, N_component, length(t_vec)))
    center_Nspks_mat = matrix(nrow=N_clus, ncol=N_component)
    for (q in 1:N_clus) {
      intensity_q = rep(0, length(t_vec))
      density_q = rep(0, length(t_vec))
      F_hat_q = 0
      if(length(clusters_list[[q]])>0){
        spks_time_q = c()
        for (id_node in clusters_list[[q]]) {
          for (id_replicate in 1:N_replicate) {
            spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
            spks_time_nodetrial = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                              spks_time_nodetrial<=max(t_vec))]
            spks_time_nodetrial = spks_time_nodetrial - v_vec[id_node]
            spks_time_nodetrial = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) &
                                                              spks_time_nodetrial<=max(t_vec))]
            spks_time_q = c(spks_time_q, spks_time_nodetrial)
          }
        }
        ### Force the median of spiking times to be the median of t_vec
        if (align_density & !fix_timeshift) {
          v_vec[clusters_list[[q]]] = v_vec[clusters_list[[q]]] + (median(spks_time_q) - median(t_vec))
          spks_time_q = spks_time_q - (median(spks_time_q) - median(t_vec))
          spks_time_q = spks_time_q[which(spks_time_q>=min(t_vec) &
                                            spks_time_q<=max(t_vec))]
        }
        
        N_spks_nodetrial_vec_q = c()
        for (id_node in clusters_list[[q]]) {
          for (id_replicate in 1:N_replicate) {
            spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
            spks_time_nodetrial = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                              spks_time_nodetrial<=max(t_vec))]
            spks_time_nodetrial = spks_time_nodetrial - v_vec[id_node]
            spks_time_nodetrial = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) &
                                                              spks_time_nodetrial<=max(t_vec))]
            N_spks_nodetrial_vec_q = c(N_spks_nodetrial_vec_q, length(spks_time_nodetrial))
          }
        }
        
        if (length(spks_time_q)>0) {
          fft_q_res = get_adaptive_fft(event_time_vec = spks_time_q, 
                                       freq_trun_max = freq_trun, 
                                       t_vec = t_vec, 
                                       bw = bw
                                       )
          fft_q = fft_q_res$fft_vec_best
        } else if(freq_trun<Inf) {
          fft_q = rep(0,2*freq_trun+1)
        } else{
          fft_q = rep(0,length(t_vec))
        }
        intensity_q = fft2f(fft = fft_q, freq_trun = freq_trun, t_vec = t_vec)
        
        N_spks_q = length(spks_time_q)
        if (N_spks_q>0) {
          density_q = intensity_q / N_spks_q
        } else{
          density_q = intensity_q*0
        }
        
        F_hat_q = sqrt( mean(N_spks_nodetrial_vec_q^2) )
        
        intensity_q = density_q*F_hat_q
      }
      
      
      
      center_intensity_array[q,1,] = intensity_q
      center_density_array[q,1,] = density_q
      center_Nspks_mat[q,1] = F_hat_q
      
    }
    
  } else if (N_component==2){
    center_intensity_array = array(0, dim=c(N_clus, N_component, length(t_vec)))
    center_density_array = array(0, dim=c(N_clus, N_component, length(t_vec)))
    center_Nspks_mat = matrix(0, nrow=N_clus, ncol=N_component)
    
    for (q in 1:N_clus) {
      intensity_q_1 = intensity_q_2 = rep(0, length(t_vec))
      density_q_1 = density_q_2 = rep(0, length(t_vec))
      F_hat_q = 0
      if (length(clusters_list[[q]])*N_replicate>=2) {
        ### Temporarily force the minimum time shift of second component to be zero
        if (!fix_timeshift) {
          v_mat_list[[2]][clusters_list[[q]], ] = v_mat_list[[2]][clusters_list[[q]], ] - quantile(v_mat_list[[2]][clusters_list[[q]], ], 0.05)
        }
        
        Y_mat_q = matrix(nrow=length(t_vec), ncol=length(clusters_list[[q]])*N_replicate )
        X_array_q = array(dim = c(length(t_vec), length(clusters_list[[q]])*N_replicate, 2))
        N_spks_nodetrial_vec_q = c()
        mid = length(t_vec)%/%2
        l_vec = c( 0:mid, (mid+1-length(t_vec)):(-1))
        for(id_node_tmp in 1:length(clusters_list[[q]])){
          id_node = clusters_list[[q]][id_node_tmp]    
          for (id_replicate in 1:N_replicate) {
            spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
            spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                              spks_time_nodetrial<=max(t_vec))]
            N_spks_nodetrial_vec_q = c(N_spks_nodetrial_vec_q, length(spks_time_vec))
            
            ### Smooth point process of id_node in id_replicate
            tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                                  freq_trun = freq_trun, 
                                  t_vec = t_vec, 
                                  bw=bw)
            intensity = tmp$intens_vec
            density = intensity/length(spks_time_vec)
            
            ### Save terms in the analytical solution of least-squares-estimator 
            Y_mat_q[ , (id_node_tmp-1)*N_replicate+id_replicate] = fft(density) / length(t_vec)
            X_array_q[ , (id_node_tmp-1)*N_replicate+id_replicate, 1] = exp(-1i*2*pi*l_vec*v_mat_list[[1]][id_node, id_replicate]/(max(t_vec)-min(t_vec)))
            X_array_q[ , (id_node_tmp-1)*N_replicate+id_replicate, 2] = exp(-1i*2*pi*l_vec*v_mat_list[[2]][id_node, id_replicate]/(max(t_vec)-min(t_vec)))
          }
        
        }
        
        ### Get the least-squares estimator
        if(length(N_spks_nodetrial_vec_q)>1){
          W_mat_q = diag(N_spks_nodetrial_vec_q)
          theta_list = lapply( 2:length(l_vec), function(l){
            Xt_W_X = t(Conj(X_array_q[l, , ])) %*% W_mat_q %*% X_array_q[l, , ]
            inv_Xt_W_X = tryCatch(solve(Xt_W_X), error=function(Xt_W_X){return(NULL)})
            if(is.null(inv_Xt_W_X)){
              tmp_2 = 0
              tmp_1 = (t(Conj(X_array_q[l, , ])) %*% W_mat_q %*% X_array_q[l, , ])[1,1]^(-1) *
                (t(Conj(X_array_q[l, , ])) %*% W_mat_q %*% Y_mat_q[l, ])[1]
              return(c(tmp_1, tmp_2))
            } else {
              return( inv_Xt_W_X %*% 
                        (t(Conj(X_array_q[l, , ])) %*% W_mat_q %*% Y_mat_q[l, ]) )
            }
             } )
        } else{
          W_mat_q = N_spks_nodetrial_vec_q
          theta_list = lapply( 2:length(l_vec), function(l){
            c((X_array_q[l, , 1])^(-1)*Y_mat_q[l, ], 0) } )
        }
        alpha_vec = sapply(theta_list, "[", 1)
        beta_vec = sapply(theta_list, "[", 2)
        alpha_beta_0 = sum(N_spks_nodetrial_vec_q*Y_mat_q[1,]) / sum(N_spks_nodetrial_vec_q)
        alpha_vec = c(alpha_beta_0, alpha_vec)
        beta_vec = c(0, beta_vec)
        
        density_q_1 = Re(fft(alpha_vec, inverse = TRUE))
        density_q_2 = Re(fft(beta_vec, inverse = TRUE))
        
        ### Flip two components
        if(FALSE & (which.max(density_q_1)[1] > which.max(density_q_2)[1]) ){
          alpha_vec = sapply(theta_list, "[", 2)
          beta_vec = sapply(theta_list, "[", 1)
          alpha_beta_0 = sum(N_spks_nodetrial_vec_q*Y_mat_q[1,]) / sum(N_spks_nodetrial_vec_q)
          alpha_vec = c(alpha_beta_0, alpha_vec)
          beta_vec = c(0, beta_vec)
          density_q_1 = Re(fft(alpha_vec, inverse = TRUE))
          density_q_2 = Re(fft(beta_vec, inverse = TRUE))
          
          v_mat_list[[1]][clusters_list[[q]], ] -> tmp
          v_mat_list[[1]][clusters_list[[q]], ] = v_mat_list[[2]][clusters_list[[q]], ]
          v_mat_list[[2]][clusters_list[[q]], ] = tmp
        }
        
        ### Force second density to be zero when t<=0
        density_q_2_before0 = mean(density_q_2[which(t_vec<=0)])
        density_q_2 = density_q_2-density_q_2_before0
        density_q_1 = density_q_1+density_q_2_before0
        density_q_2[which(t_vec<=0)] = 0
        
        ### Force the second density to be non-zero right after t=0
        if(abs(density_q_2[which(t_vec>0)[1]])<max(abs(density_q_2))*0.05){
          length_rmv = min(t_vec[which(abs(density_q_2)>=max(abs(density_q_2))*0.05)]) / t_unit
          density_q_2 = c( tail(density_q_2,length(t_vec)-length_rmv), rep(tail(density_q_2,1),length_rmv) )
          if (!fix_timeshift) {
            v_mat_list[[2]][clusters_list[[q]], ] = v_mat_list[[2]][clusters_list[[q]], ] + length_rmv*t_unit
          }
        }
        
        ### Force the tails of densities to be zero
        u_1 = v0; u_0 = v1
        density_q_1 = density_q_1 * I(t_vec<=u_1-max(v_mat_list[[1]]))
        density_q_2 = density_q_2 * I(t_vec<=u_1-max(v_mat_list[[2]]))
        
        
        F_hat_q = sqrt( mean(N_spks_nodetrial_vec_q^2) )
        
        intensity_q_1 = density_q_1*F_hat_q
        intensity_q_2 = density_q_2*F_hat_q
        
      } else if (length(clusters_list[[q]])*N_replicate==1){
        id_node_tmp = 1
        id_node = clusters_list[[q]][id_node_tmp]    
        id_replicate = 1
        
        ### The the only one point process
        spks_time_nodetrial = unlist(spks_time_mlist[id_node,id_replicate]) - stim_onset_vec[id_replicate]
        spks_time_vec = spks_time_nodetrial[which(spks_time_nodetrial>=min(t_vec) & 
                                                          spks_time_nodetrial<=max(t_vec))]
        N_spks_nodetrial_vec_q = length(spks_time_nodetrial)
        
        
        ### Smooth the point process 
        tmp = get_smoothed_pp(event_time_vec = spks_time_vec, 
                              freq_trun = freq_trun, 
                              t_vec = t_vec, 
                              bw=bw)
        intensity = tmp$intens_vec
        density = intensity/length(spks_time_vec)
        
        
        ### Let the first density to be the smoothed point process, the second density to be zero
        density_q_1 = density
        density_q_2 = 0*density
        ### Force the tails of densities to be zero
        u_1 = v0; u_0 = v1
        density_q_1 = density_q_1 * I(t_vec<=u_1-max(v_mat_list[[1]]))
        density_q_2 = density_q_2 * I(t_vec<=u_1-max(v_mat_list[[2]]))
        
        F_hat_q = sqrt( mean(N_spks_nodetrial_vec_q^2) )
        
        intensity_q_1 = density_q_1*F_hat_q
        intensity_q_2 = density_q_2*F_hat_q
      }
      
      
      center_intensity_array[q,1,] = intensity_q_1
      center_density_array[q,1,] = density_q_1
      center_Nspks_mat[q,1] = F_hat_q*sum(density_q_1*t_unit)
      
      center_intensity_array[q,2,] = intensity_q_2
      center_density_array[q,2,] = density_q_2
      center_Nspks_mat[q,2] = F_hat_q*sum(density_q_2*t_unit)
    }
    
  }
  
  return(list(center_intensity_array=center_intensity_array,
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              v_vec=v_vec,
              v_mat_list=v_mat_list))
}


