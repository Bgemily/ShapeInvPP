
### Obtain truncated fourier series (smoothed point process) for each cluster and intensity component
get_center_intensity_array = function(subjtrial_density_unsmooth_array,
                                      fft_subjtrial_density_unsmooth_array,
                                      N_spks_mat,
                                      v_trialwise_vec_list = NULL,
                                      clusters_list, 
                                      v_mat_list=NULL,
                                      N_component=1,
                                      freq_trun=5,
                                      t_vec=seq(0, 1, by=0.01),
                                      key_times_vec = c(min(t_vec), 0, max(t_vec)),
                                      bw=0,
                                      fix_timeshift=FALSE,
                                      alpha=0)
{  
  t_unit = t_vec[2]-t_vec[1]
  N_clus = length(clusters_list)
  N_subj = dim(subjtrial_density_unsmooth_array)[1]
  N_trial = dim(subjtrial_density_unsmooth_array)[2]
  
  center_intensity_baseline_vec = rep(NA,N_clus)
  center_density_baseline_vec = rep(NA,N_clus)
  center_intensity_array = array(0, dim=c(N_clus, N_component, length(t_vec)))
  center_density_array = array(0, dim=c(N_clus, N_component, length(t_vec)))
  center_Nspks_mat = matrix(0, nrow=N_clus, ncol=N_component)
  for (q in 1:N_clus) {
    intensity_q_mat = matrix(0, nrow = N_component, ncol = length(t_vec))
    density_q_mat = matrix(0, nrow = N_component, ncol = length(t_vec))
    F_hat_q = 0
    intensity_baseline = 0
    density_baseline = 0
    if (length(clusters_list[[q]])*N_trial>=2) {
      ### Temporarily force the minimum time shift of second component to be zero
      if (!fix_timeshift) {
        for (id_component in 1:N_component){
          if (id_component > 1) {
            id_trial = 1
            v_subjwise_vec = v_mat_list[[id_component]][clusters_list[[q]], id_trial] - v_trialwise_vec_list[[id_component]][id_trial]
            v_subjwise_vec = v_subjwise_vec - min(v_subjwise_vec)
            v_trialwise_vec = v_trialwise_vec_list[[id_component]]
            v_mat_list[[id_component]][clusters_list[[q]], ] = matrix(v_subjwise_vec, nrow = length(clusters_list[[q]]), ncol = N_trial) + 
              matrix(v_trialwise_vec, byrow = TRUE, nrow = length(clusters_list[[q]]), ncol = N_trial)
          }
        }
      }
      
      ### Calculate terms in the analytical solution of least-squares-estimator 
      Y_mat_q = matrix(nrow = length(t_vec), ncol = length(clusters_list[[q]])*N_trial )
      fft_curr_clus_density_unsmooth_array = fft_subjtrial_density_unsmooth_array[clusters_list[[q]], 1:N_trial, , drop = FALSE]
      Y_mat_q = t(apply(fft_curr_clus_density_unsmooth_array, MARGIN = 3, FUN = c))
      
      N_spks_subjtrial_vec_q = rep(0, length(clusters_list[[q]])*N_trial)
      N_spks_subjtrial_vec_q = c(N_spks_mat[clusters_list[[q]], ])
      
      X_array_q = array(dim = c(length(t_vec), length(clusters_list[[q]])*N_trial, N_component))
      mid = length(t_vec) %/% 2
      l_vec = c( 0:mid, (mid+1-length(t_vec)):(-1))
      for (id_component in 1:N_component) {
        v_curr_clus_mat = v_mat_list[[id_component]][clusters_list[[q]], , drop = FALSE]
        exp_timeshift_array = exp(-1i*2*pi*outer(v_curr_clus_mat, l_vec)/(max(t_vec)-min(t_vec)))
        exp_timeshift_mat = t(apply(exp_timeshift_array, MARGIN = 3, FUN = c))
        X_array_q[ , , id_component] = exp_timeshift_mat
      }
      
      ### Get the least-squares estimator
      if(length(N_spks_subjtrial_vec_q)>1){
        if (alpha==0) {
          W_mat_q = diag(N_spks_subjtrial_vec_q)
          theta_list = list()
          for (l in 2:length(l_vec)) {
            W_X = matrix(diag(W_mat_q), nrow = nrow(W_mat_q), ncol = N_component) * as.matrix(X_array_q[l, , ])
            Xt_W_X = t(Conj(as.matrix(X_array_q[l, , ]))) %*% W_X
            inv_Xt_W_X = tryCatch(solve(Xt_W_X), error=function(Xt_W_X){return(NULL)})
            W_Y = as.matrix(diag(W_mat_q) * Y_mat_q[l, ])
            if(is.null(inv_Xt_W_X)){
              tmp_2 = 0
              if (abs((t(Conj(X_array_q[l, , ])) %*% W_X)[1,1]) == 0) {
                tmp_1 = 0
              } else {
                tmp_1 = (t(Conj(X_array_q[l, , ])) %*% W_X)[1,1]^(-1) *
                  (t(Conj(X_array_q[l, , ])) %*% W_Y)[1]
              }
              res = c(tmp_1, rep(0,N_component-1))
              theta_list = c(theta_list, list(res))
            } else {
              res =  inv_Xt_W_X %*% (t(Conj(X_array_q[l, , ])) %*% W_Y) 
              theta_list = c(theta_list, list(res))
            }
          }
        } else if (alpha>0) {
          W_mat_q = diag(N_spks_subjtrial_vec_q)
          Ct_Xt_W_X_C_sum = 0
          Ct_Xt_W_Y_sum = 0
          for (l in 2:length(l_vec)) {
            W_X = matrix(diag(W_mat_q), nrow = nrow(W_mat_q), ncol = N_component) * as.matrix(X_array_q[l, , ])
            Xt_W_X = t(Conj(as.matrix(X_array_q[l, , ]))) %*% W_X
            W_Y = as.matrix(diag(W_mat_q) * Y_mat_q[l, ])
            Xt_W_Y = t(Conj(X_array_q[l, , ])) %*% W_Y
            C_l = kronecker(matrix(1*((2:length(l_vec)==l)), nrow = 1), diag(x = 1, nrow = N_component) )
            Ct_Xt_W_X_C_sum = Ct_Xt_W_X_C_sum + t(Conj(C_l)) %*% Xt_W_X %*% C_l
            Ct_Xt_W_Y_sum = Ct_Xt_W_Y_sum + t(Conj(C_l)) %*% Xt_W_Y
          }
          D_Dt_sum = 0
          for (id_component in 1:N_component) {
            D_m = kronecker( rep(1, length(l_vec)-1), 1*((1:N_component)==id_component) )
            D_Dt_sum = D_Dt_sum + D_m %*% t(D_m)
          }
          res = solve(Ct_Xt_W_X_C_sum + alpha*(diag(x = 1, nrow = nrow(D_Dt_sum)) + D_Dt_sum) ) %*% Ct_Xt_W_Y_sum
          theta_list = lapply( 2:length(l_vec), function(l){ res[((l-2)*N_component+1):((l-1)*N_component)] } )
          
        } else {
          stop("The value of alpha should be non-negative.")
        }
        
      } else{
        W_mat_q = N_spks_subjtrial_vec_q
        theta_list = lapply( 2:length(l_vec), function(l){ c((X_array_q[l, , 1])^(-1)*Y_mat_q[l, ], rep(0, N_component-1)) } )
      }
      for (id_component in 1:N_component) {
        fft_vec_tmp = sapply(theta_list, "[", id_component)
        if (TRUE | (id_component == 1)) {
          fft_l_eq_0 = -sum(fft_vec_tmp[which(abs(l_vec[-1])<=freq_trun)])
          fft_vec_tmp = c(fft_l_eq_0, fft_vec_tmp)
        } else {
          fft_vec_tmp = c(0, fft_vec_tmp)
        }
        density_q_mat[id_component, ] = Re(fft(fft_vec_tmp, inverse = TRUE))
      }
      
      ### Add constants to densities to force the tails to be zero
      if (FALSE & N_component >= 2) {
        for (id_component in 1:N_component) {
          index_tail = which(t_vec >= max(key_times_vec))
          density_q_tail = density_q_mat[id_component, index_tail]
          density_baseline = density_baseline + mean(density_q_tail)
          density_q_mat[id_component, ] = density_q_mat[id_component, ] - mean(density_q_tail)
        }
      }
      
      ### Force densities to be zero before their starting time to fix non-identifiability issue due to small variance of time shifts
      if (N_component >= 2) {
        for (id_component in 2:N_component) {
          index_before_start = which(t_vec <= key_times_vec[id_component])
          density_q_before_start = density_q_mat[id_component, index_before_start]
          density_q_mat[1, index_before_start] = density_q_mat[1, index_before_start] + density_q_before_start
          density_q_mat[id_component, index_before_start] = 0
        }
      }
      
      ### Force the second density to be non-zero right after their starting time
      if (FALSE & N_component >= 2) {
        for (id_component in 2:N_component) {
          density_value_after_key_time = density_q_mat[id_component, which(t_vec > key_times_vec[id_component])[1]]
          density_value_peak = max(abs(density_q_mat[id_component, ]))
          density_value_non_zero = density_value_peak*0.05
          if (abs(density_value_after_key_time) < density_value_non_zero) {
            length_rmv = (min(t_vec[which(abs(density_q_mat[id_component, ]) >= density_value_non_zero)]) - key_times_vec[id_component]) / t_unit
            density_q_mat[id_component, ] = c( tail(density_q_mat[id_component, ], length(t_vec)-length_rmv), 
                                               rep(tail(density_q_mat[id_component, ], 1), length_rmv) )
            if (!fix_timeshift) {
              v_mat_list[[id_component]][clusters_list[[q]], ] = v_mat_list[[id_component]][clusters_list[[q]], ] + length_rmv*t_unit
            }
          }
        }
        
      }
      
      
      ### Force the tails of densities to be zero
      if (FALSE) {
        for (id_component in 1:N_component) {
          density_q_mat[id_component, ] = density_q_mat[id_component, ] * I(t_vec <= max(key_times_vec) )
        }
      }
      
      ### Calculate intensity components
      F_hat_q = mean(N_spks_subjtrial_vec_q)
      for (id_component in 1:N_component) {
        intensity_q_mat[id_component, ] = density_q_mat[id_component, ] * F_hat_q
      }
      ### Calculate baseline density and baseline intensity 
      density_baseline = (1 - sum(density_q_mat*t_unit)) / (max(t_vec)-min(t_vec))
      intensity_baseline = density_baseline * F_hat_q
      
    } else if (length(clusters_list[[q]])*N_trial==1){
      id_subj_tmp = 1
      id_subj = clusters_list[[q]][id_subj_tmp]    
      id_trial = 1
       
      ### Let the first density to be the smoothed point process, the second density to be zero
      N_spks_subjtrial_vec_q = N_spks_mat[id_subj,id_trial]
      density = subjtrial_density_unsmooth_array[id_subj, id_trial, ]
      
      density_baseline = mean(density[t_vec>=max(key_times_vec)])
      density_q_mat[1, ] = density - density_baseline
      if (N_component >= 2) {
        density_q_mat[2:N_component, ] = 0
      }
      
      ### Force the tails of densities to be zero
      for (id_component in 1:N_component) {
        density_q_mat[id_component, ] = density_q_mat[id_component, ] * I(t_vec <= max(t_vec) - max(v_mat_list[[id_component]]))
      }
      
      ### Calculate intensity components
      F_hat_q = mean(N_spks_subjtrial_vec_q)
      for (id_component in 1:N_component) {
        intensity_q_mat[id_component, ] = density_q_mat[id_component, ] * F_hat_q
      }
      intensity_baseline = density_baseline * F_hat_q
    }
    
    # Smooth density
    if (freq_trun < Inf) {
      for (id_component in 1:N_component) {
        fft_density_tmp = fft(density_q_mat[id_component, ]) / length(t_vec)
        fft_density_tmp_smooth = c(head(fft_density_tmp, freq_trun+1), 
                                   rep(0, length(t_vec)-2*freq_trun-1),
                                   tail(fft_density_tmp, freq_trun))
        density_tmp_smooth = Re(fft(fft_density_tmp_smooth, inverse = TRUE))
        
        N_spks_tmp = sum(intensity_q_mat[id_component, ]) / sum(density_q_mat[id_component, ])
        intensity_q_mat[id_component, ] = N_spks_tmp * density_tmp_smooth
        density_q_mat[id_component, ] = density_tmp_smooth
      }
    }
    
     
    center_intensity_baseline_vec[q] = intensity_baseline
    center_density_baseline_vec[q] = density_baseline
    center_intensity_array[q, , ] = intensity_q_mat
    center_density_array[q, , ] = density_q_mat
    center_Nspks_mat[q, ] = F_hat_q * rowSums(density_q_mat * t_unit)
  }
  
  
  
  return(list(center_intensity_baseline_vec=center_intensity_baseline_vec,
              center_density_baseline_vec=center_density_baseline_vec,
              center_intensity_array=center_intensity_array,
              center_density_array=center_density_array,
              center_Nspks_mat=center_Nspks_mat,
              v_mat_list=v_mat_list))
}


