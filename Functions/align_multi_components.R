align_multi_components = function(f_target,
                                  f_origin_mat,
                                  step_size = 0.02,
                                  t_unit = 0.05, 
                                  n0_vec = c(0,0),
                                  n0_min_vec = 0,
                                  n0_max_vec = length(f_target), 
                                  pad = NULL,
                                  periodic = FALSE,
                                  MaxIter=1000, 
                                  stopping_redu=0.01, 
                                  weights=NULL)
{
  
  ### Extend f_target and f_origin from [0,T] to [-T,2T]
  if(!is.null(pad)){
    extend = function(f) {return(c( rep(pad,length(f)-1), f, rep(pad,2*round(length(f)/2)) ))} # make N:=length_of_func odd
  } else{
    extend = function(f) {return(c( rep(head(f,1),length(f)-1), f, rep(tail(f,1),2*round(length(f)/2)) ))} # make N:=length_of_func odd
  }
  f_target = extend(f_target)
  N_component = nrow(f_origin_mat)
  f_origin_mat_extend = matrix(nrow = N_component, ncol = length(f_target))
  for (id_component in 1:N_component) {
    f_origin_mat_extend[id_component, ] = extend(f_origin_mat[id_component, ])
  }
  f_origin_mat = f_origin_mat_extend
  
  ### Compute terms needed in gradients
  fft_f_target = fft(f_target) #/ length(f_target)
  fft_f_origin_mat = 0 * f_origin_mat
  for (id_component in 1:N_component) {
    fft_f_origin_mat[id_component, ] = fft(f_origin_mat[id_component, ]) #/ length(f_origin_mat[id_component, ])
  }
  
  ### Gradient descent
  iter_count = 0
  dist_redu = Inf
  dist_curr = Inf
  converge = FALSE
  while (!converge && iter_count<MaxIter) {
    iter_count = iter_count + 1
    gd_vec = gradient_multi_component(fft_f_target = fft_f_target, 
                                      fft_f_origin_mat = fft_f_origin_mat,
                                      n0_vec = n0_vec)
    hessian_mat = hessian_multi_component(fft_f_target = fft_f_target, 
                                          fft_f_origin_mat = fft_f_origin_mat,
                                          n0_vec = n0_vec)
    inv_hessian_mat = try(solve(hessian_mat), silent = TRUE)
    if (identical(class(inv_hessian_mat), "try-error")) {
      n0_vec = n0_vec - (diag(hessian_mat) + .Machine$double.eps)^(-1) * gd_vec
    } else {
      n0_vec = n0_vec - inv_hessian_mat %*% gd_vec
    }
    
    if (FALSE) {
      print((step_size)*gd_vec)
    }
    
    n0_vec[which(n0_vec<n0_min_vec)] = n0_min_vec[which(n0_vec<n0_min_vec)]
    n0_vec[which(n0_vec>n0_max_vec)] = n0_max_vec[which(n0_vec>n0_max_vec)]

    
    N = length(fft_f_target)
    l_vec = 0:(N-1)
    l_vec = c( head(l_vec, N-(N-1)%/%2),
               tail(l_vec, (N-1)%/%2) - N )
    fft_f_origin_shifted = 0
    for (id_component in 1:N_component) {
      fft_curr_comp_shifted = exp(1i*2*pi*l_vec*(-n0_vec[id_component])/N) * fft_f_origin_mat[id_component, ]
      fft_f_origin_shifted = fft_f_origin_shifted + fft_curr_comp_shifted
    }
    d = sum(abs( fft_f_origin_shifted - fft_f_target )^2)  
    d = sqrt(t_unit / N * d)
    dist_upd = d
    dist_redu = (dist_curr - dist_upd) / dist_upd
    if (is.na(dist_redu)) dist_redu = 0
    
    ### DEBUG
    # if(dist_redu<0){
      # print(step_size*gd_vec)
    # }
    
    dist_curr = dist_upd
    converge = dist_redu < stopping_redu
  
  }
  n0_vec = round(n0_vec)
  
  if (iter_count == MaxIter) {
    warning("Reached max iteration number when estimating a time shift. Consider adjusting the step size.")
  }
  
  return(list(n0_vec = n0_vec))
  
  
}