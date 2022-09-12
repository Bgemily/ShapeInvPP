
### Find the optimal permutation
find_permn = function(center_density_array_from, 
                      center_density_array_to)
{
  if (!identical(dim(center_density_array_from), dim(center_density_array_to))) {
    stop("dim(center_density_array_from) and dim(center_density_array_to) should be the same.")
  }
  
  N_clus = dim(center_density_array_from)[1]
  permn_list = combinat::permn(1:N_clus)
  
  min_dist = Inf
  for (permn in permn_list) {
    
    ### Use cross-correlation to measure similarity
    dist = 0
    for (q in 1:N_clus) {
      center_density_from_tmp_1 = center_density_array_from[permn, , ,drop=FALSE][q,1,]
      center_density_to_tmp_1 = center_density_array_to[q,1,]
      if (var(center_density_from_tmp_1)==0) {
        center_density_from_tmp_1 = c(center_density_from_tmp_1[-1],center_density_from_tmp_1[1]+1e-10)
      }
      if (var(center_density_to_tmp_1)==0) {
        center_density_to_tmp_1 = c(center_density_to_tmp_1[-1],center_density_to_tmp_1[1]+1e-10)
      }
      
      
      center_density_from_tmp_2 = center_density_array_from[permn, , ,drop=FALSE][q,2,]
      center_density_to_tmp_2 = center_density_array_to[q,2,]
      if (var(center_density_from_tmp_2)==0) {
        center_density_from_tmp_2 = c(center_density_from_tmp_2[-1],center_density_from_tmp_2[1]+1e-10)
      }
      if (var(center_density_to_tmp_2)==0) {
        center_density_to_tmp_2 = c(center_density_to_tmp_2[-1],center_density_to_tmp_2[1]+1e-10)
      }
      
      center_density_from_tmp = center_density_from_tmp_1 + center_density_from_tmp_2
      center_density_to_tmp = center_density_to_tmp_1 + center_density_to_tmp_2
      dist_tmp = 1-max(ccf(x = center_density_from_tmp, 
                           y = center_density_to_tmp, plot=FALSE)$acf)
      
      dist = dist + dist_tmp
    }
  
    if(dist < min_dist){
      the_permn = permn
      min_dist = dist
    }
    
  }
  
  return(list(permn = the_permn, dist = min_dist))
}
