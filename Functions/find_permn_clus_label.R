### Find the best permutation of cluster labels
find_permn_clus_label = function(clusters_list_true, clusters_list_est)
{
  memb_true_vec = clus2mem(clusters_list_true)
  N_clus = length(clusters_list_true)
  permn_list = combinat::permn(1:N_clus)  
  accuracy_max = 0
  the_permn = c()
  for (permn in permn_list) {
    memb_est_vec_permn = clus2mem(clusters_list_est[permn])
    accuracy_tmp = sum(memb_est_vec_permn == memb_true_vec) / length(memb_true_vec)
    if (accuracy_tmp > accuracy_max){
      the_permn = permn
      accuracy_max = accuracy_tmp
    }
  }  
  
  return(the_permn)
}


