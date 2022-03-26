
plot_intensity_array = function(center_intensity_array,
                                center_Nspks_vec,
                                clusters_list,
                                t_vec,
                                v0 = 0.2, v1 = 0.1,
                                N_component = 1
)
{
  
  t_unit = t_vec[2] - t_vec[1]
  N_clus = nrow(center_intensity_array)
  
  big.df = data.frame()
  for (id_clus in 1:N_clus) {
    for (id_component in 1:N_component) {
      intensity_tmp = center_intensity_array[id_clus,id_component, ]
      tmp.df = data.frame(intensity_val=intensity_tmp)
      tmp.df$cluster = id_clus
      tmp.df$component = id_component
      tmp.df$t = switch(id_component, 
                        `1`=t_vec,
                        `2`=t_vec-max(t_vec)+v0)
      
      big.df = rbind(big.df, tmp.df)
    }
  }
  
  
  clus_size_vec = sapply(clusters_list, length)
  
  ### Draw plots
  g_list = list()
  for (id_clus in 1:N_clus) {
    for (id_component in 1:N_component) {
      g <- big.df %>%
        mutate(cluster = as.factor(cluster), 
               component=as.factor(component)) %>%
        filter(cluster==id_clus & component==id_component) %>%
        ggplot(aes(x=t, y=intensity_val, 
                   group=interaction(cluster,component))) +
        geom_line(alpha=1)+
        geom_vline(xintercept = 0, 
                   color=switch(id_component,
                                `1`='red',
                                `2`='orange'),
                   linetype=switch(id_component,
                                   `1`='solid',
                                   `2`='dashed')) +
        annotate(geom = 'text', label = paste0('Cluster: ',id_clus,
                                               ', ','N_spks: ', round(center_Nspks_vec[id_clus],1),
                                               ', ','size: ',clus_size_vec[id_clus]), 
                 x = Inf, y = Inf, hjust = 1, vjust = 1) +
        xlab(NULL) + ylab(NULL) +
        theme_bw() +
        theme(legend.position = c(.95, .95),
              legend.background = element_blank(),
              legend.key.size = unit(-1, 'cm'), 
              legend.justification = c("right", "top")) +
        theme(strip.text.x = element_blank(),
              strip.background = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))  + 
        coord_cartesian(ylim=range(center_intensity_array))
      
      g_list = c(g_list, list(g))
    }
  }
  
  layout_mat = matrix(1:(N_clus*N_component),N_clus,N_component, byrow=TRUE)
  g = arrangeGrob(grobs=g_list, layout_matrix=layout_mat)
  
  return(list(g=g, big.df=big.df, g_list=g_list))
  
  
}
