#%%
# Clear workspace and source all R scripts in the Functions directory
rm(list=ls())
file_path = "./Functions"
file.sources = list.files(path = file_path, pattern = "*.R$", full.names = TRUE)
sapply(file.sources, source)

#%%
# Set up parameters for data generation
data_param = list(SEED=1,
                    N_subj=40,
                    N_trial=5,
                    N_clus=4, 
                    t_vec=seq(0,2.5,by=0.01),
                    key_times_vec = c(-1,0-0.2,1.5)+1,
                    N_spks_total = 150,
                    timeshift_subj_max_vec = c(1/32/4, 1/32)*2,
                    timeshift_trial_max = 0.1,
                    clus_sep = 0.5 )
  
# Define output directory for plots
output_dir <- "/Users/zitongzhang/Documents/Academic/SC/ShapeInvPP/Results/Plots/SHAP-175"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#%%
# Generate data with time-varying baseline and save true intensity curves
data_generated = do.call(what = generate_data_timevarying_baseline, args = data_param)
png(file = file.path(output_dir, "true_intensity_timevarying.png"), width = 1800, height = 2400, res = 200)
par(mfrow = c(4, 3), mar = c(2, 2, 2, 1))
for (i in 1:4) {
    for (j in 1:2) {
        plot(data_param$t_vec, data_generated$center_intensity_array_true[i, j, ],
             type = "l", main = paste("GP bsl. Cluster", i, "Component", j),
             xlab = "Time", ylab = "Intensity", ylim=c(-10,max(data_generated$center_intensity_array_true)))
    }
    sum_curve = data_generated$center_intensity_array_true[i, 1, ] + data_generated$center_intensity_array_true[i, 2, ]
    plot(data_param$t_vec, sum_curve,
         type = "l", main = paste("Sum. Cluster", i),
         xlab = "Time", ylab = "Sum Intensity", ylim=c(-10,max(data_generated$center_intensity_array_true)))
    abline(h=0, col='red')
}
par(mfrow = c(1, 1))
dev.off()

#%%
# Generate data with time-varying baseline and save true density curves
data_generated = do.call(what = generate_data_timevarying_baseline, args = data_param)
png(file = file.path(output_dir, "true_density_timevarying.png"), width = 1800, height = 2400, res = 200)
par(mfrow = c(4, 3), mar = c(2, 2, 2, 1))
for (i in 1:4) {
    for (j in 1:2) {
        plot(data_param$t_vec, data_generated$center_density_array_true[i, j, ],
             type = "l", main = paste("GP bsl. Cluster", i, "Component", j),
             xlab = "Time", ylab = "Density", ylim=c(-1, max(data_generated$center_density_array_true)))
        abline(h=0, col='red')
    }
    sum_curve = data_generated$center_density_array_true[i, 1, ] + data_generated$center_density_array_true[i, 2, ]
    print(sum(sum_curve*0.01))
    plot(data_param$t_vec, sum_curve,
         type = "l", main = paste("Sum. Cluster", i),
         xlab = "Time", ylab = "Density", ylim=c(-1, max(data_generated$center_density_array_true)))
    abline(h=0, col='red')    
}
par(mfrow = c(1, 1))
dev.off()

#%%
# Generate data with constant baseline and save true intensity curves
data_generated = do.call(what = generate_data, args = data_param)
png(file = file.path(output_dir, "true_intensity_constant.png"), width = 1800, height = 2400, res = 200)
par(mfrow = c(4, 3), mar = c(2, 2, 2, 1))
for (i in 1:4) {
    for (j in 1:2) {
        plot(data_param$t_vec, data_generated$center_intensity_array_true[i, j, ],
             type = "l", main = paste("Const bsl. Cluster", i, "Component", j),
             xlab = "Time", ylab = "Intensity", ylim=c(-10,max(data_generated$center_intensity_array_true)))
        abline(h=0, col='red')      
    }
    sum_curve = data_generated$center_intensity_array_true[i, 1, ] + data_generated$center_intensity_array_true[i, 2, ]
    plot(data_param$t_vec, sum_curve,
         type = "l", main = paste("Sum. Cluster", i),
         xlab = "Time", ylab = "Sum Intensity", ylim=c(-10,max(data_generated$center_intensity_array_true)))
    abline(h=0, col='red')
}
par(mfrow = c(1, 1))
dev.off()

#%%
# Generate data with constant baseline and save true density curves
data_generated = do.call(what = generate_data, args = data_param)
png(file = file.path(output_dir, "true_density_constant.png"), width = 1800, height = 2400, res = 200)
par(mfrow = c(4, 3), mar = c(2, 2, 2, 1))
for (i in 1:4) {
    for (j in 1:2) {
        plot(data_param$t_vec, data_generated$center_density_array_true[i, j, ],
             type = "l", main = paste("Const bsl. Cluster", i, "Component", j),
             xlab = "Time", ylab = "Density", ylim=c(-1, max(data_generated$center_density_array_true)))
        abline(h=0, col='red')      
    }
    sum_curve = data_generated$center_density_array_true[i, 1, ] + data_generated$center_density_array_true[i, 2, ]
    print(sum(sum_curve*0.01))
    plot(data_param$t_vec, sum_curve,
         type = "l", main = paste("Sum. Cluster", i),
         xlab = "Time", ylab = "Density", ylim=c(-1, max(data_generated$center_density_array_true)))
    abline(h=0, col='red')    
}
par(mfrow = c(1, 1))
dev.off()
