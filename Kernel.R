###################################################
# Creates a function to evaluate a proxy at a given age  
# 
# Marco A. Aquino-López
###################################################
# rm(list = ls()) 	# clean R environment

setwd('~/MEGA_A/Cambridge/Pulque implementation/')
# Load required packages
library(data.table)

# Read the data from the file
file_path <- './Target/EASM_PC1_1kensemble.txt'
data <- read.table(file_path,header = T,row.names = 1)

data <- as.data.frame(data )


# Initialize an empty list to store the density objects for each depth
kde_list <- list()

# Loop through each column (depth) to calculate the kernel density
for (col in names(data)) {
  values_at_depth <- data[[col]]
  
  # Compute kernel density and add to the list
  kde <- density(values_at_depth, kernel = 'gaussian', bw = 0.5)
  kde_list[[col]] <- kde
}

# Define a function to get density for a given depth and value y
depth_kernel <- function(depth, y) {
  kde <- kde_list[[depth]]
  
  # Interpolate density at value y
  density_at_y <- approx(kde$x, kde$y, xout = y)$y
  
  return(density_at_y)
}





