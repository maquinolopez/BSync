# BSynch 
#
# Synchronizes the Input and Target signals using Bayesian Stats
#
# Args:
#   Input (string): The input data that needs to be synchronized.
#   Target (string): The target data with which the input data will be synchronized.
#   folder (string, default = '~/Documents/BSynch/'): Path where output/results will be saved.
#   ta (numeric, default = 3): Parameter a of the t distribution.
#   tb (numeric, default = 4): Parameter b of the t distribution.
#   shape_acc (numeric, default = 10): Shape parameter for the compression/expansion parameter.
#   meanM (numeric, default = 0.7): Mean of the memory parameter.
#   alpha (numeric, default = 4): scale of the memory parameter.
#   iters (numeric, default = 2e+3): Number of iterations for the MCMC.
#   burn (numeric, default = 1e+3): Burn-in period for the MCMC.
#   thin (numeric, default = 150): Thinning parameter for the MCMC.
#

# Scales the data in the range of [-1,1]
# range <- function(x){
#   2*(x-min(x))/(max(x)-min(x))-1
# }

# range_quantile <- function(x,symetric=TRUE,quantile = .95,lquantile=.025,uquantile = .95) {
#   if (symetric){
#     q1 <- quantile(x, 1- (1- q)/2 )
#     q2 <- quantile(x, (1- q)/2)
#   } else{
#     q1 <- quantile(x, lquantile )
#     q2 <- quantile(x, uquantile)
#   }
#   2 * (x - q2) / (q1 - q2) - 1
# }
# 
# range <- range_quantile
# 
# 
# range2 <- function(x){
#   # Calculate mean and standard deviation
#   mu <- mean(x)
#   sigma <- sd(x)
#   
#   # Normalize by standard deviation
#   return((x - mu)/sigma) 
# }

# Student t distribution
tdistro <- function(X, Mu, sigma, a, b){
  sigma = sigma^2
  -1 * sum(( ((2*a+1.)/2.) * log(b + ((X - Mu)^2.)/(2.*sigma)) + .5 * log(sigma) ),na.rm = TRUE)
}   

t_dis <- function(X, Mu, sigma, a, b){
  (  (b + ((X - Mu)^2.)/(2.*sigma^2))^(-(2*a+1.)/2.) ) *  (sigma) 
} 

# Data loader
load_file_from_folder <- function(file_name, folder) {
  message(paste0("Getting data from '", file_name, "' file ...\n"))
  # Check for .csv extension
  csv_path <- paste0(folder,'/', file_name, ".csv")
  if (file.exists(csv_path)) {
    
    fil <- read.csv(csv_path)
    
    # Validate columns and headers
    if(ncol(fil) != 3 ) {
      stop("CSV file does not have 3 columns")
    } 
    if(!identical(names(fil), c("Depth", "Age", "ProxyValue"))) {
      names(fil) <- c("Depth", "Age", "ProxyValue") 
      # warning("CSV file headers modified to match required format")
    }
    
    return(fil)
    
  }else{
    # Check for .txt extension
    txt_path <- paste0(folder,'/', file_name, ".txt") #paste0(folder, file_name, ".txt")
    if (file.exists(txt_path)) {
      
      fil <- read.table(txt_path, header = TRUE, sep = "\t")
      
      # Validate columns and headers
      if(ncol(fil) != 3) {
        stop("TXT file does not have 3 columns")
      }
      if(!identical(names(fil), c("Depth", "Age", "ProxyValue"))) {
        names(fil) <- c("Depth", "Age", "ProxyValue")
        # warning("TXT file headers modified to match required format")
      }
    }
    return(fil)
  }
  
  # Return NULL with warning if no matching file
  warning("No matching file found for the given input in the specified folder.")
  return(NULL)
  
}

# Load twalk
source_twalk <- function(folder="~/Github/BSync/") {
  # Construct the path to the twalk.R file in the specified folder
  # twalk_path <- paste0(folder, "~/GitHub/BSynch/Rctwalk.R") # change this 
  twalk_path <- paste0(folder, "Rcpp_twalk.cpp") # change this
  # message("Loading 'twalk.R' from ", folder, ".\n")
  
  # Check if the twalk.R file exists in the specified folder
  if (file.exists(twalk_path)) {
    
    # source(twalk_path)
    Rcpp::sourceCpp(twalk_path)
    message("Successfully loaded 'twalk.R' from ", folder, " directory.\n")
  } else {
    warning("File twalk.R was not found in the specified folder.")
    break
  }
}

# Load or install dependencies
load_or_install_dependencies <- function() {
  if (!require(KernSmooth, quietly = TRUE)) {
    message("KernSmooth not found. Installing...\n")
    install.packages("KernSmooth", dependencies = TRUE)
    
    # Load after installing
    library(KernSmooth)
    message("KernSmooth installed and loaded successfully!\n")
  } else {
    message("KernSmooth loaded successfully!\n")
  }
  if (!require(coda, quietly = TRUE)) {
    message("coda not found. Installing...\n")
    install.packages("coda", dependencies = TRUE)
    
    # Load after installing
    library(coda)
    message("coda installed and loaded successfully!\n")
  } else {
    message("coda loaded successfully!\n")
  }
}

# target function as preparation
target_density <-function(tar_ages,tar){
  
  # Initialize an empty list to store the density objects for each depth
  kde_list <- list()
  
  # Loop through each column (depth) to calculate the kernel density
  for (col in 1:length(tar_ages)) {
    values_at_depth <- tar[,col]
    # Compute kernel density and add to the list
    kde <- density(values_at_depth, kernel = 'gaussian' )
    kde_list[[col]] <- kde
  }

  return(kde_list)
}

# Plot the density
density_plot <- function(tar_ages, tar_mt,xlabel,ylabel = "proxy units",
                         add = FALSE,axis=TRUE,flip=FALSE,
                         ylim=TRUE,xlim=TRUE){
  tar_ages_bw <- tar_ages[2] - tar_ages[1] 
  
  if (!add){
    if(axis){
      if (flip){
        if(ylim[1]==TRUE){
          plot(colMeans(tar_mt),tar_ages , type='l', xlab=ylabel, ylab=xlabel, ylim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1))    
        }else{
          plot(colMeans(tar_mt),tar_ages , type='l', xlab=ylabel, ylab=xlabel, ylim=ylim, xlim=xlim,col=rgb(0,0,0,1))  
        }
      }else{
        if(ylim[1]==TRUE){
          plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1))  
        }else{
          plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=xlim,ylim=ylim,col=rgb(0,0,0,1))  
        }
      }
      
    }else{
      if (flip){
        if(ylim[1]==TRUE){
          plot(colMeans(tar_mt),tar_ages , type='l', ylab=xlabel, xlab=ylabel, ylim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1),xaxt = 'n')
        }   else{
          plot(colMeans(tar_mt),tar_ages , type='l', ylab=xlabel, xlab=ylabel, ylim = ylim, xlim = xlim ,col=rgb(0,0,0,1),xaxt = 'n')
        }
      }else{
        if(ylim[1]==TRUE){
          plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1),xaxt = 'n')  
        }else{
          plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=c( tar_ages[1],tail(tar_ages,1)),ylim=ylim,col=rgb(0,0,0,1),xaxt = 'n')   
        }
        
      }
    }
  }else{
    if (flip){
      lines(colMeans(tar_mt),tar_ages ,col=rgb(0,0,0,1)) 
    }else{
      lines(tar_ages,colMeans(tar_mt) ,col=rgb(0,0,0,1))  
    }
  }
  
  for(i in 1:ncol(tar_mt) ){
    h <- hist(tar_mt[,i], plot = FALSE,breaks=150)
    cols <- gray(1-h$counts/max(h$counts),alpha = .4)
    # Plot non-zero rects 
    if (flip){
      rect(ybottom = tar_ages[i], ytop = tar_ages[i] + tar_ages_bw,
           xleft = h$breaks[-1], #head(h$breaks, -1),
           xright = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
    }else{
      rect(xleft = tar_ages[i], xright = tar_ages[i] + tar_ages_bw,
           ybottom = h$breaks[-1], #head(h$breaks, -1),
           ytop = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
    }
  }
}

# bypass filter
create_bypass_filter <- function(x, y, thresholds){
    x_t = x

    x_by_year <- seq(x_t[1,1],tail(x_t[,1],1),1 )
    x = approx(x_t[,1],x_t[,2],x_by_year )$y
    
    x_org_y <- approx(x_t[,1],x_t[,2],y[,1] )$y

    # Initialize variables
    max_corr <- -Inf
    best_cutoff <- NA
    
    # Loop through cutoff frequencies
    for(cutoff in thresholds){
      # Bandpass filter x
      x_filtered <- dplR::pass.filt(x, W=cutoff, type="low")
      
      x_filtered_t <- approx(x_by_year,x_filtered,y[,1])$y
    
      # Calculate correlation with y 
      corr <- cor(x_org_y-x_filtered_t, y[,2] )

      # Track if new max correlation
      if(corr > max_corr){
        max_corr <- corr
        best_cutoff <- cutoff
      }
    }
    
    # Return best cutoff 
    return(best_cutoff)
  }
  
# Extract slopes from input
input_slopes <- function(inp, breaks) {
  # Extract x and y values from the first and second columns of inp
  x_values <- inp[[1]]
  y_values <- inp[[2]]
  
  # Initialize a vector to store slopes
  slopes <- numeric(length(breaks) - 1)
  
  # Loop through the breaks and calculate the slopes between adjacent break points
  for (i in 1:(length(breaks) - 1)) {
    # Find the closest x-values and corresponding y-values at the breaks
    x1 <- breaks[i]
    x2 <- breaks[i + 1]
    
    # Use linear interpolation to find the corresponding y-values
    y1 <- approx(x_values, y_values, x1)$y
    y2 <- approx(x_values, y_values, x2)$y
    
    # Calculate the slope between these two points
    slopes[i] <- (y2 - y1) / (x2 - x1)
  }
  
  # Return the vector of slopes
  return(slopes)
}

## BSynch function

BSynch <- function(Input,Target,folder = '~/Documents/BSync/',
                   shape_acc = 1.5,  mean_acc = 50,
                   strength_mem = 10, mean_mem = .5,
                   n_sections = 50, 
                   sd_shape = 1.5 ,sd_scale = .01,
                   tao_mean=TRUE , tao_sd=TRUE ,
                   t_last_mean = TRUE, t_last_sd = 1000,
                   gw_z=.1,
                   depth_cm = TRUE, age_kyr=TRUE,T_age_kyr =FALSE,
                   ta = 3,tb = 4,
                   uq =FALSE,
                   double_target = FALSE,
                   depth_to_age = TRUE,
                   range_f =1, normalizar = TRUE, symmetric=TRUE,
                   quantile_val = .95,lquantile=.025,uquantile = .95,
                   section_mean_location = FALSE,
                   iters = 3e+3, burn = 5e+4,thin = 100,
                   continue_run = TRUE, verify_orientation = TRUE,
                   flip_orientation = FALSE,
                   savefig =FALSE,
                   bw_filter_hi_corr = FALSE, filter_cuyoff = TRUE,
                   last_tiepoint = TRUE,
                   twalk_fol = "~/Github/BSync/",
                   run = TRUE,
                   prior_acc_mean=TRUE,
                   plot_orig = TRUE, 
                   col.interv = rgb(30/255, 144/255, 255/255,.55),
                   plot.tar_lines = FALSE,
                   col.tar = rgb(255/255, 127/255, 14/255,.55), 
                   col.tar1 = rgb(0/255, 100/255, 0/255,.15),
                   col.tar2 = rgb(214/255, 39/255, 40/255,.15), 
                   col.ini_model = rgb(128/255, 0/255, 128/255,.55),
                   grid_col = "lightgray"
){ 
  #### Define the re-scaling function ####
  
  range <- function(x) {
    if (symmetric){
      q1 <- quantile(x, 1- (1- quantile_val)/2 )
      q2 <- quantile(x, (1- quantile_val)/2)
    } else{
      q1 <- quantile(x, uquantile )
      q2 <- quantile(x, lquantile)
    }
    2 * (x - q2) / (q1 - q2) - 1
  }
  
  # normalize the records 
  if (range_f==1){
    range_T <- range
  }else{
    range_T <- scale
  }
  
  #### Load packages ####
  load_or_install_dependencies()
  
  #### Load data and twalk ####
  source_twalk(twalk_fol)

  #### Load input ####
  inp <- load_file_from_folder(Input, folder) 
  # Remove NAs
  # inp <- inp[complete.cases(inp), ]
  
  # Get the prior accumulations from older age model if prior_acc_mean is true
  
  if (prior_acc_mean && !any(is.na(inp[[2]])) && !any(is.na(inp[[1]])) ){
    breaks <- seq(inp[[1]][1],tail(inp[[1]],1), length.out = n_sections+1)
    mean_acc <- input_slopes(inp, breaks)  
      # If depth is given in meters, convert it to cm by div by 100
    if (!depth_cm) {
      mean_acc <- mean_acc / 100
    }
    
    # If age is given in kyr (thousands of years), convert to years by multiplying by 1000
    if (age_kyr) {
      mean_acc <- mean_acc * 1000
    }
  
  }
  

  # Note: create the object inp which only containes the variable which will be alinge and the scale either depth or age
  if (depth_cm){
    depth <- inp$Depth  
  }else{
    inp$Depth <- inp$Depth * 100
    depth <- inp$Depth 
  }
  
  if (depth_to_age){
    org_time <- inp$Age
    inp <- data.frame(X=inp$Depth,ProxyValue = inp$ProxyValue)
    xlabel = "Depth"
    if (age_kyr){
      org_time <- org_time * 1000
    }
    
  }else{
    org_time <- inp$Age
    inp <- data.frame(X=inp$Age,ProxyValue = inp$ProxyValue)
    if (age_kyr){
      age_temp <- inp$X * 1000
      inp$X <- inp$X * 1000
      org_time <- org_time * 1000

    }else{
      age_temp <- inp$X
    }
    xlabel = "Age"
  }
  



  #### Load target ####
  # Note: create the object inp which only containes the variable which will be alinge and the scale either depth or age
  if (uq){
    tar <- read.table(paste0(folder,Target,'.txt'),header = T)
    tar <- as.matrix(tar)
    # Divide the data base into the ages and the iterations
    tar_ages <- tar[,1] # save ages
    tar <- t(tar[,-1]) # invert the matrix
    tar_mt <- tar
    tar_mt <-  range_T(tar)
    
    # localizador is the function which tell the target_densities function which density to use
    localizador <- function(x)(approx(tar_ages,seq(1,length(tar_ages)),method = 'constant',xout = x)$y)
    
    # Create a data frame with depth and mean
    # tar <- data.frame(X = as.numeric(tar_ages), ProxyValue = as.numeric(colMeans(tar_mt)) )
    tar <- data.frame(X = as.numeric(tar_ages), ProxyValue = apply(tar_mt, 2, median) )
    
  }else if (double_target){
    # Read the CSV file
    tar <- read.csv(paste0(folder,Target,'.csv'), header = T)
    # Rename the columns
    ori_col_names <- colnames(tar)
    colnames(tar) <- c('Age', 'target1', 'target2')
    # create the data frame
    tar <- data.frame(X=tar$Age,ProxyValue1 = tar$target1,ProxyValue2 = tar$target2, ProxyValue = c(.47*tar$target1 + .53*tar$target2) )
    # save target ages for future use
    tar_ages <- tar[,1]
  }else { # this is for the simple version
    tar <- load_file_from_folder(Target, folder) 
    ori_col_names <- colnames(tar)
    tar <- data.frame(X=tar$Age,ProxyValue = tar$ProxyValue)
    tar_ages <- tar[,1]
  }
  
  

  
  if (T_age_kyr){
    tar_ages <- tar$X * 1000
    tar$X <- tar$X * 1000
  }else{
    tar_ages <- tar$X
  }

  #### Check that the records are properly aligned ####
  if (flip_orientation){

      inp$ProxyValue <- -inp$ProxyValue

  }else{

    if (verify_orientation){
      inv_user = TRUE
      while(inv_user){

        par(mfrow=c(1,1))
        if (uq){
          density_plot(tar_ages ,tar_mt, xlabel=xlabel )
          axis(1, at=NULL, labels=FALSE, tick=FALSE)
          legend("bottomright", lty=rep(1,2),
                 c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
        }else{
          if (double_target){
            plot(tar$X,range_T(tar$ProxyValue1), type='l', xlab=xlabel, ylab="proxy units",
                 xlim=c( tar$X[1],max(tar$X)), col='gray1' )
            lines( tar$X,range_T(tar$ProxyValue2) ,col='gray' )
            legend("topright", lty=rep(1,1,1),
                   c("input", ori_col_names[3], #"target2",
                     ori_col_names[2]#"target1"
                   ), col=c(rgb(0,0,1,.9),'gray1','gray'), cex=0.5)
          }else{
            plot(tar, type='l', xlab=xlabel, ylab="proxy", xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))
            legend("bottomright", lty=rep(1,2),
                   c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
          }
        }

        par(new=TRUE)
        plot(inp$X,range_T(inp$ProxyValue), type='l', col=rgb(0,0,1,.9), xlab="", ylab="", axes=FALSE,
             xaxt='n', yaxt='n', xlim=c( inp$X[1],max(inp$X)))
        # plot(inp$X,inp$ProxyValue, type='l', col=rgb(0,0,1,.9), xlab="", ylab="", yaxt='n', xlim=c( inp$X[1],max(inp$X)))
        # legend("bottomright", lty=rep(1,2),
        #       c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
        par(new=TRUE)
        if (uq){
          density_plot(tar_ages ,tar_mt, xlabel=xlabel ,axis=FALSE)
        }else{
          plot(tar, type='l', xlab=xlabel, ylab="proxy units", xaxt='n',xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))
        }


        message('This are the two records side to side\n')
        y <- readline("Shall I invert the input? yes or not \n")
        # Check response
        if (y %in% c("y", "yes")) {
          inp$ProxyValue <- -inp$ProxyValue
        } else if (y %in% c("n", "no")) {
          # Default no inversion
        } else {
          print("Invalid input")
        }

        # Always re-plot with new inversion
        # par(mfrow=c(1,1))
        par(new = FALSE)
        plot(inp$X,inp$ProxyValue, type='l', col=rgb(0,0,1,.9), xlab="", 
             ylab="", yaxt='n', xlim=c( inp$X[1],max(inp$X)) )
        
        legend("bottomright", lty=rep(1,2),
              c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
        par(new=TRUE)
        if (uq){
          density_plot(tar_ages ,tar_mt, xlabel=xlabel,axis=FALSE )
        }else{
          if (double_target){
            plot(tar$X,range_T(tar$ProxyValue1), type='l', xlab=xlabel, ylab="proxy units",
                 xlim=c( tar$X[1],max(tar$X)), col='gray1' )
            lines( tar$X,range_T(tar$ProxyValue2) ,col='gray' )
            legend("topright", lty=rep(1,1,1),
                   c("input", ori_col_names[3], #"target2",
                     ori_col_names[2]#"target1"
                   ), 
                   col=c(rgb(0,0,1,.9),'gray1','gray'), cex=0.5)
          }else{
            plot(tar, type='l', xlab=xlabel, ylab="proxy", xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))
            legend("bottomright", lty=rep(1,2),
                   c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
          }
        }
        axis(4, col = 4)

        # Ask if records aligned
        inv_again <- readline("Are records aligned properly? Enter y/yes or n/no: \n")

        # Check response
        if (inv_again %in% c("y", "yes")) {
          inv_user <- FALSE
        } else if (inv_again %in% c("n", "no")) {
          inv_user <- TRUE
        } else {
          print("Invalid input")
        }
      }
    }
  
  }
  
  # obligatory plot
  if(!double_target){
    if(uq){
      density_plot(tar_ages ,tar_mt, xlabel=xlabel )
      # axis(1, at=NULL, labels=FALSE, tick=FALSE)
    }else{
      plot(tar$X,range_T(tar$ProxyValue), type='l', xlab=xlabel, ylab="proxy units", xlim=c( tar$X[1],max(tar$X)), 
           col=rgb(0,0,0,1))      
    }
    # plot(tar$X,range_T(tar$ProxyValue), type='l', xlab=xlabel, ylab="proxy units", xlim=c( tar$X[1],max(tar$X)), 
    #      col=rgb(0,0,0,1))
    par(new=TRUE)
    plot(inp$X,range_T(inp$ProxyValue), type='l', col=rgb(0,0,1,.9), xlab="", ylab="", axes=FALSE,
         xaxt='n', yaxt='n', xlim=c( inp$X[1],max(inp$X)))
    legend("bottomright", lty=rep(1,2),
           c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
    # legend("bottomright", lty=rep(1,2),
    #           c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
  }else{
    plot(tar$X,range_T(tar$ProxyValue1), type='l', xlab="", ylab="proxy units",
         xlim=c( tar$X[1],max(tar$X)), col='gray1' )
    lines( tar$X,range_T(tar$ProxyValue2) ,col='gray' )
    par(new=TRUE)
    plot(inp$X,range_T(inp$ProxyValue), type='l', col=rgb(0,0,1,.9), xlab="", ylab="", axes=FALSE,
         xaxt='n', yaxt='n', xlim=c( inp$X[1],max(inp$X)))
    legend("topright", lty=rep(1,1,1),
              c("input", ori_col_names[3], #"target2",
                ori_col_names[2]#"target1"
                ), col=c(rgb(0,0,1,.9),'gray1','gray'), cex=0.5)

  }


  #### Bypass filter modification ####
  if (bw_filter_hi_corr){
    t_freq <- seq(as.integer(org_time[1]),as.integer(tail(org_time,1)*.4 ),100 )
    inp_f <- data.frame(Age = org_time,inp$ProxyValue)
    if (filter_cuyoff == TRUE){
      message('Finding and applying the best cutoff for a low filter')
      bypass_results <- create_bypass_filter(tar, inp_f, t_freq)   
    }else{
      bypass_results <- filter_cuyoff
    }
    
    
    x_by_year <- seq(tar[1,1],tail(tar[,1],1),1 )
  
    if (uq){
      message('Applying filter to the target')
      for( i in 1:dim(tar_mt)[1]){
        test_serie_tar = approx(tar$X,tar_mt[i,],x_by_year )$y
        fil_tar <- dplR::pass.filt(test_serie_tar, W=bypass_results, type="low", method="Butterworth")
        tar_mt[i,] <-approx(x_by_year,test_serie_tar-fil_tar,tar$X )$y
      }
      test_serie_tar = approx(tar$X,tar$ProxyValue,x_by_year )$y
      fil_tar <- dplR::pass.filt(test_serie_tar, W=bypass_results, type="low", method="Butterworth")
      tar$ProxyValue <- approx(x_by_year,test_serie_tar-fil_tar,tar$X )$y
    }else{
      test_serie_tar = approx(tar$X,tar$ProxyValue,x_by_year )$y
      fil_tar <- dplR::pass.filt(test_serie_tar, W=bypass_results, type="low", method="Butterworth")
      tar$ProxyValue <-approx(x_by_year,test_serie_tar-fil_tar,tar$X )$y
  }
  message('Optimal frequency')
  message(bypass_results)
  message('Finish applying the filter to all the simulations ')
  }
  #### Checks that target is longer than the input ####
  # This in only necessary on the age to get alignent 
  if (!depth_to_age){
    # checks that target is longer than the input
    # Extract the scale columns
    inp_scale <- inp[,1]
    tar_scale <- tar[,1]
    
    # Check if inp scale range is within tar scale range
    if(min(inp_scale) < min(tar_scale) | 
       max(inp_scale) > max(tar_scale)) {
      
      # Trim inp to be within tar scale range
      org_time <- org_time[inp$X <= tail(tar[,1], 1) - 10 & inp$X >= head(tar[,1], 1) + 10]
      inp <- inp[inp$X <= tail(tar[,1], 1) - 10 & inp$X >= head(tar[,1], 1) + 10, ]

      message("Input dataset trimmed to be within target scale range")
    } else {
      message("Input dataset already contained within target scale range")
      breaks <- seq(inp$X[1],tail(inp$X,1), length.out = n_sections+1)
      b_length <- breaks[2] - breaks[1]
      tar <- tar[tar[, 1] >= inp[1, 1] - 1.5 * b_length & tar[, 1] <= tail(inp[, 1],1 ) + 1.5 * b_length, ]
      tar_ages <- tar$X
    }
  }
  
  #### Set variables ####
  # normalize the records 
  if (normalizar){
    inp[,2] <- range_T(inp[,2]) 
    tar$ProxyValue <- range_T(tar$ProxyValue) 
  }
  
  # Define the variables which are use for the alignment 
  breaks <- seq(inp$X[1],tail(inp$X,1), length.out = n_sections+1)
  b_length <- breaks[2] - breaks[1]
  window_labels <- as.integer(cut(inp$X, breaks, include.lowest = TRUE, labels = FALSE))
  mean_values <- as.numeric(tapply(inp$ProxyValue, window_labels, mean, na.rm = TRUE))

  # tau0 variables
  if (tao_mean == TRUE){
    if (depth_to_age){
      tao_mean <- min(tar_ages) 
      # uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
      tao_sd <-  ifelse(tao_mean  < 10 * 1000, 500, 5000)  
      # tao_mean <- min(tar_ages) + 1 * tao_sd

    }else{
      tao_mean <- min(inp$X) 
      # uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
      tao_sd <-  ifelse(head(inp$X,1) < 10 * 1000, 500, 5000)   
    }

    
  }else{ # Check if tao_mean  and tao_sd is numeric
    if(!is.numeric(tao_mean)) {
      tao_mean <- as.numeric(readline("Enter a numeric value for tao_mean: ")) 
    }
    if(!is.numeric(tao_sd)) {
      tao_sd <- as.numeric(readline("Enter a numeric value for tau_sd: "))
    }
  }
  
  # parameter for prior of alphas 
  
  if (length(mean_acc)==1){
    scale_acc = mean_acc / shape_acc  
  }else{
    # Standard deviation is provided as a constant
    sd_age_prior <- 500 / diff(breaks)[1]

    # Calculate the shape parameter k for each value in the mean_acc vector
    shape_acc <- (mean_acc / sd_age_prior)^2

    # Calculate the scale parameter theta for each value in the mean_acc vector
    scale_acc <- (sd_age_prior^2) / mean_acc
  }


  
  
  
  # parameters for memory 
  m_alpha = strength_mem * mean_mem 
  m_beta =  strength_mem * (1 - mean_mem) 
  
  # set lower and upper prior boundaries
  # note that when age to age alignment alphas will come from a lognormal (to do)
  if (depth_to_age){
    low <- c(tar_ages[1],#tar$X[1] - (3*tao_sd), # tao0
             0, # mem
             rep(0.0 , n_sections)) # alphas
    up <- c( tar_ages[1] + tao_sd,  #tao0
             1, # mem,
             rep(Inf, n_sections)) # alphas
  }else{
    low <- c(tar_ages[1], #tar$X[1] - (3*tao_sd), # tao0
             0, # mem
             # rep(.1, n_sections-1)) # arate
             rep(.25, n_sections)) # arate

    up <- c( inp$X[1]+3*tao_sd,#tar_ages[1] + (3*tao_sd),  #tao0
             1, # mem,
             rep(4, n_sections)) # alphas
             # rep(10, n_sections-1)) # alphas 
  }
  if (double_target){
    low <- c(low,0)
    up <- c(up,1)
  }
  
  mx_age <- max(tar_ages)

  if (last_tiepoint==TRUE){
    last_tiepoint = -Inf
  }
  
  
  #print(mx_age)
  # new sd (geometric mean of sd of target and input)
  # my_sd = .15 #sqrt(sd_input^2 +sd_target^2)

  
  if (uq){
    # create the density function.
    kde_list <- target_density(tar_ages,tar_mt) # target
    # Extract bandwidths using lapply
    bandwidths <- lapply(kde_list, function(x) x$bw)
    
    # result as a vector
    bandwidths <- unlist(bandwidths)
  }
  
  #### Functions ####
  # Loads the C++ functions
  Rcpp::sourceCpp("~/Github/BSync/targetDensity.cpp")
  
  # Function to ajust means
  # Cut ages into breaks and calculate mean per window in tar
  # Note that tar groups can be calculated before hand
  # breaks_tar = seq(tar$X[1],tail(tar$X,1),length.out= 4)#as.integer(n_sections/10) )
  # breaks_tar = c(tar$X[1], 10000, 20000,35000  ,tail(tar$X,1))
  # tar$age_group <- cut(tar$X, breaks = breaks_tar, include.lowest = TRUE)
  # mean_tar <- aggregate(ProxyValue ~ age_group, data = tar, mean)
  # 
  # adjust_means <- function(tar  , inp, inp_ages , breaks ,mean_tar  ) {
  #   # Apply the same breaks to inp
  #   inp$age_group <- cut(inp_ages, breaks = breaks, include.lowest = TRUE)
  #   mean_inp <- aggregate(ProxyValue ~ age_group, data = inp, mean)
  #   
  #   # Merge mean_tar with mean_inp to align means
  #   merged_data <- merge(mean_tar, mean_inp, by = "age_group", suffixes = c(".tar", ".inp"))
  # 
  #   # Calculate the difference between the means in tar and inp
  #   merged_data$difference <- merged_data$ProxyValue.tar - merged_data$ProxyValue.inp
  #   # Create a lookup table for the differences
  #   difference_lookup <- setNames(merged_data$difference, merged_data$age_group)
  #   # Adjust inp data points to match means in tar windows
  #   adjusted_inp <- inp$data + difference_lookup[as.character(inp$age_group)]
  #   return(as.numeric(adjusted_inp))
  # }
  
  
  # tau function (the function which gives the value in the target scale)
  tau <- function(x, param){
    tau0 <- param[1]
    # w <- param[2]
    slopes <- param[-c(1,2)]
    tau_is <- c(tau0, tau0 + cumsum(slopes * b_length))
    return( approx(breaks,tau_is,x)$y )
  }
  
  # Provides the alphas for calculation its priors
  n_slopes = n_sections 
  alphas <- function(param){
    # tau0 <- param[1]
    w <- param[2]
    slopes <- param[-c(1,2)]
    alf <- ( slopes[-n_slopes] - w * slopes[-1] ) / (1 - w)
    return(c(alf, tail(slopes,1)) )
  }
  
  # Prior density
  logprior <-  function(param){
    tao0 <- param[1]
    w <- param[2]
    # Tao0 prior
    l_prior <- dnorm(tao0, mean = tao_mean, sd = tao_sd,log = TRUE)
    # memory value
    l_prior <- l_prior + dbeta(w, m_alpha, m_beta, log = TRUE)
    # denisity for acc rates 
    ms_acc <- alphas(param)
    l_prior <- l_prior + sum( dgamma(ms_acc, shape = shape_acc, scale = scale_acc, log = TRUE) )
    return(l_prior)
  }
  
  # likelihood in case of UQ
  if (t_last_mean == TRUE){
    multi_last = 0
  }else{
    multi_last = 1
  }
  
  # Likelihood functions
  if (uq){
     n_kde <- length(tar_mt[,1])
     sd_convertor <- 1/(1.06*n_kde^(-1/5))
     
     loglikelihood_uq <- function(params){
       # Get the ages at which to calculate the target
       # t_times <- tau(inp$X,params)
       t_times <- tauC(inp$X,params,b_length,breaks )
       # this calculates the likelihood of the last point
       l_last <- multi_last * dnorm(tail(t_times,1),t_last_mean,t_last_sd,log=T)
       # get the values of the target at the require times.
       loca <- as.integer(localizador(t_times))
       
       # if (section_mean_location){
       #   inp_proxy <- adjust_means(tar,inp, t_times ,breaks_tar,mean_tar )  
       # }else{
       #   inp_proxy <- inp$ProxyValue
       # }
       
      
       ll <- loglikelihood_uqCpp(params, 
                                 # inp_proxy,
                                 inp$ProxyValue,
                                 bandwidths, tar_mt, sd_convertor ,loca)

       ll <- ll + l_last
       return(ll)
     }

     obj <- function(param){ 
       -(logprior(param)+loglikelihood_uq(param)) 
       }
     
   }else if (double_target){
       loglikelihood <- function(params,tar_prox,my_sd){
         
         t_times <- tauC(inp$X, params, b_length, breaks )
         l_last <- multi_last * dnorm(tail(t_times,1),t_last_mean,t_last_sd,log=T)
  
         ll <- loglikelihoodC(params,
                        tar_prox,inp$X,inp$ProxyValue,tar$X,
                        # tar_prox,inp$X, inp_proxy, tar$X,
                        my_sd, ta, tb, 
                        b_length, breaks)
         ll <- l_last + ll
         return(ll)
         }
       
       obj <- function(param){ 
         my_sd = tail(param,1)
         param0 =  param[-length(param)]
         # this is the 'prior' of the sd
         p_sd <- dgamma(my_sd,shape = sd_shape,scale = sd_scale,log = T) 
         omega_w <- tail(param0,1)# param[length(param)]
         # We create the new target
         ProxyValue <- range_T(  omega_w    * tar$ProxyValue1 + 
                               (1 - omega_w) * tar$ProxyValue2 )
         # Remove the parameter \omega which generates the target
         param1 <- param[-length(param0)]
         return( -(logprior(param1) + loglikelihood(param1,ProxyValue,my_sd) + p_sd)  )
       }
     }else{ # This is for the simple 
       loglikelihood <- function(params,tar_prox=tar$ProxyValue,my_sd){
         t_times <- tauC(inp$X,params,b_length,breaks )
         l_last <- multi_last * dnorm(tail(t_times,1),t_last_mean,t_last_sd,log=T)
         
         # if (section_mean_location){
         #   inp_proxy <- adjust_means(tar,inp, t_times ,breaks_tar,mean_tar )  
         # }else{
         #   inp_proxy <- inp$ProxyValue
         # }
         
         ll <- loglikelihoodC(params,
                        # tar_prox,inp$X, inp_proxy, tar$X, # this like was to test the moving mean
                        tar_prox,inp$X,inp$ProxyValue,tar$X,
                        my_sd, ta, tb, 
                        b_length, breaks)
         ll <- l_last + ll
         return(ll)
       }
       
       obj <- function(param){
         my_sd = tail(param,1)
         param =  param[-length(param)]
         p_sd <- dgamma(my_sd,shape = sd_shape,scale = sd_scale,log = T) # this gives almost 95% of prob bertween 0 and .1
         
         -( logprior(param) + loglikelihood(param,my_sd = my_sd) +p_sd)
       }
     }
  
  
  # support function
  if (double_target){
    supp <- function(params){
      my_sd = tail(params,1)
      params =  params[-length(params)]
      
      alp <- alphas(params[-length(params)])
      age_lim <- tauC(tail(inp$X,1),params[-length(params)],b_length,breaks)
    
      
      ifelse( all(params > low) & all(params < up)  & 
               all( alp > 0 ) &   my_sd >0 & my_sd < 0.3 &
               age_lim < mx_age ,
             return(TRUE), 
             return(FALSE)
      )
    }
  }else if (uq){
      supp <- function(params){
        alp <- alphas(params)
        age_lim <- tauC(tail(inp$X,1),params,b_length,breaks)
        
        ifelse(all(params > low) & all(params < up)  & 
                 all( alp > 0 ) & 
                 age_lim < mx_age & 
                 age_lim > last_tiepoint,
               return(TRUE), 
               return(FALSE)
        )
      }
    }else{
      # support of simple aliment
      supp <- function(params){
        my_sd = tail(params,1)
        params =  params[-length(params)]
        alp <- alphas(params)
        age_lim <- tauC(tail(inp$X,1),params,b_length,breaks)
        ifelse(all(params > low) & all(params < up)  & 
                 all( alp > 0 ) &   my_sd >0 & my_sd < 1 &
                 age_lim < mx_age & 
                 age_lim > last_tiepoint,
               return(TRUE), 
               return(FALSE)
        )
      }
    }
  
  
  # generates initial points.
  sampler <-  function(rerun=FALSE){
    # generate to0
    if (depth_to_age ){
      if (any(is.na(org_time)) ){
        tao0 <- runif(1, low[1], low[1]+1) 
      }else{
        if(org_time[1] >= low[1]){
          tao0 <- runif(1, org_time[1], org_time[1]+1) 
        }else{
          tao0 <- runif(1, low[1], low[1]+1) 
        }
      }
    }else{
      tao0 <- runif(1, inp$X[1],inp$X[1]+10)
    }
    # memory value
    w <-  runif(n = 1, min = 0, max = 1e-5)
    # Generate slopes
    if(rerun){
      if (depth_to_age){
        if (any(is.na(org_time)) ){
          tmp_mean <- floor( (tail(tar$X,1)- tar$X[1])/(tail(inp$X,1)-inp$X[1]) )
          ms <- rgamma(1, scale = tail(tmp_mean,1) / 1e+6, shape =1e+6 )
        }else{
          tmp_mean <- diff(approx(depth, org_time, xout = breaks, rule = 2)$y) / diff(breaks)
          ms <- rgamma(length(tmp_mean), scale = tmp_mean / 1e+3, shape =1e+3 )
        }
        if (length(tmp_mean) == 1){        
          for (i in 1:(n_slopes-1)){
            alpha = rgamma(1, scale = tmp_mean/shape_acc, shape = shape_acc )
            ms <- c(  w * ms[1] + (1-w) * alpha, ms)
            }
        }
      }else{
        tmp_mean <- floor( (tail(tar$X,1)- tar$X[1])/(tail(inp$X,1)-inp$X[1]) )
        ms <- rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )

        for (i in (n_slopes-1):1){
          alpha = rgamma(1,scale = scale_acc, shape = shape_acc )
          ms <- c(ms, w * ms[1] + (1-w) * alpha )
          }
      }
    }else{
      if (depth_to_age){
        if (any(is.na(org_time)) ){
          tmp_mean <- floor( (tail(tar$X,1)- tar$X[1])/(tail(inp$X,1)-inp$X[1]) )
          ms <- rgamma(1, scale = tail(tmp_mean,1) / 1e+3, shape =1e+3 )
        }else{
          tmp_mean <- diff(approx(depth, org_time, xout = breaks, rule = 2)$y) / diff(breaks)
          ms <- rgamma(length(tmp_mean), scale = tmp_mean / 1e+3, shape =1e+3 )
        }
      }else{
        tmp_mean <- floor( (tail(tar$X,1)- tar$X[1])/(tail(inp$X,1)-inp$X[1]) )
        ms <- rgamma(1, tail(tmp_mean,1) / 1e+3, shape =1e+3 )
      }
      # Here we generate the rest of the slopes
      if (length(tmp_mean) == 1){
        for (i in (n_slopes-1):1){
          alpha = rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )
          ms <- c(  w * ms[1] + (1-w) * alpha, ms)
        }
      }
    }
    
    if (double_target){
      return(c(tao0, w, ms,runif(1),runif(1,.05,.5)))
    }else{
      if (uq){
        return(c(tao0, w, ms))  
      }else{
        return(c(tao0, w, ms,runif(1,.05,.5)))  
      }
    }
  }
  
  initial_search <- function(){
    x1 <- sampler()
    while(!(supp(x1))){
      x1 <- sampler(rerun = TRUE)  
    }
    return(x1)
  }
  
  
  #### Run the t-walk ####
  
  if (run){
    # Check if the file exists
    if (file.exists(paste0(folder,"twalkstate.csv"))) {
      # Check if continue_run is TRUE
      if (continue_run) {
        # Read the data
        loaded_matrix <- as.matrix(read.csv(paste0(folder,"twalkstate.csv"), header = TRUE))
        # Extract the values from the matrix to x1 and x2
        x1 <- loaded_matrix[1,]
        x2 <- loaded_matrix[2,]
        message('Previous run was loaded.\n')
        
      } else {
        message("Searching for initial values...")
        x1 <- initial_search()
        x2 <- initial_search()
      }
      
    } else {
      message("Searching for initial values...")
      x1 <- initial_search()
      x2 <- initial_search()
    }
    
    message("Initiating the t-walk process...")
    burn = burn * length(x1)
    pphi_eval <-  min(length(x1), 4)/length(x1)
    # message(paste0('The Phi value is :' , pphi_eval))
    # measures time in the MCMC
    start.time <- Sys.time()
    
    # output <- Runtwalk(iters, Obj = obj, dim = length(x1), x0 = x1, xp0 = x2, Supp =supp,
    # thinning = length(x1) * thin, burnin= burn )
    
    output <- Runtwalk(length(x1), obj, supp, x1, x2, 
                       as.integer(iters),as.integer(burn),as.integer(length(x1) * thin) ,
                       pphi=pphi_eval )
    # measures time in the MCMC
    end.time <- Sys.time()
    time.taken <- as.numeric(difftime(end.time, start.time, units = "secs"))
    
    # Function to format the time taken
    format_time_taken <- function(time.taken) {
      if (time.taken < 60) {
        return(paste(round(time.taken, 2), "seconds"))
      } else if (time.taken < 3600) {
        return(paste(round(time.taken / 60, 2), "minutes"))
      } else if (time.taken < 86400) {
        return(paste(round(time.taken / 3600, 2), "hours"))
      } else {
        return(paste(round(time.taken / 86400, 2), "days"))
      }
    }
    
    formatted_time <- format_time_taken(time.taken)
    message_time <- paste("The MCMC run took", formatted_time, "to complete.")
    cat(message_time)

    
    # twalk state
    lastx = tail(output$output, 1)
    lastxp = tail(output$outputp, 1)
    # Write the matrix to a CSV file
    write.csv(matrix_to_save <- matrix(c(lastx, lastxp), nrow = 2,byrow = TRUE), 
              paste0(folder,"twalkstate.csv"), row.names=FALSE)#, col.names=FALSE)
    
    energy <- output$Us[-1]
    energy2 <- output$Ups[-1]
    
    c <- output$output[-1,] 

    write.csv(energy,paste0(folder,Input,'-',Target,'_',n_sections,"_energy.out"), row.names=FALSE)#, col.names=FALSE)
    write.csv(c, paste0(folder,Input,'-',Target,'_',n_sections,".out"), row.names=FALSE)#, col.names=FALSE)
    
    if (!uq){
      my_sd <- c[, ncol(c)]
      c <- c[,-ncol(c)]   # remove the sd column
    }
    # Extract the last column into a vector called my_sd
    
    
    
    if (double_target){
      w_tar <- c[, ncol(c)]
      output$w_tar <- w_tar
      c <- c[,-ncol(c)]  
    }
    
    # iat <- IAT(output,to=output$Tr)
    iat <- try(IAT(output, to = output$Tr), silent = TRUE)
    output$tar <- tar
    output$inp <- inp
    output$breaks <- breaks
    if (bw_filter_hi_corr){ output$cutoff <- bypass_results }
    
    
  }else{
    # load the last run data
    energy <- as.vector(read.csv(paste0(folder,Input,'-',Target,'_',n_sections,"_energy.out"),header = T))$x
    c <- read.csv(paste0(folder,Input,'-',Target,'_',n_sections,".out"),header = T)
    # Adjust the matrix
    c <- as.matrix(c)
    # Extract the last column into a vector called my_sd
    my_sd <- c[, ncol(c)]
    # Create a new matrix C without the last column
    c <- c[, -ncol(c)]
    
    output <- list(c[ ,-dim(c)[2]])
    if (double_target){
      w_tar <- c[, ncol(c)]
      output$w_tar <- w_tar
    }
    output$tar <- tar
    output$inp <- inp
    output$breaks <- breaks
    
    library(coda)
    # Convert to MCMC object
    chain <- mcmc(c)
    # Calculate effective sample size (Inverse of IAT)
    iat <- max(effectiveSize(chain))

  }  

  
  # 
  cat("\n================== IAT DIAGNOSTIC ==================\n")
  cat("IAT Value:", as.integer(iat), "\n")
  if (iat < 100) {
    cat("Interpretation: The chain exhibits low correlation among successive samples.\n")
  } else {
    cat("Interpretation: The chain exhibits high correlation among successive samples.\n")
    cat("Recommendation: Consider increasing the thinning value and rerunning the chain.\n")
  }
  
  #### Construct the age-depth or age-age assembles ####
  # Initialize matrix to store assembles

  tau_mat <- matrix(NA, nrow = nrow(c), ncol = length(breaks))

  
  # Loop through each row of c
  for(i in 1:nrow(c)) {
    # Apply tau to breaks with this row of c as params
    tau_row <- tauC(breaks, c[i,],b_length,breaks)
    # Store result in matrix
    tau_mat[i,] <- tau_row
  }
  
  
  # Initialize matrix to store assembles
  tau_mat2 <- matrix(NA, nrow = nrow(c), ncol = length(inp$X))
  
  # Loop through each row of c
  for(i in 1:nrow(c)) {
    # Apply tau to breaks with this row of c as params
    tau_row2 <- tauC(inp$X, c[i,],b_length,breaks)
    # Store result in matrix
    tau_mat2[i,] <- tau_row2
  }
  
  output$tau_mat = tau_mat
  output$tau_by_point = tau_mat2
  write.csv(output$tau_by_point,paste0(folder,Input,'-',Target,'_tau_',n_sections,".out"), row.names=FALSE)#,col.names = NA)
  
  #### Plot the results ####
  message("Commencing the plotting process...")
  if (savefig){  pdf(paste0(folder,'/aligment_',Input,'-',Target,'.pdf')) }

  if(!double_target){
    # Set up plot margins: bottom, left, top, and right
    par(mar=c(2.3,2.1,1.1,1.1) )
    # Set axis label size, title size, and axis title size
    par(cex.axis=0.9, cex.main=0.9, cex.lab=0.9)
    # Setting the axis title, axis labels, and axis line positions
    par(mgp=c(1.3, 0.5, 0))
    # # Set character size for points, lines etc.
    # par(cex=0.7)
    layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,
                    1,1,1,2,2,2,3,3,3,4,4,4,
                    1,1,1,2,2,2,3,3,3,4,4,4,
                    6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5), 12, 12, byrow = TRUE))  
  }else{
    # Set up plot margins: bottom, left, top, and right
    par(mar=c(2.3,2.1,1.1,1.1) )
    # Set axis label size, title size, and axis title size
    par(cex.axis=0.9, cex.main=0.9, cex.lab=0.9)
    # Setting the axis title, axis labels, and axis line positions
    par(mgp=c(1.2, 0.5, 0))
    # # Set character size for points, lines etc.
    # par(cex=0.7)
    layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,
                    1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,
                    1,1,1,2,2,2,3,3,3,4,4,4,7,7,7,
                    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
                    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                    5,5,5,5,5,5,5,5,5,5,5,5,5,5,5),nrow =  12,ncol =  15, byrow = TRUE))  
  }
  
  
  
  # plot energy
  {
    #par(mar=c(2.3,1.1,1.1,1.1) )
    #par(mgp=c(1.1, 0.5, 0))
    plot(energy,xlab='iterations',ylab='Log of Objective',type = 'l',
         col=rgb(169/255, 169/255, 169/255) ,axes=F)
    axis(side=1)
    axis(side=2)
    legend( 'topright',legend = paste("IAT:", as.integer(iat)),
            ,bg=NA, bty="n",cex=.75)
  }
  
  ## Plot posterior and prior of accumulations
  { 
    # Set up plot margins: bottom, left, top, and right
    par(mar=c(2.3,1.1,1.1,1.1) )
    d_m <- density(c[,-c(1,2)])
    # Get y-values of the gamma curve across the specified x-range
    y_gamma <- dgamma(seq(from=0, to=max(mean_acc) +  5*max(shape_acc *scale_acc), length.out=1000),
                      shape = shape_acc, scale = scale_acc)
    
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_gamma),max(d_m$y))) )
    # y_range <- c(5,40 )
    
    if(depth_to_age){
      low_lim_acc = 1e-20
      up_lim_acc = max(mean_acc) + 3 * max(shape_acc*scale_acc)
    }else{
      low_lim_acc = .25
      up_lim_acc = 4
    }
    
    # Create a sequence of x values within the limits
    x_values <- seq(low_lim_acc, up_lim_acc, length.out = 100)
    
    # Calculate gamma values for the sequence
    gamma_values <- dgamma(x_values, shape = max(shape_acc), scale = max(scale_acc) )

    # Calculate the maximum y value from both the gamma and the density plot to set y-axis limits
    max_y <- max(c(gamma_values, d_m$y[d_m$x > low_lim_acc & d_m$x < up_lim_acc]))
    
    if(depth_to_age){
      xlab_val="Acc Rates (yr/cm)"
    }else{
      xlab_val="Expanssion/Compression"
    }
    
    # Plot the gamma distribution as points with adjusted y-axis limits
    plot(x_values, gamma_values, type='n', ylim=c(0, max_y),
         xlab=xlab_val, ylab="", axes=F,
         main=" ", yaxt="n")
    
    # Fill area under gamma curve
    polygon(c(low_lim_acc, x_values, up_lim_acc), 
            c(0, gamma_values, 0), 
            col="lightgray", border=NA)
    
    # Fill area under density plot
    density_x <- d_m$x[d_m$x > low_lim_acc & d_m$x < up_lim_acc]
    density_y <- d_m$y[d_m$x > low_lim_acc & d_m$x < up_lim_acc]
    
    polygon(c(low_lim_acc, density_x, up_lim_acc), 
            c(0, density_y, 0), 
            col=col.interv, border=NA)
    
    # Re-draw the gamma curve on top of the filled area
    lines(x_values, gamma_values, col='lightgray')
    
    # Re-draw the density plot on top of the filled area
    lines(density_x, density_y, col=col.interv)
    
    axis(side=1)
    legend( 'topright',legend = c(
      paste0('mean_acc: ',round(mean(mean_acc) ,2) ),
      paste0('shape_acc: ',round(mean(shape_acc),2) ),
      paste0('n_sections: ',n_sections)
    ),
            ,bg=NA, bty="n",cex=.75)
    
    
  }

  ## plot memory prior vs posteriors
  {
    d_w <- density(c[,2])
    
    y_beta <- dbeta(seq(from=.01, to=.99, length.out=1000),
                    m_alpha, m_beta)
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_beta),max(d_w$y))) )

    # Assume d_m is your density data
    
    # Calculate beta values
    x_values_beta <- seq(0, 1, length.out = 100)
    beta_values <- dbeta(x_values_beta, m_alpha, m_beta)
    
    # Plot the beta distribution curve and prepare to fill under the curve
    plot(x_values_beta, beta_values, type='l', ylim=y_range, col="lightgray",
         xlab="memory", ylab="", main="", yaxt="n",axes=F)
    
    # Fill area under the beta curve
    polygon(c(x_values_beta, rev(x_values_beta)), c(rep(0, length(beta_values)), 
                        rev(beta_values)), col="lightgray", border=NA)
    
    # Fill area under the density plot for beta
    with(d_w, {
      relevant_indices <- x > 0 & x < 1
      polygon(c(x[relevant_indices], rev(x[relevant_indices])),
              c(rep(0, sum(relevant_indices)), rev(y[relevant_indices])), 
              col = col.interv, border=NA)
    })
    
    # Replot the beta curve for clarity
    lines(x_values_beta, beta_values, col='lightgray')
    
    # Replot the density plot lines for the beta distribution
    lines(d_w$x[d_w$x > 0 & d_w$x < 1],
          d_w$y[d_w$x > 0 & d_w$x < 1], col=col.interv)
    
    axis(side=1)
    legend( 'topright',legend = c(
      paste0('mean_mem: ',mean_mem),
      paste0('strength_mem: ',strength_mem),
      paste0('t_last_mean: ',t_last_mean),
      paste0('t_last_sd: ',t_last_sd)
    ),
    ,bg=NA, bty="n",cex=.75)

  }
  
  
  ## plot initial value tau0  or sd
  {
    if (uq){
      d_tao0 <- c[,1]
      xlabel_tao <-  expression(tau[0])
      from1=min(d_tao0)
      to1=max(d_tao0)
      # this is the prior
      y_tao <- dnorm(seq(from=min(d_tao0-.0001), to=max(d_tao0+.0001), length.out=1000),
                     mean = tao_mean, sd = tao_sd)
    }else{
      d_tao0 <- my_sd 
      xlabel_tao <- expression(sigma)
      from1=0
      to1=quantile(d_tao0,.99)
      # this is the prior
      y_tao <- dgamma(seq(from=min(d_tao0-.0001), to=max(d_tao0+.0001), length.out=5000),
                      shape = sd_shape, scale = sd_scale)
    }
    mintao <- min(d_tao0)

    # Determine the combined y-range for setting ylim
    d_tao0 <- density(d_tao0)
    max1 <- max( as.numeric(d_tao0$y))
    max2 <- max(y_tao)
    y_range <- c(0,max( c(max1,max2) ) )
    # Create a sequence of x values for the normal distribution
    x_seq <- seq(mintao-10, max(d_tao0$x) + 10, length.out = 1000)

    # Plot the normal distribution curve with no y-axis ticks
    if (uq){
    curve(dnorm(x, mean = tao_mean, sd = tao_sd), col='lightgray',
          xlab=xlabel_tao, 
          ylab="Density", main="", yaxt="n",
          # from=mintao-10, to=max(d_tao0$x) + 10, 
          from=from1, to=to1, 
          ylim=y_range,axes=F)

      # Calculate the normal distribution values for the sequence
      y_vals <- dnorm(x_seq, mean = tao_mean, sd = tao_sd)

      legend( 'topright',legend = c(
        paste0('tao_mean: ',tao_mean),
        paste0('tao_sd: ',tao_sd)
      ),
      ,bg=NA, bty="n",cex=.75)
    }else{
      curve(dgamma(x, shape = sd_shape, scale = sd_scale), col='lightgray',
            xlab=xlabel_tao, 
            ylab="Density", main="", yaxt="n",
            # from=mintao-10, to=max(d_tao0$x) + 10, 
            from=from1, to=to1, 
            ylim=y_range,axes=F)

      # Calculate the normal distribution values for the sequence
      y_vals <- dgamma(x_seq, shape = sd_shape, scale = sd_scale)
      legend( 'topright',legend = c(
        paste0('sd_shape: ',sd_shape),
        paste0('sd_scale: ',sd_scale)
      ),
      ,bg=NA, bty="n",cex=.75)
    }

    
    # Use polygon to fill the area under the normal distribution curve
    polygon(c(x_seq, max(d_tao0$x) + 10, mintao-10), c(y_vals, 0, 0), 
            col="lightgray", border=NA)
    
    # Plot the density plot for d_tao0 with blue color
    lines(d_tao0$x[d_tao0$x > mintao-10 & d_tao0$x < max(d_tao0$x) + 10],
          d_tao0$y[d_tao0$x > mintao-10 & d_tao0$x < max(d_tao0$x) + 10],
          col=col.interv)
    
    # Fill the area under the density plot for d_tao0
    polygon(c(d_tao0$x[d_tao0$x > mintao-10 & d_tao0$x < max(d_tao0$x) + 10], max(d_tao0$x) + 10, mintao-10), 
            c(d_tao0$y[d_tao0$x > mintao-10 & d_tao0$x < max(d_tao0$x) + 10], 0, 0), 
            col=col.interv, border=NA)
    
    axis(side=1)
    # axis(side=2)
  }
  
  
  # plot the age-depth models 
  {
    quants <- apply(tau_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
    common_x_lim = c(min(quants[1,]) ,max(quants[3,]) )#+ .1 * b_length  )
    # Set up plot margins: bottom, left, top, and right
    par(mar=c(3.2,2.5,.01,1.5) )
    par(mgp=c(1.5, 0.5, 0))
    
    # if (uq){
    #   plot(quants[2,] , breaks , type='l',col = rgb(0,0,0,.0),
    #        xlim = common_x_lim,
    #        ylim = c(breaks[1]-10,tail(breaks,1)+10),
    #        ylab = xlabel, xlab = "Age",axes=FALSE
    #   )
    #   # This plots the intervals 
    #   lines( quants[1,], breaks, col = col.interv)
    #   lines( quants[2,], breaks, col = col.interv)
    #   lines( quants[3,], breaks, col = col.interv)
    #   polygon(c(quants[1,], rev(quants[3,])), c(breaks, rev(breaks)), 
    #           col = col.interv,border = NA)
    #   
    # }else{
      plot( tau_mat[1,],breaks, type = "l", col = rgb(0,0,0,.0),# this plots a single iteration
           xlim = common_x_lim,
           ylim = c(breaks[1]-10,tail(breaks,1)+10),
           ylab = xlabel, xlab = "Age",axes=FALSE)

            # This plots the intervals 
      lines( quants[1,], breaks, col = col.interv)
      lines( quants[2,], breaks, col = col.interv,lwd = 1.3)
      lines( quants[3,], breaks, col = col.interv)
      polygon(c(quants[1,], rev(quants[3,])), c(breaks, rev(breaks)), 
              col = col.interv,border = NA)
      
      # This will plot the original  model if there is any.
      # In the case of age to age it plots the y=x line. 
      if(depth_to_age & plot_orig){
        lines( org_time,inp$X , col = col.ini_model , lwd = 1.3)
      }else{
        abline(0,1, col = col.ini_model)
      }
    
    # this adds the legends
      {
    leg_loc = "topleft"
    if (!any(is.na(org_time))) {
        if(double_target){
          legend(leg_loc, lty=rep(1,1,1,1,1), 
                 c("Posterior Target", 
                   ori_col_names[2],#'target 1', 
                   ori_col_names[3],#'target 2',
                   Input,
                   paste(Input , ' on original time scale') ), 
                 col=c(col=col.tar,
                       col.tar1,
                       col.tar2,
                       col.interv,
                       col.ini_model), 
                 cex=0.8,lwd=1.3,bg=NA,bty='n')
        }else{
          legend(leg_loc, lty=rep(1,1,1), 
                 c(Target, Input,paste(Input , ' on original time scale') ), 
                 col=c(col.tar,col.interv,col.ini_model), cex=0.8, lwd=1.3,bg=NA, bty="n")
        }
        
      }else{
        # legenda de age-depth without original chronology
        legend(leg_loc, lty=rep(1,2), 
               c(Target , Input), 
               col=c(col.tar,col.interv), bg=NA, bty="n",lwd=1.3)
      }
      
      }
    grid(nx = NULL, ny = NULL, col = grid_col , lty = "dotted") 
    axis(side=1)
    axis(side=2)
  }
  
  ### Plot the alignment 
  {
    # Set up plot margins: bottom, left, top, and right
    par(mar=c(.1,2.5,1.1,1.5) )
    # Calculate quantiles
    quants2 <- apply(tau_mat2, 2, quantile, probs = c(0.025, 0.5, 0.975),na.rm=T)
    mean_tau <- apply(tau_mat2, 2, mean,na.rm=T)
    output$quantiles <- quants2
    
    if (range_f==1){
      y_lim = c(-1.2,1.2)
    }else{
      y_lim = c(-3.2,3.2)
    }
    
    # Plot the target 
    if (uq){
      density_plot(tar_ages ,tar_mt, xlabel='',flip = F , ylim=y_lim,xlim=common_x_lim)
      }else if(double_target){
        
        # Here we calculate the mean mixture
        # w_t = mean(w_tar)
        w_t = median(w_tar)
        # We create the mixed target
        tar$ProxyValue <- range_T(w_t * tar$ProxyValue1 + (1 - w_t) * tar$ProxyValue2)
        # We plot the mix target
        plot(tar$X,tar$ProxyValue, type='l', xlab="", ylab="Rescaled proxy",  xaxt='n',
             xlim = common_x_lim,
             axes=FALSE,lwd = 1.3,
             col=col.tar,cex=.2)

        if(plot.tar_lines){
          lines(tar$X,range_T(tar$ProxyValue1), type='l',lty=1 ,col=col.tar1,cex=.15)
          lines(tar$X,range_T(tar$ProxyValue2), type='l', lty=1 ,col=col.tar2,cex=.15)
        }else{
          points(tar$X,range_T(tar$ProxyValue1)  ,col=col.tar1,cex=.6,pch=16)
          points(tar$X,range_T(tar$ProxyValue2) ,col=col.tar2,cex=.6,pch=16)
        }
        

        
      }else{
        plot(tar$X,tar$ProxyValue, type='l', xlab='', ylab="Rescaled proxy", xaxt='n',
             xlim = common_x_lim,
             axes=FALSE,
             col=col.tar,cex=.15)
      }
     
    
    axis(side=2)
    # Plot the aligment using same color as the age-model
    lines(mean_tau,inp$ProxyValue, col = col.interv, lty=1, type = 'l', cex=.2,lwd=1.3)
    quantiles_975 <- apply(tau_mat2, 2, quantile, probs = 0.975)
    quantiles_025 <- apply(tau_mat2, 2, quantile, probs = 0.025)
    
    # Plot lines for quantiles
    segments(y0=inp$ProxyValue, x0=quantiles_975, 
            y1=inp$ProxyValue, x1=quantiles_025,
            col=col.interv,cex=.15)
    
    if (!any(is.na(org_time))) {

      # Aqui agregar un if para cuando se hace age-depth model y cuando se hace age to age
      lines(org_time, inp$ProxyValue ,col = col.ini_model,lty=1,type='l',cex=.2)
    }

    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted") 
    usr <- par("usr")
    
    if(!double_target){
      text(x = 0.95 * usr[2] , y = 0.9 * usr[4] , 
           labels = paste0("Aligment between ",Input,' and ',Target), pos = 2,
           cex = 1.2, col = "darkblue", font = 2)
    }else{
      text(x = 0.95 * usr[2] , y = 0.9 * usr[4] , 
           labels = paste0("Aligment between ",Input,' and a mixture of ',
                           ori_col_names[2],' and ', 
                           ori_col_names[3]), pos = 2,
           cex = 1.2, col = "darkblue", font = 2)  
    }

    
    
  }
  
  
  # plot the mixture
  if(double_target){
    par(mar=c(2.3,1.1,1.1,1.1) )
    # layout(matrix(c(1)))
    # if (savefig){  pdf(paste0(folder,'/aligment_',Input,'-',Target,'_proportion.pdf')) }
    kern <- density(output$w_tar)
    plot(kern,main = '',xlab = 'mixture',ylab ='', xlim=c(0,1.),
         axes=F,col=rgb(0,0,0,0),yaxt="n")

    # plot the prior
    polygon(c(0,1,1,0),c(0,0,.3*max(kern$y),.3*max(kern$y)), col="lightgray", border=NA) 
    
    # Fill area under the density plot for beta
    with(kern, {
      relevant_indices <- x > 0 & x < 1
      polygon(c(x[relevant_indices], rev(x[relevant_indices])),
              c(rep(0, sum(relevant_indices)), rev(y[relevant_indices])), 
              col = col.interv, border=NA)
    })
    
    points(mean(output$w_tar),0,col= col.tar ,pch=18)
    # points(median(output$w_tar),0,col='blue',pch=19)
    # points(kern$x[which(kern$y == max(kern$y)) ][1],0, col='blue',pch=4)
    text(x = -.005, y = 0, labels = ori_col_names[3], pos = 4,srt = 90, adj = .5,cex=.8) # pos = 1 for below
    text(x = .91, y = 0, labels = ori_col_names[2], pos = 4,srt = 90, adj = .5,cex=.8 ) # pos = 1 for below
    
    legend( 'topright',legend = paste0('mean: ',round(mean(output$w_tar),2)),cex=.75,
            col= col.tar ,pch = c(18),bg=NA, bty="n")


    axis(side=1)
  }
  
  if(savefig){dev.off()}

  #### Bye Bye message #### 
  cat('Might sound crazy but it ain\'t no lie\n
          Bye bye bye\n')  
   # if (run){
    return (output)
   # }
  
}



#####
#' @name IAT
#' @title calculate the Integrated Autocorrelation Time
#' @description Calculate the Tntegrated Autocorrelation Time, which gives the proposed value for thinning. E.g., if the IAT is 80, it is good to thin the MCMC run by storing only every 80 iterations. This method is slower than GetAutoCov, but much better.
#' @param set This option reads the 'info' variable, which contains the data and the model output.
#' @param par The parameter to test. Defaults to 0.
#' @param from The first of the iterations of the MCMC run to check. Defaults to the first one.
#' @param to The last of the iterations of the MCMC run to check. Defaults to the last one.
#' @return The IAT value
#' @author Andres Christen
#' @export
IAT <- function(set, par=0, from=1, to) {
  ## -lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to))
  ## we get the desired time series, the parameter and the from - to selection
  if(par>0) {
    if(set$dim > 1)
      dt <- set$output[from:to, par] else
        dt <- set$output[from:to]
  } else
    dt <-  set$Us[from:to]
  
  n <- to-from
  mu <- mean(dt)  ### with its mean and variance
  s2 <- var(dt)
  
  ### The maximum lag is half the sample size
  maxlag <- max( 3, floor(n/2))
  
  #### The gammas are sums of two consecutive autocovariances
  Ga <- rep(0,2)  ## two consecutive gammas
  
  lg <- 0
  Ga[1] <- s2  #sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
  lg <- 1
  Ga[1] <- Ga[1] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
  
  m <- 1
  lg <- 2*m
  Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg)
  lg <- 2*m+1
  Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/(n-lg)
  
  IAT <- Ga[1]/s2  ### Add the autocorrelations
  
  ### RULE: while Gamma stays positive and decreasing
  while((Ga[2] > 0.0) && (Ga[2] < Ga[1])) {
    m <- m+1
    if(2*m+1 > maxlag) {
      message("Not enough data, maxlag=", maxlag, "\n")
      break
    }
    Ga[1] <- Ga[2]
    
    lg <- 2*m
    Ga[2] <- sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
    lg <- 2*m+1
    Ga[2] <- Ga[2] + sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
    
    IAT <- IAT + Ga[1]/s2
  }
  
  IAT <- -1 + 2*IAT   ##Calculates the IAT from the gammas
  return(IAT)
}
