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


BSynch <- function(Input,Target,folder = '~/Documents/BSynch/',
                   shape_acc = 2,  mean_acc = 20,
                   strength_mem = 10, mean_mem = .5,
                   n_sections = 100, 
                   sd_input = .1 ,sd_target = .1,
                   tao_mean=TRUE , tao_sd=TRUE ,
                   gw_z=.1, depth_cm = TRUE, age_kyr=TRUE,
                   ta = 3,tb = 4,
                   uq =TRUE, depth_to_age = TRUE,
                   iters = 2.5e+3, burn = 2e+3 ,thin = 150,
                   continue_run = TRUE, 
                   verify_orientation = TRUE, flip_orientation = FALSE,
                   cc_limit = FALSE
                   ){ 
  #### Load packages ####
  load_or_install_dependencies()
  
  #### Load data and twalk ####
  source_twalk(folder)
  
  #### load cc and data to check if required ####
  if (cc_limit){
    cc <- read.table(paste0(folder,'6Col_hulu_updated.14C.txt'),header = TRUE)
    test_data <- read.table(paste0(folder, 'Cariaco_14C.txt'),header=FALSE )
    # Create interpolation function 
    # up_lim <- rep(NA,nrow(test_data))
    # 
    # for (i in 1:nrow(test_data)){
    #   # Create interpolation function 
    #   f <- approxfun(cc[,1]  , cc[,2] + 3*cc[,3]  - (test_data[i,2] - 3 * test_data[i,3])  )
    #   # Find roots of interpolation function
    #   ints <- uniroot(f, lower=min(cc[,1]), upper=max(cc[,1]))$root
    #   # Find max intersection
    #   if(length(ints) > 0){
    #     max_x <- max(ints)
    #   } else {
    #     max_x <- NA 
    #   }
    #   # save intersections to up_lim 
    #   up_lim[i] <- max_x
    # }
    
    up_lim <- test_data[,2] - 3 * test_data[,3]
    depth_to_check <- test_data[,1] * 100
    
    # tmp_slope <- lm(formula = up_lim ~ depth_to_check )
    # tmp_slope <- as.numeric(tmp_slope$coefficients[2])
  
  }

  #### Load input ####
  inp <- load_file_from_folder(Input, folder) 
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
    
  }else{
    org_time <- inp$Age
    inp <- data.frame(X=inp$Age,ProxyValue = inp$ProxyValue)
    xlabel = "Age"
  }

  
  if (age_kyr){
    age_temp <- inp$X * 1000
  }else{
    age_temp <- inp$X   
  }
  
  # quick test (remove)
  # ordered <- sort(sample(1:dim(inp)[1],as.integer(.5*dim(inp)[1])))
  # 
  # inp <- inp[ordered,]
  # org_time <- org_time[ordered]
  # 
  # print(inp)
  
  #### Load target####
  # Note: create the object inp which only containes the variable which will be alinge and the scale either depth or age
  if (!uq){ # Load target when only age estimates are provided
    tar <- load_file_from_folder(Target, folder) 
    tar <- data.frame(X=inp$Age,ProxyValue = inp$ProxyValue)
  }else{# Load target when iterations of an age-depth model are provided
    tar <- read.table(paste0(folder,Target,'.txt'),header = T)
    tar <- as.matrix(tar)
    # Divide the data base into the ages and the iterations
    tar_ages <- tar[,1]
    tar <- t(tar[,-1]) # invert the matrix
    tar_mt <- tar
    # create the density function.
    kde_list <- target_density(tar_ages,tar) # target

    # localizador is the function which tell the target_densities function which density to use
    localizador <- function(x)(approx(tar_ages,seq(1,length(tar_ages)),method = 'constant',xout = x)$y)
    # Create a data frame with depth and mean
    tar <- data.frame(X = as.numeric(tar_ages), ProxyValue =as.numeric(colMeans(tar)) )
  }  
  #### Check that the records are properly aligned ####
  if (flip_orientation){inp$ProxyValue <- -inp$ProxyValue}
  
  if (verify_orientation){
    inv_user = TRUE
    while(inv_user){
      par(mfrow=c(1,1))
      if (uq){
        density_plot(tar_ages ,tar_mt, xlabel=xlabel )
        axis(1, at=NULL, labels=FALSE, tick=FALSE)
      }else{
        plot(tar, type='l', xlab=xlabel, ylab="proxy units", xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))  
      }
      par(new=TRUE)
      plot(inp$X,inp$ProxyValue, type='l', col=rgb(0,0,1,.9), xlab="", ylab="", xaxt='n', yaxt='n', xlim=c( inp$X[1],max(inp$X)))
      legend("bottomright", lty=rep(1,2),
             c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
      plot(inp$X,inp$ProxyValue, type='l', col=rgb(0,0,1,.9), xlab="", ylab="", yaxt='n', xlim=c( inp$X[1],max(inp$X)))
      legend("bottomright", lty=rep(1,2),
             c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
      par(new=TRUE)
      if (uq){
        density_plot(tar_ages ,tar_mt, xlabel=xlabel ,axis=FALSE)
      }else{
        plot(tar, type='l', xlab=xlabel, ylab="proxy units", xaxt='n',xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))  
      }
      
      
      message('This are the two records side to side\n')
      y <- readline("Shall I invert the input for a proper aligment? yes or not \n")
      # Check response
      if (y %in% c("y", "yes")) {
        inp$ProxyValue <- -inp$ProxyValue
      } else if (y %in% c("n", "no")) {
        # Default no inversion
      } else {
        print("Invalid input")
      }
      
      # Always re-plot with new inversion
      par(mfrow=c(1,1))
      plot(inp$X,inp$ProxyValue, type='l', col=rgb(0,0,1,.9), xlab="", ylab="", yaxt='n', xlim=c( inp$X[1],max(inp$X)))
      legend("bottomright", lty=rep(1,2),
             c("input", "target"), col=c(rgb(0,0,1,.9),1), cex=0.5)
      par(new=TRUE)
      if (uq){
        density_plot(tar_ages ,tar_mt, xlabel=xlabel,axis=FALSE )
      }else{
        plot(tar, type='l', xlab=xlabel, ylab="proxy units", xaxt='n',xlim=c( tar$X[1],max(tar$X)), col=rgb(0,0,0,1))  
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

  ##### Checks that target is longer than the input ####
  # This in only necessary on the age to get alignment 
  if (!depth_to_age){
    # checks that target is longer than the input
    # Extract the scale columns
    inp_scale <- inp[,1] 
    tar_scale <- tar[,1]
    
    # Check if inp scale range is within tar scale range
    if(min(inp_scale) < min(tar_scale) | 
       max(inp_scale) > max(tar_scale)) {
      
      # Trim inp to be within tar scale range
      inp <- inp[inp$X <= tail(tar[,1],1)-10, ] 
      inp <- inp[inp$X >= head(tar[,1],1)-10, ]
      
      message("input dataset trimmed to be within target scale range")
    } else {
      message("input dataset already contained within target scale range") 
    }
  }

  #### Set variables ####
  # normalize the records 
  inp[,2] <- range(inp[,2]) 
  tar[,2] <- range(tar[,2])
  
  # Define the variables which are use for the alignment 
  
  breaks <- seq(inp$X[1],tail(inp$X,1), length.out = n_sections)
  b_length <- breaks[2] - breaks[1]
  
  # tau0 variables
  if (tao_mean == TRUE){
    if (depth_to_age){
      tao_mean <- min(tar_ages) 
      # uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
      tao_sd <-  ifelse(tao_mean  < 10 * 1000, 500, 5000)  
      tao_mean <- min(tar_ages) + 1 * tao_sd
    }else{
      tao_mean <- min(inp$X) 
      # uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
      tao_sd <-  ifelse(head(inp$X,1) < 10 * 1000, 500, 5000)   
    }
    
    
  }else{ # Check if tao_mean  and tao_sd is numeric
    if(!is.numeric(tao_mean)) {
      tao_mean <- as.numeric(readline("Enter a numeric value for tao_mean: ")) 
    }
    if(!is.numeric(tau_sd)) {
      tau_sd <- as.numeric(readline("Enter a numeric value for tau_sd: "))
    }
  }
  
  # parameter for prior of alphas 
  scale_acc = mean_acc / shape_acc
  
  # parameters for memory 
  m_alpha = strength_mem * mean_mem 
  m_beta =  strength_mem * (1 - mean_mem) 
  
  # set lower and upper prior boundaries
  # note that when age to age alignment alphas will come from a lognormal (to do)
  if (depth_to_age){
    low <- c(tar_ages[1],#tar$X[1] - (3*tao_sd), # tao0
             0, # mem
             rep(0, n_sections-1)) # m_s
    up <- c( tar_ages[1] + (2*tao_sd),  #tao0
             1, # mem,
             rep(50, n_sections-1)) # alphas
             # rep(Inf, n_sections-1)) # alphas
  }else{
    low <- c(tar_ages[1], #tar$X[1] - (3*tao_sd), # tao0
             0, # mem
             rep(0.25, n_sections-1)) # arate
    up <- c( tar_ages[1] + (2*tao_sd),  #tao0
             1, # mem,
             rep(4, n_sections-1)) # alphas 
  }
  
  mx_age <- max(tar_ages)
  #print(mx_age)
  # new sd (geometric mean of sd of target and input)
  my_sd = sqrt(sd_input^2 +sd_target^2)
  
  # get bandwidths
  if (uq){
    # Extract bandwidths using lapply
    bandwidths <- lapply(kde_list, function(x) x$bw)
    
    # result as a vector
    bandwidths <- unlist(bandwidths)
  }
  
  
  #### Functions ####
  
  # tau function (the function which gives the value in the target scale)
  tau <- function(x, param){
    tau0 <- param[1]
    # w <- param[2]
    slopes <- param[-c(1,2)]
    tau_is <- c(tau0, tau0 + cumsum(slopes * b_length))
    return( approx(breaks,tau_is,x)$y )
  }
  
  # Provides the alphas for calculation its priors
  n_slopes = n_sections -1 
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
    # To0 prior
    l_prior <- dnorm(tao0, mean = tao_mean, sd = tao_sd,log = TRUE)
    # memory value
    l_prior <- l_prior + dbeta(w, m_alpha, m_beta, log = TRUE)
    # denisity for acc rates 
    l_prior <- l_prior + sum( dgamma(alphas(param), shape = shape_acc, scale = scale_acc, log = TRUE) )
    return(l_prior)
  }
  
  n_dat <- length(inp$ProxyValue)
  # loglikelihood
  loglikelihood <- function(params){
    # Get the ages at which to calculate the target 
    t_times <- tau(inp$X,params)
    # get the values of the target at the required times. 
    new_target <- approx(tar$X,tar$ProxyValue,t_times)$y
    # predicted values should be close to observed values. Use t distribution
    likelihood <- tdistro(X = inp$ProxyValue, Mu = new_target, sigma = my_sd, a = ta, b = tb) # change 0.1 to the reported sd
    
    # this is a experiment  (remove or comment)
    # subsample <- sample(1:n_dat,as.integer(n_dat*.25))
    # likelihood <- tdistro(X = inp$ProxyValue[subsample], Mu = new_target[subsample], sigma = my_sd, a = ta, b = tb) # change 0.1 to the reported sd
    
    # sum the likelihood
    sumll = sum(likelihood)
    return(sumll)
  }
  
  # likelihood in case of UQ
  if (uq){
    Rcpp::sourceCpp("~/GitHub/BSynch/targetDensity.cpp")
    loglikelihood_uqC <- function(params){
      # Get the ages at which to calculate the target 
      t_times <- tau(inp$X,params)
      # get the values of the target at the require times. 
      loca <- localizador(t_times)
      loca <- as.integer(loca)
      
      ll <-loglikelihood_uqCpp(params, inp$ProxyValue, bandwidths, tar_mt, sd_convertor ,loca)
      
      return(ll)
    }
    
    n_kde <- length(tar_mt[,1])
    sd_convertor <- 1/(1.06*n_kde^(-1/5))
    
    # Define a function to get the log density for a given depth and value y
    l_target_kernel <- function(d, new_x,bws) {
      # d is the location in the age vector
      ll_vector <- t_dis(X = new_x ,  Mu =  t( tar_mt[,d] ),  sigma =  sd_convertor * bws[d], 
                         a = ta, b = tb) 
      
      uq_l <-   sum( log( colSums(t(ll_vector))  ))
      # return(sum( dnorm(new_x , tar_mt[,d],sd_convertor * kde$bw)) / n_kde )
      return (uq_l)
    }
    
    loglikelihood_uq <- function(params){
      # Get the ages at which to calculate the target 
      t_times <- tau(inp$X,params)
      # get the values of the target at the require times. 
      loca <- localizador(t_times)
      ll <- l_target_kernel(d = loca,new_x = inp$ProxyValue,bws = bandwidths) 
      return(ll)
    }
  }
  
  

  # objective distritro
  if (uq){
    if(cc_limit){
      invcal <- approxfun(cc[,1],cc[,2])
      
      obj <- function(param){ 
        cal_yr <- tau(depth_to_check,param)
        cc_yr <- invcal(cal_yr)
        newlogll <- sum( dnorm((cc_yr - test_data[,2] ) , 30,50,log = T ) )
        return( -(logprior(param) + loglikelihood_uqC(param) + newlogll))
        }
      
    }else{
      obj <- function(param){ -(logprior(param)+loglikelihood_uqC(param)) }
    }
    
  }else{
    obj <- function(param)(-(logprior(param)+loglikelihood(param)))
  }
  
  
  # support function
  if (cc_limit){
    invcal_supp <- approx(cc[,1]  , cc[,2] + 3*cc[,3])
    supp <- function(params){
      alp <- alphas(params)
      age_lim <- tauC(tail(inp$X,1),params) 
      
      # ages_to_check <- tau(depth_to_check,params)
      # ages_to_check <-invcal_supp( ages_to_check )
      # 
      # print(ages_to_check > up_lim  )
      # print('______')
      # print(all(alp > 0 ))
      # tau0 <- param[1]
      # w <- param[2]

      ifelse(all( c(params > low, params < up,alp > 0,#ages_to_check > up_lim ,
                    age_lim < mx_age ) ) ,
             return(TRUE), 
             return(FALSE)
      )
    }
    # generates initial points.
    sampler <-  function(rerun=FALSE){
      # generate to0
      tao0 <- runif(1, low[1],low[1]+1) 
      # memory value
      indx <- which(breaks >= depth_to_check[1])[1]-1
      indx2 <-  tail(which(breaks <= tail(depth_to_check,1)),1)+2
      # get the first slopes
      
      f <- approxfun(cc[,1]  , cc[,2] + 3 * cc[,3] - max(up_lim[1:5] )  )
      # Find roots of interpolation function
      xx_lim <- uniroot(f, lower=min(cc[,1]), upper=max(cc[,1]))$root
       
      m0 <-  ( max(xx_lim + 1000) - tao0 ) / ( depth_to_check[1] - inp$X[1] )
      ms <- rnorm(indx, m0, .01 )  
      #get the slopes for the place where there is cc
      xx_lim <- approx( cc[,2] + 3 * cc[,3] ,cc[,1]  ,up_lim)$y
      slopes <- lm(formula = (xx_lim +500 ) ~ depth_to_check  )
      slopes <- ceiling(as.numeric(slopes$coefficients[2]))
      ns <- indx2 - indx
      slopes <- rnorm(ns,slopes,.01)
      ms <- c(ms,slopes)

      #finisg the ms with the remaning 
      tau_is <- c(tao0, tao0 + cumsum(ms * b_length))
      last_ms = ((mx_age-5) - tail(tau_is,1))/(tail(breaks,1)-breaks[indx2])
      ns <- length(breaks) - indx2 - 1
      last_mss <- rnorm(ns,last_ms,.01)
      ms <- c(ms,last_mss)
      w_lim <- min (ms[-length(ms)]/ms[-1])

      w <-  runif(n = 1, min = 0, max = w_lim)
      ini_points <- c(tao0, w, ms)

      # print(length(breaks))
      # print(length(ms))
      # print(supp(c(tao0, w, ms)) )
      return(ini_points)
    }
    
    
    
  }else{
      supp <- function(params){
        # tau0 <- param[1]
        # w <- param[2]
        alp <- alphas(params)
        age_lim <- tau(tail(inp$X,1),params)
        ifelse(all( c(params > low, params < up) ) & 
                 all( alp > 0 ) &   
                 age_lim < mx_age ,
               return(TRUE), 
               return(FALSE)
        )
      }
      # generates initial points.
      sampler <-  function(rerun=FALSE){
        # generate to0
        tao0 <- runif(1, low[1],low[1]+10) 
        # memory value
        w <-  runif(n = 1, min = .6, max = 1)
        #  accumulation rates 
        # m_{k-1} = \omega m_k + (1-\omega) alpha_{k-1}
        # alpha_{k-1} \sim gamma(a,b)
        if(rerun){
          tmp_mean <- floor( (tail(tar$X,1)- tar$X[1])/(tail(inp$X,1)-inp$X[1]) )
          if (cc_limit){
            tao0 <- rnorm(1,org_time[1]+500,1)
            print(breaks)
            print(test_data)

            print( diff(approx(inp$X,org_time,breaks)$x ) )
            tmp_mean <- diff( approx(inp$X,org_time*1000,breaks)$y ) / diff( approx(inp$X,org_time*1000,breaks)$x ) + 5
            # print(tmp_mean)
          }
          ms <- rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )  
          for (i in 1:(n_slopes-1)){
            alpha = rgamma(1, scale = tmp_mean/ shape_acc, shape = shape_acc )  
            ms <- c(ms , w * ms[1] + (1-w) * alpha)
          }
        }else{
          ms <- rgamma(1, scale = scale_acc , shape = shape_acc )
          for (i in 1:(n_slopes-1)){
            alpha = rgamma(1,scale = scale_acc, shape = shape_acc )
            ms <- c(ms, w * ms[1] + (1-w) * alpha )
          }
        }
        
        return(c(tao0, w, ms))
      }
  }
  

  
  #### Run the twalk ####

  initial_search <- function(){
    x1 <- sampler()
    while(!(supp(x1))){
      x1 <- sampler(rerun = TRUE)  
    }
    return(x1)
  }
  

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
    x1 <- initial_search()
    x2 <- initial_search()
  }
  
  message("Initiating the t-walk process...")
  # burn = burn * length(x1)
  
  output <- Runtwalk(iters,Obj = obj,dim = length(x1),x0 = x1,xp0 = x2,Supp =supp,
                     thinning = length(x1)*thin, burnin= length(x1)*burn )
                     # thinning = thin,burnin= burn )
  
  # twalk state
  lastx = tail(output$output, 1)
  lastxp = tail(output$outputp, 1)
  # Write the matrix to a CSV file
  write.csv(matrix_to_save <- matrix(c(lastx, lastxp), nrow = 2,byrow = TRUE), 
            paste0(folder,"twalkstate.csv"), row.names=FALSE)#, col.names=FALSE)

  # Check convergance
  # plot(output$Us,type='l')
  # test <- geweke.diag(as.mcmc(as.numeric(output$Us)),frac1=0.5, frac2=0.5)
  # cat("\n================== Geweke DIAGNOSTIC ==================\n")
  # cat('Geweke value: ',round(test$z,3),"\n")
  # cat('It appers that the chain has converge,',"\n") 
  # cat('please refer to the energy plot to validate the convergance\n')
  # cat("====================================================\n\n")  
  # 
  
  c <- output$output[-1,]
  energy <- output$Us[-1]
  energy2 <- output$Ups[-1]
  write.csv(c,paste0(folder,Input,'-',Target,'_',n_sections,".out"), row.names=FALSE)#, col.names=FALSE)

  iat <- IAT(output,to=output$Tr)
  
  # 
  cat("\n================== IAT DIAGNOSTIC ==================\n")
  cat("IAT Value:", iat, "\n")
  if (iat < 5) {
    cat("Interpretation: The chain exhibits low correlation among successive samples.\n")
    cat("Recommendation: Current settings appear satisfactory.\n")
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
    tau_row <- tau(breaks, c[i,])
    # Store result in matrix
    tau_mat[i,] <- tau_row
  }
  
  
  # Initialize matrix to store assembles
  tau_mat2 <- matrix(NA, nrow = nrow(c), ncol = length(inp$X))
  
  # Loop through each row of c
  for(i in 1:nrow(c)) {
    # Apply tau to breaks with this row of c as params
    tau_row2 <- tau(inp$X, c[i,])
    # Store result in matrix
    tau_mat2[i,] <- tau_row2
  }
  
  

  #### Plot the results ####
  message("Commencing the plotting process...")
  # pdf(paste0(folder,'/aligment.pdf'))
layout(matrix(c(1,1,1,1,5,5,5,5,5,5,
                1,1,1,1,5,5,5,5,5,5,
                2,2,2,2,5,5,5,5,5,5,
                2,2,2,2,5,5,5,5,5,5,
                3,3,3,3,6,6,6,6,6,6,
                3,3,3,3,6,6,6,6,6,6,
                4,4,4,4,6,6,6,6,6,6,
                4,4,4,4,6,6,6,6,6,6), 8, 10, byrow = TRUE))
  # plot energy
  plot(energy,type = 'l',col="grey40")
  lines(energy2,col="grey60")
  

  ## Plot posterior and prior of accumulations
  { 
    d_m <- density(output$output[,-c(1,2)])
    # Get y-values of the gamma curve across the specified x-range
    y_gamma <- dgamma(seq(from=0, to=mean_acc +  2*(shape_acc*scale_acc), length.out=1000),
                      shape = shape_acc, scale = scale_acc)
    
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_gamma),max(d_m$y))) )
    
    # Plot gamma curve with adjusted y-axis limits
    curve(dgamma(x, shape = shape_acc, scale = scale_acc), 
          from=0, to=mean_acc +  3*(shape_acc*scale_acc),  yaxt="n",
          xlab="Accumulation Rates (yr/cm)", ylab="Density", 
          main=" ", ylim=y_range)
    
    # Add density plot
    lines(d_m, col='blue')
  }
  
  ## plot memory prior vs posteriors
  {
    d_w <- density(output$output[,2])
    
    y_beta <- dbeta(seq(from=0, to=1, length.out=1000),
                    m_alpha, m_beta)
    # Determine the combined y-range for setting ylim
    y_range <- c(0,max(c(max(y_beta),max(d_w$y))) )
    
    curve(dbeta(x, m_alpha, m_beta), 
          from=0, to=1, ylim=y_range, yaxt="n",
          xlab="memory", ylab="Density", main="")
    lines(d_w,col='blue')
  }
  
  
  ## plot initial value
  {
    d_tao0 <- output$output[,1]
    mintao <- min(d_tao0)
    
    y_tao <- dnorm(seq(from=min(d_tao0-.0001), to=max(d_tao0+.0001), length.out=1000),
                   mean = tao_mean, sd = tao_sd)
    # Determine the combined y-range for setting ylim
    d_tao0 <- density(output$output[,1])
    max1 <- max( as.numeric(d_tao0$y))
    max2 <- max(y_tao)
    y_range <- c(0,max( c(max1,max2) ) )
    curve(dnorm(x, mean = tao_mean, sd = tao_sd), 
          xlab="Tao0 ", ylab="Density", main="",yaxt="n",
          from=mintao-10, to=max(d_tao0$x+10), ylim=y_range)

    lines(d_tao0,col='blue')

  }
  

  # plot the age-depth models 
  {
    if (uq){
      density_plot(breaks,tar_mt = tau_mat ,
                   xlabel = 'Depth',ylabel ='Age',flip = T) 
        # Get column quarantines
      quants <- apply(tau_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
      lines( quants[1,], breaks ,col = "red")
      lines( quants[2,], breaks ,col = "red")
      lines( quants[3,], breaks ,col = "red")
      lines( org_time*1000, inp$X, col=rgb(0,1,1,.5))
      
    }else{
      plot(breaks, tau_mat[1,], type = "l", col = rgb(0,0,0,.01),
           ylim=c(min(tau_mat[,1]),max(tau_mat)),
           xlab = xlabel, ylab = "Age")
    # Get column quarantines
    quants <- apply(tau_mat, 2, quantile, probs = c(0.025, 0.5, 0.975))
    lines( breaks, quants[1,], col = "red")
    lines( breaks, quants[2,], col = "red")
    lines( breaks, quants[3,], col = "red")
    lines( rev(org_time)*1000, inp$X, col=rgb(0,1,1,.5))

    }
    
  }

  # Plot the alignment 
  {
    quants2 <- apply(tau_mat2, 2, quantile, probs = c(0.25, 0.5, 0.75))
    if (uq){
      density_plot(tar_ages ,tar_mt, xlabel='Age',flip = F )
    }else{
      plot(tar, type='l', xlab=xlabel, ylab="proxy units", xlim=c( tar$X[1],max(tar$X)),col=rgb(1,0,1,.7))
    }  
    
    lines(quants2[2,],inp$ProxyValue, col=rgb(0,0,1,.5))
    if (!is.na(org_time[1])){
      lines(org_time*1000,inp$ProxyValue ,col=rgb(1,0,1,.4))
      legend("bottomright", lty=rep(1,2,1), c("target", "input",'input orig time scale'), col=c(1,rgb(0,0,1,.5),rgb(1,0,1,.4)), cex=0.5)
    }else{
      legend("bottomright", lty=rep(1,2), c("input", "target"), col=c(1,4), cex=0.5)
    }
  
  }
  
  if (cc_limit){
    layout(matrix(c(1), 1, 1, byrow = TRUE))
    {
      # test_data[,2] - 3 * test_data[,3]
      cal_yr <- approx(breaks, quants[2,],depth_to_check)$y
      plot(cal_yr,test_data[,2],
           xlim = c(min(cal_yr)-100, max(cal_yr)+100),
           ylim = c(min(test_data[,2])-100,max(test_data[,2])+100),
           xlab = 'cal yr', ylab = 'Radiocarbon age',
           pch=16,col=rgb(1,0,0,.6))
      lines(cc[,1],cc[,2],col='black')
      lines(cc[,1],cc[,2] - 3 * cc[,3],col='black')
      lines(cc[,1],cc[,2] + 3 * cc[,3],col='black')
      points(approx(inp$X,org_time*1000,depth_to_check)$y,test_data[,2],col = rgb(0,0,1,.3),pch = 16)
      legend("topleft",legend = c('Cariaco','Calibration Curve'), pch = c(16,NA),
             lty = c(NA,1),bg=NA,col = c(rgb(1,0,0,.5),1),bty = 'n' )
      
    }    
  }


  #### finishing message ####
  cat('Might sound crazy but it ain\'t no lie\n
          Bye bye bye\n')  
  return (output)
  }




# Scales the data in the range of [-1,1]
range <- function(x){
  2*(x-min(x))/(max(x)-min(x))-1
  }

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
  message(paste0("Getting data from '", file_name, "'...\n"))
  # Check for .csv extension
  csv_path <- paste0(folder, file_name, ".csv")
  if (file.exists(csv_path)) {
    
    fil <- read.csv(csv_path)
    
    # Validate columns and headers
    if(ncol(fil) != 3) {
      stop("CSV file does not have 3 columns")
    } 
    if(!identical(names(fil), c("Depth", "Age", "ProxyValue"))) {
      names(fil) <- c("Depth", "Age", "ProxyValue") 
      warning("CSV file headers modified to match required format")
    }
    
    return(fil)
    
  }else{
    # Check for .txt extension
    txt_path <- paste0(folder, file_name, ".txt")
    if (file.exists(txt_path)) {
      
      fil <- read.table(txt_path, header = TRUE, sep = "\t")
      
      # Validate columns and headers
      if(ncol(fil) != 3) {
        stop("TXT file does not have 3 columns")
      }
      if(!identical(names(fil), c("Depth", "Age", "ProxyValue"))) {
        names(fil) <- c("Depth", "Age", "ProxyValue")
        warning("TXT file headers modified to match required format")
      }
    }
    return(fil)
  }
  
  # Return NULL with warning if no matching file
  warning("No matching file found for the given input in the specified folder.")
  return(NULL)
  
}

# Load twalk
source_twalk <- function(folder) {
  # Construct the path to the twalk.R file in the specified folder
  twalk_path <- paste0(folder, "..", "/twalk.R")
  # Check if the twalk.R file exists in the specified folder
  if (file.exists(twalk_path)) {
    source(twalk_path)
    message("Successfully loaded 'twalk.R' from", folder, "directory.\n")
  } else {
    warning("File twalk.R was not found in the specified folder.")
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
target_density <-function(tar_ages,tar,bw = 0.05){

  # Initialize an empty list to store the density objects for each depth
  kde_list <- list()
  
  # Loop through each column (depth) to calculate the kernel density
  for (col in 1:length(tar_ages)) {
    values_at_depth <- tar[,col]
    # Compute kernel density and add to the list
    kde <- density(values_at_depth, kernel = 'gaussian',bw =bw )
    kde$x <- c(-.Machine$double.xmax,kde$x,.Machine$double.xmax) 
    kde$y <- c(.Machine$double.xmin,kde$y+1e-200,.Machine$double.xmin  )
    kde$y <- log(kde$y)
    kde_list[[col]] <- kde
    # print(kde$bw)
  }

  return(kde_list)
}

# Plot the density
density_plot <- function(tar_ages, tar_mt,xlabel,ylabel = "proxy units",add = FALSE,axis=TRUE,flip=FALSE){
  tar_ages_bw <- tar_ages[2] - tar_ages[1] 
  if (!add){
    if(axis){
      if (flip){
        plot(colMeans(tar_mt),tar_ages , type='l', xlab=ylabel, ylab=xlabel, ylim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1))  
      }else{
        plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1))  
      }
       
    }else{
      if (flip){
        plot(colMeans(tar_mt),tar_ages , type='l', ylab=xlabel, xlab=ylabel, ylim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1),xaxt = 'n')   
      }else{
        plot(tar_ages,colMeans(tar_mt) , type='l', xlab=xlabel, ylab=ylabel, xlim=c( tar_ages[1],tail(tar_ages,1)),col=rgb(0,0,0,1),xaxt = 'n')    
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
      rect(ybottom = tar_ages[i]-.5*tar_ages_bw, 
           ytop = tar_ages[i]+.5*tar_ages_bw,
           xleft = h$breaks[-1], #head(h$breaks, -1),
           xright = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
    }else{
      rect(xleft = tar_ages[i], 
           xright = tar_ages[i]+tar_ages_bw,
           ybottom = h$breaks[-1], #head(h$breaks, -1),
           ytop = head(h$breaks, -1),# h$breaks[-1],
           col = cols, border = NA)
    }
  }
}




##### Run the model for testing (using Uq)

output = BSynch(Input='input',Target='target',folder = '~/Documents/BSynch/Uq_cc/',
                n_sections = 40,
                thin=10,burn=3e+3,iters=1.5e+3,
                shape_acc = 100,  mean_acc = 20,
                strength_mem = 10, mean_mem = .5,
                depth_cm = FALSE,depth_to_age = TRUE,
                age_kyr = TRUE,continue_run = T,
                cc_limit = TRUE,
                verify_orientation = FALSE,flip_orientation = TRUE  )




