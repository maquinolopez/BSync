# BSynch 
#
# Synchronizes the Input and Target signals using Bayesian Stats
#
# Args:
#   Input (data.frame): The input data that needs to be synchronized.
#   Target (data.frame): The target data with which the input data will be synchronized.
#   folder (string, default = '~/Documents/BSynch/'): Path where output/results will be saved.
#   ta (numeric, default = 3): Parameter a of the t distribution.
#   tb (numeric, default = 4): Parameter b of the t distribution.
#   shape_acc (numeric, default = 10): Shape parameter for the expansion/compression parameter.
#   meanM (numeric, default = 0.7): Mean of the memory parameter.
#   alpha (numeric, default = 4): scale of the memory parameter.
#   iters (numeric, default = 2e+3): Number of iterations for the MCMC.
#   burn (numeric, default = 1e+3): Burn-in period for the MCMC.
#   thin (numeric, default = 150): Thinning parameter for the MCMC.
#


BSynch <- function(Input,Target,folder = '~/Documents/BSynch/',ta = 3,tb = 4,
                   shape_acc = 10,meanM = 0.7, alpha = 4 ,
                   iters = 2e+3, burn = 1e+3 ,thin = 50){ 

  # load data and twalk
  source_twalk(folder)
  inp <- load_file_from_folder(Input, folder) #read.csv(paste0(folder,Input), header=TRUE)
  tar <- load_file_from_folder(Target, folder) #read.csv(paste0(folder,Target), header=TRUE)
  load_or_install_dependencies()
  
  
  # parameters of the t-distribution
  # ta <- 3
  # tb <- 4
  # ## Estimate the shape and scale parameters using this formulation
  # shape_acc = 10 # NB: you can introduce autocorrelation in the gamma distribution by specifying a 
  # # memory parameters
  # meanM <- 0.7 # mean memory
  # alpha <- 4 # strength (higher values result in more peaked shapes)
  # # number of MCMC steps
  # iters <- 2e+3
  # burn <- 1e+3 
  # thin <- 150
  ########################
  # set.seed(123)
  # inflate the stdev of the target if necessary
  inf <- 1 
  # error associated with the input percentage fraction
  mySD <- 0.1
  
  
  #################################################################
  ## useful functions
  # let's scale the data (-1:1) accounting for the standard deviation
  range <- function(x){2*(x-min(x))/(max(x)-min(x))-1}
  #
  # density estimate of Student t distribution
  tdistro <- function(X, Mu, sigma, a, b){
    sigma = sigma^2
    -1* sum(( ((2*a+1.)/2.) * log(b + ((X - Mu)^2.)/(2.*sigma)) + .5 * log(sigma) ),na.rm = TRUE)
  }                  
  #################################################################
  
  ######################################
  # select input data
  
  # inp <- read.csv("~/Desktop/Francesco/LIG compilation/inputs/ODP1082.csv", header=TRUE)
  # remove data older than 140 ka (based on published age model)
  inp <- inp[inp$Age <= 140, ]
  # let's store the depth and age scale
  depth <- inp$MBSF 
  age_temp <- inp$Age
  # let's store the input on its published age scale
  inp <- cbind(inp$Age, inp$ProxyValue)
  #####################################
  # load target
  
  # tar <- read.table(file="~/Desktop/Francesco/LIG compilation/targets/EDC.txt", header=TRUE)
  #####################################
  
  ######################################
  # plot data
  par(mfrow=c(1,1))
  plot(tar, type='l', xlab="Age ka BP", ylab="proxy units", xlim=c(0,150))
  par(new=TRUE)
  plot(inp, type='l', col=4, xlab="", ylab="", xaxt='n', yaxt='n', xlim=c(0,150))
  legend("bottomright", lty=rep(1,2),
         c("input", "target"), col=c(1,4), cex=0.5)
  axis(4, col = 4)
  ######################################
  
  ######################################
  # subset both data sets so that target is longer than input
  # NB this setup assumes that the target is always longer than the input
  start1 <- head(inp[,1])[1]; end1 <- tail(inp[,1])[6] 
  start2 <- head(tar[,1])[1]; end2 <- tail(tar[,1])[6] 
  step1 <- round(mean(diff(inp[,1])),1); step2 <- round(mean(diff(tar[,1])),2) 
  # subset input data before synchronization
  # Note that the input doesn't need to be evenly interpolated for the synchronization to work
  input <- inp[inp[,1] >= start1 & inp[,1] <= end1,]
  target <- as.data.frame(approx(tar[,1], tar[,2], seq(start2, end2, by = step2)))
  ######################################
  
  ######################################
  # scale data
  input[,2] <- range(input[,2]) # flip the data here if necessary
  target[,2] <- range(target[,2])
  
  
  ##########################################################
  ##  Let's set up some parameters for our synchronization model
  # guesstimates of start and end age for input
  edge1 <- head(input[,1])[1]; edge2 <- tail(input[,1])[6] 
  # uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
  edge1_sd <- ifelse(head(input[,1])[1] < 10, 0.5, 5); edge2_sd <- 5 
  
  
  ######################################
  
  ######################################
  # plot scaled data
  # par(mfrow=c(1,1))
  # plot(input, type="l", col=4, ylim=c(-1,1),
  #      xlim=c(0,max(target[,1])), xlab="Age (years BP)", ylab="norm")
  # par(new=TRUE)
  # plot(target, type="l", xaxt='n', yaxt='n', xlab="", ylab="",
  #      col=1, ylim=c(-1,1), xlim=c(0,max(target[,1])))
  # legend("bottomleft", lty=rep(1,2),
  #        c("input", "target"), col=c(4,1), cex=0.5)
  # axis(4, col = 4)
  # abline(v=seq(0,200,10), lty=3, lwd=0.5)
  # #
  # hist(target[,2], breaks=50, xlab="norm", main="")
  # hist(input[,2], breaks=50, col=4, add=TRUE)
  ##########################################################
  
  
  ##########################################################
  # Calculate the interval size and the distance between consecutive nodes
  myLength <- tail(input[,1])[6] - head(input[,1])[1] # length of input
  N <- round(myLength,0) # This is the number of segments (one segment every ~1kyr)
  interval <- (myLength)/N # in meters
  nodes <- seq(head(input[,1])[1], tail(input[,1])[6], by = interval) # position of each node along stratigraphic depth
  # generate a random realization of the mean sed rates (in m/yr) of the core
  # sr <- (nodes[length(nodes)] - nodes[1])/(rnorm(1, mean = edge2, sd = edge2_sd) - rnorm(1, mean = edge1, sd = edge1_sd))
  sr <- (nodes[length(nodes)] - nodes[1])/(edge2 - edge1)
  ##########################################################
  # Set the parameters for the gamma distribution of sedimentation rates 
  mean <- sr
  # large value for the shape parameter. This results in a distribution with a small variance and 
  # a peak close to the mean.
  scale_acc = mean/shape_acc
  ##########################################################
  # Set the parameters for the beta distribution of the ory
  # Set the desired mean and alpha parameter
  
  # Calculate the beta parameter
  beta <- (alpha / meanM) - alpha
  ##########################################################
  
  ##########################################################
  # Set LIKELIHOOD
  ##########################################################
  likelihood <- function(params){
    ## Extract parameters from sampler and calculate acc rates
    d <- numeric(N)  # Initialize a numeric vector of length N
    # Acc rate initial segment
    d[1] <- params[3]  # Initialize the first element
    # Innovations and downcore accumulation rates 
    arate <- params[4:length(params)] 
    for (i in 2:N) {
      d[i] <- (params[2] * d[i-1]) + ((1 - params[2])*arate[i-1])
    } 
    # Calculate the sequence of points based on the counts and spacing
    points <- c(params[1], params[1] + cumsum(interval/d))
   # apply new age scale to input data
    new_input <- cbind(approx(nodes, points, input[,1])$y, input[,2])
   # interpolate target and its stdev to res new input
    new_target <- approx(target[,1], target[,2], new_input[,1])$y
   # predicted values should be close to observed values. Use t distribution
    likelihood <- tdistro(X = new_input[,2], Mu = new_target, sigma = inf*mySD, a = ta, b = tb)
    # sum the likelihood
    sumll = sum(likelihood)
    return(sumll)
  }
  
  ##########################################################
  # PRIORS: create a sequence of positive sed rates values
  # where each value depend on the previous one
  ##########################################################
  # uninformative proposal distribution
  sampler <-  function(n=1){
    # density of the first and last age point
    d1 <- runif(n, min=head(target[,1])[1],max=edge1 + (2*edge1_sd))
    # memory value
    mem <- runif(n, min=0.5, max=1)
    # downcore accumulation rates 
    arate <- runif(N, min = 0.5*sr, max = 2*sr) # Initialize a numeric vector of innovations of length N-1
    # arate <- runif(N, min = 0.9, max = 1.1) # Initialize a numeric vector of innovations of length N-1
    #
    d <- c(d1, mem, arate)
    return(cbind(d))
  }
  ##########################################################
  
  ##########################################################
  # Prior density
  density <-  function(params){
    d1 <- dnorm(params[1], mean = edge1, sd = edge1_sd,log = TRUE)
    # memory value
    dmem <- dbeta(params[2], alpha, beta, log = TRUE)
    # denisity for acc rates 
    darate <- dgamma(params[3:length(params)], shape = shape_acc, scale = scale_acc, log = TRUE)
    #
    d <- c(d1, dmem, darate) 
    return(sum(d))
  }
  ##########################################################
  
  ##########################################################
  # objective distritro
  obj <- function(param)(-(density(param)+likelihood(param)))
  
  ##########################################################
  
  ##########################################################
  # support function
  supp <- function(params){
    d <- numeric(N)  # Initialize a numeric vector of length N
    # Acc rate initial segment
    d[1] <- params[3]  # Initialize the first element
    # Innovations and downcore accumulation rates 
    arate <- params[4:length(params)] 
    for (i in 2:N) {
      d[i] <- (params[2] * d[i-1]) + ((1 - params[2])*arate[i-1])
    } 
    # Calculate the sequence of points based on the counts and spacing
    points <- c(params[1], params[1] + cumsum(interval/d))
    # ensure that the sequence is not longer than the target. 
    # note: here I change the or (|) to an and (&), double check with Francesco
    ifelse(all(c(params> low, params<up)) & (tail(points,1) < tail(target[,1],1) & 
                                               (points[length(points)] < (edge2 + (2*edge2_sd)) & 
                                                  points[length(points)] > (edge2 - (2*edge2_sd)))), 
           return(TRUE), 
           return(FALSE)
           )
  }
  
  
  # set lower and upper prior boundaries
  low <- c(head(target[,1])[1], # d1
           0, # mem
           rep(0.5*sr, N)) # arate
  up <- c(edge1+(2*edge1_sd), # d1
          1, # mem,
          rep(2*sr, N)) # arate
  
  # Run the twalk
  #############################
  
  
  x1 <- sampler()
  while(!(supp(x1))){
    x1 <- sampler()  
  }
  x2 <- sampler()
  while(!(supp(x2))){
    x2 <- sampler()  
  }
  
  tester <- TRUE
  while(tester){
    if(!file.exists(paste0(folder,'twalk_state.csv'))){
      output <- Runtwalk(iters,Obj = obj,dim = length(x1),x0 = x1,xp0 = x2,Supp =supp,
                         thinning = length(x1)*thin,burnin= length(x1)*burn   )
      # save the twalk state
      last.points <- matrix(c(tail(output$output,1), tail(output$outputp,1)),nrow = 2,byrow = T)
      write.table(x = last.points,file = paste0(folder,"twalk_state.csv"), row.names = F,col.names = F,sep=',')
      # Check convergance
      mcmc_samples <- as.mcmc(as.numeric(output$Us))
      plot(output$Us,type='l')
      test <- geweke.diag(mcmc_samples,frac1=0.5, frac2=0.5)
      if(abs(test$z) < .1){tester=FALSE}
      else{print(paste0('The geweke z-score is ',test$z,', we will run another set of iterations'))}
    }else{
      ini = as.matrix(read.csv(paste0(folder,'twalk_state.csv'),sep=',',header = F))
      # run for iterations
      print('There is a previus run. We will continue the chain')
      output <- Runtwalk(iters,Obj = obj,dim = length(x1) , x0 = as.numeric(ini[1,]),xp0 = as.numeric(ini[2,]),Supp =supp,
                         thinning = length(x1)*thin,burnin= length(x1)*burn   )
      # save the twalk state
      last.points <- matrix(c(tail(output$output,1), tail(output$outputp,1)),nrow = 2,byrow = T)
      write.table(x = last.points,file = paste0(folder,"twalk_state.csv"), row.names = F,col.names = F,sep=',')
      # Check convergance
      mcmc_samples <- as.mcmc(as.numeric(output$Us))
      plot(output$Us,type='l')
      test <- geweke.diag(mcmc_samples,frac1=0.5, frac2=0.5)
      if(abs(test$z) < .1){tester=FALSE}
      else{print(paste0('The geweke z-score is ',test$z,', we will run another set of iterations'))}
    }
  }
  cat("\n================== Geweke DIAGNOSTIC ==================\n")
  cat('Geweke value: ',round(test$z,3),"\n")
  cat('It appers that the chain has converge,',"\n") 
  cat('please refer to the energy plot to validate the convergance\n')
  cat("====================================================\n\n")  
  
  
  c <- output$output[-1,]
  energy <- output$Us
  energy2 <- output$Ups
  
  iat <- IAT(output,to=output$Tr)
  
  cat("\n================== IAT DIAGNOSTIC ==================\n")
  
  if (iat < 2) {
    cat("IAT Value:", iat, "\n")
    cat("Interpretation: The chain exhibits low correlation among successive samples.\n")
    cat("Recommendation: Current settings appear satisfactory.\n")
  } else {
    cat("IAT Value:", iat, "\n")
    cat("Interpretation: The chain exhibits high correlation among successive samples.\n")
    cat("Recommendation: Consider increasing the thinning value and rerunning the chain.\n")
  }
  
  cat("====================================================\n\n")  
  
  
  ####################################################
  ########### Plots                  #################
  ####################################################
  
  layout(matrix(c(1,1,5,5,5,5,5,5,
                  1,1,5,5,5,5,5,5,
                  2,2,6,6,6,6,6,6,
                  2,2,6,6,6,6,6,6,
                  3,3,6,6,6,6,6,6,
                  3,3,6,6,6,6,6,6,
                  4,4,6,6,6,6,6,6,
                  4,4,6,6,6,6,6,6), 8, 8, byrow = TRUE))
  
  
  #############################
  #############################

  plot(energy,type = 'l',col="grey40")
  lines(energy2,col="grey60")
  ####################################################
  # extract MCMC results and convert acc rates to ages
  accrate <- sapply(1:nrow(c), function(i){
    # extract relevant parameters and calculate acc rates
    mem <- c[i,2]
    d <- numeric(N)
    d[1] <- c[i,3]
    arate <- c[i,4:ncol(c)]
    for (j in 2:N) {
      d[j] <- (mem * d[j-1]) + ((1 - mem) * arate[j-1])
    } 
    result <- d
  })
  # estimate age of nodes
  ages <- sapply(1:nrow(c), function(i){
    points <- c(c[i,1], c[i,1] + cumsum(interval/accrate[,i]))
  })
  # estimate posterior memory
  memory <- quantile(c[,2], probs = c(0.5, 0.05, 0.95), na.rm =TRUE)
  # estimate quantiles
  quant <- apply(ages, 1, quantile, probs = c(0.5, 0.05, 0.95, 0.32, 0.68), na.rm =TRUE)
  # synchronized input
  align <- cbind(approx(nodes, quant[1,], input[,1])$y, input[,2])
  ####################################################
  
  ####################################################
  # compute the kernel density estimates for the priors and posteriors
  # of acc rates and memory, respectively
  {
    # generate the data for plotting ppriors
    x1 <- rgamma(100000, shape = shape_acc, scale = scale_acc)
    x2 <- rbeta(100000, alpha, beta)
    # estimate kernels
    h1_pr <- bkde(x1)
    h2_pr <- bkde(x2)
    h1_po <- bkde(interval/accrate[1:length(nodes)-1,], bandwidth=0.1)
    h2_po <- bkde(c[,2], bandwidth=0.1)
    # plot graphs
    if (max(h1_po$y)>max(h1_pr$y)) {hgt = max(h1_po$y)} else {hgt = max(h1_pr$y)}
    plot(h1_pr$x, h1_pr$y, type='l', col=1, ylim=c(0,hgt),xlim=c(tail(low,1),tail(up,1)),
         xlab="", ylab="density")
    title(xlab = "acc rate ratios", line = 2, cex.lab = 1)
    polygon(c(h1_pr$x, rev(h1_pr$x)), c(h1_pr$y, rep(0, length(h1_pr$y))),
            col = "grey", border = NA)
    lines(h1_po$x, h1_po$y, col="darkblue")
    polygon(c(h1_po$x, rev(h1_po$x)), c(h1_po$y, rep(0, length(h1_po$y))),
            col = adjustcolor("dodgerblue", alpha=0.25), border = NA)
    legend("topright", c("prior", "posterior"), 
           fill = c("grey", "dodgerblue"), bg = "white", cex=0.5)
  #######
    if (max(h2_po$y)>max(h2_pr$y)) {hgt = max(h2_po$y)} else {hgt = max(h2_pr$y)}
    plot(h2_pr$x, h2_pr$y, type='l', col=1, xlim=c(0,1), ylim=c(0,hgt),
         xlab="", ylab="")
    title(xlab = "memory", line = 2, cex.lab = 1)
    polygon(c(h2_pr$x, rev(h2_pr$x)), c(h2_pr$y, rep(0, length(h2_pr$y))),
            col = "grey", border = NA)
    lines(h2_po$x, h2_po$y, col="darkblue")
    polygon(c(h2_po$x, rev(h2_po$x)), c(h2_po$y, rep(0, length(h2_po$y))),
            col = adjustcolor("dodgerblue", alpha=0.25), border = NA)
  }
  ####################################################
  
  {
    plot(quant[2,]-quant[1,], quant[1,], yaxt='n',
         ylim=c(0,max(align[,1])), 
         xlim=c(-(max(c(quant[3,]-quant[1,], quant[2,]-quant[1,]))), max(c(quant[3,]-quant[1,], quant[2,]-quant[1,]))),
         lty=3, type='l', ylab="", xlab="Posterior credibility (ka)")
    lines(quant[3,]-quant[1,], quant[1,], lty=3)
    polygon(c(quant[2,]-quant[1,], rev(quant[3,]-quant[1,])), 
            c(quant[1,], rev(quant[1,])),
            col = adjustcolor("grey", alpha=0.5), border = NA)
    axis(2, labels = FALSE)
    if (max(c(quant[3,]-quant[2,]))>1) {
      grid = seq(-10,10,by=1.25)
    } else {
      grid = seq(-10,10,by=0.5)
    }
    abline(v=grid,lty=3, lwd=1.5)
    abline(h=seq(0,500,by=10), lty=3, lwd=0.5)
    # par(mar=c(5.1,4.1,4.1,5.1))
    
    mydepth <- approx(age_temp, depth, nodes)$y
    # par(mar=c(5.1,4.1,4.1,23.1))
    plot(quant[1,],mydepth,  type='l', xlim=c(0,max(align[,1])),
         ylab="Input depth (MBSF)", xlab="Target age (ka BP)")
    lines(quant[2,],mydepth,  lty=3, lwd=0.5)
    lines(quant[3,],mydepth,  lty=3, lwd=0.5)
    lines(age_temp,depth,  col="orange", lty=5)
    polygon(
            c(quant[2,], rev(quant[3,])),
            c(mydepth, rev(mydepth)),
            col = adjustcolor("grey", alpha=0.5), border = NA)
    abline(h=seq(0,500,by=1),lty=3, lwd=0.5)
    abline(v=seq(0,500,by=10), lty=3, lwd=0.5)
    legend("topleft", c("posterior","median","published age model", "95% Cl"),
           col = c("grey", 1, "orange",1), pch=c(15,NA,NA,NA), lty=c(NA,1,5,3), cex=0.5, bg = "white")
  }
  
  
  
  
  ####################################################
  # plot results
  {
    plot(align, type='o', xlim=c(0,150), cex=0.35, pch=19,
         ylim=c(-1,1), xlab="", ylab="norm", col=4)
    title(xlab = "AICC2012 ka BP", line = 2, cex.lab = 1)
    lines(input, col="orange", lwd=0.5, lty=5)
    lines(target, col=1, lwd=0.5)
    abline(v=seq(0,500,10), lty=3, lwd=0.5)
    legend("topleft", c("input (synchronized)", "input (published age model)", "target"),
           lty=c(1,5,1), col=c(4, "orange", 1), cex=0.5)
  }
  ####################################################
  

  # Save output
  write.csv(as.matrix(c),paste0(folder,"twalk_output.csv"),row.names = F)
  # Specify the file path and name for the PDF
  pdf_file_path <- paste0(folder,'/aligment.pdf')
  # Open the PDF device
  pdf(pdf_file_path)
  # Close the PDF device
  dev.off()


}



# loader function
load_file_from_folder <- function(Input, folder) {
  
  # Check for .csv extension
  csv_path <- paste0(folder, Input, ".csv")
  if (file.exists(csv_path)) {
    fil <- read.csv(csv_path)
    return(fil)
  }
  
  # Check for .txt extension
  txt_path <- paste0(folder, Input, ".txt")
  if (file.exists(txt_path)) {
    fil <- read.table(txt_path, header = TRUE, sep = "\t") # Assuming tab-separated for .txt files
    return(fil)
  }
  
  # Return NULL with a warning if no matching file is found
  warning("No matching .csv or .txt file found for the given input in the specified folder.")
  return(NULL)
}




source_twalk <- function(folder) {
  # Construct the path to the twalk.R file in the specified folder
  twalk_path <- file.path(folder, "..", "twalk.R")
  # Check if the twalk.R file exists in the specified folder
  if (file.exists(twalk_path)) {
    source(twalk_path)
    cat("Successfully sourced twalk.R from", folder, "\n")
  } else {
    warning("File twalk.R was not found in the specified folder.")
  }
}


load_or_install_dependencies <- function() {
  if (!require(KernSmooth, quietly = TRUE)) {
    cat("KernSmooth not found. Installing...\n")
    install.packages("KernSmooth", dependencies = TRUE)
    
    # Load after installing
    library(KernSmooth)
    cat("KernSmooth installed and loaded successfully!\n")
  } else {
    cat("KernSmooth loaded successfully!\n")
  }
  if (!require(coda, quietly = TRUE)) {
    cat("coda not found. Installing...\n")
    install.packages("coda", dependencies = TRUE)
    
    # Load after installing
    library(coda)
    cat("coda installed and loaded successfully!\n")
  } else {
    cat("coda loaded successfully!\n")
  }
}



# BSynch(Input='MD01_2444',Target='GreenStack',folder = '~/Documents/BSynch/MD01_2444-GreenStack/',thin=50,burn=1e+3,iters=2e+3)



