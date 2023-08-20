# load libraries
require(BayesianTools)
require(truncnorm)
require(coda)
require(KernSmooth)
require(fields)

# set working directory and source the data pre-treatment code
# setwd("~/Desktop/Francesco/LIG compilation")
setwd("~/MEGA_A/Cambridge/alignment MCMC/")

#### 
# iniciate profiling
# profvis(
  # {
##  Let's set up some parameters for our synchronization model

# inflate the stdev of the target if necessary
inf <- 1 
# error associated with the input percentage fraction
mySD <- 0.1
# parameters of the t-distribution
ta <- 3
tb <- 4
# number of MCMC steps
iters <- 4e+5
burn <- 1e+5

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
inp <- read.csv("~/MEGA_A/Cambridge/alignment MCMC//MD01-2444.csv", header=TRUE)
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
tar <- read.table(file="~/MEGA_A/Cambridge/alignment MCMC//GreenStack.txt", header=TRUE)
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

# guesstimates of start and end age for input
edge1 <- head(input[,1])[1]; edge2 <- tail(input[,1])[6] 
# uncertainty of start and end age for input (assume an error of 0.25ka for start <10ka BP and 5ka for >10ka BP)
edge1_sd <- ifelse(head(input[,1])[1] < 10, 0.5, 5); edge2_sd <- 5 

######################################

######################################
# plot scaled data
par(mfrow=c(1,1))
plot(input, type="l", col=4, ylim=c(-1,1),
     xlim=c(0,max(target[,1])), xlab="Age (years BP)", ylab="norm")
par(new=TRUE)
plot(target, type="l", xaxt='n', yaxt='n', xlab="", ylab="",
     col=1, ylim=c(-1,1), xlim=c(0,max(target[,1])))
legend("bottomleft", lty=rep(1,2),
       c("input", "target"), col=c(4,1), cex=0.5)
axis(4, col = 4)
abline(v=seq(0,200,10), lty=3, lwd=0.5)
#
hist(target[,2], breaks=50, xlab="norm", main="")
hist(input[,2], breaks=50, col=4, add=TRUE)
##########################################################

##########################################################

##########################################################
# Calculate the interval size and the distance between consecutive nodes
myLength <- tail(input[,1])[6] - head(input[,1])[1] # length of input
N <- round(myLength,0) # This is the number of segments (one segment every ~1kyr)
interval <- (myLength)/N # in meters
nodes <- seq(head(input[,1])[1], tail(input[,1])[6], by = interval) # position of each node along stratigraphic depth
# generate a random realization of the mean sed rates (in m/yr) of the core
sr <- (nodes[length(nodes)] - nodes[1])/(rnorm(1, mean = edge2, sd = edge2_sd) - rnorm(1, mean = edge1, sd = edge1_sd))
##########################################################
# Set the parameters for the gamma distribution of sedimentation rates 
mean <- sr
## Estimate the shape and scale parameters using this formulation
shape_acc = 1.5 # NB: you can introduce autocorrelation in the gamma distribution by specifying a 
# large value for the shape parameter. This results in a distribution with a small variance and 
# a peak close to the mean.
scale_acc = mean/shape_acc
##########################################################
# Set the parameters for the beta distribution of the memory
# Set the desired mean and alpha parameter
meanM <- 0.5 # mean memory
alpha <- 4 # strength (higher values result in more peaked shapes)
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
  # ensure that the sequence is not longer than the target. If so, yield -Inf 
  ifelse(points[length(points)] > tail(target[,1])[6] | points[length(points)] < (edge2 - (2*edge2_sd)), 
         {likelihood = -Inf},
         {
           # apply new age scale to input data
           new_input <- cbind(approx(nodes, points, input[,1])$y, input[,2])
           # interpolate target and its stdev to res new input
           new_target <- approx(target[,1], target[,2], new_input[,1])$y
           # predicted values should be close to observed values. Use t distribution
           likelihood <- tdistro(X = new_input[,2], Mu = new_target, sigma = inf*mySD, a = ta, b = tb)
         })
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
  d1 <- runif(n, min=head(target[,1])[1],
              max=edge1 + (2*edge1_sd))
  # memory value
  mem <- runif(n, min=0, max=1)
  # downcore accumulation rates 
  arate <- runif(N, min = 0.5*sr, max = 2*sr) # Initialize a numeric vector of innovations of length N-1
  #
  d <- c(d1, mem, arate)
  return(cbind(d))
}
##########################################################

##########################################################
# Prior density
density <-  function(params){
  d1 <- log(dtruncnorm(params[1], mean = edge1, sd = edge1_sd,
                       a=head(target[,1])[1], b=edge1 + (2*edge1_sd)))
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
# set lower and upper prior boundaries
low <- c(head(target[,1])[1], # d1
         0, # mem
         rep(0.5*sr, N)) # arate
up <- c(edge1+(2*edge1_sd), # d1
        1, # mem
        rep(2*sr, N)) # arate
# create the prior setup
prior <- createPrior(sampler = sampler, density = density, lower = low, upper = up, best = NULL)

# collect prior, likelihood and posterior functions
bayesianSetup <- createBayesianSetup(likelihood = likelihood, prior = prior,
                                     catchDuplicates = FALSE)
# setup MCMC simulation
settings = list(iterations = 3*(iters+burn), burnin = 0)
# run MCMC chains

start.time <- Sys.time()
chain <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs",settings = settings)
end.time <- Sys.time()
time.taken.BT <- end.time - start.time
time.taken.BT


# chain <- runMCMC(bayesianSetup = chain, sampler = "DEzs",settings = settings) # continue a previous MCMC run
#############################

#############################

gelmanDiagnostics(chain)
layout(matrix(c(1,1,2,2,3,3,4,4,
                1,1,2,2,3,3,4,4,
                5,5,6,6,6,6,6,6,
                5,5,6,6,6,6,6,6,
                5,5,6,6,6,6,6,6,
                5,5,6,6,6,6,6,6,
                5,5,6,6,6,6,6,6), 7, 8, byrow = TRUE))

# Plot chain

energy = as.numeric(chain$chain[[1]][,length(low)+1])
energy <- energy[!is.infinite(energy)]
plot(energy,
          ylim=c(min(energy[!is.na(energy)]),max(energy[!is.na(energy)])),
          type='l', ylab="log Likelihood",col="grey40")
for(i in 2:3){
  energy = chain$chain[[i]][,length(low)+1]
  energy <- energy[!is.infinite(energy)]
  col = c( "grey60", "grey80")
  lines(energy, col=col[i])
}
# extract posterior samples
c <- getSample(chain, start = 0)
#############################
# additional burn-in
burnIn = iters/2
c <- c[-(1:burnIn),] # remove burn in steps
# thinning MCMC output
thin <- 10
c <- c[seq(1, nrow(c), thin),] # extract only every nth output
#############################
# check convergence (should be below 1.05 for all parameters to demonstrate convergence)

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
  # par(mfrow=c(1,1))
  # par(mar=c(23.1,4.1,1.1,25.1))
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
  plot(h1_pr$x, h1_pr$y, type='l', col=1, ylim=c(0,hgt),
       xlab="", ylab="density")
  title(xlab = "acc rate ratios", line = 2, cex.lab = 1)
  polygon(c(h1_pr$x, rev(h1_pr$x)), c(h1_pr$y, rep(0, length(h1_pr$y))),
          col = "grey", border = NA)
  lines(h1_po$x, h1_po$y, col="darkblue")
  polygon(c(h1_po$x, rev(h1_po$x)), c(h1_po$y, rep(0, length(h1_po$y))),
          col = adjustcolor("dodgerblue", alpha=0.25), border = NA)
  legend("topright", c("prior", "posterior"), 
         fill = c("grey", "dodgerblue"), bg = "white", cex=0.5)
  # par(new=TRUE)
  # par(mar=c(23.1,27.1,1.1,2.1))
  if (max(h2_po$y)>max(h2_pr$y)) {hgt = max(h2_po$y)} else {hgt = max(h2_pr$y)}
  plot(h2_pr$x, h2_pr$y, type='l', col=1, xlim=c(0,1), ylim=c(0,hgt),
       xlab="", ylab="")
  title(xlab = "memory", line = 2, cex.lab = 1)
  polygon(c(h2_pr$x, rev(h2_pr$x)), c(h2_pr$y, rep(0, length(h2_pr$y))),
          col = "grey", border = NA)
  lines(h2_po$x, h2_po$y, col="darkblue")
  polygon(c(h2_po$x, rev(h2_po$x)), c(h2_po$y, rep(0, length(h2_po$y))),
          col = adjustcolor("dodgerblue", alpha=0.25), border = NA)
  # par(mar=c(5.1,4.1,4.1,5.1))
}
####################################################

{
  
  # par(new=TRUE)
  # par(mar=c(5.1,26.1,4.1,5.1))
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
    grid = seq(-10,10,by=0.25)
  } else {
    grid = seq(-10,10,by=0.1)
  }
  abline(v=grid,lty=3, lwd=0.5)
  abline(h=seq(0,500,by=10), lty=3, lwd=0.5)
  # par(mar=c(5.1,4.1,4.1,5.1))
  
  mydepth <- approx(age_temp, depth, nodes)$y
  # par(mar=c(5.1,4.1,4.1,23.1))
  plot(mydepth, quant[1,], type='l', ylim=c(0,max(align[,1])),
       xlab="Input depth (MBSF)", ylab="Target age (ka BP)")
  lines(mydepth, quant[2,], lty=3, lwd=0.5)
  lines(mydepth, quant[3,], lty=3, lwd=0.5)
  lines(depth, age_temp, col="orange", lty=5)
  polygon(c(mydepth, rev(mydepth)),
          c(quant[2,], rev(quant[3,])),
          col = adjustcolor("grey", alpha=0.5), border = NA)
  abline(v=seq(0,500,by=1),lty=3, lwd=0.5)
  abline(h=seq(0,500,by=10), lty=3, lwd=0.5)
  legend("topright", c("posterior","median","published age model", "95% Cl"),
         col = c("grey", 1, "orange",1), pch=c(15,NA,NA,NA), lty=c(NA,1,5,3), cex=0.5, bg = "white")
}




####################################################
# plot results
{
  # par(new=TRUE)
  # par(mar=c(3.1,4.1,9.1,2.1))
  # par(mfrow=c(1,1))
  plot(align, type='o', xlim=c(0,150), cex=0.35, pch=19,
       ylim=c(-1,1), xlab="", ylab="norm", col=4)
  title(xlab = "AICC2012 ka BP", line = 2, cex.lab = 1)
  lines(input, col="orange", lwd=0.5, lty=5)
  # polygon(c(align[,1], rev(align[,1])),
  #         c(align[,2]+mySD, rev(align[,2]-mySD)),
  #         col = adjustcolor("dodgerblue", alpha=0.25), border = NA)
  lines(target, col=1, lwd=0.5)
  abline(v=seq(0,500,10), lty=3, lwd=0.5)
  legend("topleft", c("input (synchronized)", "input (published age model)", "target"),
         lty=c(1,5,1), col=c(4, "orange", 1), cex=0.5)
  # par(mar=c(5.1,4.1,4.1,5.1))
}
####################################################

####################################################

# ##########
# # write up
# write.table(cbind(depth, round(align[,1],3), 
#                   round(inp[,2],3)),
#                   file="~/Desktop/Francesco/LIG compilation/output/MD01-2444.txt", 
#                   row.names=FALSE,
#                   col.names=c("Depth_mbsl","Median_ka", "proxy"))
# ##########

print(time.taken.BT)
energybt <- energy
samplebt <- c