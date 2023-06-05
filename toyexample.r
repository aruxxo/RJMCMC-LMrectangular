############# 
## EXAMPLE ##
#############



# Generata data from a Latent Markov model
# with Ti=4 and true sequence 4,4,3,3
# Observations are in the form of 200 multivariate responses, of dimensionality r=3,
# over 4 time occasions. The groups have a degree of separation equal to 4.

# We set default value b=1.5 for all the coefficients modulating
# initial and transition probs, except reference categories and intercepts.


source("_functionsrlm.r")
load("datarlm.RData")


# Model-based clustering for initial values
initials=getInits(data$y,data$Z,5,b=1.5)

# Calibration for Centroids
mxi=0 # prior mean
vxi=1 # prior variance

# Misspecify the sequence of latent states 
ki = c(4,3,4,3)

# Run the algorithm for 1000 iterations performing reversible jumps only in firts 100
set.seed(321)
#out = rlmrjmcmc(y=data$y, w=1, Z=data$Z, n.iter = 1000, ki=initials$k,kmax=5,inits = TRUE, initials = initials, verbose = T,cutOffRj=100)
# uncomment the line above if you wish to actually run the algorithm
# the following loads object out: 
load("outrlm.RData")

apply(out$xi[,3,1:3,],2:3,mean)
apply(out$xi[,4,1:4,],2:3,mean)

apply(out$sigma[,3,1:3,],2:3,mean)
apply(out$sigma[,4,1:4,],2:3,mean)


## Estimated coefficients for initial and transitions probs ##

apply(out$be[,,4,1:4],2:3,mean)

# Transition 44 #

apply(out$Be[,,4,4,1:4,1:4],2:4,mean)

# Transition 43 #

apply(out$Be[,,4,3,1:4,1:3],2:4,mean)

# Transition 33 #

apply(out$Be[,,3,3,1:3,1:3],2:4,mean)

