This file describes the main "rlmrjmcmc" function to perform Reversible Jump MCMC for a rectangular Latent Markov model with unknown number of regime profiles and latent parameters modulated by covariates.

It accompanyes the paper "Covariate-modulated rectangular Latent Markov models with an unknown number of regime profiles" by Alfonso Russo, Alessio Farcomeni, Maria Grazia Pittau, and Roberto Zelli

## INPUT ARGUMENTS ##

y : an N x Ti x r array of response variables
w : an NxTi matrix of weights, default is 1 in each entry 
Z : an N x Ti x p array of covariates
n.iter : number of iterations to be performed
ki : a vector of length Ti specifying initial values for the sequence of latent states
kmax : an integer number giving the maximum number of states admitted a priori at each time occasion
inits : a logical value indicating whether initial values are given to the function or should be generated
initials : a list of initial values for all the parameters involved
verbose : a logical value indicating whether or not to print progress in the algorithm
msi : prior parameter for the ads of the components, default is 2
vsi : prior parameter for the ads of the components, default is 2
mxi : prior location parameter for the centroids of the components, default is 0
vxi : prior scale parameter for the centroids of the components, default is 1
Tauxi : calibration parameter for auxiliary variables in MH steps for the centroids, default is rep(0.01,r)
mbe : prior location parameter for coefficients modulating initial probs, default is 0
vbe : prior scale parameter for coefficients modulating initial probs, default is 1
mBe : prior location parameter for coefficients modulating transition probs, default is 0
vBe : prior scale parameter for coefficients modulating transition probs, default is 1.5
mu : calibration parameter for auxiliary variables in split/combine steps for the centroids, default is 0
vu : calibration parameter for auxiliary variables in split/combine steps for the centroids, default is 1
mw : calibration parameter for auxiliary variables in split/combine steps for the SDs, default is 2
vw : calibration parameter for auxiliary variables in split/combine steps for the SDs, default is 2
mepsb : calibration parameter for auxiliary variables in split/combine steps for the coeffs initial probs, default is 0
vepsb : calibration parameter for auxiliary variables in split/combine steps for the coeffs initial probs, default is 1
mepsB : calibration parameter for auxiliary variables in split/combine steps for the coeffs transit. probs, default is 0
vepsB : calibration parameter for auxiliary variables in split/combine steps for the coeffs transit. probs, default is 1
thB : threshold for initial and transition probabilities, default is 10^-5
cutOffRj : indicates the iteration at which transdimensional steps should stop and continue only with MH, default is Inf

## STRUCTURE ##
After initialisation the function performs three different types of steps.
Birth/Death and Split/Combine are meant for transdimensional sampling to update the sequence of latent states
Metropolis-Hasting steps update all the parameters involved, conditionally on the actual sequence

## OUTPUT ##
The function returns a list with the following objects

liks : a matrix of suitable dimensions storing the values of the likelihood at each iteration of the algorithm
k : the output matrix storing the values for the sequence of latent states as estimated at each iteration of the alghorithm
be : the output array storing the values for the coefficients modulating initial probs at each iteration of the algorithm, conditionally on the sequence of latent groups
Be : the output array storing the values for the coefficients modulating transition probs at each iteration of the algorithm, conditionally on the sequence of latent groups
xi : the output array storing the values for the centroids at each iteration of the algorithm, conditionally on the sequence of latent groups
sigma : the output array storing the values for the SDs at each iteration of the algorithm, conditionally on the sequence of latent groups



