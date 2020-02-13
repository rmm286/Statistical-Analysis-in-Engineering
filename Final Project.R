##### Final Project ######
### MAE 207
### Rayne Milner

#DRE: r = L*9.81*W/(2*pi*r^3)

#Data Input
#trial 1
r1 <- 0.5*(1/1000)*c(9.4,9.71,9.63,9.64,9.62,9.72,9.63,9.54,9.57,9.97,9.7) #radius of chaulk in m
L1 <- 46.91*(1/1000) #length of gap in m
force1 <-  c(2.375,2.005,2.47,2.325,2.48,2.145,2.04,2.39,2.115,2.27,2.515) #Breaking force in N

#trial 2
r2 <- 0.5*(1/1000)*c(9.73,9.7,9.58,9.67,9.55,9.7,9.6,9.58,9.65,9.57,9.61)#radius of chaulk in m
L2 <- 6.81*(1/1000)#length of gap in m
force2 <-  c(1.455,1.39,1.04,1.175,1.175,1.385,1.075,1.21,1.46,1.375,0.84)#Breaking force in N

#trial 3
r3 <- 0.5*(1/1000)*c(14.7,14.68,14.71,14.78,14.72,14.87,14.78,14.7,14.75,14.71,14.69)#radius of chaulk in m
L3 <- 45.94*(1/1000)#length of gap in m
force3 <-  c(4.95,4.92,4.7,6.515,5.86,2.835,3.855,6.005,4.53,5.64,5.605)#Breaking force in N

#trial 4
r4 <- 0.5*(1/1000)*c(14.68,14.67,14.98,14.66,15.03,14.7,14.81,15.09,14.67,15.03,14.97)#radius of chaulk in m
L4 <- 35.62*(1/1000)#length of gap in m
force4 <-  c(7.135,2.385,5.825,6.22,6.025,5.28,7.345,5.17,5.32,6.255,5.54)#Breaking force in N


#Outlier Identification
#Here we identify the outliers using the modified thompson tau technique. 
#we choose to keep the outliers because we already have very little data.
P <- 0.95

#trial 1
#for diameter
x <- r1
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there is one

#for breaking force
x <- force1
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there are none

#trial 2
#for diameter
x <- r2
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there are none
#for breaking force
x <- force2
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there is one

#trial 3
#for diameter
x <- r3
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there are none

#for breaking force
x <- force3
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there is one

#trial 4
#for diameter
x <- r4
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there are none

#for breaking force
x <- force4
N <- length(x)
t <- qt(p = (1+P)/2, df = N-1)
tau <- t*(N-1)/(sqrt(N)*sqrt(N-2+t^2))
deltaMax <- tau*sd(x)
outliers <- x[which(abs(x-mean(x))>deltaMax)]
outliers #there is one

#Uncertainty Analysis using Taylor Series Approach
#Here we aim to estimate the mean and uncertianty of the rupture strength
#of the chalk using the Taylor Series approach.
#note that we are assuming with 11 data points that the value of K is approximately 2
#and we are not using the welch satterwise eqn.
#DRE: r = L*F/(2*pi*r^3)

betas <- c(0.01*1E-3,0.01*1E-3,10*1E-3)
# uncertaintites, described below
# betas[1] => accuracy error associated with length measurments, assumed to be normal
# betas[2] => digitization error associated with length measuments, assumed to be normal
# betas[3] => accuracy uncertanty associated with force measuremnts(F), assumed ot be normal


#calculate B values
#b11 associated with Length measurment
b11 <- (betas[1]^2 + betas[2]^2)/1.96^2

#b22 associated with force measurement
b22 <- (betas[3]^2)/1.96^2

#b33 associated with radius measurement
b33 <- b11

#common errors: 1 is common with 3 but 2 is not common with either
b12 <- 0
b23 <- 0
b13 <- b11

#trial 1
#calculate mean and sd values
N <- length(r1)
radbar <- mean(r1)
sdrad <- sd(r1)
sdradbar <-sdrad/sqrt(N)
Fbar <- mean(force1)
sdF <- sd(force1)
sdFbar <- sdF/sqrt(N)
Lbar <- L1

#calculate r, rbar and sd r
r <- L1*force1/(2*pi*(r1)^3)
rbar <- mean(r)
rbar
sdr <- sd(r)
sdrbar <- sdr/sqrt(N) #standard deviation of rbar using direct approach, this includes covariance effects

#calculate Theta values

#Theta1 = dr/dL = g*w/(2*pi*r^3)
Theta1 <- 9.81*Fbar/(2*pi*(radbar)^3)

#Theta2 = dr/dF = L/(2*pi*r^3)
Theta2 <- 9.81*Lbar/(2*pi*(radbar)^3)

#Theta3 = dr/drad = 3*F*L/(2*pi*r^4)
Theta3 <- 9.81*3*Fbar*Lbar/(2*pi*(radbar)^4)

#calculate br
brsqr <- Theta1^2*b11 + Theta2^2*b22 + Theta3^2*33 + 2*Theta1*Theta2*b12 + 2*Theta1*Theta3*b13 + 2*Theta2*Theta3*b23

#calculate urbar and Urbar
urbar <- (sdrbar^2+brsqr)^(1/2)
Urbar <- 2*urbar

rbar
Urbar
rbar+Urbar
rbar-Urbar

#trial 2
#calculate mean and sd values
N <- length(r2)
radbar <- mean(r2)
sdrad <- sd(r2)
sdradbar <-sdrad/sqrt(N)
Fbar <- mean(force2)
sdF <- sd(force2)
sdFbar <- sdF/sqrt(N)
Lbar <- L2

#calculate r, rbar and sd r
r <- L2*force2/(2*pi*(r2)^3)
rbar <- mean(r)
rbar
sdr <- sd(r)
sdrbar <- sdr/sqrt(N) #standard deviation of rbar using direct approach, this includes covariance effects

#calculate Theta values

#Theta1 = dr/dL = F/(2*pi*r^3)
Theta1 <- 9.81*Fbar/(2*pi*(radbar)^3)

#Theta2 = dr/dF = L/(2*pi*r^3)
Theta2 <- 9.81*Lbar/(2*pi*(radbar)^3)

#Theta3 = dr/drad = 3*F*L/(2*pi*r^4)
Theta3 <- 9.81*3*Fbar*Lbar/(2*pi*(radbar)^4)

#calculate br
brsqr <- Theta1^2*b11 + Theta2^2*b22 + Theta3^3*b33 + 2*Theta1*Theta2*b12 + 2*Theta1*Theta3*b13 + 2*Theta2*Theta3*b23

#calculate urbar and Urbar
urbar <- (sdrbar^2+brsqr)^(1/2)
Urbar <- 2*urbar

rbar
Urbar
rbar+Urbar
rbar-Urbar

#trial 3
#calculate mean and sd values
N <- length(r3)
radbar <- mean(r3)
sdrad <- sd(r3)
sdradbar <-sdrad/sqrt(N)
Fbar <- mean(force3)
sdF <- sd(force3)
sdFbar <- sdF/sqrt(N)
Lbar <- L3

#calculate r, rbar and sd r
r <- L3*force3/(2*pi*(r3)^3)
rbar <- mean(r)
rbar
sdr <- sd(r)
sdrbar <- sdr/sqrt(N) #standard deviation of rbar using direct approach, this includes covariance effects

#calculate Theta values

#Theta1 = dr/dL = g*F/(2*pi*r^3)
Theta1 <- 9.81*Fbar/(2*pi*(radbar)^3)

#Theta2 = dr/dF = g*L/(2*pi*r^3)
Theta2 <- 9.81*Lbar/(2*pi*(radbar)^3)

#Theta3 = dr/drad = g*3*F*L/(2*pi*r^4)
Theta3 <- 9.81*3*Fbar*Lbar/(2*pi*(radbar)^4)

#calculate br

brsqr <- Theta1^2*b11 + Theta2^2*b22 + Theta3^3*b33 + 2*Theta1*Theta2*b12 + 2*Theta1*Theta3*b13 + 2*Theta2*Theta3*b23

urbar <- (sdrbar^2+brsqr)^(1/2)

Urbar <- 2*urbar

rbar
Urbar
rbar+Urbar
rbar-Urbar

#trial 4
#calculate mean and sd values
N <- length(r4)
radbar <- mean(r4)
sdrad <- sd(r4)
sdradbar <-sdrad/sqrt(N)
Fbar <- mean(force4)
sdF <- sd(force4)
sdFbar <- sdF/sqrt(N)
Lbar <- L4

#calculate r, rbar and sd r
r <- L4*force4/(2*pi*(r4)^3)
rbar <- mean(r)
rbar
sdr <- sd(r)
sdrbar <- sdr/sqrt(N) #standard deviation of rbar using direct approach, this includes covariance effects

#calculate Theta values

#Theta1 = dr/dL = g*F/(2*pi*r^3)
Theta1 <- 9.81*Fbar/(2*pi*(radbar)^3)

#Theta2 = dr/dF = g*L/(2*pi*r^3)
Theta2 <- 9.81*Lbar/(2*pi*(radbar)^3)

#Theta3 = dr/drad = g*3*F*L/(2*pi*r^4)
Theta3 <- 9.81*3*Fbar*Lbar/(2*pi*(radbar)^4)

#calculate br
brsqr <- Theta1^2*b11 + Theta2^2*b22 + Theta3^3*b33 + 2*Theta1*Theta2*b12 + 2*Theta1*Theta3*b13 + 2*Theta2*Theta3*b23

#caculate urbar and Urbar
urbar <- (sdrbar^2+brsqr)^(1/2)
Urbar <- 2*urbar

rbar
Urbar
rbar+Urbar
rbar-Urbar

#Uncertainty Analysis using Monte Carlo Approach
#Here we aim to estimate the mean and uncertianty of the rupture strength
#We are using a parametric approach
M <- 2*10^5 # nuber of replicates
P <- 0.95

#uncertainties
b1 <- betas[1]/2 #accuracy associated with length
b2 <- betas[2]/2 #digitization associated with length
b3 <- betas[3]/2 #systematic associated with force

#trial 1
#calculate mean and sd values

rad <- r1
L <- L1
force <- force1

boot.r <- numeric(M)
for (i in 1:M){
  boot.sample.rad <- sample(rad,size=length(rad), replace=T)
  boot.sample.force <- sample(force,size=length(force),replace=T)
  beta1rad <- rnorm(n=1,mean=0,sd=b1)
  beta2rad <- rnorm(n=1,mean=0,sd=b2)
  beta1L <- rnorm(n=1,mean=0,sd=b1)
  beta2L <- rnorm(n=1,mean=0,sd=b2)
  beta1force <- rnorm(n=1,mean=0,sd=b3)
  
  rad_s <- mean(boot.sample.rad) + beta1rad+beta2rad
  L_s <- L + beta1L + beta2L
  force_s <- mean(boot.sample.force)+beta1force
  
  boot.r[i] <- 9.81*force_s*L_s/(2*pi*(rad_s)^3)
}

hist(boot.r, main = "Histogram of MC calculation", xlab = "Result")
mean(boot.r)
quantile(boot.r,probs = c((1-P)/2,(1+P)/2))

#trial 2
#calculate mean and sd values

rad <- r2
L <- L2
force <- force2

boot.r <- numeric(M)
for (i in 1:M){
  boot.sample.rad <- sample(rad,size=length(rad), replace=T)
  boot.sample.force <- sample(force,size=length(force),replace=T)
  beta1rad <- rnorm(n=1,mean=0,sd=b1)
  beta2rad <- rnorm(n=1,mean=0,sd=b2)
  beta1L <- rnorm(n=1,mean=0,sd=b1)
  beta2L <- rnorm(n=1,mean=0,sd=b2)
  beta1force <- rnorm(n=1,mean=0,sd=b3)
  
  rad_s <- mean(boot.sample.rad) + beta1rad+beta2rad
  L_s <- L + beta1L + beta2L
  force_s <- mean(boot.sample.force)+beta1force
  
  boot.r[i] <- 9.81*force_s*L_s/(2*pi*(rad_s)^3)
}

hist(boot.r, main = "Histogram of MC calculation", xlab = "Result")
mean(boot.r)
quantile(boot.r,probs = c((1-P)/2,(1+P)/2))

#trial 3
#calculate mean and sd values

rad <- r3
L <- L3
force <- force3

boot.r <- numeric(M)
for (i in 1:M){
  boot.sample.rad <- sample(rad,size=length(rad), replace=T)
  boot.sample.force <- sample(force,size=length(force),replace=T)
  beta1rad <- rnorm(n=1,mean=0,sd=b1)
  beta2rad <- rnorm(n=1,mean=0,sd=b2)
  beta1L <- rnorm(n=1,mean=0,sd=b1)
  beta2L <- rnorm(n=1,mean=0,sd=b2)
  beta1force <- rnorm(n=1,mean=0,sd=b3)
  
  rad_s <- mean(boot.sample.rad) + beta1rad+beta2rad
  L_s <- L + beta1L + beta2L
  force_s <- mean(boot.sample.force)+beta1force
  
  boot.r[i] <- 9.81*force_s*L_s/(2*pi*(rad_s)^3)
}

hist(boot.r, main = "Histogram of MC calculation", xlab = "Result")
mean(boot.r)
quantile(boot.r,probs = c((1-P)/2,(1+P)/2))

#trial 4
#calculate mean and sd values

rad <- r4
L <- L4
force <- force4

boot.r <- numeric(M)
for (i in 1:M){
  boot.sample.rad <- sample(rad,size=length(rad), replace=T)
  boot.sample.force <- sample(force,size=length(force),replace=T)
  beta1rad <- rnorm(n=1,mean=0,sd=b1)
  beta2rad <- rnorm(n=1,mean=0,sd=b2)
  beta1L <- rnorm(n=1,mean=0,sd=b1)
  beta2L <- rnorm(n=1,mean=0,sd=b2)
  beta1force <- rnorm(n=1,mean=0,sd=b3)
  
  rad_s <- mean(boot.sample.rad) + beta1rad+beta2rad
  L_s <- L + beta1L + beta2L
  force_s <- mean(boot.sample.force)+beta1force
  
  boot.r[i] <- 9.81*force_s*L_s/(2*pi*(rad_s)^3)
}

hist(boot.r, main = "Histogram of MC calculation", xlab = "Result")
mean(boot.r)
quantile(boot.r,probs = c((1-P)/2,(1+P)/2))
