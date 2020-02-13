# Rayne Milner
# MAE 2017 Final Exam





########################## Q1 #######################################
P1_data <- read.csv("mae_207_final_exam_f_2019_problem_1_data.csv") #imports the data for the first problem
# Problem 3.7
x <- P1_data$T1
# Chauvenet's Criterion
N <- length(x) 							# number of data points
z <- (x-mean(x))/sd(x) 					# calculate the data z values
p_chauv <- 1-1/(2*N) 					# calculate the Chauvenet probability band value
z_max <- qnorm(0.5+0.5*p_chauv) 		# find the maximum allowable z value
z_min <- qnorm(0.5-0.5*p_chauv) 		# find the minimum allowable z value
x_max <- mean(x) + z_max*sd(x) 			# find the maximum allowable data value
x_min <- mean(x) + z_min*sd(x) 			# find the minimum allowable data value
plot(x, ylim = c(min(x)-1,max(x)+1), main = "Outlier Identification", ylab = "T1 Data")			# plot the data
abline(h=x_min,col="red",lw=2)			# mark the minimum allowable data value
abline(h=x_max,col="red",lw=2)			# mark the maximum allowable data value

outliers <- x[which(abs(z)>z_max)] #reports the T1 values found to be outliers
#data rows 10(T1 = 312.1) and 32(T2 = 308.2) are the rows found to contain outliers

# (a) Outlier Identification
x <- P1_data$T2
# Chauvenet's Criterion
N <- length(x) 							# number of data points
z1 <- (x-mean(x))/sd(x) 					# calculate the data z values
p_chauv <- 1-1/(2*N) 					# calculate the Chauvenet probability band value
z_max <- qnorm(0.5+0.5*p_chauv) 		# find the maximum allowable z value
z_min <- qnorm(0.5-0.5*p_chauv) 		# find the minimum allowable z value
x_max <- mean(x) + z_max*sd(x) 			# find the maximum allowable data value
x_min <- mean(x) + z_min*sd(x) 			# find the minimum allowable data value
plot(x, ylim = c(min(x)-1,max(x)+1), main = "Outlier Identification", ylab = "T2 Data")			# plot the data
abline(h=x_min,col="red",lw=2)			# mark the minimum allowable data value
abline(h=x_max,col="red",lw=2)			# mark the maximum allowable data value

outliers <- x[which(abs(z1)>z_max)] #reports the T2 values found to be outliers
#data point 10 and 32 are the points found to be outliers

P1_data <- P1_data[which((abs(z)<z_max) & (abs(z1)<z_max)),] #removes all rows which contained outliers in both T1 and T2 data
# Rows 10, 18, and 32 were deleted.

# (b) Estimate Mean and SD of Theta using MC approach
T1 <- P1_data$T1
T2 <- P1_data$T2
Q <- P1_data$Q

N <- length(T1)			# number of data points
M <- 1e5				# number of MC replicates
b2 <- 0.1

z.dis <- numeric(length(M))

for (i in 1:M) {
  beta1 <- rnorm(n=1,mean=0,sd=0.5) #systematic error for T1
  betadigT1 <- rnorm(n=1, mean=0, sd = 0.01) #systematic digitization error for T1
  betadigT2 <- rnorm(n=1, mean=0, sd = 0.01) #systematic digitization error for T2
  betadigQ <- rnorm(n=1, mean = 0, sd = 0.02) #systematic digitization error for Q
  beta2 <- rnorm(n=1,mean=0,sd=1) #systematic error for T2
  beta3 <- rnorm(n=1,mean=0,sd=0.5) #systematic for Q
  
  T1.dis <- sample(T1, size = N, replace = TRUE) + beta1 + betadigT1
  T2.dis <- sample(T2, size = N, replace = TRUE) + beta2 + betadigT2
  Q.dis <- sample(Q, size = N, replace = TRUE) + beta3 + betadigQ
  
  z.dis[i] <- mean((T1.dis - T2.dis)/(Q.dis))
}
mean(z.dis) #= 0.490
sd(z.dis) #= 0.0579

# (c) Combined PDF's from multiple experiments
dx <- 0.000001
x <- seq(0,1, by = dx)
mua <- 0.551 #K/W, mean of pdf a
mub <- 0.49 #K/W, mean of pdf b
siga <- 0.074 #K/W, sig of pdf a
sigb <- 0.059 #K/W, sig b of pdf b

pdfa <- dnorm(x, mean = 0.551, sd = 0.074) #pdf for the distribution a
pdfb <- dnorm(x, mean = 0.490, sd = 0.059) #pdf for the distribution b

pdfc_unnorm <- pdfa*pdfb #unnormalized pdf 

#normailzation constant for pdf c integral
c <- 1/sum(pdfc_unnorm*dx) # = 0.292

pdfc <- c*pdfc_unnorm #normalized pdf

#plots the three pdfs
plot(x, pdfc, main = "Combined PDFs", xlab = "Theta", ylab = "PDF", type = "l")
lines(x,pdfb, col = "blue")
lines(x,pdfa, col = "red")
legend("topleft",
       c("pdf C","pdf B", "pdf A"),
       fill = c("black","blue", "red"))

# Calculate Mean and SD of combined pdf.
mean_pdfc <- sum(pdfc*x*dx) # mean of this combined pdf= 0.514
sd_pdfc <- sqrt(sum((x-mean_pdfc)^2*pdfc*dx)) #sd of this combined pdf = 0.0461

########################## Q2 #######################################
P2_data <- read.csv("mae_207_final_exam_f_2019_problem_2_data.csv") #imports the data for the secondproblem
#Balance Check
x1 <- P2_data$Pipe_1_Flow_Rate #kg/s, flow in
x2 <- P2_data$Pipe_2_Flow_Rate #kg/s, flow out 1
x3 <- P2_data$Pipe_3_Flow_Rate #kg/s, flow out 2

bins <- 0.02 #kg/s, error from installation
bcal <- 0.04 #kg/s, error from calibration, 

r <- x1 - (x2 + x3) #DRE Sum(flow in) - Sum(flow out)

rbar <- mean(r) #mean value of r, = 0.06474226

N <- length(r)

sr <- 1/(N-1)*sum((r-rbar)^2) #calculating sr via "direct approach"

#sensativity values
Theta1 <- 1 #dr/dx1
Theta2 <- -1 #dr/dx2
Theta3 <- -1 #dr/dx3

#b values
#linearity is not-common, calibration is common
b11 <- (bins/1.96)^2 + (bcal/1.96)^2 #= 0.0005206164
b22 <- b11
b33 <- b11

b12 <- (bcal/1.96)^2 #=0.0004164931
b13 <- b12
b23 <- b12

#calculate br
brsqr <- Theta1^2*b11 + Theta2^2*b22 + Theta3^2*b33 + 2*Theta1*Theta2*b12 + 2*Theta1*Theta3*b13 + 2*Theta2*Theta3*b23 # =0.000728863

#calculate expanded uncertainty assuming coverage factor is two
Urbar <- 2*sqrt((sr^2)/N + brsqr) # =0.05399539

if (abs(rbar) > Urbar){ #0.06474226 > 0.05399539, TRUE
  print("Balance Check has failed")
}else{
  print("Balance Check has passed, DRE o.k.")}

 #The balance check fails
########################## Q3 #######################################
# pt a completed in Latex

#### b
set.seed(10)
x <- rexp(n=20, rate = 0.2) 

#### c
M <- 1e5
N <- length(x)

l.dis <- numeric(length(M))

for (i in 1:M) {
  x <- sample(x, size = N, replace = TRUE)
  
  l.dis[i] <- 1/mean(x) #lambda = 1/mean(x)
}

mu_lambda1 <- mean(l.dis) #mu_lambda=0.152

### d

xbar <- mean(x) #=4.96

mu_lambda2 <- N^2*xbar #2633.1


