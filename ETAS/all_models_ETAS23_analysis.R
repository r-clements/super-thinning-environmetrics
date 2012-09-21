#AUTHOR: Robert Clements
#DATE: 12/11/2011
#SUMMARY: Perform super-thinning on ETAS model example Helm.ETAS in the paper "Evaluation
#         of space-time point process models using super-thinning"
library(splancs); library(spatstat)

#the conditional intensity function
lamb <- function(x, y, t, hm, mu, theta)
{
  K0 = theta[1]; c = theta[2];
  p = theta[3]; a = theta[4]; d = theta[5]; q = theta[6]
  
  lam1 = rep(mu[1], length(x))
  if(length(lam1) == 1) 
    return(lam1)
  for(i in 2:length(x)){
    r2  = (x[i]-x[1:(i-1)])^2+(y[i]-y[1:(i-1)])^2 # not using Great Circle Distance
    lam1[i] = mu[i] + K0 * sum((t[i]-t[1:(i-1)]+c)^(-p)*(r2*exp(-a*hm[1:(i-1)])+d)^-q)
  }
  return(lam1)
}

super.thin <- function(z, bg.rate, mu, theta2, k) {
  mu <- c(mu[,1])
  mu.total <- bg.rate$rate
  
  n = nrow(z)
  t = z$dt
  x = z$Longitude
  y = z$Latitude
  hm = z$Magnitude - 3.05  #SUBTRACT M0
  
  lambda <- lamb(x, y, t, hm, mu, theta2)
  
  #PERFORM SUPER-THINNING
  #######################
  
  #create the observation region
  temp <- bg.rate
  rows <- length(unique(temp[, 4]))
  
  x.left <- c(); x.right <- c(); y.left <- c(); y.right <- c()
  for(i in 1:rows) {
    mtop <- max(temp[,4])
    row <- temp[which(temp[ ,4] == mtop), ]
    xl <- min(row[,1])
    xr <- max(row[,2])
    x.left <- c(x.left, rep(xl, 2))
    x.right <- c(x.right, rep(xr, 2))
    y.left <- c(y.left, row[1,4], row[1,3])
    y.right <- y.left
    temp <- temp[which(temp[ ,4] != mtop), ]
  }
  
  x.right <- rev(x.right)
  y.right <- rev(y.right)
  
  top <- cbind(x.left, y.left)
  bottom <- cbind(x.right, y.right)
  
  space <- rbind(top, bottom)
  space <- unique(space)	
  window<-owin(poly=list(x=space[,1], y=space[,2]))
  
  ###################
  ###################
  
  #save as z2 
  z2 <- z[,c(1,2,3,4)]
  z2 <- cbind(z2, lambda)
  
  #thinning
  keep <- z2[which(z2$lambda <= k), ]
  thin.data <- z2[which(z2$lambda > k), ]
  prob <- runif(nrow(thin.data))
  retain <- (prob < k/thin.data$lambda)
  keep2 <- thin.data[retain, ]
  
  #superposition
  num <- k*(365*5-12)
  sim.pts <- rpoispp(num, win = window)
  x.s <- sim.pts$x; y.s <- sim.pts$y; t.s <- runif(length(x.s), 0, (365*5-12)); hm.s <- runif(length(x.s), 0, 4.95)
  S <- data.frame(cbind(x.s, y.s, t.s, hm.s))
  S <- S[order(t.s), ]
  
  #extract mu for each point in S
  mu.s <- c(); error <- c()
  for(i in 1:nrow(S)){
    place <- which((bg.rate$minlon <= S[i,1])&(bg.rate$maxlon > S[i,1])&(bg.rate$minlat <= S[i,2])&(bg.rate$maxlat > S[i,2]))
    mu.s <- c(mu.s, bg.rate[place, 9])
    if(length(mu.s) < i) 
    {
      mu.s <- c(mu.s, i)
      error <- c(error, i)}
  }
  
  if (length(error) > 0){
    S <- data.frame(cbind(S[-error, ], mu.s[-error]))
  } else S <- data.frame(cbind(S, mu.s))
  names(S) <- names(z2)
  lambda.s <- c()
  
  for(i in 1:nrow(S)){
    temp <- rbind(z2, S[i,])
    temp <- temp[which(temp$dt <= S[i, 3]), ]
    x.s <- temp[,1]; y.s <- temp[,2]; t.s <- temp[,3]; hm.s <- temp[,4]; mu.s <- temp[,5]
    t.l <- lamb(x.s, y.s, t.s, hm.s, mu.s, theta2)
    lambda.s <- c(lambda.s, tail(t.l, 1))
  }
  probs <- runif(nrow(S))
  retain <- (probs < (k - lambda.s)/lambda.s)
  sim.pts <- S[retain, ]
  
  type <- c(rep(1, nrow(keep)+nrow(keep2)), rep(2, nrow(sim.pts)))
  
  #final residuals
  resids <- data.frame(cbind(rbind(keep, keep2, sim.pts), type))
  return(resids)
}

#THE SUPER-THINNING RATES, k, FOLLOW. PICK ONE VALUE OF k, AND THEN PROCEED TO THE PROCEDURE BELOW
#--- vol IS SPECIFIC TO EACH 
#DATASET, SO FIRST DEFINE vol BELOW, AND THEN CHOOSE k
###################
###################
#USE THIS FOR THINNING AND SUPERPOSING SAME NUMBER OF POINTS 
k <- nrow(z)/2/vol 
#USE THIS FOR PURE THINNING
k <- min(lambda)
#USE THIS FOR 800 POINTS
k <- 800/vol
#USE THIS FOR 400 POINTS
k <- 400/vol
#USE THIS FOR 1500 POINTS
k <- 1500/vol

#READ IN THE DATA AND RUN SUPER-THINNING
#################
#HELM.ETAS
z <- read.table("ETAS/data/eqdata305_revised.dat", h=T) #the earthquake catalog
bg.rate <- read.table("ETAS/data/bg_rate.dat", h=T) #the background rate 
mu <- read.table("ETAS/data/eqcat_helm_mu_summed.txt", h=F) #the Helm.ETAS background rate
#estimates of theta obtained from EM-type algorithm
theta2 <- c(3.047139e-05, 3.737948e-02, 1.304083e+00, 9.097675e-01, 1.375795e-04, 1.643300e+00)
#volume of spatial-temporal region
vol <- nrow(bg.rate)*.1*.1*(5*365-12)
helm <- super.thin(z, bg.rate, mu, theta2, k)


#SHEN.ETAS
z <- read.table("ETAS/data/eqdata_shen_revised.dat", h=T) #the earthquake catalog
bg.rate <- read.table("ETAS/data/bg_rate_shen.dat", h=T) #the background rate 
mu <- read.table("ETAS/data/eqcat_shen_mu_summed.txt", h=F) #the Helm.ETAS background rate
#estimates of theta obtained from EM-type algorithm
theta2 <- c(2.795276e-05, 3.931234e-02, 1.308017e+00, 9.382453e-01, 2.240144e-04, 1.722722e+00)
#volume of spatial-temporal region
vol <- nrow(bg.rate)*.1*.1*(5*365-12)
shen <- super.thin(z, bg.rate, mu, theta2, k)


#KAGAN.ETAS
z <- read.table("ETAS/data/eqdata_kag_revised.dat", h=T) #the earthquake catalog
bg.rate <- read.table("ETAS/data/bg_rate_kagan.dat", h=T) #the background rate 
mu <- read.table("ETAS/data/eqcat_kagan_mu_summed.txt", h=F) #the Helm.ETAS background rate
#estimates of theta obtained from EM-type algorithm
theta2 <- c(2.247393e-05, 3.746132e-02, 1.299134e+00, 9.397220e-01, 2.212194e-04, 1.750720e+00)
#volume of spatial-temporal region
vol <- nrow(bg.rate)*.1*.1*(5*365-12)
kagan <- super.thin(z, bg.rate, mu, theta2, k)


##################
##################