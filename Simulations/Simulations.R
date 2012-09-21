#AUTHOR: Robert Clements
#DATE: 12/11/2011
#SUMMARY: Simulates several point processes and runs super-thinned residual analysis using the 'true' models, i.e. 
#         the models that were used to simulate the points. These are the simulations used in the paper "Evaluation
#         of space-time point process models using super-thinning"

#required packages
library(spatstat); library(stppResid)


#SIMULATION 1 - corresponding to Figure 1.
#########################################

x1<-rpoispp(80, win=owin(c(0,1),c(0,1)))
x2<-rpoispp(80, win=owin(c(1,2),c(0,1)))
x3<-rpoispp(20, win=owin(c(0,1), c(1,2)))
x4<-rpoispp(80, win=owin(c(1,2),c(1,2)))

x <- c(x1$x, x2$x, x3$x, x4$x)
y <- c(x1$y, x2$y, x3$y, x4$y)
t <- runif(length(x))

X <- stpp(x, y, t, stw = stwin(xcoord = c(0,2), ycoord = c(0,2), tcoord = c(0,1)))

#the conditional intensity function definition
hom1 <- function(X) {
	x <- X$x; y <- X$y; t <- X$t
	q1 <- which(y <= 1)
	q2 <- which((x > 1) & (y > 1))
	q3 <- which((x <= 1) & (y > 1))
	lam <- rep(0, length(x))
	lam[c(q1, q2)] <- 80
	lam[q3] <- 20
	lam
} 


#SIMULATION 2 - corresponding to Figure 2.
##########################################

y1<-rpoispp(20, win=owin(c(0,1),c(0,1)))
y2<-rpoispp(20, win=owin(c(1,2),c(0,1)))
y3<-rpoispp(80, win=owin(c(0,1), c(1,2)))
y4<-rpoispp(20, win=owin(c(1,2),c(1,2)))

x <- c(y1$x, y2$x, y3$x, y4$x)
y <- c(y1$y, y2$y, y3$y, y4$y)
t <- runif(length(x))

X2 <- stpp(x, y, t, stw = stwin(xcoord = c(0,2), ycoord = c(0,2), tcoord = c(0,1)))

#the conditional intensity function
hom2 <- function(X) {
	x <- X$x; y <- X$y; t <- X$t
	q1 <- which(y <= 1)
	q2 <- which((x > 1) & (y > 1))
	q3 <- which((x <= 1) & (y > 1))
	lam <- rep(0, length(x))
	lam[c(q1, q2)] <- 20
	lam[q3] <- 80
	lam
} 


#SIMULATION 3 - corresponding to Figures 3 and 4.
#################################################

#the function that we will be simulating from
fx<-function(x,y){3000*exp(-3*x - 4*y)}
fxpp<-rpoispp(fx)

x <- fxpp$x
y <- fxpp$y
t <- runif(length(x))

X3 <- stpp(x, y, t, stw = stwin(xcoord = c(0,1), ycoord = c(0,1), tcoord = c(0,1)))

#the conditional intensity function properly defined
inhom1 <- function(X) {
	x <- X$x; y <- X$y; t <- X$t
	lam <- 3000*exp(-3*x - 4*y)
	lam
}

#SUPER-THINNED RESIDUALS
########################

stX1 <- superthin(X, hom1, k = 80) #pure superposition 
stX2 <- superthin(X2, hom2, k=20)  #pure thinning 
stX3.1 <- superthin(X3, inhom1, k=2.735646) #thinning 
stX3.2 <- superthin(X3, inhom1, k=3000)     #superposition 
stX3.3 <- superthin(X3, inhom1, k=233.2023) #super-thinning 

plot(stX1) # --- Figure 2
plot(stX2) # --- Figure 1
plot(stX3.1) # --- Figure 3c
plot(stX3.2) # --- Figure 3d
plot(stX3.3) # --- Figure 4a

#INHOMOGENEOUS CENTERED L-FUNCTIONS FOR SIMULATIONS
###################################################

#true conditional intensities for each simulation
lamb1 <- hom1(X)
lamb2 <- hom2(X2)
lamb3 <- inhom1(X3)

#set-up spatial windows
W1 <- owin(c(0,2), c(0,2))
W2 <- owin(c(0,1), c(0,1))

#create ppp objects using 'spatstat' package
pp1 <- as.ppp(data.frame(cbind(X$x, X$y)), W1)
pp2 <- as.ppp(data.frame(cbind(X2$x, X2$y)), W1)
pp3 <- as.ppp(data.frame(cbind(X3$x, X3$y)), W2)

#estimate the centered inhomogeneous L-function for each observed point pattern
l1 <- Lest(pp1, lambda=lamb1, correction = "Ripley")
l2 <- Lest(pp2, lambda=lamb2, correction = "Ripley")
l3 <- Lest(pp3, lambda=lamb3, correction = "Ripley")

r1 <- l1$r; r2 <- l2$r; r3 <- l3$r

l1f <- l1$iso - r1
l2f <- l2$iso - r2
l3f <- l3$iso - r3

#INHOMOGENEOUS CENTERED L-FUNCTIONS FOR RESIDUALS
#################################################

#conditional intensities of the residual point process
lamb1 <- rep(80, nrow(stX1$residuals))
lamb2 <- rep(20, nrow(stX2$residuals))
lamb3 <- rep(233.2023, nrow(stX3.3$residuals))

W1 <- owin(c(0,2), c(0,2))
W2 <- owin(c(0,1), c(0,1))

pp1 <- as.ppp(stX1$residuals, W1)
pp2 <- as.ppp(stX2$residuals, W1)
pp3 <- as.ppp(stX3.3$residuals, W2)

l1 <- Linhom(pp1, lambda=lamb1, correction = "Ripley")
l2 <- Linhom(pp2, lambda=lamb2, correction = "Ripley")
l3 <- Linhom(pp3, lambda=lamb3, correction = "Ripley")

r1 <- l1$r; r2 <- l2$r; r3 <- l3$r

l1fres <- l1$iso - r1
l2fres <- l2$iso - r2
l3fres <- l3$iso - r3

#INHOM CENTERED L-FUNCTION CONFIDENCE BOUNDS
############################################

#simulated 95% confidence bounds based on 5000 simulated homogeneous Poisson processes
lhf1 <- lhf2 <- lhf3 <- lgf1 <- lgf2 <- data.frame()
for(i in 1:5000) {
	h1 <- rpoispp(80, win = W1)
	g1 <- rpoispp(65, win = W1)
	h2 <- rpoispp(20, win = W1)
	g2 <- rpoispp(35, win = W1)
	h3 <- rpoispp(233.2023, win = W2)
	
	
	lambdah.1 <- rep(80, h1$n)
	lambdag.1 <- rep(65, g1$n)
	lambdah.2 <- rep(20, h2$n)
	lambdag.2 <- rep(35, g2$n)
	lambdah.3 <- rep(233.2023, h3$n)
	
	lh1 <- Linhom(h1, lambdah.1, correction = "Ripley")
	lg1 <- Linhom(g1, lambdag.1, correction = "Ripley")
	lh2 <- Linhom(h2, lambdah.2, correction = "Ripley")
	lg2 <- Linhom(g2, lambdag.2, correction = "Ripley")
	lh3 <- Linhom(h3, lambdah.3, correction = "Ripley")
	
	lhf1 <- rbind(lhf1, lh1$iso - r1)
	lgf1 <- rbind(lgf1, lg1$iso - r1)
	lhf2 <- rbind(lhf2, lh2$iso - r2)
	lgf2 <- rbind(lgf2, lg2$iso - r2)
	lhf3 <- rbind(lhf3, lh3$iso - r3)
}

Quantile<-function(x, probs=c(.025, .975), ...){
   	       quantile(x, probs, na.rm=TRUE, ...)
   	       }
   	     
ci1 <- apply(as.matrix(lgf1), 2, Quantile)
ci.res1 <- apply(as.matrix(lhf1), 2, Quantile)
ci2 <- apply(as.matrix(lgf2), 2, Quantile)
ci.res2 <- apply(as.matrix(lhf2), 2, Quantile)
ci3 <- apply(as.matrix(lhf3), 2, Quantile)
ci.res3 <- apply(as.matrix(lhf3), 2, Quantile)

#PLOT L-FUNCTIONS FOR SIMULATIONS
#################################

ylim1 <- c(min(min(ci1[1,]), min(l1f)), max(max(ci1[2,]), max(l1f)))
ylim2 <- c(min(min(ci2[1,]), min(l2f)), max(max(ci2[2,]), max(l2f)))
ylim3 <- c(min(min(ci3[1,]), min(l3f)), max(max(ci3[2,]), max(l3f)))

plot(r1, l1f, ylim = ylim1, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r1, ci1[1,], lty=2)
lines(r1, ci1[2,], lty=2)

plot(r2, l2f, ylim = ylim2, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r2, ci2[1,], lty=2)
lines(r2, ci2[2,], lty=2)

plot(r3, l3f, ylim = ylim3, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r3, ci3[1,], lty=2)
lines(r3, ci3[2,], lty=2)


#PLOT L-FUNCTIONS FOR RESIDUALS
###############################

ylim4 <- c(min(min(ci.res1[1,]), min(l1fres)), max(max(ci.res1[2,]), max(l1fres)))
ylim5 <- c(min(min(ci.res2[1,]), min(l2fres)), max(max(ci.res2[2,]), max(l2fres)))
ylim6 <- c(min(min(ci.res3[1,]), min(l3fres)), max(max(ci.res3[2,]), max(l3fres)))

plot(r1, l1fres, ylim = ylim4, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r1, ci.res1[1,], lty=2)
lines(r1, ci.res1[2,], lty=2)

plot(r2, l2fres, ylim = ylim5, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r2, ci.res2[1,], lty=2)
lines(r2, ci.res2[2,], lty=2)

plot(r3, l3fres, ylim = ylim6, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r3, ci.res3[1,], lty=2)
lines(r3, ci.res3[2,], lty=2)

#PLOTS FOR PAPER
################

# --- Figure 1
par(mfrow=c(2,2))
plot(X)
plot(r1, l1f, ylim = ylim1, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r1, ci1[1,], lty=2)
lines(r1, ci1[2,], lty=2)
plot(stX1)
plot(r1, l1fres, ylim = ylim4, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r1, ci.res1[1,], lty=2)
lines(r1, ci.res1[2,], lty=2)

# --- Figure 2
par(mfrow=c(2,2))
plot(X2)
plot(r2, l2f, ylim = ylim2, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r2, ci2[1,], lty=2)
lines(r2, ci2[2,], lty=2)
plot(stX2)
plot(r2, l2fres, ylim = ylim5, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r2, ci.res2[1,], lty=2)
lines(r2, ci.res2[2,], lty=2)

# --- Figure 3
par(mfrow=c(2,2))
plot(X3)
plot(r3, l3f, ylim = ylim3, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r3, ci3[1,], lty=2)
lines(r3, ci3[2,], lty=2)
plot(stX3.1)
plot(stX3.2)

# --- Figure 4
par(mfrow=c(1,2))
plot(stX3.3)
plot(r3, l3fres, ylim = ylim6, type="l", xlab = "r", ylab = "Lw(r)-r")
lines(r3, ci.res3[1,], lty=2)
lines(r3, ci.res3[2,], lty=2)