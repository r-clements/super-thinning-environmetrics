#AUTHOR: Robert Clements
#DATE: 12/11/2011
#SUMMARY: Plots for the ETAS examples in the paper "Evaluation
#         of space-time point process models using super-thinning"

library(maps); library(spatstat)

#shen defines the spatial window that will be plotted for all 3 ETAS models
shen <- read.table("ETAS/data/bg_rate_shen.dat", h=T)

plot.function <- function(resids, data, l.data, r.data, ci1) {
  data[(data$rate >= .16), 9] <- .16
  mx1 <- min(data[,1]); mx2 <- max(data[,1])
  my1 <- min(data[,3]); my2 <- max(data[,3])
  xlines <- seq(mx1, mx2, by=.1)
  ylines <- seq(my1, my2, by=.1)
  
  x.g <- (data$maxlon + data$minlon)/2
  y.g <- (data$maxlat + data$minlat)/2
  nxc <- seq(min(x.g), max(x.g), .1)
  nyc <- seq(min(y.g), max(y.g), .1)
  nxg <- rep(nxc, length(nyc))
  nyg <- c()
  for(i in 1:length(nyc)){
    nyg <- c(nyg, rep(nyc[i], length(nxc)))
  }
  r2 <- rep(0, length(nxg))
  space.n <- data.frame(cbind(x.g, y.g, data$rate))
  space2 <- data.frame(cbind(nxg, nyg, r2))
  space.n[,4] <- paste(space.n[,1], space.n[,2], sep=" ")
  space2[,4] <- paste(space2[,1], space2[,2], sep=" ")
  space2[match(space.n[,4], space2[,4]),3] <- space.n[,3]
  
  immat.data <- matrix(space2$r2, byrow=F, nrow=length(nxc))
  
  temp <- data
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
  
  cuts <- seq(0, .16, length.out=66)
  
  #PLOT WITH THE L-FUNCTIONS
  ci1 <- ci1[, -c(128:513)]
  l.data <- l.data[-c(128:513), ]
  r.data <- r.data[-c(128:513), ]
  r.data <- r.data*cos(pi*7/36)*110.85
  ylim1 <- c(min(min(ci1[1,]), min(l.data)), max(max(ci1[2,]), max(l.data)))
  
  layout(matrix(c(1,2,3), ncol=3), widths=c(4,.3,3))
  par(mar=c(4,4,2,.1), bty="n")
  image(nxc, nyc, immat.data, breaks=cuts, col=rev(gray((0:64)/64)), 
        xlim = c(min(shen$minlon)-.2, max(shen$maxlon)+1), 
        ylim = c(min(shen$minlat)-.2, max(shen$maxlat)+.2), xlab = "Longitude", ylab = "Latitude")
  map("state", "ca", add=T, col=gray(0.45))
  points(resids[which(resids$type==1),1], resids[which(resids$type==1),2])
  points(resids[which(resids$type==2),1], resids[which(resids$type==2),2], pch=3, col = grey(.4))
  lines(space[,1], space[,2])
  lines(space[c(1,nrow(space)),1], space[c(1,1),2])
  
  par(mar=c(4,.1,2,1))
  key	<-	(0:64)/64	
  plot(NULL, ylim=c(-.2,1.2), xlim=c(-2, 5), type="n", axes=F, xlab="", ylab="", main="")
  image(-2:0, key, t(matrix(rep(key,2), nrow=65, byrow=F)), add=T, col=rev(gray((0:64)/64)))
  text(c(.3,.3,.3,.3,.3), c(0,.25,.5,.75,1), c(0,.04,.08,.12,paste(.16, "+")), adj=0, cex=.75)
  
  par(mar=c(9,5,7,1))
  plot(r.data, l.data, ylim = ylim1, type="l", xlab = "km", ylab = "Lw(r)-r")
  lines(r.data, ci1[1,], lty=2)
  lines(r.data, ci1[2,], lty=2)
}

#PLOT HELMSTETTER RESIDS WITH BG RATE
#####################################
#####################################
resids <- read.table("ETAS/data/helm_resids_1500.txt", h=T)
helm <- read.table("ETAS/data/bg_rate.dat", h=T)
l.helm <- read.table("ETAS/data/l.helm.txt", h=F)
r.helm <- read.table("ETAS/data/r.helm.txt", h=F)
ci1 <- read.table("ETAS/data/ci1_helm.txt", h=F)

plot.function(resids, helm, l.helm, r.helm, ci1)
#############################
#############################

#PLOT KAGAN RESIDS WITH BG RATE
###############################
###############################

resids <- read.table("ETAS/data/kagan_resids_1000.txt", h=T)
kagan <- read.table("ETAS/data/bg_rate_kagan.dat", h=T)
l.kagan <- read.table("ETAS/data/l.kagan.txt", h=F)
r.kagan <- read.table("ETAS/data/r.kagan.txt", h=F)
ci1 <- read.table("ETAS/data/ci1_kagan.txt", h=F)

plot.function(resids, kagan, l.kagan, r.kagan, ci1)
##############################
##############################

#PLOT SHEN RESIDS WITH BG RATE
##############################
##############################

resids.shen <- read.table("ETAS/data/shen_resids_1000.txt", h=T)
shen <- read.table("ETAS/data/bg_rate_shen.dat", h=T)
l.shen <- read.table("ETAS/data/l.shen.txt", h=F)
r.shen <- read.table("ETAS/data/r.shen.txt", h=F)
ci1 <- read.table("ETAS/data/ci1_shen.txt", h=F)

plot.function(resids, shen, l.shen, r.shen, ci1)
#############################
#############################