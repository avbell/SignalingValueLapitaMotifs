######################
# Supplementary Information: Founding Motifs Maintain Signaling Roles Three Thousand Years Among the Descendent Population
# Adrian Viliami Bell and Gracyn Hinojosa
######################

rm(list = ls())

# Set working directory housing required files
# CodeBell2019C.R
# TongaUS-TriadData-ForPub.rdata
# images m1.jpeg through m8.jpeg
# AnalysisProgress.rdata (for saved estimations)

wkdir <- "mydirectory"
setwd(wkdir)
# Packages used in plotting
library(jpeg)
library(igraph)
library(visNetwork)

## ---- Load dataset ----
load( "TongaUS-TriadData-ForPub.rdata" )

# number of objects (motifs)
nObj <- max(d.my[[1]]$motif_tricomb[[1]], na.rm = TRUE)

## ---- Similarity matrix builder (dynamic N) ----
# triadCount2(): counts how many times each object was selected within each triad type
triadCount2 <- function(y, triad) {
  # triad: nTri x 3 matrix of motif IDs
  # y:     nTri x >=4 matrix; columns 2:4 are one-hot
  stopifnot(is.matrix(triad), ncol(triad) == 3)
  stopifnot(is.matrix(y), ncol(y) >= 4)
  
  combin <- t(combn(1:nObj, 3)) # all possible triads for nObj motifs
  ntri   <- nrow(combin)
  
  count <- matrix(0, nrow = ntri, ncol = 3)
  
  for (i in 1:nrow(triad)) {
    pr <- triad[i, ]
    
    # Find which row in combin matches this triad (ignoring order)
    wrow <- which(apply(combin, 1, function(r) setequal(r, pr)))
    if (length(wrow) != 1) next
    
    # Identify which object was "chosen" by the one-hot outcome
    chosen_obj <- triad[i, which(y[i, 2:4] == 1)]
    add.index  <- which(combin[wrow, ] == chosen_obj)
    
    if (length(add.index) == 1) {
      count[wrow, add.index] <- count[wrow, add.index] + 1
    }
  }
  
  count
}

# simatrix2(): converts triad choice counts into an nObj x nObj similarity matrix
simatrix2 <- function(d) {
	  
  combin <- t(combn(1:nObj, 3))
  ans <- matrix(0, nrow = nrow(combin), ncol = 3)
  
  # Sum triad counts across participants
  for (i in 1:length(d)) {
    if (!is.null(d[[i]]$motif_tricomb[[1]])) {
      triad <- d[[i]]$motif_tricomb[[1]]
      y     <- d[[i]]$motif_triad_data[[1]]
    } else {
      triad <- d[[i]]$tricomb
      y     <- d[[i]]$triad_data
    }
    ans <- ans + triadCount2(y = y, triad = triad)
  }
  
  # Convert counts (per triad) into a symmetric similarity matrix
  ans2 <- matrix(0, nrow = nObj, ncol = nObj)
  
  for (i in 1:nrow(combin)) {
    ind <- combin[i, ]
    ch  <- ans[i, ]  # counts corresponding to which item "won" in that triad
    
    # If object ind[k] was chosen, increment similarity of the *other two*.
    ans2[ind[1], ind[2]] <- ans2[ind[1], ind[2]] + ch[3]
    ans2[ind[2], ind[1]] <- ans2[ind[2], ind[1]] + ch[3]
    
    ans2[ind[1], ind[3]] <- ans2[ind[1], ind[3]] + ch[2]
    ans2[ind[3], ind[1]] <- ans2[ind[3], ind[1]] + ch[2]
    
    ans2[ind[2], ind[3]] <- ans2[ind[2], ind[3]] + ch[1]
    ans2[ind[3], ind[2]] <- ans2[ind[3], ind[2]] + ch[1]
  }
  
  ans2
}

# Compute similarity matrices, normalized
# number of possible times two motifs can be paired is 6
# combn(8,3) # check here
d.my.sim   <- simatrix2(d.my)   / (6 * length(d.my))
d.t23.sim  <- simatrix2(d.t23m) / (6 * length(d.t23m))   
# verify matrix dimensions
dim(d.my.sim); dim(d.t23.sim)

## ---- Dendrograms ----
motif.dendro_icons <- function(m, title=NULL, type="m", show_icons=TRUE) {
  
  imageDir <- file.path(wkdir)
  stopifnot(nObj == ncol(m))
  
  # scale to [0,1] and build distance
  m2 <- m / max(m, na.rm = TRUE)
  rownames(m2) <- as.character(1:nObj)
  colnames(m2) <- as.character(1:nObj)
  
  hc <- hclust(as.dist(1 - m2))
    
  dev.new()
  par(mar = c(0,0,0,4), xpd = NA)
  plot(as.dendrogram(hc), horiz = TRUE, axes = FALSE, dLeaf = -0.25)
  if (!is.null(title)) mtext(text = title, side = 3, line = -1, xpd = TRUE)
  
  if (show_icons) {
    hts <- order.dendrogram(as.dendrogram(hc))
    for (i in 1:nObj) {
      f <- file.path(imageDir, paste0(type, hts[i], ".jpeg"))
      if (file.exists(f)) {
        img <- jpeg::readJPEG(f)
        rasterImage(img, xleft=-0.01, xright=-0.11, ybottom=i-0.3, ytop=i+0.3)
      }
    }
  }
 }

# Plot dendrograms
motif.dendro_icons(d.t23.sim, title="Tonga 2023")
motif.dendro_icons(d.my.sim,  title="Utah 2025")

## ---- Network view ----
motif.network_icons <- function(m, title=NULL, type="m") {
  imageDir <- paste0(wkdir,"/")  
  the.icons <- function(z) paste0(imageDir, type, z, ".jpeg")
  
  N <- nrow(m)
  stopifnot(N == ncol(m))
  
  from <- combn(1:N, 2)[1,]
  to   <- combn(1:N, 2)[2,]
  
  edgeweights <- sapply(1:length(from), function(z) m[from[z], to[z]])
  edgelabs <- format(edgeweights, digits = 2)
  
  edges <- data.frame(
    from = from,
    to = to,
    width = edgeweights,
    length = (1 - edgeweights),
    color = "gray",
    label = edgelabs,
    font.size = 6
  )
  
  nodes <- data.frame(
    id = 1:N,
    shape = "image",
    image = the.icons(1:N),
    label = as.character(1:N)  # fallback label if image fails
  )
  
  visNetwork(nodes, edges, width = "100%") %>%
    visNodes(size = 20,        # larger size helps images render, change if too large
             shape = "image",
             brokenImage = the.icons(1)) %>%  # fallback if image breaks
    visIgraphLayout(layout = "layout_with_mds", dist = 1 - m, smooth = TRUE) %>%
    visOptions(highlightNearest = TRUE)
}


# Draw networks
motif.network_icons(d.my.sim, title="Utah 2025 network")
motif.network_icons(d.t23.sim, title="Tonga 2023 network")

# -----------------------------
# EVOLUTIONARY MODEL OF SOCIAL COORDINATION
# TO BE USED IN STATISTICAL MODEL
# -----------------------------
# define inverse logit function
logit <- function(x) exp(x)/(1+exp(x)) 

# define function with recursions
xtl <- function( x12,x13,d1,d2,d3,Bs=1 ){
  x23 <- 1 - x12 - x13
  
  x12p <- 2*x12*x23*(1/2)*logit(Bs*d2) + 2*x12*x13*(1/2)*logit(Bs*d1) + 2*x23*x13*(1-logit(Bs*d3)) + x12^2*(1 - (1/2)*logit(Bs*(d1-d2)) - (1/2)*logit(Bs*(d2-d1)) ) + x23^2*(1/2)*logit(Bs*(d2-d3)) + x13^2*(1/2)*logit(Bs*(d1-d3))  
  
  x13p <-  2*x12*x23*(1 - logit(Bs*d2)) + 2*x12*x13*(1/2)*logit(Bs*d1) + 2*x23*x13*(1/2)*logit(Bs*d3) + x12^2*(1/2)*logit(Bs*(d1-d2)) + x23^2*(1/2)*logit(Bs*(d3-d2)) + x13^2*(1 - (1/2)*logit(Bs*(d1-d3))-(1/2)*logit(Bs*(d3-d1))) 
  
  x23p <- 1 - x12p - x13p
  ans <- list( x12p, x13p, x23p )
  names(ans) <- c("x12", "x13", "x23")
  ans
}

xtl(x12=0.3,x13=0.3,d1=1,d2=1,d3=0.3,Bs=0.1)

# function to run and/or plot evolutionary
# dynamics of association
assoc.sim.l <- function(time=100, x12i=1/3,x13i=1/3,d1=0.2,d2=0.3,d3=0.4, Bs=1/3, plotit=FALSE ){
  
  x12.tm <- rep(NA, time)
  x13.tm <- rep(NA, time)
  x23.tm <- rep(NA, time)
  
  # starting frequencies
  x12.tm[1] <- x12i
  x13.tm[1] <- x13i
  x23.tm[1] <- 1-x12i-x13i	
  
  for ( t in 1:(time-1) ){
    xtp <- xtl( x12=x12.tm[t], x13=x13.tm[t],d1=d1,d2=d2,d3=d3,Bs=Bs )
    x12.tm[t+1] <- xtp$x12
    x13.tm[t+1] <- xtp$x13
    x23.tm[t+1] <- xtp$x23		
  }
  x <- list( x12.tm, x13.tm, x23.tm )
  names(x) <- c("x12", "x13","x23")	
  
  if( plotit==TRUE ){
    par( mar=c(2,3,0.5,0))
    plot(1:time, x12.tm, type="l", axes=FALSE, ylim=c(0,1), ylab="", xlab="", col="blue")
    lines( 1:time, x13.tm, lty=2, col="red" )
    lines( 1:time, x23.tm, lty=3, col="green" )
    legend( "topright", legend=c("x12", "x13", "x23"), lty=c(1,2,3), box.lwd=0, col=c("blue", "red", "green") )
    axis(side=1)
    axis(side=2)		
  }
  x
}

# check model simulation
y <- assoc.sim.l(time=100, x12i=1/3,x13i=1/3,d1=2,d2=8, d3=10, Bs=0.3, plotit=TRUE)
yt <- data.frame(y$x12,y$x13,y$x23)
rowSums(yt) # should be a vector of ones
yt[100,]


# -----------------------------
# Model of inherent coordination from reference population
# -----------------------------
# Train model with data from non-Tongans

# use function triadCount defined above
# to count responses to each triad

inherent.coordination <- function(d){
  combin <- t(combn(1:nObj,3))
  ans <- matrix( rep(0,dim(combin)[1]*dim(combin)[2] ), nrow=dim(combin)[1] )
  for ( i in 1:length(d) ){
    triad <- d[[i]]$motif_tricomb[[1]]
    y <- d[[i]]$motif_triad_data[[1]]
    ans <- ans + triadCount2( y = y, triad = triad ) 
  }
  ans	
}

# the rows of this nTri x 3 matrix gives the 
# frequency of an object chosen to be most 
# different from the other two objects
# based on object IDs found in combin
inher.sim.freq <- inherent.coordination( d.my ) / length( d.my )

################################
# Likelihood model

mtriadProb <- function( y=y, triad=triad, delta, Bs.par, u.par, inher.sim )
{
  # first check conditional probabilities
  # make sense, if not then return log likelihood
  
  # needs data triad (number combination)
  # and y, which is a matrix of zeros and ones
  # x are the parameters to the model
  combin.motif <- t(combn(1:8,3))
  ntri <- dim(combin.motif)[1] # no. of triads
  
  pr12 <- pr13 <- pr23 <- rep(NA, ntri)	
  
  llik <- 0 # set up likelihood	
  for ( i in 1:ntri ){
    pr <- triad[i,]
    # what should the starting frequencies be of the
    # evolutionary model? from simulations above, the 
    # equilibrium is not sensitive to initial frequencies
    # specify near-equal frequencies of the three strategies
    em <- 	assoc.sim.l(time=100, x12=0.333,x13=0.333,d1=delta[pr[1]],d2=delta[pr[2]],d3=delta[pr[3]],Bs=Bs.par, plotit=FALSE)
    emt <- data.frame(em$x12,em$x13,em$x23)
    # check results
    # rowSums(emt)
    
    # -----------
    # trained model for inherent coordination
    isr <- which( sapply( 1:ntri, function(z) any(combin.motif[z,1:3]==pr[1]) & any(combin.motif[z,1:3]==pr[2]) & any(combin.motif[z,1:3]==pr[3] ) )	)	
    
    # prob of pair 1,2
    pr1 <- u.par * emt[100,1] + ( 1 - u.par ) * inher.sim[isr,3]
    # prob of pair 1,3
    pr2 <- u.par * emt[100,2] + ( 1 - u.par ) * inher.sim[isr,2]
    # prob of pair 2,3
    pr3 <- u.par * emt[100,3] + ( 1 - u.par ) * inher.sim[isr,1]
    
    # final log probabilities
    # first and second object grouped
    pr12[i] <- log( pr1 )# * ( 1 - pr2 ) * ( 1 - pr3 ) )
    # first and third object grouped		
    pr13[i] <- log( pr2 )#* ( 1 - pr1 ) * ( 1 - pr3 ) )
    # second and third object grouped		
    pr23[i] <- log( pr3 )#* ( 1 - pr1 ) * ( 1 - pr2 ) )
    
    obj <- which( y[i,2:4] == 1 )
    # probability
    prob <- ifelse( obj==1, pr23[i], ifelse( obj==2, pr13[i], ifelse( obj==3, pr12[i], 0 ) ) )
    llik <- llik + prob # add the log probabilities
  }
  -llik	
}


m1loglik <- function( x ){
  # all values are >zero with
  # u from zero to one
  # so use logit transformation
  
  nObj <- 8 
  int <- t(combn(1:nObj,3)) # triad combinations
  ntri <- dim(int)[1] # no. of triads
  
  # second set of parameters in x, delta,
  # benefits to social coordination on one
  # signal is passed to the evolutionary model
  delta <- exp(x[1:nObj]) 
  
  # parameter scaling payoff differentials
  # to strength of selection
  Bs.par <- exp(x[(nObj+1)])
  
  # parameter weighing the relative weight of
  # social coordination explaining data versus
  # inherent coordination
  u.par <- plogis(x[(nObj+2)])
  
  # number of participants
  # Tongan data: d.t23m
  nObs <- length(d.t23m)
  
  neglogLik <- rep(0,nObs)
  for( i in 1:nObs ){
    neglogLik[i] <- mtriadProb( y=d.t23m[[i]]$triad_data, triad=d.t23m[[i]]$tricomb, delta=delta, Bs.par=Bs.par, u.par=u.par, inher.sim=inher.sim.freq  )
  }
  sum(neglogLik)		
}

# test likelihood function
# initial parameter values into vector  x

delta.inits <- runif( nObj )
Bs.init <- runif( 1 )
u.init <- runif( 1 )
x.inits <- c( delta.inits, Bs.init, u.init )
m1loglik (x.inits)
# very slow calculation

# load likelihood coded in C to speed things up
require(Rcpp)
source("CodeBell2019C.r") 

# -- test code for consistencies --

# show simulation endpoints are equivalent for R and C code
assoc.sim.l(x12i=1/3,x13i=1/3,d1=2,d2=8, d3=10, Bs=0.3, plotit=FALSE)$x12[100]
evoSim(x12eq=numeric(1), x13eq=numeric(1), d1=2, d2=8, d3 = 10, Bs=0.3)$x12eq
# YES!!

# -------
#likelihood
# -------
# prep data for C function of negLogLikelihood

# y is a nTri x ( nObs * 3 ) matrix of observations data
y <- d.t23m[[1]]$triad_data[,2:4]
# combine observations into one matrix
for( i in 2:length(d.t23m) ) y <- cbind( y,d.t23m[[i]]$triad_data[,2:4] )
# stack data into one vector by column
y <- as.integer(y)

# standard triad indices corresponding to inher.sim.freq
combin <- t(combn(1:nObj,3))

# randomized object indices from focal population
tricomb <- d.t23m[[1]]$tricomb
# combine indices across participants into one matrix
for( i in 2:length(d.t23m) ) tricomb <- cbind( tricomb,d.t23m[[i]]$tricomb )
# stack index matrix into one vector by column
tricomb <- as.integer(tricomb)

# test function
funLL( 	x = x.inits, 
		y=y, 
		nObj=as.integer(nObj), 
		combin=as.integer(combin), 
		tricomb=tricomb, 
		inhersim=as.numeric(inher.sim.freq), 
		candLL=numeric(1), 
		N=as.integer(length(d.t23m)), 
		nTri=as.integer(dim(combin)[1]), 
		debug=numeric(1) )$candLL


################################
# INVESTIGATE THE SURACE OF THE LIKELIHOOD
# FOR TWO PARAMETERS

# loop over values of u and Bs
d.init <- c(14,1.7e-06,15,19,18,20,23,18)
zseq <- seq( 0, 1, .01 ) # range for u
yseq <- seq( 0, 1, .01 ) # range for Bs
mB.surf <- matrix( rep(NA,length(zseq)*length(yseq)), nrow=length(zseq) )
for( i in 1:length(zseq) ){
  for(j in 1:length(yseq) ){
    x.pars <- c( d.init, c( yseq[i],zseq[j] ) )
    	mB.surf[i,j] <- -funLL( 
    						x = x.pars, 
    						y=y, 
    						nObj=as.integer(nObj), 
    						combin=as.integer(combin), 
    						tricomb=tricomb, 
    						inhersim=as.numeric(inher.sim.freq), 
    						candLL=numeric(1), 
    						N=as.integer(length(d.t23m)), 
    						nTri=as.integer(dim(combin)[1]), 
    						debug=numeric(1) )$candLL
  				}		
			}	

sum(is.na(mB.surf))
head(mB.surf)

library(graphics)
# find point where the maximum likelihood is
minval <- which( mB.surf==min(mB.surf,na.rm=TRUE), arr.ind=TRUE )
#####
zseq[minval[1]]
yseq[minval[2]]
# plot the contour with a red point showing the MLEs
dev.new()
mB.surf2 <- ifelse( mB.surf>10000, NA, mB.surf )
lvls <- 30
filled.contour(x=zseq, y=yseq, z=mB.surf2, xlab="u", ylab="Bs", nlevels=lvls, col=gray(0:lvls/lvls), plot.axis={ 	points(zseq[minval[1]],yseq[minval[2]],pch=3, col="red", bg="black", lwd=2)	} )

# likelihood looks smooth

################################
# use an optimization algorithm
# to find the estimates at the 
# minimum negative log-likelihood

# wrapper for the negative log-likelihood
nll.tri.wrp <- function(xinits) {
	candpar <-  as.numeric( c( exp(xinits[1:(nObj+1)]), plogis(xinits[nObj+2]) ) )
	-funLL( x = candpar, y=y, nObj=as.integer(nObj), combin=as.integer(combin), tricomb=tricomb, inhersim=as.numeric(inher.sim.freq), candLL=numeric(1), N=as.integer(length(d.t23m)), nTri=as.integer(dim(combin)[1]), debug=numeric(1) )$candLL
}

# run below multiple times to find global minimum
mfit <- optim( x.inits, nll.tri.wrp, method="SANN", hessian=TRUE, control=list(trace=TRUE, maxit=10000))
for ( i in 1:5 ){
	mfit <- optim( mfit$par, nll.tri.wrp, method="SANN", hessian=TRUE, control=list(trace=TRUE, maxit=10000))
}
# check estimates
exp(mfit$par[1:nObj]) # signaling values
exp(mfit$par[nObj+1]) # beta, scaling value
plogis(mfit$par[(nObj+2)]) # u
sqrt( diag( solve( mfit$hessian ) ) )

# save objects to file when ready
#save.image( "AnalysisProgress.rdata")

# point precision is hard to estimate 

# try Monte Carlo methods,
# set up the input for the 
# Metropolis-Hastings sampler

library("BayesianTools")

# liklihood wrapper (not negative log-likelihood)
ll.tri.wrp <- function(xinits) {
	candpar <-  as.numeric( c( exp(xinits[1:(nObj+1)]), plogis(xinits[nObj+2]) ) )
	funLL( x = candpar, y=y, nObj=as.integer(nObj), combin=as.integer(combin), tricomb=tricomb, inhersim=as.numeric(inher.sim.freq), candLL=numeric(1), N=as.integer(length(d.t23m)), nTri=as.integer(dim(combin)[1]), debug=numeric(1) )$candLL
}

# parameter labels
par.labs <- c(paste0("d",1:nObj), "Bs", "u" )

# priors informed by optimization process above (mfit)
# constrain parameter Bs since it may be correlated
prior_density <- function(par) {
  sum(dnorm(par, mean = c(mfit$par), sd = c(rep(5,nObj), 0.01, 5), log = TRUE))
}

prior_sampler <- function(n = 1) {
  matrix(rnorm(n * length(mfit$par), mfit$par, c(rep(5,nObj), 0.01, 5)),
         ncol = length(mfit$par), byrow = TRUE)
}

prior <- createPrior(
  density = prior_density,
  sampler = prior_sampler,
  lower = rep(-Inf,10),
  upper = rep(Inf,10)
)

bayesianSetup = createBayesianSetup(likelihood = ll.tri.wrp, prior=prior, parallel = F, names=par.labs )
settings = list(iterations = 50000, nrChains= 3, message = TRUE, adapt=TRUE, optimize=F, startValue = mfit$par)

start <- proc.time()
out2 <- runMCMC(bayesianSetup = bayesianSetup, sampler = "Metropolis", settings = settings)
end <- proc.time()

(end - start)["elapsed"] 
# running across cores only gives a two-second improvement
# for a five and a half minute process

# MCMC samples plot for 
samp1 <- getSample( out2, start=25000, coda=TRUE, whichParameters=1:4)
dev.new(); plot(samp1)
samp2 <- getSample( out2, start=25000, coda=TRUE, whichParameters=5:7)
dev.new(); plot(samp2)
samp3 <- getSample( out2, start=25000, coda=TRUE, whichParameters=8:10)
dev.new(); plot(samp3)
sumMCMC <- summary(out2)

samp <- getSample( out2, start=10000, coda=FALSE, whichParameters=1:10)
dev.new(); plot(samp)

# save objects to file when ready
#save.image( "AnalysisProgress.rdata")

# <<<<<<<<<<<<<<<<<<<<<<<<
# FIGURE 5
niter <- dim(samp)[1]
transf <- function(x) c(exp(x[1:9]),plogis(x[10]))
#est <- transf( colMeans( samp[1:niter,]) )
est <- transf( sapply( 1:10, function(z) quantile( samp[1:niter,z], probs=0.5) ) )
est.lb <- transf( sapply( 1:10, function(z) quantile( samp[1:niter,z], probs=0.025) ) )
est.ub <- transf( sapply( 1:10, function(z) quantile( samp[1:niter,z], probs=0.975) ) )
data.frame( est, est.lb, est.ub )

imageDir <- paste0(wkdir,"/")


dev.new(width=9, height=5)

par(mfrow=c(1,1))
delt <- est[1:nObj]
delt.lb <- est.lb[1:nObj]
delt.ub <- est.ub[1:nObj]
b <- 0.25 # image spacer
delt.vals <- delt
delt.or <- order(delt.vals)
par(las=1, mar=c(4,12,1.5,17), cex=0.8)
mv <- 0.2
ymax <- ceiling(max(delt.ub))
plot( c(0.5,(nObj+0.5)), c(0,ymax), type="n", axes=FALSE, ylab="",xlab="" )
points( 1:nObj-mv, delt.vals[delt.or], pch=16, cex=2)
axis(1, at=1:nObj, labels=rep("",nObj))
axis(2)

library(jpeg)
for ( i in 1:nObj ){
whr <- delt.or[i]
ybottom <- -0.6; ytop=-0.3 # adjust to correct for stretching
img<-readJPEG(paste0(imageDir,type="m",whr,".jpeg"))
rasterImage(img, xleft=i-b,ybottom=ybottom,xright=i+b,ytop=ytop, xpd=TRUE)
lines( c(i-mv,i-mv), c(delt.lb[whr],delt.ub[whr]) )
lines( c(i-0.1-mv,i+0.1-mv), c(delt.lb[whr],delt.lb[whr])  )
lines( c(i-0.1-mv,i+0.1-mv), c(delt.ub[whr],delt.ub[whr])  )
}

mtext(side=2, "Signaling value\nof motif", line=10, las=1, adj=0)
utext <- bquote(italic(u) == .(format(est[nObj+2],digits=2 )) )	
utext <- bquote(italic(u) == .(paste(format(est[nObj+2],digits=3 )," (",format(est.lb[nObj+2],digits=3 ),", ",format(est.ub[nObj+2],digits=3 ),")",sep="") ))	
text( x=1, y=0.99*ymax, utext, adj=0, cex=1 )

# coefficient of selection Beta_s
Bs.set <- est[nObj+1]
Bs.set * delt.vals
invlgt.Bs <- plogis( Bs.set * delt.vals )
invlgt.Bs.lb <- plogis( Bs.set * delt.lb )
invlgt.Bs.ub <- plogis( Bs.set * delt.ub )
#plot(logit(Bs.set * delt.vals[delt.or]))

# add right axis
minaxis <- plogis(Bs.set*0)
maxaxis <- plogis(Bs.set*ymax)
intvls <- ( maxaxis - minaxis ) / ymax
lbsl <- seq( minaxis, maxaxis, intvls)
axis(4, at=0:ymax, labels=format(lbsl,digits=2), las=1)
mtext(side=4, "Association probability\nwith motif", line=4.5, las=1, cex=1, adj=0 )
scl <- function( x ) ( x - minaxis ) / intvls # scale points
for ( i in 1:nObj ){
whr <- delt.or[i]
points( i+mv, scl( invlgt.Bs[whr] ), pch=2, cex=1 )
points( i+mv, scl( invlgt.Bs[whr] ), pch=2, cex=2 )

lines( c(i+mv,i+mv), c(scl(invlgt.Bs.lb[whr]),scl(invlgt.Bs.ub[whr]) ) )
lines( c(i+mv-0.1,i+mv+0.1), c(scl(invlgt.Bs.lb[whr]),scl(invlgt.Bs.lb[whr]) )  )
lines( c(i+mv-0.1,i+mv+0.1), c(scl(invlgt.Bs.ub[whr]),scl(invlgt.Bs.ub[whr]) )  )
}


