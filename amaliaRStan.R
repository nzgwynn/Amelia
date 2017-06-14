
##
##  Assumes no missing data (other than due to diagnosis)
## Fixes number of time points at 5, uses same lambda for each time
##

## install.packages("rootSolve")
library("rootSolve")

Gamma<-function(to,do) phi^(to-1)*phi^(1-do)*(1-phi)^do
Delta<-function(ti, to,do) phi^(ti-1)*(1-theta)^(to-ti)*(1-theta)^(1-do)*theta^do

like<-function(to,do,x,beta,gamma){
	like0(to,do,x,beta,gamma)+likek(to,do,x,beta,gamma,k=1)+
		likek(to,do,x,beta,gamma,k=2)+likek(to,do,x,beta,gamma,k=3)+
		likek(to,do,x,beta,gamma,k=4)+likek(to,do,x,beta,gamma,k=5)
}

like0<-function(to,do,x,beta,gamma){
	risk<-exp(x*beta)
	haz<-1+exp(gamma)
	(haz^-risk)^to*Gamma(to,do)
}

likek<-function(to,do,x,beta,gamma,k){
	risk<-exp(x*beta)
	haz<-1+exp(gamma)
	ifelse(k<=to,((haz^-risk)^(k-1))*(1-haz^-risk)*Delta(k,to,do),0)
}

makedev <- function(d){
	function(pars) {-sum(log(like(d$to,d$do,d$x,pars[1],pars[2])))}
}

logit<-function(lambda) log(lambda/(1-lambda))
expit<-function(gamma) exp(gamma)/(1+exp(gamma))

make.data<-function(halfn, beta, gamma, phi, theta){
	n<-halfn*2
	x<-rep(0:1,halfn)
	risk<-exp(x*beta)
	lambda<-exp(gamma)/(1+exp(gamma))
	p<- (1- (1-lambda)^exp(x*beta))
	ti<-rgeom(n,p)+1
	to<-rep(Inf,n)
	do<-rep(0, n)
	for(s in 1:5){
		to<-ifelse(do==0, 
				ifelse(ifelse(ti<=s, rbinom(n,1,theta)==1,
                                              rbinom(n,1,1-phi)==1 ),s,to),to)
		do[to==s]<-1
		
	}
	to[to==Inf]<-5
	data.frame(to=to,do=do,x=x)
}

theta<-0.8;phi<-0.9

## one run
the.data = list(
  TDX = as.array(as.matrix(make.data(400, beta = 1.3, gamma = logit(0.05), 
                                                  theta = theta, phi = phi))),
  N = 800)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/gwynn/Documents/SurvivalAnalysis/Amelia")

fit <- stan(file="model.stan",
            data=the.data, init=myinits_f,
            iter=1000, chains=4)

# optim(c(1.3,logit(0.05)), the.minusloglike,method="BFGS")
# 
# 
# ## many runs
# aph<-replicate(100, optim(c(1.3,logit(0.05)), makedev(make.data(400,beta=1.3,gamma=logit(0.05),theta=0.8,phi=.9)), method="BFGS")$par)
# apply(aph,1,mean)-c(1.3,logit(0.05))
# apply(aph,1,sd)
# 
# 
# ## or with the improved optimx package, and comparing two optimisers to make sure they agree
# library(optimx)
# 
# aph<-replicate(100,coef(optimx(c(1.3,logit(0.05)),makedev(make.data(400,beta=1.3,gamma=logit(0.05),theta=0.8,phi=.9)))))
# 
# apply(aph,1:2.,mean)
# apply(aph,1:2.,sd)

## To work in RStan we need to invert the function and sample from this hazard rate. The
## following function can be used to find the inverse of a function. Our function should
## involve inserting a beta and gamma (lambda) then returning "T_real", putting it in an 
## optimiser, which will then spit out approximately the gammas and beta (lambda) we 
## used earlier.

## The makedev from above is for the.minusloglike. We need the likelihood.
# makedev.like <- function(d){
# 	function(pars) {sum(like(d$to,d$do,d$x,pars[1],pars[2]))}
# }
# 
# the.like = makedev.like(d = the.data)




# ## Now we need to invert the.like wrt to - the observed values of t
# 
# ## It might be easier to do this piece by piece.
# Gamma<-function(to,do) phi^(to-1)*phi^(1-do)*(1-phi)^do
# Gamma.Inv <- function(do,to) log((to*phi)/(phi^(1-do)*(1-phi)^do))/log(phi)
# 
# ## Using uniroot to force R to find the inverse
# Delta<-function(ti, to, do) phi^(ti-1)*(1-theta)^(to-ti)*(1-theta)^(1-do)*theta^do
# make.Delta.Inv <- function(ti, do){
#   function(to) Delta(ti = ti, to, do = do)
# }
# 
# find.Delta.Inv = make.Delta.Inv(ti = 1, do = 1) ## A function of to only
# Delta.Inv = inverse(find.Delta.Inv) ## works
# 
# try = Delta(ti = 1, to=2, do=1)
# Delta.Inv(try) ## should be 2
# 
# ## This does not work as must take log(log(1-theta)) and 0<theta<1
# ## Delta.Inv <- function(ti, to, do) log(to*(ti*log(1-theta))/ (phi^(ti-1)*(1-theta)^(1-do)*theta^do))/log(1-theta)
# 
# ## Using uniroot.all to force the R to find the inverse function
# inverse.all = function (f, lower = 0, upper = 16) {
#    function (y) uniroot.all((function (x) f(x) - y), lower = lower, upper = upper, tol = .Machine$double.eps^0.4, maxiter = 10000, n = 10000)
# }
# 
# ## Making the function
# like<-function(to,do,x,beta,gamma){
# 	like0(to,do,x,beta,gamma)+likek(to,do,x,beta,gamma,k=1)+
# 		likek(to,do,x,beta,gamma,k=2)+likek(to,do,x,beta,gamma,k=3)+
# 		likek(to,do,x,beta,gamma,k=4)+likek(to,do,x,beta,gamma,k=5)
# }
# 
# make.like.Inv <- function(do, x, beta, gamma){
#   function(to) like(to, do = do, x = x, beta = beta, gamma = gamma)
# }
# 
# 
# 
# find.like.Inv = make.like.Inv(do = 1, x = 1, beta = 1.3, gamma = logit(0.05)) ## A function of to only
# like.Inv = inverse.all(find.like.Inv) 

# try = like(to = 3, do = 1, x = 1, beta = 1.3, gamma = logit(0.05))
# like.Inv(try)[1] ## It works!! This should be to from the input above!


