
##
##  Assumes no missing data (other than due to diagnosis)
## Fixes number of time points at 5, uses same lambda for each time
##

## functions to be used in the likelihood below which account 
## for mismeasurement in the data
Gamma<-function(to,do) phi^(to-1)*phi^(1-do)*(1-phi)^do
Delta<-function(ti, to,do) phi^(ti-1)*(1-theta)^(to-ti)*(1-theta)^(1-do)*theta^do

## the full likelihood which combines subsequent functions below
like<-function(to,do,x,beta,gamma){
	like0(to,do,x,beta,gamma)+likek(to,do,x,beta,gamma,k=1)+
		likek(to,do,x,beta,gamma,k=2)+likek(to,do,x,beta,gamma,k=3)+
		likek(to,do,x,beta,gamma,k=4)+likek(to,do,x,beta,gamma,k=5)
	
}

## The first product from equation 4 in Amalia's paper
like0<-function(to,do,x,beta,gamma){
	risk<-exp(x*beta)
	haz<-1+exp(gamma)
	(haz^-risk)^to*Gamma(to,do)
}

## The summation of the horrible product thing from equation 4 in Amalia's paper
likek<-function(to,do,x,beta,gamma,k){
	risk<-exp(x*beta)
	haz<-1+exp(gamma)
	ifelse(k<=to,((haz^-risk)^(k-1))*(1-haz^-risk)*Delta(k,to,do),0)
}

## This function makes the likelihood by putting all of the data into it
## it returns a function of pars[1] and pars[2] that we then optimize below
makedev <- function(d){
	function(pars) {-sum(log(like(d$to,d$do,d$x,pars[1],pars[2])))}
}

## Random functions that we need to use to make it easier to optimize.
## In section 3.2.3 we reparameterize lambda to avoid bounds 
logit<-function(lambda) log(lambda/(1-lambda))
expit<-function(gamma) exp(gamma)/(1+exp(gamma))

## A function that makes data to be input into the likelihood made above. It 
## makes survival analysis data then adds some error according to phi and theta 
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
				ifelse( ifelse(ti<=s, rbinom(n,1,theta)==1,rbinom(n,1,1-phi)==1 ),
					 s,to),
				to)
		do[to==s]<-1
		
	}
	to[to==Inf]<-5
	data.frame(to=to,do=do,x=x)
}

## giving the computer phi and theta
theta<-0.8;phi<-0.9

## one run
## making the data once then making the likelihood below and optimizing
the.data<-make.data(400,beta=1.3,gamma=logit(0.05),theta=theta,phi=phi)
the.minusloglike<-makedev(the.data)

optim(c(1.3,logit(0.05)), the.minusloglike,method="BFGS")


## many runs
## making the data 100 times to get 100 estimates to see if the estimates are
## biased. Here theta and phi are the same as they are above (ie accurate)
aph<-replicate(100,optim(c(1.3,logit(0.05)),
                         makedev(make.data(400,beta=1.3,gamma=logit(0.05),
                                           theta=0.8,phi=.9)),method="BFGS")$par)
apply(aph,1,mean)-c(1.3,logit(0.05)) ## beta and gamma are 1.3 and logit(0.05)
apply(aph,1,sd)


##or with the improved optimx package, and comparing two optimisers to make sure 
## they agree
## doing above with a different optimization method.
## making the data 100 times to get 100 estimates to see if the estimates are
## biased. Here theta and phi are the same as they are above (ie accurate)
install.packages("optimx")
library(optimx)

aph<-replicate(100,coef(optimx(c(1.3,logit(0.05)),
                               makedev(make.data(400,beta=1.3,gamma=logit(0.05),
                                                 theta=0.8,phi=.9)))))

apply(aph,1:2.,mean) ## you see that the mean are close to 1.3 and logit(0.05)
apply(aph,1:2.,sd)

## Here I slightly misspecify phi and theta and the whole thing breaks.
## Which is where I'd like to use Bayesian Stats to put a prior on for 
## beta and gamma to see if I can force the optimiser to find something 
## reasonable.
aph.bias<-replicate(100,coef(
  optimx(c(1.3,logit(0.05)), makedev(make.data(400,beta=1.3,gamma=logit(0.05),
                                               theta=0.75,phi=.85)))))

apply(aph.bias,1:2.,mean) ## Nowhere near 1.3 logit(0.05)
apply(aph.bias,1:2.,sd)