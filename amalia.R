
##
##  Assumes no missing data (other than due to diagnosis)
## Fixes number of time points at 5, uses same lambda for each time
##

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
				ifelse( ifelse(ti<=s, rbinom(n,1,theta)==1,rbinom(n,1,1-phi)==1 ),
					 s,to),
				to)
		do[to==s]<-1
		
	}
	to[to==Inf]<-5
	data.frame(to=to,do=do,x=x)
}

theta<-0.8;phi<-0.9

## one run
the.data<-make.data(400,beta=1.3,gamma=logit(0.05),theta=theta,phi=phi)
the.minusloglike<-makedev(the.data)

optim(c(1.3,logit(0.05)), the.minusloglike,method="BFGS")


## many runs

aph<-replicate(100,optim(c(1.3,logit(0.05)),makedev(make.data(400,beta=1.3,gamma=logit(0.05),theta=0.8,phi=.9)),method="BFGS")$par)
apply(aph,1,mean)-c(1.3,logit(0.05))
apply(aph,1,sd)


##or with the improved optimx package, and comparing two optimisers to make sure they agree
library(optimx)

aph<-replicate(100,coef(optimx(c(1.3,logit(0.05)),makedev(make.data(400,beta=1.3,gamma=logit(0.05),theta=0.8,phi=.9)))))

apply(aph,1:2.,mean)
apply(aph,1:2.,sd)