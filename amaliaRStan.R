## Extraneous functions used to not put bounds on lambda

logit<-function(lambda) log(lambda/(1-lambda))
expit<-function(gamma) exp(gamma)/(1+exp(gamma))

## Function used to make data here

make.Adata<-function(halfn, beta, gamma, phi, theta){
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

## Inputs for function above to make data
theta<-0.8;phi<-0.9

S = 1000
## one run
Adata = list(
  TDX = as.array(as.matrix(make.Adata(S, beta = 1.3, gamma = logit(0.05), 
                                                  theta = theta, phi = phi))),
  N = 2*S)


library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/gwynn/Documents/SurvivalAnalysis/Amelia")

myinits_Af <- function(chainnum){
  return(list(beta = 1.3, gamma = logit(0.05)))
}

Afit <- stan(file = "model.stan",
            data = Adata, init = myinits_Af,
            iter = 2000, chains = 4)

traceplot(Afit, warmup = "FALSE")
pairs(Afit)
hist(extract(Afit)$beta)
hist(extract(Afit)$gamma)


