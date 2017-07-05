make.one.BP<-function(startrange, trend, error.sd, treatment, times, end.of.treatment, carryover){
  untreated <- runif(1,startrange[1], startrange[2]) + trend*times
  if (carryover == 0)
  {treated <- untreated + treatment*(times<=end.of.treatment)}
  else {treated <- untreated + treatment*(times <= end.of.treatment) -
          treatment*(times > end.of.treatment &
                       times < end.of.treatment + carryover)*
          ((times - end.of.treatment - carryover)/carryover)}
  measured<- treated + rnorm(length(times), m=0, s=error.sd)
  return(c(measured, treated))
}

rule1over <- function(person){
  min(which(person>140))
} 

rule2consec <- function(person){
   min(which(c(FALSE,person>140) & c(person>140,FALSE)))
}

rule3over <- function(person){
     w<-which(person>140)
     if (length(w)>2) w[3] else Inf
}

ruleavg2consec <- function(person) {
        w<-which(.5*(person+c(NA,person[-length(person)]))>140)
        if (length(w)==0)
          Inf
        else
          w[1]
}

ruleavg3consec <- function(person) {
     w<-which((person+c(NA, person[-length(person)])+c(NA, NA, person[-(length(person)-(0:1))]))/3>140)
     if (length(w)==0)
       Inf
     else
       w[1]
}

logit<-function(lambda) log(lambda/(1-lambda))
expit<-function(gamma) exp(gamma)/(1+exp(gamma))

nreps = 400
nreps.phi.theat = 10^6
## startrange = c(125, 140)
## THOMAS I HAD TO CHANGE THE START RANGE B/C NO ONE WAS BEING DIAGNOSED! 
## For use with less error and fewer measurement times.
startrange = c(135,145)
par = rev(expand.grid(rev(
   list(treatment= c(-5,-10),
        trend =  c(1, 2),  ## No one is truly diagnosed when trend = 0
        error.sd = c(3, 5, 7),
        times = list( seq (.5, 4, by = .5), seq ( 1, 4, by = 1), seq ( .25, 4, by = .25) ),
        end.of.treatment = c(1, 1.5, 2, 2.5, 3),
        carryover = seq(0, 2, by =.5),
        rule = c("rule1over", "rule2consec", "ruleavg2consec", "rule3over", "ruleavg3consec"))),
   stringsAsFactors = FALSE))

########### INPUT BP MEASUREMENTS AND OUTPUTS SURVIVAL DATA ###########################
## trt = 1, for treatment grp, true = 1 for Long Term Avg, true = 0 for Instantaneous
makeSurvData <- function(data, times, rule, nreps, trt, true){
  x = rep(trt, nreps)

  ## For adjusted Parameters with std.error smaller and 2 times use this
  ##times.length = length(times) also need to change T below.
  times.length = length(times[[1]])
  if(true == 0){BP = data[1:times.length,]
                Diag = suppressWarnings(apply(BP, 2, rule))
              }else{
                BP = data[(times.length + 1):(times.length*2),]
                Diag = suppressWarnings(apply(BP,2, rule1over))
              }

  lengths = as.vector(table(c(Diag, 1:times.length, Inf)) - 1)  
  D = c(rep(0, lengths[times.length+ 1]), rep(1, nreps - (lengths[times.length+ 1])))

  ## a = No. of times measurements occur.
  a = ifelse(times.length == 8, 2, ifelse(times.length == 16, 4, 1))

  ## Changing measurement times as a portion of a year (.25, .5, .75, ...) to an integer (1, 2, 3, ...) 
  T = c(rep(max(a*times[[1]]), lengths[max(times.length)+ 1]),  rep(a*times[[1]], lengths[1:times.length]))
  
  ## For adjusted Parameters with std.error smaller and 2 times
  ## T = c(rep(max(a*times), lengths[times.length+1]),  rep(a*times, lengths[1:times.length])) 

 matrix(c(to=T, do=D, x=x), ncol = 3)
}

make.data <- function(P){
  C = replicate(nreps, with(P, make.one.BP(startrange = startrange,
                                                            trend = trend,
                                                            error.sd = error.sd,
                                                            treatment = 0,
                                                            times = as.numeric(unlist(times)),
                                                            end.of.treatment = end.of.treatment,
                                                            carryover = as.numeric(carryover))))
  
  Tt = replicate(nreps, with(P, make.one.BP(startrange = startrange,
                                                            trend = trend,
                                                            error.sd = error.sd,
                                                            treatment = 0,
                                                            times = as.numeric(unlist(times)),
                                                            end.of.treatment = end.of.treatment,
                                                            carryover = as.numeric(carryover))))
  
  
  D = rbind(makeSurvData(data = Tt, times = with(P, times), 
                         rule = get(P$rule), nreps = 400, trt = 1, true = 0),
            makeSurvData(data = C, times = with(P, times), 
                         rule = get(P$rule), nreps = 400, trt = 0, true = 0))
  
  rbind(D)
}


## Inputs for function above to make data
theta<-0.8;phi<-0.9

i = round(runif(1,1,4500),0)
## one run
data = list(
  TDX = as.array(make.data(P = par[i,])),
  N = 800)



library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd("/Users/gwynn/Documents/SurvivalAnalysis/Amelia")

myinits_f <- function(chainnum){
  return(list(beta = 1, gamma = 1))
}

fit <- stan(file = "MyBrokenModel.stan",
            data = data, init = myinits_f,
            iter = 1000, chains = 4)





## One way to plot things.
for(i in 1:length(estimates)){
  True.Gamma[i,] = estimates[[i]]$'True HZD'
  gamma.hats[i,] = exp(estimates[[i]]$Method2.NLM)
}

pdf("Bias.pdf")
plot(y=gamma.hats, x = True.Gamma,  main = "True.HZD vs. exp(Optim)", 
     ylab = "Gamma.Hats", xlab = "True")
abline(a=0,b=1)
dev.off()


