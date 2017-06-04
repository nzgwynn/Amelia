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

colcumprod <- function(x) {
  y = x
  for(j in 2:ncol(x))
    y[,j] = y[,j] * y[,j-1]
  y
}

makeIndic <-function(t0, ncols) {
  matrix(
         as.numeric(rep(1:ncols, each = length(t0)) <= t0),
         nc = ncols)
}

logit<-function(lambda) log(lambda/(1-lambda))
expit<-function(gamma) exp(gamma)/(1+exp(gamma))

##install.packages("optimx")
library("optimx")

nreps = 400
nreps.phi.theat = 10^6
## startrange = c(125, 140)
## THOMAS I HAD TO CHANGE THE START RANGE B/C NO ONE WAS BEING DIAGNOSED! For use with less error and fewer measurement times.
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
  times.length = length(times)
  ##times.length = length(times[[1]])
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
  ##T = c(rep(max(a*times[[1]]), lengths[max(times.length)+ 1] + 1),  rep(a*times[[1]], lengths[1:times.length]))
  
  ## For adjusted Parameters with std.error smaller and 2 times
  T = c(rep(max(a*times), lengths[times.length+1]),  rep(a*times, lengths[1:times.length])) 

 data.frame(to=T, do=D, x=x)
}

findPhiTheta <- function(data, times, rule){
  ## For less error and 2 measurement times use:
  times.length = length(times)
  ## This confuses the poor thing:
  ## times.length = length(times[[1]])
  Instant.BP = data[1:times.length,]
  Long.Term.BP = data[(times.length + 1):(times.length*2),]
  
  Truth = c(suppressWarnings(apply(Long.Term.BP, 2, rule1over)),1:times.length, Inf)
  Diagnosis = c(suppressWarnings(apply(Instant.BP, 2, rule)),1:times.length, Inf)

  D = cbind(Truth, Diagnosis)
  T = table(Diagnosis, Truth)
  
  num = diag(T)[1:times.length]-1
  
  theta.denom.part1 = table(Truth)[1:times.length]-1
  phi.denom.part1 = rev(as.numeric(cumsum(rev(table(Diagnosis)-1))))
  theta.denom.part2 = phi.part2 = 0
  
  ## The minus 1 to theta.num is added since above we add 1:times.length to both Truth and Diagnosis so that T has 1 extra along the diagonal. Note that phi.num & phi.denom are for 1-phi
  for(k in 1:times.length){
    theta.denom.part2[k+1] = sum(T[1:k,(k+1)])
    phi.part2[k] =  sum(T[k,1:k])-1
  }
  
  phi.num = as.numeric(table(Diagnosis))[1:times.length] - phi.part2 - 1
  ## The plus 1 to phi.denom one is added since above we add 1:times.length to both Truth and Diagnosis so that T has 1 extra along the diagonal. 

  fill.phi = rep(1, times.length - length(theta.denom.part1))
  fill.theta = rep(0, times.length - length(phi.num))
  #theta.denom = theta.denom.part1 - theta.denom.part2[1:times.length]
  #The minus 1:times.length to theta.denom is there as we add 1:times.length to both Truth. 
  theta.denom = table(Truth)[1:times.length]-theta.denom.part2[1:times.length]-1
  phi.denom = phi.denom.part1[1:times.length] - phi.part2[1:times.length]
  
  return(c(fill.phi, 1-phi.num/phi.denom, fill.theta, (num + 1)/(theta.denom+1)))
}


## Vectorised functions that agree with Thomas' functions 
Gamma.Vec<-function(to,do,phi){
  lots.of.phi<-matrix(phi,nrow=length(to),ncol=length(phi),byrow=TRUE)
  lots.of.i<-matrix(1:length(phi), nrow=length(to),ncol=length(phi),byrow=TRUE)
  lots.of.phi[lots.of.i>=to]<-1
  apply(lots.of.phi,1,prod)*phi[to]^(1-do)*(1-phi[to])^do
}

Delta.Vec<- function(ti,to,do, phi, theta){
  lots.of.phi <- matrix(phi,nrow=length(to),ncol=length(phi), byrow=TRUE)
  lots.of.i<-matrix(1:length(phi), nrow=length(to),ncol=length(phi),byrow=TRUE)
  lots.of.phi[lots.of.i>=ti]<-1
  part.one = apply(lots.of.phi,1,prod)

  lots.of.theta = 1 - matrix(theta, nrow = length(to), ncol = length(theta), byrow=TRUE)
  lots.of.theta[lots.of.i<ti | lots.of.i>=to]<-1

  part.one*apply(lots.of.theta,1,prod)*(1-theta[to])^(1-do)*theta[to]^do
}

###### THIS FUNCTION INPUTS DATA AND OUTPUTS GAMMA'S ##################################
max.optimx <- function(Phi, Theta, truth, data, times.length, initial.vals, method){
  if(method == 1) makedev = makedev.log.like1
  if(method == 2) makedev = make.log.like.Vec
  if(method == 3) makedev = makedev.log.like2
  if(method == 4) makedev = makedev.log.like3
  
  est = list(2)
  the.minusloglike<-makedev(data, phi=Phi, theta=Theta)
  
  est[[2]] = deriv = delta = NA
  
  est[[1]] = optimx(initial.vals, the.minusloglike, method = "L-BFGS-B", control = list(factr = 1e-10))

  optim.gamma = as.numeric(coef(est[[1]]))
  
  if(truth == 1){
    phi = theta = rep(1, times.length)
    for(i in 1:times.length){
      delta = rep(0,times.length)
      delta[i] = 1e-5
      deriv[i] = (the.minusloglike(pars = optim.gamma + delta) -
                  the.minusloglike(pars = optim.gamma))*1e5
    }
    est[[2]] = deriv
  }
  
  est
}

baseline.hazard <- function(SurvData, times.length){
  if(sum(SurvData[,"do"]) > 0){
    num = as.vector(table(SurvData[which(SurvData[,"do"]==1),]))
    denom = nreps - cumsum(as.vector(table(SurvData[which(SurvData[,"do"]==
      1),])))
  }else{
    num = rep(0, times.length)
    denom = rep(1, times.length)
  }
  
  HZD = num/(num+denom)
  HZD.length = length(HZD)

  if(HZD.length != times.length){
    HZD = c(rep(0, times.length-HZD.length),HZD)
  }

  HZD
}
    
find.one.baseline.hazard <- function(Parameters){
  data.ctrl = replicate(nreps, with(Parameters, make.one.BP(startrange = startrange,
    trend = trend,
    error.sd = error.sd,
    treatment = 0,
    times = as.numeric(unlist(times)),
    end.of.treatment = end.of.treatment,
    carryover = as.numeric(carryover))))

  mis.data = makeSurvData(data = data.ctrl, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 0, true = 0)

  mis.HZD = baseline.hazard(mis.data, times.length = length(unlist(with(Parameters, times))))
  
  true.data = makeSurvData(data = data.ctrl, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 0, true = 1)

  true.HZD = baseline.hazard(true.data, times.length = length(unlist(with(Parameters, times))))
  
  c(mis.HZD,true.HZD)
}

find.one.bias.NB <- function(Parameters,error,PT.Tr){
  Parameters = c(unlist(Parameters[c(1:2,5:7)]), times = list(1:4), error.sd = error)

  Parameters = c(unlist(Parameters[c(1:4)]), times = list(1:4), error.sd = error, rule = "rule1over")
  
  data.trt = replicate(nreps, with(Parameters, make.one.BP(startrange = startrange,
    trend = as.numeric(trend),
    error.sd = error.sd,
    treatment = 0,
    times = as.numeric(unlist(times)),
    end.of.treatment = as.numeric(end.of.treatment),
    carryover = as.numeric(carryover))))
  
  data.ctrl = replicate(nreps, with(Parameters, make.one.BP(startrange = startrange,
    trend = as.numeric(trend),
    error.sd = error.sd,
    treatment = 0,
    times = as.numeric(unlist(times)),
    end.of.treatment = as.numeric(end.of.treatment),
    carryover = as.numeric(carryover))))

  PT.Est = findPhiTheta(data = cbind(data.trt, data.ctrl), times = with(Parameters,times), rule = get(Parameters$rule))

  mis.data = rbind(makeSurvData(data = data.ctrl, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 0, true = 0), makeSurvData(data = data.trt, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 1, true = 0))

  true.data = rbind(makeSurvData(data = data.ctrl, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 0, true = 1), makeSurvData(data = data.trt, times = with(Parameters,times), rule = get(Parameters$rule), nreps = nreps, trt = 1, true = 1))

  solution = list()

  ## When using less error and shorter measurement times use:
  times.length = length(with(Parameters,times))
  ## Not this, it confuses the poor thing
  ## times.length = length(with(Parameters,times)[[1]])

  Phi.Est = PT.Est[1:times.length]
  Theta.Est = PT.Est[(times.length+1):(2*times.length)]

  Phi.Tr = PT.Tr[1:times.length]
  Theta.Tr = PT.Tr[(times.length+1):(2*times.length)]

  IntVals = log(baseline.hazard.NB(SurvData = true.data, times.length = times.length))

  if(sum(is.finite(IntVals)) == times.length){
    ##method1.optimx = max.optimx(Phi = Phi, Theta = Theta, truth = 0, data= mis.data, times.length = times.length,initial.vals = IntVals, method = 1)

    ##method2.optimx = max.optimx(Phi = Phi, Theta = Theta, truth = 0, data = mis.data, times.length = times.length, initial.vals = IntVals, method = 2)

    ##method3.optimx = max.optimx(Phi = Phi, Theta = Theta, truth = 0, data= mis.data, times.length = times.length,initial.vals = IntVals, method = 3)

    ##method4.optimx = max.optimx(Phi = Phi, Theta = Theta, truth = 0, data = mis.data, times.length = times.length, initial.vals = IntVals, method = 4)

    ##method1.nlm = method1(phi = Phi, theta = Theta, data = mis.data, times.length = times.length, initial.vals = IntVals)

    method2.nlm.Tr = method2(phi = Phi.Tr, theta = Theta.Tr, data = mis.data, times.length = times.length, initial.vals = IntVals/2)
    
    method2.nlm.Est = method2(phi = Phi.Est, theta = Theta.Est, data = mis.data, times.length = times.length, initial.vals = IntVals/2)

    ##method3.nlm = method3(phi = Phi, theta = Theta, data = mis.data, times.length = times.length, initial.vals = IntVals)

    ##method4.nlm = method4(phi = Phi, theta = Theta, data = mis.data, times.length = times.length, initial.vals = IntVals)
  
    solution[[1]] = Parameters
    solution[[2]] = Phi.Tr
    solution[[3]] = Theta.Tr
    solution[[4]] = as.numeric(method2.nlm.Tr$estimate)
    solution[[5]] = c(method2.nlm.Tr$code)
    
    solution[[6]] = Phi.Est
    solution[[7]] = Theta.Est
    solution[[8]] = as.numeric(method2.nlm.Est$estimate)
    solution[[9]] = c(method2.nlm.Est$code)
   
    solution[[10]] = IntVals
    
    ## solution[[5]] = as.numeric(coef(method1.optimx[[1]]))
    ## solution[[6]] = as.numeric(coef(method2.optimx[[1]]))
    ## solution[[7]] = as.numeric(coef(method3.optimx[[1]]))
    ## solution[[8]] = as.numeric(coef(method4.optimx[[1]]))
    
    ## solution[[9]] = c(method1.optimx[[1]]$kkt1,method1.optimx[[1]]$kkt2)
    ## solution[[10]] = c(method2.optimx[[1]]$kkt1,method2.optimx[[1]]$kkt2)
    ## solution[[11]] = c(method3.optimx[[1]]$kkt1,method3.optimx[[1]]$kkt2)
    ## solution[[12]] = c(method4.optimx[[1]]$kkt1,method4.optimx[[1]]$kkt2)
    
    ## solution[[13]] = as.numeric(method1.nlm$estimate)
    ## solution[[14]] = as.numeric(method2.nlm$estimate)
    ## solution[[15]] = as.numeric(method3.nlm$estimate)
    ## solution[[16]] = as.numeric(method4.nlm$estimate)
    
    ## solution[[17]] = c(method1.nlm$code)
    ## solution[[18]] = c(method2.nlm$code)
    ## solution[[19]] = c(method3.nlm$code)
    ## solution[[20]] = c(method4.nlm$code)

    names.optimx = names.nlm = code.optimx = code.nlm = character(4)
    for(m in 1:4) {names.optimx[m] = paste0("Method",m,".Optimx")
                   names.nlm[m] = paste0("Method",m,".NLM")
                   code.optimx[m] = paste0("Code",m,".Optimx")
                   code.nlm[m] =  paste0("Code",m,".NLM")}
    names(solution) = c("Pars", "Phi.LargeN", "Theta.LargeN","Est.LargeN","Code.LargeN", "Phi.SmallN", "Theta.SmallN","Est.SmallN","Code.SmallN", "True HZD")##, names.optimx, code.optimx, names.nlm, code.nlm)
    
  }else{
    solution[[1]] = Parameters
    solution[[2]] = Phi
    solution[[3]] = Theta

    names(solution) = c("Pars", "Phi", "Theta")
  }
  
  solution
}

############################################################################################## THIS IS THE COMBINED BASELINE HZD FOR BOTH TREATMENT AND CONTROL GROUPS ##################################################################################################### IT ASSUMES THAT BOTH TRT & CTRL GROUPS HAVE FIRST DIAGNOSED AT SAME TIME #############################################################################################
####### THIS FINDS LAMBDA/(1-LAMBDA) #########################################################################################################################################
baseline.hazard.NB <- function(SurvData, times.length){
  if(sum(SurvData[,"do"]) > 0){
    index = length(as.numeric(row.names(table(SurvData[which(SurvData[,"do"]==1),]))))
    num = as.vector(table(SurvData[which(SurvData[,"do"]==1),]),)[1:index] +
      as.vector(table(SurvData[which(SurvData[,"do"]==1),]))[(index + 1):(index*2)]
    denom = nreps*2 - (cumsum(as.vector(table(SurvData[which(SurvData[,"do"]==1),])))[1:index] + cumsum(as.vector(table(SurvData[which(SurvData[,"do"]==1),]))[(index + 1):(index*2)])) 
  }else{
    num = rep(0, times.length)
    denom = rep(1, times.length)
  }
  
  HZD = num/denom
  HZD.length = length(HZD)

  if(HZD.length != times.length){
    HZD = c(rep(0, times.length-HZD.length),HZD)
  }
  
  HZD
}

############################################################################################################## CODING THE LOG(LIKE) - NEW METHOD ####################################################################################################################

# k = 1, ... no of time points.  Not times themselves.
log.like1<-function(to,do,x,gamma,phi,theta){
  product = matrix(NA, nrow = length(x), ncol = length(gamma))
  summation = matrix(NA, nrow = length(x), ncol = length(gamma))
  solution = numeric(length(to))
  
  for(j in 1:length(gamma)){
      product[,j] = like0(x,gamma[j])
      summation[,j] = likek.nonconstant1(to,do,x,gamma,k=j,phi,theta)
    }

  for(i in 1:length(to)){
    solution[i] = Gamma.Vec(to[i],do[i],phi)*prod(product[i,][1:to[i]])+sum(summation[i,])
  }
  
  sum(log(solution))
}

like1<-function(to,do,x,gamma,phi,theta){
  product = matrix(NA, nrow = length(x), ncol = length(gamma))
  summation = matrix(NA, nrow = length(x), ncol = length(gamma))
  solution = numeric(length(to))
  
  for(j in 1:length(gamma)){
      product[,j] = like0(x,gamma[j])
      summation[,j] = likek.nonconstant1(to,do,x,gamma,k=j,phi,theta)
    }

  for(i in 1:length(to)){
    solution[i] = Gamma.Vec(to[i],do[i],phi)*prod(product[i,][1:to[i]])+
      sum(summation[i,])
  }

  solution
}

likek.nonconstant1<-function(to,do,x,gamma,k,phi,theta){
  risk<-exp(x)
  haz<-1+exp(gamma)
  out = outer(rep(1,2*nreps),haz)
  solution = numeric(length(to))

  for(j in 1:length(gamma)){
    out[,j] = out[,j]^-risk
  }

  for(l in 1:length(to)){
    if(k<=to[l]){
      if(k>1){
        solution[l] = prod(out[l,1:(k-1)])*(1-out[l,k])*
          Delta.Vec(k,to[l],do[l],phi,theta)
      }else{
        solution[l] = (1-out[l,k])*Delta.Vec(k,to[l],do[l],phi,theta)
    }}else{
      solution[l] = 0
    }
  }
  solution
}

like0<-function(x,gamma){
	risk<-exp(x)
	haz<-1+exp(gamma)
	haz^-risk
}

makedev.log.like1 <- function(d, phi, theta){
	function(pars){-sum(log.like1(d$to,d$do,d$x,pars,phi,theta))}
}

makedev.like1 <- function(d, phi, theta){
	function(pars){-sum(like1(d$to,d$do,d$x,pars,phi,theta))}
}
######################################################################################## NO BETA'S (NB) INCLUDED HERE - CODING FOR THE LIKELIHOOD PG 949 AMALIA  VECTORISED #########################################################################################
log.like.Vec<-function(to,do,x,gamma,phi,theta){
  n = length(to)
  p = length(gamma)
  D = matrix(NA, nrow = n, ncol = p)
  lots.of.i = matrix(1:p, nrow = n, ncol =p, byrow = TRUE)
  
  risk<-exp(x*0)
  haz<-1+exp(gamma)
  B = B1 = outer(risk, haz, function(x,y) y^-x)
  Ind = makeIndic(t0 = to, ncols = p)
  B1[lots.of.i>to]=1

  for(i in 1:p){D[,i] = Delta.Vec(ti=i,to=to,do=do, phi=phi,theta=theta)
              }
  
  sum(log(apply((cbind(1,colcumprod(B)[,-p])*(1-B)*D*Ind),1,sum) +
    Gamma.Vec(to,do,phi)*colcumprod(B1)[,p]))
}

like.Vec<-function(to,do,x,gamma,phi,theta){
  n = length(to)
  p = length(gamma)
  D = matrix(NA, nrow = n, ncol = p)
  lots.of.i = matrix(1:p, nrow = n, ncol =p, byrow = TRUE)
  
  risk<-exp(x*0)
  haz<-1+exp(gamma)
  B = B1 = outer(risk, haz, function(x,y) y^-x)
  Ind = makeIndic(t0 = to, ncols = p)
  B1[lots.of.i>to]=1

  for(i in 1:p){D[,i] = Delta.Vec(ti=i,to=to,do=do, phi=phi,theta=theta)
              }
  
  apply((cbind(1,colcumprod(B)[,-p])*(1-B)*D*Ind),1,sum) +
    Gamma.Vec(to,do,phi)*colcumprod(B1)[,p]
}

make.log.like.Vec <- function(d,phi,theta){
  function(pars){-sum(log.like.Vec(to=d$to,do=d$do,x=d$x, gamma=pars,phi=phi,theta=theta))}
}

make.like.Vec <- function(d,phi,theta){
  function(pars){-sum(like.Vec(to=d$to,do=d$do,x=d$x, gamma=pars,phi=phi,theta=theta))}
}

############################################################################################## NO BETA'S (NB) INCLUDED HERE - CODING FOR THE LIKELIHOOD PG 949 AMALIA ###############################################################################################
# k = 1, ... no of time points.  Not times themselves.
log.like2<-function(to,do,x,gamma,phi,theta){
  product = matrix(NA, nrow = length(x), ncol = length(gamma))
  summation = matrix(NA, nrow = length(x), ncol = length(gamma))
  solution = numeric(length(to))
  
  for(j in 1:length(gamma)){
      product[,j] = like0.NB(to,x,gamma[j])
      summation[,j] = likek.nonconstant2(to,do,x,gamma,k=j,phi,theta)
    }

  for(i in 1:length(to)){
    solution[i] = Gamma.Vec(to[i],do[i],phi)*prod(product[i,1:to[i]]) + sum(summation[i,])
  }

  sum(log(solution))
}

like2<-function(to,do,x,gamma,phi,theta){
  product = matrix(NA, nrow = length(x), ncol = length(gamma))
  summation = matrix(NA, nrow = length(x), ncol = length(gamma))
  solution = numeric(length(to))
  
  for(j in 1:length(gamma)){
      product[,j] = like0.NB(to,x,gamma[j])
      summation[,j] = likek.nonconstant2(to,do,x,gamma,k=j,phi,theta)
    }

  for(i in 1:length(to)){
    solution[i] = Gamma.Vec(to[i],do[i],phi)*prod(product[i,1:to[i]]) + sum(summation[i,])
  }

  solution
}

likek.nonconstant2<-function(to,do,x,gamma,k,phi,theta){
  risk<-exp(x)
  haz<-1+exp(gamma)
  out = outer(rep(1,2*nreps),haz)
  solution = numeric(length(to))

  for(j in 1:length(gamma)){
    out[,j] = out[,j]^-risk
  }

  for(l in 1:length(to)){
    if(k<=to[l]){
      if(k>1){
        solution[l] = prod(out[l,1:(k-1)])*(1-out[l,k])*
          Delta.Vec(ti=k,to=to[l],do=do[l],phi=phi,theta=theta)
      }else{
        solution[l] = (1-out[l,k])*Delta.Vec(ti=k,to=to[l],do=do[l],phi=phi,theta=theta)
    }}else{
      solution[l] = 0
    }
  }
  solution
}

like0.NB<-function(to,x,gamma){
  risk<-exp(x)
  haz<-1+exp(gamma)
  haz^-risk
}

makedev.log.like2 <- function(d, phi, theta){
	function(pars){-sum(log.like2(d$to,d$do,d$x,pars,phi=phi,theta=theta))}
}

makedev.like2 <- function(d, phi, theta){
	function(pars){sum(like2(d$to,d$do,d$x,pars,phi=phi,theta=theta))}
}

## MAKES THE LIKELIHOOD (NEW METHOD) ##################################################
log.like3 <- function(to,do,x,gamma,phi,theta){
  sum(log(like.part1.NB3(to=to,do=do,x=x,phi=phi,gamma=gamma)+like.part2(to=to,do=do,x=x,theta=theta,phi=phi,gamma=gamma)))
}

like3 <- function(to,do,x,gamma,phi,theta){
  like.part1.NB3(to=to,do=do,x=x,phi=phi,gamma=gamma)+like.part2(to=to,do=do,x=x,theta=theta,phi=phi,gamma=gamma)
}

## THE FIRST PART OF THE LIKELIHOOD VECTORISED IN to ##################################
like.part1.NB3 <- function(to,do,x,phi,gamma){
  matrix.prod = outer(exp(x), (1+exp(gamma)), function(x,y) y^-x)
  matrix.i = matrix(1:length(gamma), nrow=length(to), ncol=length(gamma), byrow=TRUE)
  matrix.prod[matrix.i>to]=1
  ToSum = apply(matrix.prod,1,prod)*Gamma.Vec(to=to,do=do,phi=phi)
  ToSum
}

## THE SECOND PART OF THE LIKELIHOOD ##################################################
like.part2 <- function(to,do,x,theta,phi,gamma){
  one.k.like = one.to.like = numeric()
  for(m in 1:length(to)){
    for(k in 1:to[m]){
      one.k.like[k] = one.like.part2(k=k,to=to[m],do=do[m],x=x[m],theta=theta,phi=phi,gamma=gamma)
    }
    one.to.like[m] = sum(one.k.like)
    one.k.like = numeric()
  }
  one.to.like
}

## USED FOR EACH k IN 1:to ############################################################
one.like.part2 <- function(k,to,do,x,theta,phi,gamma){
  product = (1+exp(gamma))^(-exp(x))
  product2 = 1 - product[k]

  ifelse(k != 1, prod(product[1:(k - 1)]) * (product2) * Delta.Vec(ti = k, 
           to = to, do = do, phi = phi, theta = theta), (product2) * 
         Delta.Vec(ti = k, to = to, do = do, phi = phi, theta = theta))
}

makedev.log.like3 <- function(d, phi, theta){
	function(pars){-sum(log.like3(d$to,d$do,d$x,pars,phi=phi,theta=theta))}
}

makedev.like3 <- function(d, phi, theta){
	function(pars){sum(like3(d$to,d$do,d$x,pars,phi=phi,theta=theta))}
}

makePhiTheta <- function(Parameters,error){
  Parameters = c(unlist(Parameters[c(1:2,5:7)]), times = list(1:16/4), error.sd = error)
  
  data.trt = replicate(nreps.phi.theta, with(Parameters, make.one.BP(startrange = startrange,
    trend = as.numeric(trend),
    error.sd = error.sd,
    treatment = 0,
    times = as.numeric(unlist(times)),
    end.of.treatment = as.numeric(end.of.treatment),
    carryover = as.numeric(carryover))))
  
  data.ctrl = replicate(nreps.phi.theta, with(Parameters, make.one.BP(startrange = startrange,
    trend = as.numeric(trend),
    error.sd = error.sd,
    treatment = 0,
    times = as.numeric(unlist(times)),
    end.of.treatment = as.numeric(end.of.treatment),
    carryover = as.numeric(carryover))))
  
  PT = findPhiTheta(data = cbind(data.trt, data.ctrl), times = with(Parameters,times), rule = get(Parameters$rule))

  PT
}


############################################################################################ CODING TO USE NLM MAXIMISER INPUT:THETA OUTPUT:FUNCTION VALUE & DERIVATIVE ################################ METHOD 1 ###################################################
make.deriv.method1 <- function(d,phi,theta,times.length){
  first.deriv = list()
  for(i in 1:times.length){
    first.deriv[[i]] = function(pars){deriv.log.like.NB1(to=d$to,do=d$do,gamma=pars,phi=phi,theta=theta,which.gamma=i)}
  }
  function(pars){
     answer = numeric(times.length)
     for(i in 1:times.length) answer[i] = sum(first.deriv[[i]](pars))
     answer
  }
}

method1 <- function(data,phi,theta,times.length,initial.vals){
 the.minusloglike = makedev.log.like1(d=data,phi=phi,theta=theta)
 
 to.min <- function(pars){
   value = the.minusloglike(pars) 
   value
 }
 nlm(f = to.min, p = initial.vals)
}

############################################################################################ CODING TO USE NLM MAXIMISER INPUT:THETA OUTPUT:FUNCTION VALUE & DERIVATIVE ################################ METHOD 2 ###################################################
make.deriv.method2 <- function(d,phi,theta,times.length){
  first.deriv = list()
  for(i in 1:times.length){
    first.deriv[[i]] = function(pars){deriv.log.like.NB2(to=d$to,do=d$do,gamma=pars,phi=phi,theta=theta,which.gamma=i)}
  }
   function(pars){
     answer = numeric(times.length)
     for(i in 1:times.length) answer[i] = sum(first.deriv[[i]](pars))
     answer
  }
}

method2 <- function(data,phi,theta,times.length,initial.vals){
 the.minusloglike = make.log.like.Vec(d=data,phi=phi,theta=theta)

 to.min <- function(pars){
   value = the.minusloglike(pars)
   value
 }

 nlm(f = to.min, p = initial.vals)
}


method3 <- function(data,phi,theta,times.length,initial.vals){
 the.minusloglike = makedev.log.like2(d=data,phi=phi,theta=theta)

 to.min <- function(pars){
   value = the.minusloglike(pars)
   value
 }

 nlm(f = to.min, p = initial.vals)
}

method4 <- function(data,phi,theta,times.length,initial.vals){
 the.minusloglike = makedev.log.like3(d=data,phi=phi,theta=theta)

 to.min <- function(pars){
   value = the.minusloglike(pars)
   value
 }

 nlm(f = to.min, p = initial.vals)
}
#################################################################################################### CODING THE DERIVATIVE LOG(LIKE) WRT GAMMMA - RIGHT #############################################################################################################

deriv.log.like.1.new <- function(x,to,do,phi,theta,gamma,which.gamma){
  recip.like = (1/like.nonconstant(x=x,to=to,do=do,gamma=gamma,phi=phi,theta=theta))
  constant.prod = -exp(x + gamma[which.gamma])*(1+exp(gamma[which.gamma]))^(-(exp(x)+1))

  product = matrix(NA, nrow = length(x), ncol = length(gamma))
  solution = numeric(length(to))

  for(j in 1:length(gamma)){
    product[,j] = like0.deriv(x=x,gamma=gamma[j])
  }
  ## Removing the column where we are taking the derivative
  product[,which.gamma] = 1
  ## If you are diagnosed before than you add nothing to the derivative
  product[which(to<which.gamma),]=0

  To.Sum = numeric(length(to))
  for(i in 1:length(to)){
    Sum.for.i = numeric(to[i])
    for(k in 1:to[i])
      Sum.for.i[k] = deriv.likek.nonconstant.new(to=to[i],do=do[i],x=x[i],gamma=gamma,k=k,phi=phi,theta=theta,which.gamma=which.gamma)*Delta.Vec(ti=k,to=to[i],do=do[i],phi=phi,theta=theta)
    To.Sum[i] = sum(Sum.for.i)
  }
  
  for(i in 1:length(to)){
    solution[i] = Gamma.Vec(to=to[i],do=do[i],phi=phi)*prod(product[i,1:to[i]])*constant.prod[i] + To.Sum[i]
  }
  solution*recip.like   
}

deriv.likek.nonconstant.new <- function(to,do,x,gamma,k,phi,theta, which.gamma){
  risk<-exp(x)
  haz<-1+exp(gamma)
  To.Prod = haz^-risk
  deriv = -exp(x + gamma[which.gamma])*
    (1 + exp(gamma[which.gamma]))^(-(exp(x) + 1))

  if(k < which.gamma) solution = 0
  if(k == which.gamma){
    if(k == 1){solution = -deriv
             }else{
               solution = -deriv*prod(To.Prod[1:(k-1)])}
  }
  if(k > which.gamma)
    solution = deriv*prod(To.Prod[1:(k-1)][-which.gamma])*(1-To.Prod[k])
  
  solution
}

like0.deriv<-function(x,gamma){
  risk<-exp(x)
  haz<-1+exp(gamma)
  haz^-risk
}

## Phi.Theta = list()
## Phi.Theta = makePhiTheta(par[r,])
## Phi = Phi.Theta$Phi
## Theta = Phi.Theta$Theta


########################################################################################################## To run on the cluster ####################################################################################################################################

jobid = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
set.seed(jobid)
##par[round(runif(1,1,nrow(par))),]
r = jobid

estimates = list()
error = seq(0,7,by=.5)

for(j in 1:length(error)){
  PT = makePhiTheta(Parameters=par[r,], error=error[j])
  estimates[[j]] = find.one.bias.NB(Parameters=par[r,], error=error[j], PT.Tr=PT)
  print(j)
}

filename = paste0("5June",jobid,".rda")
##filename = "1June.rda"

save(estimates,file=filename)

gamma.optim = c(-3.001855, -3.204178, -1.892595, -2.421256)
phi= rep(.9,4)
theta=rep(.8,4)

################################################################################################## CODING TO CHECK THE NUMERIC DERIVATIVE AGAINST THE ANALYTIC DERIV ########################################### WRT GAMMA ##########################################
to = rep(4,5)
x = rep(0,5)
do = rep(1,5)
d = data.frame(to,x,do)

to = c(rep(1,10),rep(2,10),rep(3,10),rep(4,10))
x = c(rep(1,10),rep(0,10),rep(1,10),rep(0,10))
do = c(rep(1,10),rep(0,10),rep(1,10),rep(0,10))
d = data.frame(to,x,do)

table(like1(to=d$to,do=d$do,x=d$x,gamma=gamma.optim,phi=phi,theta=theta))
table(like2(to=d$to,do=d$do,x=d$x,gamma=gamma.optim,phi=phi,theta=theta))
table(like3(to=d$to,do=d$do,x=d$x,gamma=gamma.optim,phi=phi,theta=theta))
table(like.Vec(to=d$to,do=d$do,x=d$x,gamma=gamma.optim,phi=phi,theta=theta))

## The numeric derivative 
minus.log.like1 = make.log.like.Vec(d=d, phi=phi, theta=theta)
-(minus.log.like1(gamma.optim + c(.00001,0,0,0)) - minus.log.like1(gamma.optim))/.00001

minus.log.like2 = makedev.log.like1(d=d, phi=phi, theta=theta)
-(minus.log.like2(gamma.optim + c(.00001,0,0,0)) - minus.log.like2(gamma.optim))/.00001

minus.log.like3 = makedev.log.like2(d=d, phi=phi, theta=theta)
-(minus.log.like3(gamma.optim + c(.00001,0,0,0)) - minus.log.like3(gamma.optim))/.00001

minus.log.like4 = makedev.log.like3(d=d, phi=phi, theta=theta)
-(minus.log.like4(gamma.optim + c(.00001,0,0,0)) - minus.log.like4(gamma.optim))/.00001

################################################################################################################################################ CODING TO ANAYLSE RESULTS ########################################################################################################################################################################

True = gamma.hats = matrix(NA, nrow = length(estimates), ncol = length(estimates[[1]]$Phi))


for(i in 1:length(estimates)){
  True.Gamma[i,] = estimates[[i]]$'True HZD'
  gamma.hats[i,] = exp(estimates[[i]]$Method2.NLM)
}

pdf("Bias.pdf")
plot(y=gamma.hats, x = True.Gamma,  main = "True.HZD vs. exp(Optim)", ylab = "Gamma.Hats", xlab = "True")
abline(a=0,b=1)
dev.off()

################################################################################################################ Using RSTAN to put a prior on the hazards #################################################################################################################################################################################################

library("rootSolve")
## this can be used to find the inverse of a function
inverse.all = function (f, lower = 0, upper = 16) {
  function (y) uniroot.all((function (x) f(x) - y), lower = lower,
    upper = upper, tol = .Machine$double.eps^0.4,maxiter = 10000, n = 10000)
}

make.like.Inv <- function(do, x, gamma, phi, theta){
  function(to) like3(to, do = do, x = x, gamma = gamma, 
                        phi = phi, theta = theta)
}

find.like.Inv = make.like.Inv(do = 1, x = 1, gamma = gamma.optim, phi = phi, theta = theta) ## A function of to only
like.Inv = inverse.all(find.like.Inv) 


try = like3(to = 3, do = 1, x = 1, gamma = gamma.optim, phi = phi, 
           theta = theta)
like.Inv(try) ## It works!! This should be to from the input above!


