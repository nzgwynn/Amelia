functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(real to, real dO, real phi);

real Delta(real ti,  real to, real dO, real phi, real theta);

real like_lpdf(vector TDX, real beta, real gamma, real phi, real theta);

real like0(vector TDX, real beta, real gamma, real phi);

real likek(vector TDX, real k, real beta, real gamma, real phi, real theta);

##################
## Definitions 
##################
real Delta(real ti, real to, real dO, real phi, real theta){
    return phi^(ti-1)*(1-theta)^(to-ti)*(1-theta)^(1-dO)*theta^dO;
    }

real Gamma(real to, real dO, real phi){
    return phi^(to-1)*phi^(1-dO)*(1-phi)^dO;
    }

real like0(vector TDX, real beta, real gamma, real phi){
    real to;
    real dO;
    real x;
    real risk;
    real haz;

    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    return (1+exp(gamma))^(-exp(x*beta))*Gamma(to, dO, phi);
    }
    
real likek(vector TDX, real k, real beta, real gamma, real phi, real theta){
    real to;
    real dO;
    real x;
    real risk;
    real haz;

    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    if (k <= to) 
    {
      risk = exp(x*beta);
      haz = 1+exp(gamma);
      return pow(pow(haz, -risk), k-1) * (1-pow(haz, -risk)) * Delta(k, to, dO, phi, theta);
     } 
     else 
     {
      return 0;
     }
}

real like_lpdf(vector TDX, real beta, real gamma, real phi, real theta){
    real to;
    real dO;
    real x;
    vector [6] S; // Upate this here too!
    
    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    S[1] = like0(TDX, beta, gamma, phi); //TL tried to do log(...)
    
    for (j  in 1:5) //TL not sure if this is correct, just one TDX at at time.
    {
      S[j+1]= likek(TDX, j, beta, gamma, phi, theta); //TL tried to do log(...) didn't like it at all.
    }
    
    return sum(S); //TL this keeps giving me log(0) and gets fussy.
    }
}


data {

  int N; // TL tried sample size 1000 in trt & ctl still no
  vector[3] TDX[N];

}

parameters {

  real gamma;
  real beta;

}

model {
  
  for(i in 1:N){
    
    TDX[i] ~ like_lpdf(beta, gamma, 0.9, 0.8);
    
  }
  
  beta ~ normal(1.3, 10);
  gamma ~ normal(-2.9, 10);

}
