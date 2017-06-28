functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(vector TDX, real phi);

real Delta(real ti,  vector TDX, real phi, real theta);

real like(vector TDX, real beta, real gamma, real phi, real theta);

real like0(vector TDX, real beta, real gamma, real phi);

real likek(vector TDX, real k, real beta, real gamma, real phi, real theta);

##################
## Definitions 
##################
real Delta(real ti, vector TDX, real phi, real theta){
    real to;
    real dO;
    
    to = TDX[1];
    dO = TDX[2];
  
    return phi^(ti-1)*(1-theta)^(to-ti)*(1-theta)^(1-dO)*theta^dO;
    }

real Gamma(vector TDX, real phi){
    real to;
    real dO;
 
    to = TDX[1];
    dO = TDX[2];
  
    return phi^(to-1)*phi^(1-dO)*(1-phi)^dO;
    }

real like0(vector TDX, real beta, real gamma, real phi){
    real to;
    real dO;
    real x;

    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    return (1+exp(gamma))^(-exp(x*beta))*Gamma(TDX, phi);
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
      return pow(pow(haz, -risk), k-1) * (1-pow(haz, -risk)) * Delta(k, TDX, phi, theta);
     } 
     else 
     {
      return 0;
     }
}

real like(vector TDX, real beta, real gamma, real phi, real theta){
    real to;
    real dO;
    real x;
    real S;
    
    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    S = like0(TDX, beta, gamma, phi);
    
    //for(int i = 5; i >= 1; --i){
      //S += likek(TDX, i, beta, gamma, phi, theta);
    //}
    
    return S;
    }
}


data {

  int<lower=0> N;
  vector[2] TDX[N];

}

parameters {

  real<lower=0> gamma;
  real<lower=0> beta;

}

model {
  
  for(i in 1:N){
    
    TDX[i] ~ like(beta, gamma, 0.9, 0.8);
    
  }
  
  beta ~ normal(0,1);
  gamma ~ normal(0,1);

}
