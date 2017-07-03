functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(real to, real dO, vector phi);

real Delta(int ti,  real to, real dO, vector phi, vector theta);

real like_lpdf(vector TDX, real beta, vector gamma, vector phi, vector theta);

real like0(vector TDX, real beta, vector gamma, vector phi);

real like.part2(vector TDX, real beta, vector gamma, vector phi, vector theta);

real one.like.part2(int ti, real to, real dO, real x, real beta, vector gamma, vector phi, vector theta);



##################
## Definitions 
##################
real Delta(int ti, real to, real dO, vector phi, vector theta){
  real T;
  real P1;
  int i;
  
  P1 = 1;
  i = ti;
   
  while(ti <= i <= (to-1)){
    P1 = P1*(1-theta[i]);
    i = i+1;
    }
  
  return prod(phi[1:(ti-1)])*P1*(1-theta[static_cast<int>(to)])^(1-dO)*theta[static_cast<int>(to)]^dO;
    }

real Gamma(real to, real dO, vector phi){
    return prod(phi[1:(to-1)])*phi[to]^(1-dO)*(1-phi[to])^dO;
    }

real like0(vector TDX, real beta, vector gamma, vector phi){
    real to;
    real dO;
    real x;
    real P;
    int i;

    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];
    P = 1;
    i = 1;
    
    while(i <= to){
     P = P * (1 + exp(gamma[i]))^(-exp(x*beta));
     i = i + 1;
    }

    return P*Gamma(to, dO, phi);
    }
    
// THE SECOND PART OF THE LIKELIHOOD ###########################################
real like.part2(vector TDX, real beta, vector gamma, vector phi, vector theta){
  real to;
  real dO;
  real x;
  real k;
  vector [16] klike; // number of measurement times
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  k = 1;
  
    while (k <= to){
      klike[k] = one.like.part2(k, to, dO, x, beta, gamma, phi, theta);
      k = k + 1;
    }
  }
  return sum(klike);
}

// USED FOR EACH k IN 1:to #####################################################
real one.like.part2(int ti, real to, real dO, real x, real beta, vector gamma, vector phi, vector theta){
  real to;
  real dO;
  real x;
  real P2;
  real i;
  vector [to] P;
  
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  P = 1;
  i = 1;
    
    while(i <= to){
      P = P * (1 + exp(gamma[i]))^(-exp(x*beta))
      i = i + 1
    }
     
  P2 = 1 - P[ti];

  if(ti != 1) 
  {
    return prod(P[1:(ti - 1)]) * (P2) * Delta(ti = ti, 
           to = to, dO = dO, phi = phi, theta = theta);
  }
  else
  {
    return (P2) * Delta(ti = ti, to = to, dO = dO, phi = phi, theta = theta));
  }
}


real like_lpdf(vector TDX, real beta, vector gamma, vector phi, vector theta){
    real to;
    real dO;
    real x;
    vector [6] S;
    
    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    S[1] = like0(TDX, beta, gamma, phi);
    
    for (j  in 1:5)
    {
      S[j+1]= likek(TDX, j, beta, gamma, phi, theta);
    }
    
    return sum(S);
    }
}


data {

  int<lower=0> N;
  vector[2] TDX[N];

}

parameters {

  vector [16] <lower=-20,upper=20> gamma;
  real<lower=-20,upper=20> beta;

}

model {
  
  for(i in 1:N){
    
    TDX[i] ~ like_lpdf(beta, gamma, 0.9, 0.8);
    
  }
  
  beta ~ normal(0,10);
  gamma ~ normal(0,10);

}
