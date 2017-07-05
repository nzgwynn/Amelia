functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(vector phi, int to, int dO);

real Delta(vector phi, vector theta, int ti,  int to, int dO);

real like_lpdf(real beta, vector gamma, vector phi, vector theta, int TDX);

real like0(real beta, vector gamma, vector phi, int TDX);

real like.part2(real beta, vector gamma, vector phi, vector theta, int TDX);

real one.like.part2(real beta, vector gamma, int ti, int to, int dO, int x, vector phi, vector theta);



##################
## Definitions 
##################
real  Delta(vector phi, vector theta, int ti,  int to, int dO){
  real T;
  real P1;
  
  P1 = 1;
   
  for(i in ti:(to-1)){
    P1 = P1*(1-theta[i]);
    }
  
  return prod(phi[1:(ti-1)])*P1*(1-theta[to])^(1-dO)*theta[to]^dO;
    }

real Gamma(vector phi, int to, int dO){
    return prod(phi[1:(to-1)])*phi[to]^(1-dO)*(1-phi[to])^dO;
    }

real like0(real beta, vector gamma, vector phi, int TDX){
    int to;
    int dO;
    int x;
    real P;

    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];
    P = 1;
    i = 1;
    
    for(i in 1:to){
     P = P * (1 + exp(gamma[i]))^(-exp(x*beta));
    }

    return P*Gamma(phi, to, dO);
    }
    
// THE SECOND PART OF THE LIKELIHOOD ###########################################
real like.part2(real beta, vector gamma, vector phi, vector theta, int TDX){
  int to;
  int dO;
  int x;
  int k;
  vector [4] klike; // number of measurement times
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  
    for(k in 1:to)
      klike[k] = one.like.part2(beta, gamma, k, to, dO, x, phi, theta);
      
  return sum(klike);
}

// USED FOR EACH k IN 1:to #####################################################
real one.like.part2(real beta, vector gamma, int ti, int to, int dO, int x, vector phi, vector theta){
  int to;
  int dO;
  int x;
  real P2;
  real i;
  vector[to] P;
  
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  P = 1;
    
    for(i in 1:to)
      P = P * (1 + exp(gamma[i]))^(-exp(x*beta))
     
  P2 = 1 - (1 + exp(gamma[ti]))^(-exp(x*beta));

  if(ti != 1) 
  {
    return prod(P[1:(ti - 1)]) * (P2) * Delta(phi, theta, ti, to, dO);
  }
  else
  {
    return (P2) * Delta(phi, theta, ti, to, dO));
  }
}


real like_lpdf(real beta, vector gamma, vector phi, vector theta, int TDX){
    int to;
    int dO;
    int x;
    vector [6] S;
    
    to = TDX[1];
    dO = TDX[2];
    x = TDX[3];

    S[1] = like0(beta, gamma, phi, TDX);
    
    for (j  in 1:5)
      S[j+1]= likek(TDX, j, beta, gamma, phi, theta);
    
    return sum(S);
}


data {

  int<lower=0> N;
  int TDX[3];

}

parameters {

  vector [4] <lower=-20,upper=20> gamma;
  real<lower=-20,upper=20> beta;

}

model {
  
  for(i in 1:N){
    
    TDX[i] ~ like_lpdf(beta, gamma, 0.9, 0.8);
    
  }
  
  beta ~ normal(0,10);
  gamma ~ normal(0,10);

}
