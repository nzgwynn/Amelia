functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(int to, real dO, vector phi);

real Delta(int ti,  int to, real dO, vector phi, vector theta);

real litie(vector TDX, real beta, vector gamma, vector phi, vector theta);

real like0(vector TDX, real beta, vector gamma, vector phi);

real like.part2(vector TDX, real beta, vector gamma,  vector phi, vector theta);

real one.like.part2(int ti, int to, real dO, real x, real beta, vector gamma, vector phi, vector theta);



##################
## Definitions 
##################
real Delta(int ti, int to, real dO, vector phi, vector theta){
  real T;
  real P1;
  
    P1 = 1;
   
    for(i in ti:(to-1))
      P1 = P1*(1-theta[i]);
  
    return prod(phi[1:(ti-1)])*P1*(1-theta[to])^(1-dO)*theta[to]^dO;
    }

real Gamma(int to, real dO, vector phi){
    return prod(phi[1:(to-1)])*phi[to]^(1-dO)*(1-phi[to])^dO;
    }

real like0(vector TDX, real beta, vector gamma, vector phi){
    int to;
    real dO;
    real x;
    real P;

    to = floor(TDX[1]);
    dO = TDX[2];
    x = TDX[3];
    P = 1
    
    for(i in 1:to)
     P = P * (1 + exp(gamma[i]))^(-exp(x*beta));

    return P*Gamma(to, dO, phi);
    }
    
// THE SECOND PART OF THE LIKELIHOOD ###########################################
real like.part2(vector TDX, real beta, vector gamma, vector phi, vector theta){
  int to;
  real dO;
  real x;
  vector one.k.like; 
  vector one.to.like;
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  
  for(m in 1:length(to)){
    for(k in 1:to[m]){
      one.k.like[k] = one.like.part2(k, to, dO, x, beta, gamma, phi, theta)
    }
    one.to.like[m] = sum(one.k.like)
    one.k.like = numeric()
  }
  return one.to.like;
}

// USED FOR EACH k IN 1:to #####################################################
real one.like.part2(int ti, int to, real dO, real x, real beta, vector gamma, vector phi, vector theta){
  int to;
  real dO;
  real x;
  vector product;
  vector product 2;
  
  to = TDX[1];
  dO = TDX[2];
  x = TDX[3];
  
  product = (1+exp(gamma))^(-exp(x*beta));
  product2 = 1 - product[ti];

  if(ti != 1) 
  {
    return prod(product[1:(ti - 1)]) * (product2) * Delta(ti = ti, 
           to = to, dO = dO, phi = phi, theta = theta);
  }
  else
  {
    return (product2) * Delta(ti = ti, to = to, dO = dO, phi = phi, theta = theta));
  }
}


real like(vector TDX, real beta, real gamma, real phi, real theta){
    int to;
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
    
    TDX[i] ~ like(beta, gamma, 0.9, 0.8);
    
  }
  
  beta ~ normal(0,10);
  gamma ~ normal(0,10);

}
