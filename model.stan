functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real Gamma(real to,  real dO);

real Delta(real ti,  real to, real dO);

real like(real to,  real dO, real x, real beta, real gamma);

real like0(real to,  real dO, real x, real beta, real gamma);

real likek(real to,  real dO, real x, real beta, real gamma, real k);

##################
## Definitions 
##################

real like(vector TDX, real beta, real gamma){
    real to;
    real dO;
    real x;
    real beta;
    real gamma;

    to <- TDX[1];
    dO <- TDX[2];
    x <- TDX[3];

    return like0(to,dO,x,beta,gamma) + likek(to,dO,x,beta,gamma,k=1) +
		likek(to,dO,x,beta,gamma,k=2) + likek(to,dO,x,beta,gamma,k=3) + 
		likek(to,dO,x,beta,gamma,k=4) + likek(to,dO,x,beta,gamma,k=5);

    }

real like0(vector TDX, real beta, real gamma){
    real to;
    real dO;
    real x;
    real beta;
    real gamma;

    to <- TDX[1];
    dO <- TDX[2];
    x <- TDX[3];

    return (1+exp(gamma))^(-exp(x*beta))*Gamma(to, dO);
    }
    
real likek(vector TDX, real beta, real gamma, real k){
    real to;
    real dO;
    real x;
    real beta;
    real gamma;
    real k;

    to     <- TDX[1];
    dO <- TDX[2];
    x <- TDX[3];

    return
    REWRITE IN ~R LANG
    ifelse(k<=to,((haz^-risk)^(k-1))*(1-haz^-risk)*Delta(k,to,dO),0);

    }
}

data {

 

}
parameters {




}
model {

  

}
