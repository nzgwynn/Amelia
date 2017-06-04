functions {

##################
## Declarations 
##################

# Have to declare the functions before I can use them


real log_h(real T,  real B, real theta);

real H_t(real T,  real B, real theta);

##################
## Definitions 
##################

real surv_dens_log(vector T_and_delta,real B, real theta){
    real T;
    real delta;

    T     <- T_and_delta[1];
    delta <- T_and_delta[2];

    return (delta * log_h(T, B, theta))
                 -    H_t(T, B, theta) ;

    }

real log_h(real T,  real B, real theta){

      return log(B) + theta*T;

}
real H_t(real T,  real B, real theta){

      return (B/theta)*( exp(theta*T)-1 ) ;

}

}

data {

 

}
parameters {




}
model {

  

}
