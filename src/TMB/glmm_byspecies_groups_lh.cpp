#include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_INTEGER(n);
    DATA_INTEGER(ngrp);
    DATA_INTEGER(nspc);
    DATA_VECTOR(y);
    DATA_MATRIX(U);
    DATA_MATRIX(X);
    DATA_MATRIX(Xint);
    DATA_IVECTOR(Xspc);  // index of species IDs
    DATA_MATRIX(A);
    DATA_IVECTOR(LH); // index of which life history strategy for each species
    DATA_IVECTOR(taxgroup); // index of taxonomic group for each species
    DATA_MATRIX(Z); // fishIDs
    
    PARAMETER_VECTOR(bo); // intercepts for fish / parasite combos
    PARAMETER_VECTOR(b); // fixed effects
    PARAMETER_VECTOR(g); // random effects of time by species
    PARAMETER_VECTOR(araw); // random effects of time by species group
    PARAMETER_VECTOR(abar); // group level mean effects
    PARAMETER(logsigma_g); // variance of random effects within groups
    PARAMETER(logsigma_a); // variance of random effects across groups
    PARAMETER_VECTOR(z); // random effects of fish ID
    PARAMETER(logsigma_z); // variance of fish ID random effects
    PARAMETER_VECTOR(logphi); // negative binomial dispersion parameters by species
   
    
    Type jnll = 0.0;
    Type sigma_g = exp(logsigma_g);
    Type sigma_a = exp(logsigma_a);
    Type sigma_z = exp(logsigma_z);
    
    jnll -= sum(dnorm(z, 0, sigma_z, true)); 
    
    vector<Type> phi(nspc);
    for(int i = 0; i < nspc; i++){ 
      phi[i]= exp(logphi[i]);
    }
    ADREPORT(phi);
    ADREPORT(sigma_g);
    ADREPORT(sigma_a);
    ADREPORT(sigma_z);
    
    
    vector<Type> a(ngrp);
    // get species group effects on slopes
    for(int i = 0; i < ngrp; i++){ 
      a[i] = abar[LH[i]] + araw[i] * sigma_a;
    }
    ADREPORT(a);
    // get the overall time effect for each species
    vector<Type> spc_slopes(nspc);
    for(int i = 0; i < nspc; i++){ 
      spc_slopes[i] = g[i] + a[taxgroup[i]];
    }
    ADREPORT(spc_slopes);
    
    vector<Type>  mu(n);
    mu = exp(X*b + U * g + Xint * bo + A * a + Z * z) ;
    
    // for integrating over species level effects
    jnll -= sum(dnorm(g, 0, sigma_g, true)); 
    
    // for integrating over group level effects
    jnll -= sum(dnorm(araw, 0, 1, true)); 
    
    for(int i = 0; i < n; i++){ 
      Type var = mu[i] + pow(mu[i],2) * pow (phi[Xspc[i]], -1);
      jnll -=  	dnbinom2(y[i], mu[i],var, true);
    }
    
    return jnll;
  }
  
