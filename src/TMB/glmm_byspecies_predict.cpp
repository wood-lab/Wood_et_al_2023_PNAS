#include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_INTEGER(n);
    DATA_INTEGER(nspc);
    DATA_VECTOR(y);
    DATA_MATRIX(U);
    DATA_MATRIX(Xlat);
    DATA_MATRIX(Xlength);
    DATA_MATRIX(Xint);
    DATA_IVECTOR(Xspc);  // index of species IDs
    DATA_MATRIX(Z); // fishIDs
    DATA_VECTOR(Upredict);
    DATA_INTEGER(npredict);
    
    PARAMETER_VECTOR(bo); // intercepts for fish / parasite combos
    PARAMETER(bobar);
    PARAMETER(logsigma_bo);
    PARAMETER_VECTOR(blat); // fixed effects
    PARAMETER(blatbar);
    PARAMETER(logsigma_blat);
    PARAMETER_VECTOR(blength); // fixed effects
    PARAMETER(blengthbar);
    PARAMETER(logsigma_blength);
    PARAMETER_VECTOR(g); // random effects of species
    PARAMETER(gbar); // mean of time effects
    PARAMETER(logsigma_g); // variance of species level random effects
    PARAMETER_VECTOR(z); // fish ID random effects
    PARAMETER(logsigma_z);
    PARAMETER_VECTOR(logphi); // negative binomial dispersion parameters
    PARAMETER(logphibar);
    PARAMETER(loglogsigma_phi);
    
    Type sigma_g = exp(logsigma_g);
    Type sigma_z= exp(logsigma_z);
    Type sigma_bo = exp(logsigma_bo);
    Type sigma_blat = exp(logsigma_blat);
    Type sigma_blength = exp(logsigma_blength);
    Type logsigma_phi= exp(loglogsigma_phi);
    
    Type jnll = 0;
    ADREPORT(sigma_g);
    ADREPORT(sigma_z);
    ADREPORT(sigma_bo);
    ADREPORT(sigma_blat);
    ADREPORT(sigma_blength);
    ADREPORT(logsigma_phi);
    
    // random effects likelihood
    jnll -= sum(dnorm(z,0,sigma_z, true));
    jnll -= sum(dnorm(g, gbar, sigma_g, true)); 
    jnll -= sum(dnorm(bo, bobar, sigma_bo, true)); 
    jnll -= sum(dnorm(blat, blatbar, sigma_blat, true)); 
    jnll -= sum(dnorm(blength, blengthbar, sigma_blength, true)); 
    jnll -= sum(dnorm(logphi, logphibar, logsigma_phi, true)); 
    
    vector<Type> phi(nspc);
    for(int i = 0; i < nspc; i++){ 
      phi[i]= exp(logphi[i]);
    }
    ADREPORT(phi);

    vector<Type>  mu(n);
    mu = exp( Xint * bo + Xlat*blat + Xlength *blength + U * g + Xint * bo + Z * z) ;
    
    
    for(int i = 0; i < n; i++){ 
      Type var = mu[i] + pow(mu[i],2) * pow (phi[Xspc[i]], -1);
      jnll -=  	dnbinom2(y[i], mu[i],var, true);
    }
    
    // Now get predicted values
    vector<Type>  logmu(npredict);
    for(int i = 0; i < npredict; i++){ 
      //logmu[i] =   Upredict[i] * g[parindex] + bo[hostparindex];
      logmu[i] =   Upredict[i] * gbar + bobar;
    }
    ADREPORT(logmu);
    
    return jnll;
  }
  
