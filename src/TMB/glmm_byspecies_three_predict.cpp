#include <TMB.hpp>
  
  template<class Type>
  Type objective_function<Type>::operator() ()
  {
    DATA_INTEGER(n);
    DATA_INTEGER(nspc);
    DATA_VECTOR(y);
    DATA_MATRIX(U1);
    DATA_MATRIX(U2);
    DATA_MATRIX(U3);
    DATA_MATRIX(Xlat);
    DATA_MATRIX(Xlength);
    DATA_MATRIX(Xint);
    DATA_IVECTOR(spcindex);  // index of species IDs
    DATA_MATRIX(Z); // fishIDs
    
    PARAMETER_VECTOR(bo); // intercepts for fish / parasite combos
    PARAMETER(bobar);
    PARAMETER(logsigma_bo);
    PARAMETER_VECTOR(blat); // fixed effects
    PARAMETER(blatbar);
    PARAMETER(logsigma_blat);
    PARAMETER_VECTOR(blength); // fixed effects
    PARAMETER(blengthbar);
    PARAMETER(logsigma_blength);
    PARAMETER_VECTOR(g1); // random effects of species
    PARAMETER(gbar1); // mean of time effects
    PARAMETER(logsigma_g1); // variance of species level random effects
    PARAMETER_VECTOR(g2); // random effects of species
    PARAMETER(gbar2); // mean of time effects
    PARAMETER(logsigma_g2); // variance of species level random effects
    PARAMETER_VECTOR(g3); // random effects of species
    PARAMETER(gbar3); // mean of time effects
    PARAMETER(logsigma_g3); // variance of species level random effects
    PARAMETER_VECTOR(z); // fish ID random effects
    PARAMETER(logsigma_z);
    PARAMETER_VECTOR(logphi); // negative binomial dispersion parameters
    PARAMETER(logphibar);
    PARAMETER(loglogsigma_phi);  
    
    Type jnll = 0.0;
    Type sigma_g1 = exp(logsigma_g1);
    Type sigma_g2 = exp(logsigma_g2);
    Type sigma_g3 = exp(logsigma_g3);
    Type sigma_z = exp(logsigma_z);
    Type sigma_bo = exp(logsigma_bo);
    Type sigma_blat = exp(logsigma_blat);
    Type sigma_blength = exp(logsigma_blength);
    Type logsigma_phi= exp(loglogsigma_phi);
    
    ADREPORT(sigma_g1);
    ADREPORT(sigma_g2);
    ADREPORT(sigma_g3);
    ADREPORT(sigma_z);
    ADREPORT(sigma_bo);
    ADREPORT(sigma_blat);
    ADREPORT(sigma_blength);
    ADREPORT(logsigma_phi);
    
    jnll -= sum(dnorm(z,0,sigma_z, true));
    jnll -= sum(dnorm(g1, gbar1, sigma_g1, true)); 
    jnll -= sum(dnorm(g2, gbar2, sigma_g2, true)); 
    jnll -= sum(dnorm(g3, gbar3, sigma_g3, true)); 
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
    mu = exp(Xint * bo + Xlat*blat + Xlength * blength + U1 * g1 + U2 * g2 + U3 * g3 + Z * z) ;
    
    Type nll_data = 0;
    for(int i = 0; i < n; i++){ 
      Type var = mu[i] + pow(mu[i],2) * pow (phi[spcindex[i]], -1);
      jnll -=  	dnbinom2(y[i], mu[i],var, true);
      nll_data -= dnbinom2(y[i], mu[i],var, true);
    }
    REPORT(nll_data)
    
    return jnll;
  }
  
