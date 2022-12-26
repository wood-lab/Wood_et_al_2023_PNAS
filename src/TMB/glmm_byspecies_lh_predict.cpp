#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(n);
  DATA_INTEGER(nspc);
  DATA_VECTOR(y);
  DATA_MATRIX(U);
  DATA_MATRIX(Xlat)
  DATA_MATRIX(Xlength);
  DATA_MATRIX(Xg);
  DATA_MATRIX(Xint);
  DATA_IVECTOR(LH); // index of which life history strategy for each species
  DATA_MATRIX(Z); // fishIDs
  DATA_IVECTOR(Xspc);  // index of species IDs
  DATA_VECTOR(Upredict); // years to predict
  DATA_INTEGER(npredict); // how many predictions to make
  
  
  PARAMETER(bobar);
  PARAMETER_VECTOR(bo); // intercepts for fish / parasite combos
  PARAMETER(logsigma_bo);
  PARAMETER_VECTOR(bog);
  PARAMETER(logsigma_bog);
  PARAMETER_VECTOR(blat); // fixed effects
  PARAMETER(blatbar);
  PARAMETER(logsigma_blat);
  PARAMETER_VECTOR(blength); // fixed effects
  PARAMETER(blengthbar);
  PARAMETER(logsigma_blength);    
  PARAMETER_VECTOR(graw); // random effects of time by species
  PARAMETER_VECTOR(gbar); // group level mean effects
  PARAMETER(logsigma_g); // variance of random effects within groups
  PARAMETER_VECTOR(z); // random effects of fish ID
  PARAMETER(logsigma_z); // variance of fish ID random effects
  PARAMETER_VECTOR(logphi); // negative binomial dispersion parameters by species
  PARAMETER(logphibar);
  PARAMETER(loglogsigma_phi);
  
  Type jnll = 0.0;
  Type sigma_g = exp(logsigma_g);
  Type sigma_z= exp(logsigma_z);
  Type sigma_bo = exp(logsigma_bo);
  Type sigma_bog = exp(logsigma_bog);
  Type sigma_blat = exp(logsigma_blat);
  Type sigma_blength = exp(logsigma_blength);
  Type logsigma_phi= exp(loglogsigma_phi);
  
  
  
  ADREPORT(sigma_g);
  ADREPORT(sigma_z);
  ADREPORT(sigma_bo);
  ADREPORT(sigma_bog);
  ADREPORT(sigma_blat);
  ADREPORT(sigma_blength);
  ADREPORT(logsigma_phi);
  
  jnll -= sum(dnorm(z, 0, sigma_z, true));
  jnll -= sum(dnorm(bo, 0, sigma_bo, true)); 
  jnll -= sum(dnorm(bog,bobar, sigma_bog, true)); 
  jnll -= sum(dnorm(blat, blatbar, sigma_blat, true)); 
  jnll -= sum(dnorm(blength, blengthbar, sigma_blength, true)); 
  jnll -= sum(dnorm(logphi, logphibar, logsigma_phi, true)); 
  jnll -= sum(dnorm(graw, 0, 1, true)); 
  

  vector<Type> g(nspc);
  
  // get species group effects on slopes
  for(int i = 0; i < nspc; i++){ 
    g[i] = gbar[LH[i]] + graw[i] * sigma_g;
  }
  
  
  vector<Type> phi(nspc);
  for(int i = 0; i < nspc; i++){ 
    phi[i]= exp(logphi[i]);
  }
  ADREPORT(phi);
  
  vector<Type>  mu(n);
  mu = exp(Xint * bo + Xg * bog +  Xlat*blat + Xlength * blength + U * g +  Z * z) ;
  

  
  for(int i = 0; i < n; i++){ 
    Type var = mu[i] + pow(mu[i],2) * pow (phi[Xspc[i]], -1);
    jnll -=  	dnbinom2(y[i], mu[i],var, true);
  }
    
    // Now get predicted values
    vector<Type>  mu_predict1(npredict);
  vector<Type>  mu_predict2(npredict);
  vector<Type>  mu_predict3(npredict);
  
    for(int i = 0; i < npredict; i++){ 
      mu_predict1[i] =  exp(bog[0] + Upredict[i] * gbar[0]);
      mu_predict2[i] =  exp(bog[1] + Upredict[i] * gbar[1]);
      mu_predict3[i] =  exp(bog[2] + Upredict[i] * gbar[2]);
    }
    ADREPORT(mu_predict1);
    ADREPORT(mu_predict2);
    ADREPORT(mu_predict3);
    
    return jnll;
  }
  
