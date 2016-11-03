#include <TMB.hpp>                                
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER(mu);
  PARAMETER(log_resid_sd);
  Type resid_sd=exp(log_resid_sd);
  ADREPORT(resid_sd);
  
  //CALCULATE NEGATIVE LOG LIKELIHOOD
  Type nll;       
  for(int i=0; i<y.size(); i++)
  {
    nll -= dnorm(y(i), mu, resid_sd, true);
    SIMULATE {
      y(i) = rnorm(mu, resid_sd);
    }
  }
  
  SIMULATE {
    REPORT(y);
  }
  return nll;
}
