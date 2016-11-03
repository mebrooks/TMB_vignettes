#include <TMB.hpp>                                
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  DATA_MATRIX(X);                             
  PARAMETER(log_resid_sd);
  Type resid_sd=exp(log_resid_sd);
  ADREPORT(resid_sd);
  PARAMETER_VECTOR(beta);
  
  //CALCULATE NEGATIVE LOG LIKELIHOOD
  Type nll;       
  vector<Type> Xbeta= X*beta;
  for(int i=0; i<y.size(); i++)
  {
    nll -= dnorm(y(i), Xbeta(i), resid_sd, true);
    SIMULATE {
      y(i) = rnorm(Xbeta(i), resid_sd);
    }
  }
  
  SIMULATE {
    REPORT(y);
  }
  return nll;
}
