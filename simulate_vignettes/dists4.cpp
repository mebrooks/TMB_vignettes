#include <TMB.hpp>                                
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(B);
  DATA_VECTOR(P);
  DATA_VECTOR(NB);
  DATA_VECTOR(G);
  
  PARAMETER(logit_prob); //for binomial
  Type prob=invlogit(logit_prob);
  ADREPORT(prob);
  PARAMETER(log_lambda); // for Poisson
  Type lambda=exp(log_lambda);
  ADREPORT(lambda);
  PARAMETER(log_mu); // for negative binomial
  Type mu=exp(log_mu);
  ADREPORT(mu);
  PARAMETER(log_var); // for negative binomial
  Type var=exp(log_var);
  ADREPORT(var);
  PARAMETER(log_shape); // for gamma
  Type shape=exp(log_shape);
  ADREPORT(shape);
  PARAMETER(log_scale); // for gamma
  Type scale=exp(log_scale);
  ADREPORT(scale);
  
  Type nll=0;
  nll -= sum(dbinom(B, Type(10), prob, true));
  nll -= sum(dpois(P, lambda, true));
  nll -= sum(dnbinom2(NB, mu, var, true));
  nll -= sum(dgamma(G, shape, scale, true));
  
  SIMULATE
  {
    for(int i=0; i<B.size(); i++)
    {
      B(i) = rbinom(Type(10), prob);
      P(i) = rpois(lambda);
      NB(i) = rnbinom2(mu, var);
      G(i) = rgamma(shape, scale);
    }
    REPORT(B);
    REPORT(P);
    REPORT(NB);
    REPORT(G);
  }
  return nll;
}
