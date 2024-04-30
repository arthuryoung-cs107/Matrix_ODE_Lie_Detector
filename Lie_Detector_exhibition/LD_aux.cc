#include "LD_aux.hh"
#include <gsl/gsl_randist.h>

LD_linalg::LD_linalg() {}
LD_linalg::~LD_linalg() {}

double LD_rng::rand_uni(double low_,double high_)
  {return (gsl_rng_uniform(gsl_gen)-0.5)*(high_-low_) + 0.5*(high_+low_);}
double LD_rng::rand_gau(double mu_, double sigma_)
  {return mu_+gsl_ran_gaussian(gsl_gen,sigma_);}
