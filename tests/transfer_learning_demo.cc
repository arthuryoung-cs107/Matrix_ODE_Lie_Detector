#include "matrix_Lie_detector.hh"
#include "LD_ode_system.hh"
#include "LD_integrators.hh"
#include <sys/stat.h>

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

const char eqn_name[] = "Duffing"; ode_solspc_meta meta0(2,1);

const int noise_level = -1;

// number of curves and uniform number of points per curve for dataset
const int nc = 50, // number of curves
          np_min = 300; // MINIMUM points per curve for extrapolation test

// specify order of embedding function space
const int bor = 10;
// const int bor = 9;
// const int bor = 8;

// class of embedding function space
orthopolynomial_space fspace0(meta0,bor);

// orthogonal polynomial family for function space configuration
// const char fam_name[] = "Legendre";
const char fam_name[] = "Chebyshev1";
// const char fam_name[] = "Chebyshev2";

int main()
{
  mkdir(dir_name, S_IRWXU); strcpy(eqn_name,ode.name); strcpy(intgen_name,ode_integrator.name);



  return 0;
}
