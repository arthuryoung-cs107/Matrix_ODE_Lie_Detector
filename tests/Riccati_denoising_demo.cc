// #include <sys/stat.h>
// #include "named_odes.hh"
// #include "LD_trivial_flows.hh"

#include "LD_noise_aux.hh"
#include "Lie_detector.hh"

/*
  This file reads in noisy observed trajectories of the Riccati equation:

    dx u = (a*u)/x + b*(u^2)*(x^2),

    where x is the independent variable and u is the real valued dependent variable.
    Both a and b are real numbers, a*u denotes the scalar multiplication , and u^2 denotes the square of u.

  This is a quadratically nonlinear ordinary differential equation that is singular at x = 0.
  We observe trajectories and build a uniquely defined n=1 order model of the system.

  Assuming that the observations of the dependent variable, u, are corrupted by noise, we proceed to "denoise" the input data.
*/
// const char dir_name[] = "./denoise_data_directory/Gaussian_IC_perturbation"; // data directory
// const char obs_name[] = "Riccati_xrange0_noise0_DoP853gen"; // name of observations file
// const char dat_suff[] = "lddat"; // data file suffix (*.lddat)
// name of data file
const char data_name[] =
"./denoise_data_directory/Gaussian_IC_perturbation/"
"Riccati_xrange0_noise0_DoP853gen"
".lddat";

ode_curve_observations obs(data_name);
Lie_detector lie(obs); //

int main()
{
  obs.print_details();

  return 0;
}
