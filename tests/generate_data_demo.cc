#include "LD_framework.hh"
#include "LD_ode_system.hh"
#include "LD_integrators.hh"

#ifdef _OPENMP
  #include "omp.h"
  int thread_id() {return omp_get_thread_num();}
  int numthreads() {return omp_get_max_threads();}
#else
  int thread_id() {return 0;}
  int numthreads() {return 1;}
#endif

const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

Duffing_ode ode(-1.0,1.0,0.5,0.3,1.2); // chaos
ode_solspc_meta meta0(ode.eor,ode.ndep);

dop853_settings ode_integrator_settings;
dop853_integrator ode_integrator(ode,ode_integrator_settings);

const int nc = 50, // number of curves
          np = 300; // points per curve

const int bor = 10;
orthopolynomial_space fspace0(meta0,bor);

const bool  write_JFs = true,
            write_dnp1xu = true;

char name_buffer[200];

int main()
{
  generated_ode_observations inputs_gen(ode,nc,np);
  inputs_gen.set_random_ICs(LD_rng(9365),ode.get_default_IC_range());  // shay's number
  inputs_gen.generate_solution_curves(ode_integrator, ode.get_default_IC_indep_range());

  sprintf(name_buffer,"%s/%s_%s.%s",dir_name,ode.name,"true_obs",dat_suff);
  inputs_gen.write_solution_curves(name_buffer);

  if (write_JFs)
  {

  }
  if (write_dnp1xu)
  {

  }

  sprintf(name_buffer,"%s/%s_%s.%s",dir_name,ode.name,"true_obs",dat_suff);
  LD_observations_set Strue(meta0,input_ode_observations(name_buffer));

  fspace0.set_Legendre_coeffs(); Strue.configure_centered_domain(fspace0);
  fspace0.debugging_description();
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,ode.name,"Legendre",bor,dat_suff);
  fspace0.write_configuration_file(name_buffer);

  fspace0.set_Chebyshev1_coeffs(); Strue.configure_0maxmag_0pi05_domain(fspace0);
  fspace0.debugging_description();
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,ode.name,"Chebyshev1",bor,dat_suff);
  fspace0.write_configuration_file(name_buffer);

  fspace0.set_Chebyshev2_coeffs(); Strue.configure_0maxmag_0pi05_domain(fspace0);
  fspace0.debugging_description();
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,ode.name,"Chebyshev2",bor,dat_suff);
  fspace0.write_configuration_file(name_buffer);
}
