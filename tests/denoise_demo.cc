#include "LD_ode_system.hh"
#include "LD_integrators.hh"

#include <sys/stat.h>

#include "LD_cokernal_policies.hh"
// #include "LD_cokernal_vfields.hh"

#include "LD_noise_aux.hh"

// specify data directory for writing binary files
const char dir_name[] = "./denoise_data_directory";
const char dat_suff[] = "lddat";
const char addtl_prefix[] = "";

// specify subject ordinary differential equation for tests
Duffing_ode ode;
// VanDerPol_ode ode;
// Pendulum_ode ode;
// Bessel_ode ode;
// Riccati_ode ode;
// Brusselator_ode ode;

ode_solspc_meta meta0(ode.eor,ode.ndep);

const int nc = 50, // number of curves
          np_min = 300, // min number of points for extrapolation experiment
          np = np_min; // points per curve

// identifier for independent variable range of data
const int xrange = 0;

const bool  write_dnp1xu = true,
            write_JFs = true;

// level of noise applied to observational data. If <0, then unnoised
// const int noise_level = -1;
const int noise_level = 2;

LD_observations_set Sobs(meta0,nc,np,write_dnp1xu,write_JFs);

// specify runge kutta integrator for the generation of synthetic data
DoP853_settings integrator_settings; DoP853 ode_integrator(ode,integrator_settings);
// DoPri5_settings integrator_settings; DoPri5 ode_integrator(ode,integrator_settings);

// specify order of embedding function space
const int bor = 10;

orthopolynomial_space fspace0(meta0,bor); // class of embedding function space
orthopolynomial_basis fbasis0(fspace0); // prolongation workspace based on function space

// orthogonal polynomial family for function space configuration
// const char fam_name[] = "Legendre";
const char fam_name[] = "Chebyshev1";
// const char fam_name[] = "Chebyshev2";

const bool  write_gen_obs_data = true,
            write_fspace_config = true;

char  eqn_name[50],
      intgen_name[50],
      noise_name[20],
      data_name[100],
      bse_name[75];
      // intrec_name[75],
      // recon_mat_name[50];

// character buffer for writing / reading data
char  name_buffer[350];


/*

  task functions invoked by main

*/


int generate_trajectories()
{
  mkdir(dir_name, S_IRWXU);

  if (np>np_min) sprintf(eqn_name,"%s_xrange%d_extrap",ode.name,xrange);
  else sprintf(eqn_name,"%s_xrange%d",ode.name,xrange);

  strcpy(intgen_name,ode_integrator.name);

  generated_ode_observations inputs_gen(ode,nc,np);
  inputs_gen.set_random_ICs(LD_rng(9365),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));
  inputs_gen.generate_solution_curves(ode_integrator);

  if (noise_level>=0)
  {
    // inputs_gen.apply_noise_pts();
    LD_noise_aux::perturb_pts(inputs_gen,noise_level,1234,0.1);
    sprintf(noise_name,"noise%d",noise_level);
  }
  else sprintf(noise_name,"true");

  sprintf(data_name,"%s%s_%s_%sgen",addtl_prefix,eqn_name,noise_name,intgen_name);

  if (write_gen_obs_data)
  {
    sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
    inputs_gen.write_observed_solutions(name_buffer);
    if (write_dnp1xu)
    {
      inputs_gen.generate_dnp1xu();
      sprintf(name_buffer,"%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_dnp1xu(name_buffer);
    }
    if (write_JFs)
    {
      inputs_gen.generate_JFs();
      sprintf(name_buffer,"%s/%s_JFs.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_JFs(name_buffer);
    }
  }
  return 0;
}

int configure_function_space(bool check_fspaces_=false)
{
  sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer),true);

  if (write_fspace_config)
  {
    fspace0.set_Legendre_coeffs(); Sobs.configure_centered_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Legendre",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    fspace0.set_Chebyshev1_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Chebyshev1",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    // fspace0.set_Chebyshev2_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
    // fspace0.set_Chebyshev2_coeffs(); Sobs.configure_centered_domain(fspace0);
    fspace0.set_Chebyshev2_coeffs(); Sobs.configure_center_mass_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Chebyshev2",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);
  }
  sprintf(bse_name, "%s.%d",fam_name,bor);
  sprintf(name_buffer, "%s/%s_%s.domain_config.%s",dir_name,data_name,bse_name,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  return 0;
}

int main()
{

  generate_trajectories();

  int gen_check = generate_trajectories(),
      cnf_check = configure_function_space();
}
