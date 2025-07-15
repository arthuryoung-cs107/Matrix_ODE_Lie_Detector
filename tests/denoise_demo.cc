#include "LD_ode_system.hh"
#include "LD_integrators.hh"

#include <sys/stat.h>

// #include "LD_cokernal_policies.hh"
// #include "LD_cokernals.hh"
#include "LD_cokernals.hh"
// #include "LD_cokernal_vfields.hh"
#include "LD_encodings.hh"

#include "LD_noise_aux.hh"

// specify data directory for writing binary files
// const char dir_name[] = "./denoise_data_directory";
const char dir_name[] = "./denoise_data_directory/Gaussian_IC_perturbation";
// const char dir_name[] = "./dense_data_directory/Gaussian_IC_perturbation";
const char dat_suff[] = "lddat";
const char addtl_prefix[] = "";

const char * LD_noise_aux::dir_name = dir_name;

// specify subject ordinary differential equation for tests
Duffing_ode ode;
// VanDerPol_ode ode;
// Pendulum_ode ode(2.0,1.0,0.0,0.0,1.0,-9.81); // Pendulum_ode ode;
// Bessel_ode ode;
// Riccati_ode ode;
// Brusselator_ode ode;

ode_solspc_meta meta0(ode.eor,ode.ndep);

LD_noise_aux nse(meta0); // helper class for setting up and running noise experiments

const int nc = 50, // 50, number of curves
          np_min = 300, // 300, min number of points for extrapolation experiment
          np = np_min; // points per curve

// identifier for independent variable range of data
const int xrange = 0;

const bool  write_dnp1xu = true,
            write_JFs = true;

// level of noise applied to observational data. If <0, then unnoised
const int noise_level = -1;
// const int noise_level = 1;
// const int noise_level = 2;

LD_observations_set Sobs(meta0,nc,np,write_dnp1xu,write_JFs);

// specify runge kutta integrator for the generation of synthetic data
DoP853_settings integrator_settings; DoP853 ode_integrator(ode,integrator_settings);
// DoPri5_settings integrator_settings; DoPri5 ode_integrator(ode,integrator_settings);

// specify order of embedding function space
const int bor = 10;

orthopolynomial_space fspace0(meta0,bor); // class of embedding function space
// orthopolynomial_basis fbasis0(fspace0); // prolongation workspace based on function space
orthopolynomial_basis **fbases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);
orthopolynomial_basis &fbasis0 = *(fbases0[0]);

// orthogonal polynomial family for function space configuration
// const char fam_name[] = "Legendre";
const char fam_name[] = "Chebyshev1";
// const char fam_name[] = "Chebyshev2";

const bool  write_gen_obs_data = true,
            write_fspace_config = true,
            write_data_matrices = true;
const bool load_true_fspace_config = true; // override default noisy domain configuration file

char  eqn_name[50],
      intgen_name[50],
      noise_name[20],
      data_name[100],
      bse_name[75];
      // intrec_name[75],
      // recon_mat_name[50];

// character buffer for writing / reading data
char  name_buffer[350],
      name_dnp1xu_buffer[300],
      name_JFs_buffer[300];

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
  // inputs_gen.set_random_ICs(LD_rng(9365),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));
  inputs_gen.set_Gaussian_random_ICs(LD_rng(9365),ode.get_default_IC_range(),ode.get_default_IC_indep_range(xrange));
  inputs_gen.generate_solution_curves(ode_integrator);
  if (write_dnp1xu) inputs_gen.generate_dnp1xu();
  if (write_JFs) inputs_gen.generate_JFs();

  if (noise_level>=0)
  {
    sprintf(noise_name,"noise%d",noise_level);
    LD_noise_aux::perturb_pts(inputs_gen,nse.sigma_s,noise_level,1234,1.0,0.0);
    if (write_dnp1xu) LD_noise_aux::perturb_dnp1xu(inputs_gen,nse.sigma_dnp1xu,noise_level,12345,1);
    if (write_JFs) LD_noise_aux::perturb_JFs(inputs_gen,nse.sigma_JFs,noise_level,123456,1);
  }
  else sprintf(noise_name,"true");

  sprintf(data_name,"%s%s_%s_%sgen",addtl_prefix,eqn_name,noise_name,intgen_name);

  if (write_gen_obs_data)
  {
    sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
    inputs_gen.write_observed_solutions(name_buffer);
    if (write_dnp1xu)
    {
      sprintf(name_buffer,"%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_dnp1xu(name_buffer);
    }
    if (write_JFs)
    {
      sprintf(name_buffer,"%s/%s_JFs.%s",dir_name,data_name,dat_suff);
      inputs_gen.write_JFs(name_buffer);
    }
    if (noise_level>=0)
    {
      sprintf(name_buffer,"%s/%s_sigmas.%s",dir_name,data_name,dat_suff);
      nse.write_sigmas(name_buffer);
    }
  }
  return 0;
}

int configure_function_space(bool check_fspaces_=false)
{
  sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer),true);
  sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
  sprintf(name_dnp1xu_buffer, "%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_dnp1xu_buffer),false);
  sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
  sprintf(name_JFs_buffer, "%s/%s_JFs.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_JFs_buffer),false);

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
  if (load_true_fspace_config)
  {
    char true_data_name[100];
    sprintf(true_data_name,"%s%s_true_%sgen",addtl_prefix,eqn_name,intgen_name);
    sprintf(name_buffer, "%s/%s_%s.domain_config.%s",dir_name,true_data_name,bse_name,dat_suff);
  }
  else sprintf(name_buffer, "%s/%s_%s.domain_config.%s",dir_name,data_name,bse_name,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  return 0;
}

int compute_data_matrices()
{
  LD_L_encoder Lenc(meta0); printf("\n\nLmat\n");
  Jet_function_vector_space jfvs_L(Sobs,fspace0,Lenc,fbases0,false);
  jfvs_L.compute_svds(true);
  sprintf(name_buffer, "%s/%s_%s.L.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_L.Acode.write_LD_encoding_bundle(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.Lsvd.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_L.svd0.write_LD_svd_bundle(name_buffer);

  LD_G_encoder Genc(meta0); printf("\n\nGmat\n");
  Jet_function_vector_space jfvs_G(Sobs,fspace0,Genc,fbases0,false);
  jfvs_G.compute_svds(true);
  sprintf(name_buffer, "%s/%s_%s.G.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_G.Acode.write_LD_encoding_bundle(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.Gsvd.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_G.svd0.write_LD_svd_bundle(name_buffer);

  LD_R_encoder Rnenc(meta0,meta0.eor); printf("\n\nRnmat\n");
  Jet_function_vector_space jfvs_Rn(Sobs,fspace0,Rnenc,fbases0,false);
  jfvs_Rn.compute_svds(true);
  sprintf(name_buffer, "%s/%s_%s.Rn.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_Rn.Acode.write_LD_encoding_bundle(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.Rnsvd.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_Rn.svd0.write_LD_svd_bundle(name_buffer);

  LD_R_encoder Rnp1enc(meta0,meta0.eor+1); printf("\n\nRnp1mat\n");
  Jet_function_vector_space jfvs_Rnp1(Sobs,fspace0,Rnp1enc,fbases0,false);
  jfvs_Rnp1.compute_svds(true);
  sprintf(name_buffer, "%s/%s_%s.Rnp1.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_Rnp1.Acode.write_LD_encoding_bundle(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.Rnp1svd.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_Rnp1.svd0.write_LD_svd_bundle(name_buffer);

  LD_OG_encoder OGenc(meta0); printf("\n\nOGmat\n");
  Jet_function_vector_space jfvs_OG(Sobs,fspace0,OGenc,fbases0,false);
  jfvs_OG.compute_svds(true);
  sprintf(name_buffer, "%s/%s_%s.OG.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_OG.Acode.write_LD_encoding_bundle(name_buffer);
  sprintf(name_buffer, "%s/%s_%s.OGsvd.%s",dir_name,data_name,bse_name,dat_suff);
    jfvs_OG.svd0.write_LD_svd_bundle(name_buffer);

  return 0;
}

int compute_jets()
{
  LD_trivial_jet_experiment tjexp(Sobs);
  tjexp.determine_trivial_jets();

  tjexp.print_details();

  return 0;
}

int main()
{
  int gen_check = generate_trajectories(),
      cnf_check = configure_function_space(),
      mat_check = compute_data_matrices(),
      jet_check = compute_jets();
}
