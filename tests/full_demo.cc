#include "matrix_Lie_detector.hh"
#include "LD_ode_system.hh"
#include "LD_integrators.hh"
#include <sys/stat.h>

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

// specify subject ordinary differential equation for tests
Duffing_ode ode(-1.0,1.0,0.5,0.3,1.2); // chaos
ode_solspc_meta meta0(ode.eor,ode.ndep);

// number of curves and uniform number of points per curve for dataset
const int nc = 50, // number of curves
          np = 300; // points per curve

const bool  write_dnp1xu = true,
            write_JFs = true;

// level of noise applied to observational data. If <0, then unnoised
const int noise_level = -1;
// const int noise_level = 0;

LD_observations_set Sobs(meta0,nc,np,write_dnp1xu,write_JFs);

// specify runge kutta integrator for the generation of synthetic data
// DoP853_settings integrator_settings; DoP853 ode_integrator(ode,integrator_settings);
DoPri5_settings integrator_settings; DoPri5 ode_integrator(ode,integrator_settings);

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

// names of encoded data matrices
const char Rmat_name[] = "Rmat";
const char Pmat_name[] = "Pmat";
const char Qmat_name[] = "Qmat";
const char Gmat_name[] = "Gmat";

LD_R_matrix Amat(fspace0,Sobs); const char * const Amat_name = Rmat_name;
// LD_P_matrix Amat(fspace0,Sobs); const char * const Amat_name = Pmat_name;
// LD_Q_matrix Amat(fspace0,Sobs); const char * const Amat_name = Qmat_name;
// LD_G_matrix Amat(fspace0,Sobs); const char * const Amat_name = Gmat_name;

LD_matrix_svd_result Amat_svd(Amat.ncrvs_tot,Amat.ndof_full);

const bool  write_gen_obs_data = true,
            write_fspace_config = true,
            write_encoded_mats = true,
            write_decoded_mats = true;

// buffers for writing file names via meta strings
char  eqn_name[50],
      intgen_name[50],
      noise_name[20],
      data_name[100],
      bse_name[75],
      mat_name[20],
      intrec_name[50];

// character buffer for writing / reading data
char  name_buffer[200],
      name_dnp1xu_buffer[200],
      name_JFs_buffer[200];

int generate_observational_data(bool write_data_)
{
  mkdir(dir_name, S_IRWXU); strcpy(eqn_name,ode.name); strcpy(intgen_name,ode_integrator.name);

  generated_ode_observations inputs_gen(ode,nc,np);
  inputs_gen.set_random_ICs(LD_rng(9365),ode.get_default_IC_range());  // shay's number
  inputs_gen.generate_solution_curves(ode_integrator,ode.get_default_IC_indep_range());

  if (noise_level>=0)
  {
    // inputs_gen.apply_noise_pts();
    sprintf(noise_name,"noise%d",noise_level);
  }
  else sprintf(noise_name,"true");

  sprintf(data_name,"%s_%s_%sgen",eqn_name,noise_name,intgen_name);

  if (write_data_)
  {
    sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
    inputs_gen.write_solution_curves(name_buffer);
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

int configure_function_space(bool write_data_,bool check_fspaces_=false)
{
  sprintf(name_buffer,"%s/%s.%s",dir_name,data_name,dat_suff);
  Sobs.load_additional_inputs(ode_curve_observations(name_buffer),true);

  if (write_data_)
  {
    fspace0.set_Legendre_coeffs(); Sobs.configure_centered_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Legendre",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    fspace0.set_Chebyshev1_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
    if (check_fspaces_) fspace0.debugging_description();
    sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s",dir_name,data_name,"Chebyshev1",bor,dat_suff);
    fspace0.write_configuration_file(name_buffer);

    fspace0.set_Chebyshev2_coeffs(); Sobs.configure_0maxmag_0pi05_domain(fspace0);
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

template<class BSIS> int encode_data_matrices(BSIS **bases_, bool write_data_)
{
  bases_[LD_threads::thread_id()]->debugging_description();

  if (write_data_)
  {
    LD_R_matrix Rmat(fspace0,Sobs); Rmat.populate_R_matrix<orthopolynomial_basis>(bases_);
    sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Rmat_name,dat_suff);
    Rmat.write_matrix(name_buffer);

    if (write_dnp1xu)
    {
      sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
      sprintf(name_dnp1xu_buffer, "%s/%s_dnp1xu.%s",dir_name,data_name,dat_suff);
      Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_dnp1xu_buffer),false);

      LD_P_matrix Pmat(fspace0,Sobs); Pmat.populate_P_matrix<orthopolynomial_basis>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Pmat_name,dat_suff);
      Pmat.write_matrix(name_buffer);

      LD_Q_matrix Qmat(fspace0,Sobs); Qmat.populate_Q_matrix<orthopolynomial_basis>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Qmat_name,dat_suff);
      Qmat.write_matrix(name_buffer);
    }
    if (write_JFs)
    {
      sprintf(name_buffer, "%s/%s.%s",dir_name,data_name,dat_suff);
      sprintf(name_JFs_buffer, "%s/%s_JFs.%s",dir_name,data_name,dat_suff);
      Sobs.load_additional_inputs(ode_curve_observations(name_buffer,name_JFs_buffer),false);

      LD_G_matrix Gmat(fspace0,Sobs); Gmat.populate_G_matrix<orthopolynomial_basis>(bases_);
      sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Gmat_name,dat_suff);
      Gmat.write_matrix(name_buffer);
    }
  }
  strcpy(mat_name,Amat_name);
  sprintf(name_buffer, "%s/%s_%s.%s.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
  Amat.read_matrix(name_buffer);

  return 0;
}

int decode_data_matrices(bool write_data_)
{
  sprintf(name_buffer, "%s/%s_%s.%s_svd.%s", dir_name,data_name,bse_name,Amat_name,dat_suff);
  if (write_data_)
  {
    matrix_Lie_detector::compute_curve_svds(Amat,Amat_svd,Amat.min_nrow_curve());
    Amat_svd.write_svd_results(name_buffer);
  }
  else Amat_svd.read_svd_results(name_buffer);

  Amat_svd.print_details();

  return 0;
}

int main()
{
  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);

  int gen_check = generate_observational_data(write_gen_obs_data),
      cnf_check = configure_function_space(write_fspace_config),
      enc_check = encode_data_matrices<orthopolynomial_basis>(bases0,write_encoded_mats),
      dec_check = decode_data_matrices(write_decoded_mats);

  rspace_infinitesimal_generator rinfgen0(fspace0,Amat_svd.kappa_def(),Amat_svd.VTtns);

  DoP853_settings intgr_rec_settings; DoP853 intgr_rec(rinfgen0,intgr_rec_settings);
  // DoPri5_settings intgr_rec_settings; DoPri5 intgr_rec(rinfgen0,intgr_rec_settings);

  strcpy(intrec_name,intgr_rec.name);

  generated_ode_observations inputs_recon(rinfgen0,Sobs.ncrvs_tot,Sobs.min_npts_curve());
  inputs_recon.set_solcurve_ICs(Sobs.curves);

  inputs_recon.parallel_generate_solution_curves<rspace_infinitesimal_generator,DoP853>(rinfgen0,intgr_rec,Sobs.get_default_IC_indep_range());
  // inputs_recon.parallel_generate_solution_curves<rspace_infinitesimal_generator,DoPri5>(rinfgen0,intgr_rec,Sobs.get_default_IC_indep_range());

  sprintf(name_buffer, "%s/%s_%s.%s.%srec.%s", dir_name,data_name,bse_name,Amat_name,intrec_name,dat_suff);
  inputs_recon.write_solution_curves(name_buffer);

  ode_curve_observations inputs_extnd(meta0.eor,meta0.ndep,inputs_recon.nobs);
  rspace_infinitesimal_generator::init_extended_observations(inputs_extnd,inputs_recon);
  matrix_Lie_detector::extend_ode_observations<rspace_infinitesimal_generator,orthopolynomial_basis>(inputs_extnd,rinfgen0,bases0);
  sprintf(name_buffer, "%s/%s_%s.%s.%sext.%s", dir_name,data_name,bse_name,Amat_name,intrec_name,dat_suff);
  inputs_extnd.write_observed_solutions(name_buffer);

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;
}
