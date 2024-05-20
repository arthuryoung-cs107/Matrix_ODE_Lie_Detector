#include "matrix_Lie_detector.hh"
#include "LD_integrators.hh"

const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

const char eqn_name[] = "Duffing"; ode_solspc_meta meta0(2,1);

// const char bse_name[] = "Legendre";
const char bse_name[] = "Chebyshev1";
// const char bse_name[] = "Chebyshev2";

const int bor = 10;
// const int bor = 9;
// const int bor = 8;
// const int bor = 7;
// const int bor = 6;

const char exp_name[] = "DoP853_true_obs"; const char rec_name[] = "DoP853_true_rec"; const char ext_name[] = "DoP853_true_ext";
// const char exp_name[] = "DoPri5_true_obs"; const char rec_name[] = "DoPri5_true_rec"; const char ext_name[] = "DoPri5_true_ext";

const char Rmat_name[] = "Rmat";
const char Pmat_name[] = "Pmat";
const char Qmat_name[] = "Qmat";
const char Gmat_name[] = "Gmat";

orthopolynomial_space fspace0(meta0, bor);

const bool  extend_observations = true;

char  name_buffer[200],
      name_dnp1xu_buffer[200],
      name_JFs_buffer[200];

int main()
{
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s", dir_name,eqn_name,exp_name,dat_suff);
  LD_observations_set Sdat(meta0,ode_curve_observations(name_buffer));
  // LD_R_matrix Amat(fspace0,Sdat); const char * const mat_name = Rmat_name;
  // LD_P_matrix Amat(fspace0,Sdat); const char * const mat_name = Pmat_name;
  LD_Q_matrix Amat(fspace0,Sdat); const char * const mat_name = Qmat_name;
  // LD_G_matrix Amat(fspace0,Sdat); const char * const mat_name = Gmat_name;
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,mat_name,exp_name,dat_suff);
  Amat.read_matrix(name_buffer);

  LD_matrix_svd_result Amat_svd(Amat.ncrvs_tot,Amat.ndof_full);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,mat_name,exp_name,dat_suff);
  Amat_svd.read_svd_results(name_buffer);
  Amat_svd.print_details();

  rspace_infinitesimal_generator rinfgen0(fspace0,Amat_svd.kappa_def(),Amat_svd.VTtns);

  DoP853_settings integrator_settings; DoP853 infgen_integrator(rinfgen0,integrator_settings);
  // DoPri5_settings integrator_settings; DoPri5 infgen_integrator(rinfgen0,integrator_settings);

  generated_ode_observations inputs_recon(rinfgen0,Sdat.ncrvs_tot,Sdat.min_npts_curve());
  inputs_recon.set_solcurve_ICs(Sdat.curves);
  inputs_recon.generate_solution_curves(infgen_integrator,Sdat.get_default_IC_indep_range());
  sprintf(name_buffer,"%s/%s_%s.%s",dir_name,eqn_name,rec_name,dat_suff);
  inputs_recon.write_solution_curves(name_buffer);

  const int nbases0 = LD_threads::numthreads();
  orthopolynomial_basis ** bases0 = make_evaluation_bases<orthopolynomial_basis,orthopolynomial_space>(fspace0);
  bases0[LD_threads::thread_id()]->debugging_description();

  if (extend_observations)
  {
    ode_curve_observations inputs_extnd(meta0.eor,meta0.ndep,inputs_recon.nobs);
    rspace_infinitesimal_generator::init_extended_observations(inputs_extnd,inputs_recon);
    matrix_Lie_detector::extend_ode_observations<rspace_infinitesimal_generator,orthopolynomial_basis>(inputs_extnd,rinfgen0,bases0);
    sprintf(name_buffer,"%s/%s_%s.%s",dir_name,eqn_name,ext_name,dat_suff);
    inputs_extnd.write_observed_solutions(name_buffer);
  }

  free_evaluation_bases<orthopolynomial_basis>(bases0);

  return 0;
}
