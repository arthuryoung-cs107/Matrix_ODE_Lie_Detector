#include "matrix_Lie_detector.hh"

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

// const char exp_name[] = "true_obs";
// const char exp_name[] = "DoPri5_true_obs";
const char exp_name[] = "DoP853_true_obs";

const char Rmat_name[] = "Rmat";
const char Pmat_name[] = "Pmat";
const char Qmat_name[] = "Qmat";
const char Gmat_name[] = "Gmat";

orthopolynomial_space fspace0(meta0, bor);

char  name_buffer[200],
      name_dnp1xu_buffer[200],
      name_JFs_buffer[200];

int main()
{
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s", dir_name,eqn_name,exp_name,dat_suff);
  sprintf(name_dnp1xu_buffer, "%s/%s_%s_%s.%s", dir_name,eqn_name,exp_name,"dnp1xu",dat_suff);
  sprintf(name_JFs_buffer, "%s/%s_%s_%s.%s", dir_name,eqn_name,exp_name,"JFs",dat_suff);

  LD_observations_set Sdat(meta0,ode_curve_observations(name_buffer,name_dnp1xu_buffer,name_JFs_buffer));

  LD_R_matrix Rmat(fspace0,Sdat);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Rmat_name,exp_name,dat_suff);
  Rmat.read_matrix(name_buffer);
  LD_matrix_svd_result Rmat_svd(Rmat.ncrvs_tot,Rmat.ndof_full);
  matrix_Lie_detector::compute_curve_svds(Rmat,Rmat_svd,Rmat.min_nrow_curve());
  Rmat_svd.print_details();
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,Rmat_name,exp_name,dat_suff);
  Rmat_svd.write_svd_results(name_buffer);

  LD_P_matrix Pmat(fspace0,Sdat);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Pmat_name,exp_name,dat_suff);
  Pmat.read_matrix(name_buffer);
  LD_matrix_svd_result Pmat_svd(Pmat.ncrvs_tot,Pmat.ndof_full);
  matrix_Lie_detector::compute_curve_svds(Pmat,Pmat_svd,Pmat.min_nrow_curve());
  Pmat_svd.print_details();
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,Pmat_name,exp_name,dat_suff);
  Pmat_svd.write_svd_results(name_buffer);

  LD_Q_matrix Qmat(fspace0,Sdat);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Qmat_name,exp_name,dat_suff);
  Qmat.read_matrix(name_buffer);
  LD_matrix_svd_result Qmat_svd(Qmat.ncrvs_tot,Qmat.ndof_full);
  matrix_Lie_detector::compute_curve_svds(Qmat,Qmat_svd,Qmat.min_nrow_curve());
  Qmat_svd.print_details();
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,Qmat_name,exp_name,dat_suff);
  Qmat_svd.write_svd_results(name_buffer);

  LD_G_matrix Gmat(fspace0,Sdat);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Gmat_name,exp_name,dat_suff);
  Gmat.read_matrix(name_buffer);
  LD_matrix_svd_result Gmat_svd(Gmat.ncrvs_tot,Gmat.ndof_full);
  matrix_Lie_detector::compute_curve_svds(Gmat,Gmat_svd,Gmat.min_nrow_curve());
  Gmat_svd.print_details();
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,Gmat_name,exp_name,dat_suff);
  Gmat_svd.write_svd_results(name_buffer);

  return 0;
}
