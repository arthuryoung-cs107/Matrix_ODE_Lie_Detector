#include "matrix_Lie_detector.hh"

const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

const char eqn_name[] = "Duffing";
ode_solspc_meta meta0(2,1);

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

const char mat_name[] = "Rmat";

orthopolynomial_space fspace0(meta0, bor);

char name_buffer[200];

int main()
{
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s", dir_name,eqn_name,exp_name,dat_suff);
  LD_observations_set Sdat(meta0,input_ode_observations(name_buffer));
  LD_R_matrix Rmat(fspace0,Sdat);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,mat_name,exp_name,dat_suff);
  Rmat.read_matrix(name_buffer);

  LD_matrix_svd_result Rmat_svd(Rmat.ncrvs_tot,Rmat.ndof_full);
  matrix_Lie_detector::compute_curve_svds(Rmat,Rmat_svd,Rmat.min_nrow_curve());
  Rmat_svd.print_details();
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s_svd.%s", dir_name,eqn_name,bse_name,bor,mat_name,exp_name,dat_suff);
  Rmat_svd.write_svd_results(name_buffer);
}
