#include "matrix_Lie_detector.hh"

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

// const char data_name[] = "Duffing_xrange0_true_DoP853gen";

const char eqn_name[] = "Duffing"; ode_solspc_meta meta0(2,1);

LD_name_buffer name(dir_name,dat_suff,100);
int main()
{
  // load observational data
  LD_name_buffer obs_name; obs_name.name_observational_data(eqn_name,0,-1,"DoP853");
  LD_observations_set Sobs(meta0,ode_curve_observations(name.name_file(obs_name)));

  LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",8);
  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(name.name_domain_config_file(obs_name,fam_name)));
  fspace0.debugging_description();

  LD_matrix Lmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"L")));
  LD_matrix_svd_result Lmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L")));
  Lmat_svd.print_details("Lmat_svd");

  LD_matrix Rmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"R")));
  LD_matrix_svd_result Rmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R")));
  Rmat_svd.print_details("Rmat_svd");

  LD_matrix Qmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"Q")));
  LD_matrix_svd_result Qmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"Q")));
  Qmat_svd.print_details("Qmat_svd");

  LD_matrix Omat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"O")));
  LD_matrix_svd_result Omat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"O")));
  Omat_svd.print_details("Omat_svd");

  LD_matrix Gmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"G")));
  LD_matrix_svd_result Gmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"G")));
  Gmat_svd.print_details("Gmat_svd");


  return 0;
}
