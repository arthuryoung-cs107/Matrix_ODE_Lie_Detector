#include "LD_framework.hh"

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

  const int nbases0 = LD_threads::numthreads();
  orthopolynomial_basis ** bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0,nbases0);
  bases0[LD_threads::thread_id()]->debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s", dir_name,eqn_name,exp_name,dat_suff);
  sprintf(name_dnp1xu_buffer, "%s/%s_%s_%s.%s", dir_name,eqn_name,exp_name,"dnp1xu",dat_suff);
  sprintf(name_JFs_buffer, "%s/%s_%s_%s.%s", dir_name,eqn_name,exp_name,"JFs",dat_suff);

  LD_observations_set Sdat(meta0,input_ode_observations(name_buffer,name_dnp1xu_buffer,name_JFs_buffer));

  LD_R_matrix Rmat(fspace0,Sdat);
  Rmat.populate_R_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Rmat_name,exp_name,dat_suff);
  Rmat.write_matrix(name_buffer);

  LD_P_matrix Pmat(fspace0,Sdat);
  Pmat.populate_P_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Pmat_name,exp_name,dat_suff);
  Pmat.write_matrix(name_buffer);

  LD_G_matrix Gmat(fspace0,Sdat);
  Gmat.populate_G_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.%s_%s.%s", dir_name,eqn_name,bse_name,bor,Gmat_name,exp_name,dat_suff);
  Gmat.write_matrix(name_buffer);

  free_evaluation_bases<orthopolynomial_basis>(nbases0,bases0);

  return 0;
}
