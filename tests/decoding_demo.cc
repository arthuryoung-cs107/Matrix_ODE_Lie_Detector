#include "LD_framework.hh"

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

const char eqn_name[] = "Duffing";
ode_solspc_meta meta0(2,1);

const char bse_name[] = "Legendre";
// const char bse_name[] = "Chebyshev1";
// const char bse_name[] = "Chebyshev2";
const int bor = 10;

const char exp_name[] = "true_obs";

orthopolynomial_space fspace0(meta0, bor);

char name_buffer[200];

int main()
{
  sprintf(name_buffer, "%s/%s_%s.%d.domain_config.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  fspace0.configure_self(name_buffer);
  fspace0.debugging_description();

  const int nbases0 = numthreads();
  orthopolynomial_basis ** bases0 = make_evaluation_bases<orthopolynomial_basis,orthopolynomial_space>(fspace0,nbases0);
  bases0[thread_id()]->debugging_description();

  sprintf(name_buffer, "%s/%s_%s.%s", dir_name,eqn_name,exp_name,dat_suff);
  LD_observations_set Sdat(meta0,input_ode_observations(name_buffer));
  LD_R_matrix Rmat(fspace0,Sdat);
  // Rmat.populate_R_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.Rmat_%s.%s", dir_name,eqn_name,bse_name,bor,exp_name,dat_suff);
  // Rmat.write_matrix(name_buffer);
  Rmat.read_matrix(name_buffer);

  free_evaluation_bases<orthopolynomial_basis>(nbases0,bases0);
}
