#include "LD_infinitesimal_generators.hh"
#include "LD_integrators.hh"
#include "LD_encodings.hh"
#include "LD_vspace_evaluators.hh"
// #include "solution_space_cokernals.hh"
// #include "LD_cokernals.hh"
#include "LD_cokernal_policies.hh"
#include "LD_cokernal_vfields.hh"

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

// const char data_name[] = "Duffing_xrange0_true_DoP853gen";

const char eqn_name[] = "Duffing"; ode_solspc_meta meta0(2,1);

LD_name_buffer  nbuf0(dir_name,dat_suff,100),
                nbuf1(dir_name,dat_suff,100),
                nbuf2(dir_name,dat_suff,100);

int main()
{

  // load observational data
  LD_name_buffer obs_name; obs_name.name_observational_data(eqn_name,0,-1,"DoP853");
  // LD_observations_set Sobs(meta0,ode_curve_observations(nbuf0.name_file(obs_name)));
  // LD_observations_set Sobs(meta0,ode_curve_observations(  nbuf0.name_file(obs_name),
                                                          // nbuf1.name_file(obs_name,"_dnp1xu")));
                                                          // nbuf2.name_file(obs_name,"_JFs")));
  LD_observations_set Sobs(meta0,ode_curve_observations(  nbuf0.name_file(obs_name),
                                                          nbuf1.name_file(obs_name,"_dnp1xu"),
                                                          nbuf2.name_file(obs_name,"_JFs")));

  const int bor = 10;
    LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",bor);

  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(nbuf0.name_domain_config_file(obs_name,fam_name)));
    fspace0.debugging_description();

  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);
  orthopolynomial_basis &basis0 = *(bases0[0]);

  const bool normalize_flag = true;
  // const char rec_suf[] = ".DoP853rnrec";

  const bool wdistance = true;
  // const bool wdistance = false;
  const bool normalize_cokernals = normalize_flag;
  LD_OG_encoder OGenc(meta0);

  // Frobenius_vspace_measure msr(fspace0.ndof_full,Sobs.ncrvs_tot);
  indepcomp_vspace_measure msr(fspace0.ndof_full,Sobs.ncrvs_tot);

  nullspace_clst_policy pol(msr,wdistance);
  // nullspace_near_policy pol(msr,wdistance);

  Jet_function_vector_space jfvs(Sobs,fspace0,OGenc,bases0,normalize_flag);
  cokernal_refinement rfne(jfvs,pol,true);


  // LD_integration_package int_pckg(Sobs);

  printf("ending main\n"); free_evaluation_bases<orthopolynomial_basis>(bases0); return 0;
}
