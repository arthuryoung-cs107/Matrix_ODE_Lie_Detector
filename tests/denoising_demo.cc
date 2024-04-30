#include "LD_framework.hh"
// #include "LD_function_space.hh"
// #include "LD_aux.hh"
// #include <cstdio>

#ifdef _OPENMP
  #include "omp.h"
  int thread_id() {return omp_get_thread_num();}
  int numthreads() {return omp_get_max_threads();}
#else
  int thread_id() {return 0;}
  int numthreads() {return 1;}
#endif

const char dir_name[] = "./writing_dat_dir";
const char dat_suff[] = "lddat";

const char eqn_name[] = "Duffing";
ode_solspc_meta meta0(2,1);

const char bse_name[] = "Chebyshev";
// const char bse_name[] = "Legendre";
const int bor = 8;
// const int bor = 10;

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

  sprintf(name_buffer, "%s/%s_true_obs.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_true(name_buffer);
  // inputs_true.print_details();
  // LD_observations_set Strue(meta0,inputs_true);
  LD_observations_set Strue(meta0,input_ode_observations(name_buffer));
  LD_R_matrix Rmat_true(fspace0,Strue);
  Rmat_true.populate_R_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.Rmat_true.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  Rmat_true.write_matrix(name_buffer);

  sprintf(name_buffer, "%s/%s_nois_obs.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_nois(name_buffer);
  // inputs_nois.print_details();
  // LD_observations_set Snois(meta0,inputs_nois);
  LD_observations_set Snois(meta0,input_ode_observations(name_buffer));
  LD_R_matrix Rmat_nois(fspace0,Snois);
  Rmat_nois.populate_R_matrix<orthopolynomial_basis>(bases0);
  sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  Rmat_nois.write_matrix(name_buffer);

  // sprintf(name_buffer, "%s/%s_nois_obs_sm1_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm1_rr(name_buffer);
  // inputs_sm1_rr.print_details();
  // LD_observations_set Ssm1_rr(meta0,inputs_sm1_rr);
  // LD_R_matrix Rmat_sm1_rr(fspace0,Ssm1_rr);
  // Rmat_sm1_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm1_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm1_rr.write_matrix(name_buffer);
  //
  // sprintf(name_buffer, "%s/%s_nois_obs_sm2_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm2_rr(name_buffer);
  // inputs_sm2_rr.print_details();
  // LD_observations_set Ssm2_rr(meta0,inputs_sm2_rr);
  // LD_R_matrix Rmat_sm2_rr(fspace0,Ssm2_rr);
  // Rmat_sm2_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm2_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm2_rr.write_matrix(name_buffer);
  //
  // sprintf(name_buffer, "%s/%s_nois_obs_sm3_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm3_rr(name_buffer);
  // inputs_sm3_rr.print_details();
  // LD_observations_set Ssm3_rr(meta0,inputs_sm3_rr);
  // LD_R_matrix Rmat_sm3_rr(fspace0,Ssm3_rr);
  // Rmat_sm3_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm3_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm3_rr.write_matrix(name_buffer);



  // sprintf(name_buffer, "%s/%s_nois_obs_sm1_kfull_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm1_kfull_rr(name_buffer);
  // inputs_sm1_kfull_rr.print_details();
  // LD_observations_set Ssm1_kfull_rr(meta0,inputs_sm1_kfull_rr);
  // LD_R_matrix Rmat_sm1_kfull_rr(fspace0,Ssm1_kfull_rr);
  // Rmat_sm1_kfull_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm1_kfull_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm1_kfull_rr.write_matrix(name_buffer);

  // sprintf(name_buffer, "%s/%s_nois_obs_sm1_k1_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm1_k1_rr(name_buffer);
  // inputs_sm1_k1_rr.print_details();
  // LD_observations_set Ssm1_k1_rr(meta0,inputs_sm1_k1_rr);
  // LD_R_matrix Rmat_sm1_k1_rr(fspace0,Ssm1_k1_rr);
  // Rmat_sm1_k1_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm1_k1_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm1_k1_rr.write_matrix(name_buffer);
  // sprintf(name_buffer, "%s/%s_nois_obs_sm1_k2_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm1_k2_rr(name_buffer);
  // inputs_sm1_k2_rr.print_details();
  // LD_observations_set Ssm1_k2_rr(meta0,inputs_sm1_k2_rr);
  // LD_R_matrix Rmat_sm1_k2_rr(fspace0,Ssm1_k2_rr);
  // Rmat_sm1_k2_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm1_k2_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm1_k2_rr.write_matrix(name_buffer);

  // sprintf(name_buffer, "%s/%s_nois_obs_sm2_k1_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm2_k1_rr(name_buffer);
  // inputs_sm2_k1_rr.print_details();
  // LD_observations_set Ssm2_k1_rr(meta0,inputs_sm2_k1_rr);
  // LD_R_matrix Rmat_sm2_k1_rr(fspace0,Ssm2_k1_rr);
  // Rmat_sm2_k1_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm2_k1_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm2_k1_rr.write_matrix(name_buffer);
  // sprintf(name_buffer, "%s/%s_nois_obs_sm2_k2_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm2_k2_rr(name_buffer);
  // inputs_sm2_k2_rr.print_details();
  // LD_observations_set Ssm2_k2_rr(meta0,inputs_sm2_k2_rr);
  // LD_R_matrix Rmat_sm2_k2_rr(fspace0,Ssm2_k2_rr);
  // Rmat_sm2_k2_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm2_k2_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm2_k2_rr.write_matrix(name_buffer);
  //
  // sprintf(name_buffer, "%s/%s_nois_obs_sm3_k1_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm3_k1_rr(name_buffer);
  // inputs_sm3_k1_rr.print_details();
  // LD_observations_set Ssm3_k1_rr(meta0,inputs_sm3_k1_rr);
  // LD_R_matrix Rmat_sm3_k1_rr(fspace0,Ssm3_k1_rr);
  // Rmat_sm3_k1_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm3_k1_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm3_k1_rr.write_matrix(name_buffer);
  // sprintf(name_buffer, "%s/%s_nois_obs_sm3_k2_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm3_k2_rr(name_buffer);
  // inputs_sm3_k2_rr.print_details();
  // LD_observations_set Ssm3_k2_rr(meta0,inputs_sm3_k2_rr);
  // LD_R_matrix Rmat_sm3_k2_rr(fspace0,Ssm3_k2_rr);
  // Rmat_sm3_k2_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm3_k2_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm3_k2_rr.write_matrix(name_buffer);
  //
  // sprintf(name_buffer, "%s/%s_nois_obs_sm4_k1_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm4_k1_rr(name_buffer);
  // inputs_sm4_k1_rr.print_details();
  // LD_observations_set Ssm4_k1_rr(meta0,inputs_sm4_k1_rr);
  // LD_R_matrix Rmat_sm4_k1_rr(fspace0,Ssm4_k1_rr);
  // Rmat_sm4_k1_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm4_k1_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm4_k1_rr.write_matrix(name_buffer);
  // sprintf(name_buffer, "%s/%s_nois_obs_sm4_k2_raw_reg.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm4_k2_rr(name_buffer);
  // inputs_sm4_k2_rr.print_details();
  // LD_observations_set Ssm4_k2_rr(meta0,inputs_sm4_k2_rr);
  // LD_R_matrix Rmat_sm4_k2_rr(fspace0,Ssm4_k2_rr);
  // Rmat_sm4_k2_rr.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_sm4_k2_raw_reg.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm4_k2_rr.write_matrix(name_buffer);





  // input_ode_observations inputs_smooth1("./writing_dat_dir/Duffing_nois_obs_smooth1.lddat");
  // inputs_smooth1.print_details();
  // LD_observations_set Ssmooth1(meta0,inputs_smooth1);
  // LD_R_matrix Rmat_smooth1(fspace0,Ssmooth1);
  // Rmat_smooth1.populate_R_matrix<orthopolynomial_basis>(bases0);
  // Rmat_smooth1.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_smooth1.lddat");
  // Rmat_smooth1.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_smooth1.lddat");


  // sprintf(name_buffer, "%s/%s_nois_obs_smooth1_kfull.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm1_kfull(name_buffer);
  // inputs_sm1_kfull.print_details();
  // LD_observations_set Ssm1_kfull(meta0,inputs_sm1_kfull);
  // LD_R_matrix Rmat_sm1_kfull(fspace0,Ssm1_kfull);
  // Rmat_sm1_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_smooth1_kfull.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm1_kfull.write_matrix(name_buffer);


  // input_ode_observations inputs_br1_kfull("./writing_dat_dir/Duffing_nois_obs_brgman1_kfull.lddat");
  // inputs_br1_kfull.print_details();
  // LD_observations_set Sbr1_kfull(meta0,inputs_br1_kfull);
  // LD_R_matrix Rmat_br1_kfull(fspace0,Sbr1_kfull);
  // Rmat_br1_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_br1_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_brgman1_kfull.lddat");
  // Rmat_br1_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_brgman1_kfull.lddat");

  // input_ode_observations inputs_prjsm1_kfull("./writing_dat_dir/Duffing_nois_obs_prjsmooth1_kfull.lddat");
  // inputs_prjsm1_kfull.print_details();
  // LD_observations_set Sprjsm1_kfull(meta0,inputs_prjsm1_kfull);
  // LD_R_matrix Rmat_prjsm1_kfull(fspace0,Sprjsm1_kfull);
  // Rmat_prjsm1_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_prjsm1_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_prjsmooth1_kfull.lddat");
  // Rmat_prjsm1_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_prjsmooth1_kfull.lddat");


  // input_ode_observations inputs_smooth2("./writing_dat_dir/Duffing_nois_obs_smooth2.lddat");
  // inputs_smooth2.print_details();
  // LD_observations_set Ssmooth2(meta0,inputs_smooth2);
  // LD_R_matrix Rmat_smooth2(fspace0,Ssmooth2);
  // Rmat_smooth2.populate_R_matrix<orthopolynomial_basis>(bases0);
  // Rmat_smooth2.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_smooth2.lddat");
  // Rmat_smooth2.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_smooth2.lddat");


  // sprintf(name_buffer, "%s/%s_nois_obs_smooth2_kfull.%s", dir_name,eqn_name,dat_suff);
  // input_ode_observations inputs_sm2_kfull(name_buffer);
  // inputs_sm2_kfull.print_details();
  // LD_observations_set Ssm2_kfull(meta0,inputs_sm2_kfull);
  // LD_R_matrix Rmat_sm2_kfull(fspace0,Ssm2_kfull);
  // Rmat_sm2_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // sprintf(name_buffer, "%s/%s_%s.%d.Rmat_nois_smooth2_kfull.%s", dir_name,eqn_name,bse_name,bor,dat_suff);
  // Rmat_sm2_kfull.write_matrix(name_buffer);


  // input_ode_observations inputs_br2_kfull("./writing_dat_dir/Duffing_nois_obs_brgman2_kfull.lddat");
  // inputs_br2_kfull.print_details();
  // LD_observations_set Sbr2_kfull(meta0,inputs_br2_kfull);
  // LD_R_matrix Rmat_br2_kfull(fspace0,Sbr2_kfull);
  // Rmat_br2_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_br2_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_brgman2_kfull.lddat");
  // Rmat_br2_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_brgman2_kfull.lddat");

  // input_ode_observations inputs_prjsm2_kfull("./writing_dat_dir/Duffing_nois_obs_prjsmooth2_kfull.lddat");
  // inputs_prjsm2_kfull.print_details();
  // LD_observations_set Sprjsm2_kfull(meta0,inputs_prjsm2_kfull);
  // LD_R_matrix Rmat_prjsm2_kfull(fspace0,Sprjsm2_kfull);
  // Rmat_prjsm2_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_prjsm2_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_prjsmooth2_kfull.lddat");
  // Rmat_prjsm2_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_prjsmooth2_kfull.lddat");


  // input_ode_observations inputs_sm3_kfull("./writing_dat_dir/Duffing_nois_obs_smooth3_kfull.lddat");
  // inputs_sm3_kfull.print_details();
  // LD_observations_set Ssm3_kfull(meta0,inputs_sm3_kfull);
  // LD_R_matrix Rmat_sm3_kfull(fspace0,Ssm3_kfull);
  // Rmat_sm3_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_sm3_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_smooth3_kfull.lddat");
  // Rmat_sm3_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_smooth3_kfull.lddat");
  // // Rmat_sm3_kfull.write_matrix("./writing_dat_dir/Duffing_Legendre.10.Rmat_nois_smooth3_kfull.lddat");

  // input_ode_observations inputs_br3_kfull("./writing_dat_dir/Duffing_nois_obs_brgman3_kfull.lddat");
  // inputs_br3_kfull.print_details();
  // LD_observations_set Sbr3_kfull(meta0,inputs_br3_kfull);
  // LD_R_matrix Rmat_br3_kfull(fspace0,Sbr3_kfull);
  // Rmat_br3_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_br3_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_brgman3_kfull.lddat");
  // Rmat_br3_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_brgman3_kfull.lddat");

  // input_ode_observations inputs_prjsm3_kfull("./writing_dat_dir/Duffing_nois_obs_prjsmooth3_kfull.lddat");
  // inputs_prjsm3_kfull.print_details();
  // LD_observations_set Sprjsm3_kfull(meta0,inputs_prjsm3_kfull);
  // LD_R_matrix Rmat_prjsm3_kfull(fspace0,Sprjsm3_kfull);
  // Rmat_prjsm3_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // Rmat_prjsm3_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_prjsmooth3_kfull.lddat");


  // input_ode_observations inputs_sm4_kfull("./writing_dat_dir/Duffing_nois_obs_smooth4_kfull.lddat");
  // inputs_sm4_kfull.print_details();
  // LD_observations_set Ssm4_kfull(meta0,inputs_sm4_kfull);
  // LD_R_matrix Rmat_sm4_kfull(fspace0,Ssm4_kfull);
  // Rmat_sm4_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_sm4_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_smooth4_kfull.lddat");
  // Rmat_sm4_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_smooth4_kfull.lddat");
  // // Rmat_sm4_kfull.write_matrix("./writing_dat_dir/Duffing_Legendre.10.Rmat_nois_smooth4_kfull.lddat");


  // input_ode_observations inputs_sm5_kfull("./writing_dat_dir/Duffing_nois_obs_smooth5_kfull.lddat");
  // inputs_sm5_kfull.print_details();
  // LD_observations_set Ssm5_kfull(meta0,inputs_sm5_kfull);
  // LD_R_matrix Rmat_sm5_kfull(fspace0,Ssm5_kfull);
  // Rmat_sm5_kfull.populate_R_matrix<orthopolynomial_basis>(bases0);
  // // Rmat_sm4_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.8.Rmat_nois_smooth5_kfull.lddat");
  // Rmat_sm5_kfull.write_matrix("./writing_dat_dir/Duffing_Chebyshev.10.Rmat_nois_smooth5_kfull.lddat");
  // // Rmat_sm5_kfull.write_matrix("./writing_dat_dir/Duffing_Legendre.10.Rmat_nois_smooth5_kfull.lddat");


  free_evaluation_bases<orthopolynomial_basis>(nbases0,bases0);
}
