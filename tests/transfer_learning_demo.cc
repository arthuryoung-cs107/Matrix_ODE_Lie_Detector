#include "LD_matrices.hh"
#include "matrix_Lie_detector.hh"
#include "LD_infinitesimal_generators.hh"
#include "LD_integrators.hh"

// specify data directory for writing binary files
const char dir_name[] = "./data_directory";
const char dat_suff[] = "lddat";

// const char data_name[] = "Duffing_xrange0_true_DoP853gen";

const char eqn_name[] = "Duffing"; ode_solspc_meta meta0(2,1);

LD_name_buffer name(dir_name,dat_suff,100);
LD_name_buffer name_dnp1xu(dir_name,dat_suff,100);
LD_name_buffer name_JFs(dir_name,dat_suff,100);
int main()
{
  // load observational data
  LD_name_buffer obs_name; obs_name.name_observational_data(eqn_name,0,-1,"DoP853");
  // LD_observations_set Sobs(meta0,ode_curve_observations(name.name_file(obs_name)));
  // LD_observations_set Sobs(meta0,ode_curve_observations(  name.name_file(obs_name),
                                                          // name_dnp1xu.name_file(obs_name,"_dnp1xu")));
                                                          // name_JFs.name_file(obs_name,"_JFs")));
  LD_observations_set Sobs(meta0,ode_curve_observations(  name.name_file(obs_name),
                                                          name_dnp1xu.name_file(obs_name,"_dnp1xu"),
                                                          name_JFs.name_file(obs_name,"_JFs")));

  // const int bor = 8;
  // const int bor = 9;
  const int bor = 10;

  LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",bor);
  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(name.name_domain_config_file(obs_name,fam_name)));
  fspace0.debugging_description();

  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);
  orthopolynomial_basis &basis0 = *(bases0[0]);

  // bool normalize_flag = true;
  bool normalize_flag = false;

  LD_encoding_bundle Lcode(Sobs.ncrvs_tot,fspace0.perm_len,Sobs.npts_per_crv,1);
  LD_L_encoder::encode_L_bundle<orthopolynomial_basis>(Lcode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Lcode_svd(Lcode); // Lcode_svd.print_details("Lcode_svd");
  // LD_matrix_svd_result Lmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L"))); Lmat_svd_check.print_details("Lmat_svd_check");

  LD_encoding_bundle Gcode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  LD_G_encoder::encode_G_bundle<orthopolynomial_basis>(Gcode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Gcode_svd(Gcode); // Gcode_svd.print_details("Gcode_svd");
  // LD_matrix_svd_result Gmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"G"))); Gmat_svd_check.print_details("Gmat_svd_check");

  LD_encoding_bundle Ocode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  LD_R_encoder::encode_O_bundle<orthopolynomial_basis>(Ocode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Ocode_svd(Ocode); // Ocode_svd.print_details("Ocode_svd");
  // LD_matrix_svd_result Omat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"O"))); Omat_svd_check.print_details("Omat_svd_check");

  LD_encoding_bundle Rncode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep*meta0.eor);
  LD_R_encoder::encode_Rn_bundle<orthopolynomial_basis>(Rncode,Sobs,bases0,meta0.eor,normalize_flag);
  LD_svd_bundle Rncode_svd(Rncode); // Rncode_svd.print_details("Rncode_svd");
  // LD_matrix_svd_result Rmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R"))); Rmat_svd_check.print_details("Rmat_svd_check");

  LD_encoding_bundle OGcode(Ocode,Gcode);
  LD_svd_bundle OGcode_svd(OGcode); // OGcode_svd.print_details("OGcode_svd");
  // LD_matrix_svd_result OGmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"OG"))); OGmat_svd_check.print_details("OGmat_svd_check");

  double  atol_use_R = 1e-10,
          rtol_use_R = 1e-06,
          tol_use_G = 1e-08;

  LD_vspace_record  Ocode_Osat(Ocode_svd.rec),
                    OGcode_Osat(OGcode_svd.rec);

  Rk_vspace_eval::evaluate_kth_ratio_condition<orthopolynomial_basis>(Ocode_Osat,Ocode_svd.rec,Sobs,bases0,atol_use_R,rtol_use_R);
    Ocode_svd.rec.print_subspace_details(Ocode_Osat,"Ocode_Osat");
  Rk_vspace_eval::evaluate_kth_ratio_condition<orthopolynomial_basis>(OGcode_Osat,OGcode_svd.rec,Sobs,bases0,atol_use_R,rtol_use_R);
    OGcode_svd.rec.print_subspace_details(OGcode_Osat,"OGcode_Osat");

  Ocode_svd.rec.compare_subspaces(Ocode_Osat,"Ocode_Osat",OGcode_Osat,"OGcode_Osat");

  // G_vspace_eval::evaluate_infinitesimal_criterion<orthopolynomial_basis>(Gsat_rec,svd_rec,Sobs,bases0,tol_use_G);
  //   svd_rec.print_subspace_details(Gsat_rec,"Gsat_rec");
  // LD_vspace_record OGsat_rec(Osat_rec,Gsat_rec);
  //   svd_rec.print_subspace_details(OGsat_rec,"OGsat_rec");

  
  // r_xu_infgen R_Kbndl_rxu_ign(fspace0,R_Kbndl.Vspaces),
  //             RYL_T_x_rxu_ign(fspace0,RYL_Tbndl_x.Vspaces),
  //             RYL_T_xu_rxu_ign(fspace0,RYL_Tbndl_xu.Vspaces);
  //
  // DoP853_settings intgr_rec_settings; DoP853 intgr_rec(R_Kbndl_rxu_ign,intgr_rec_settings);
  //
  //   generated_ode_observations gen_RK_rxu(R_Kbndl_rxu_ign,Sobs.ncrvs_tot,Sobs.min_npts_curve());
  //   gen_RK_rxu.set_solcurve_ICs(Sobs.curves);
  //   gen_RK_rxu.parallel_generate_solution_curves<r_xu_infgen,DoP853>(R_Kbndl_rxu_ign,intgr_rec,Sobs.get_default_IC_indep_range());
  //
  //   generated_ode_observations gen_RYLT_x_rxu(RYL_T_x_rxu_ign,Sobs.ncrvs_tot,Sobs.min_npts_curve());
  //   gen_RYLT_x_rxu.set_solcurve_ICs(Sobs.curves);
  //   gen_RYLT_x_rxu.parallel_generate_solution_curves<r_xu_infgen,DoP853>(RYL_T_x_rxu_ign,intgr_rec,Sobs.get_default_IC_indep_range());
  //
  //   generated_ode_observations gen_RYLT_xu_rxu(RYL_T_xu_rxu_ign,Sobs.ncrvs_tot,Sobs.min_npts_curve());
  //   gen_RYLT_xu_rxu.set_solcurve_ICs(Sobs.curves);
  //   gen_RYLT_xu_rxu.parallel_generate_solution_curves<r_xu_infgen,DoP853>(RYL_T_xu_rxu_ign,intgr_rec,Sobs.get_default_IC_indep_range());

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;
}
