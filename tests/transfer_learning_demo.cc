#include "LD_matrices.hh"
#include "LD_infinitesimal_generators.hh"
#include "LD_integrators.hh"
#include "LD_encodings.hh"
#include "LD_vspace_evaluators.hh"

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

  // const int bor = 8;
  // const int bor = 9;
  const int bor = 10;

  LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",bor);
  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(nbuf0.name_domain_config_file(obs_name,fam_name)));
  fspace0.debugging_description();

  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);
  orthopolynomial_basis &basis0 = *(bases0[0]);

  bool normalize_flag = true;
  // bool normalize_flag = false;

  LD_encoding_bundle Lcode(Sobs.ncrvs_tot,fspace0.perm_len,Sobs.npts_per_crv,1);
  LD_L_encoder::encode_L_bundle<orthopolynomial_basis>(Lcode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Lcode_svd(Lcode); // Lcode_svd.print_details("Lcode_svd");

  LD_encoding_bundle Gcode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  LD_G_encoder::encode_G_bundle<orthopolynomial_basis>(Gcode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Gcode_svd(Gcode); // Gcode_svd.print_details("Gcode_svd");

  LD_encoding_bundle Ocode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  LD_R_encoder::encode_O_bundle<orthopolynomial_basis>(Ocode,Sobs,bases0,normalize_flag);
  LD_svd_bundle Ocode_svd(Ocode); // Ocode_svd.print_details("Ocode_svd");

  LD_encoding_bundle Rncode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep*meta0.eor);
  LD_R_encoder::encode_Rn_bundle<orthopolynomial_basis>(Rncode,Sobs,bases0,meta0.eor,normalize_flag);
  LD_svd_bundle Rncode_svd(Rncode); // Rncode_svd.print_details("Rncode_svd");

  LD_encoding_bundle OGcode(Ocode,Gcode);
  LD_svd_bundle OGcode_svd(OGcode); // OGcode_svd.print_details("OGcode_svd");

  // double  tol_use_L = 1e-10,
  //         tol_use_G = 1e-8,
  //         atol_use_R = 1e-12,
  //         rtol_use_R = 1e-08;
  double  tol_use_L = 1e-10,
          tol_use_G = 1e-06,
          atol_use_R = 1e-08,
          // atol_use_R = 1e-09,
          // rtol_use_R = 0e-00;
          rtol_use_R = 1e-08;
  // const char rec_suf[] = ".DoP853rnrec";
  const char rec_suf[] = "_1.DoP853rnrec";

  LD_vspace_record  L_rec0(Lcode_svd.rec);
  LD_vector_bundle LYbndl(L_rec0);
  L_vspace_eval::evaluate_lambda_signal_strength<LD_vector_bundle,orthopolynomial_basis>(LYbndl,L_rec0,Sobs,bases0,tol_use_L);
    L_rec0.print_subspace_details(LYbndl.rec,"LYbndl.rec");

  LD_vector_bundle  OGbndl0(OGcode_svd.rec);
    OGbndl0.rec.set_record_nullspc(OGcode_svd.rank_vec);
      OGbndl0.set_Vspaces();
        OGcode_svd.rec.print_subspace_details(OGbndl0.rec,"OGbndl0.rec");

  LD_vector_bundle  OGbndl(OGcode_svd.rec);
  LD_vspace_record  OG_rec0(OGbndl.rec),
                    OG_Osat_rec(OG_rec0),
                    OG_Gsat_rec(OG_rec0);
  Rk_vspace_eval::evaluate_kth_ratio_condition<orthopolynomial_basis>(OG_Osat_rec,OG_rec0,
    Sobs,bases0,atol_use_R,rtol_use_R);
  G_vspace_eval::evaluate_infinitesimal_criterion<orthopolynomial_basis>(OG_Gsat_rec,OG_rec0,
    Sobs,bases0,tol_use_G);
  LD_vspace_record  OG_OGsat_rec(OG_Osat_rec,OG_Gsat_rec); OGbndl.set_Vspaces(OG_OGsat_rec);
  OG_rec0.print_subspace_details(OGbndl.rec,"OGbndl.rec");

  LD_Theta_bundle OGYL_Tbndl_x(meta0,Sobs.ncrvs_tot,fspace0.ndof_full,LYbndl,0);
  LD_svd_bundle::project_Theta_bundle(OGYL_Tbndl_x,OGcode); // LD_svd_bundle OGYL_x_svd(OGcode,OGYL_Tbndl_x,true);
  LD_vspace_record  OGYL_x_rec0(OGYL_Tbndl_x.rec),
                    OGYL_x_Osat_rec(OGYL_x_rec0),
                    OGYL_x_Gsat_rec(OGYL_x_rec0);
  Rk_vspace_eval::evaluate_kth_ratio_condition<orthopolynomial_basis>(OGYL_x_Osat_rec,OGYL_x_rec0,
    Sobs,bases0,atol_use_R,rtol_use_R);
  G_vspace_eval::evaluate_infinitesimal_criterion<orthopolynomial_basis>(OGYL_x_Gsat_rec,OGYL_x_rec0,
    Sobs,bases0,tol_use_G);
  LD_vspace_record  OGYL_x_OGsat_rec(OGYL_x_Osat_rec,OGYL_x_Gsat_rec); OGYL_Tbndl_x.set_Vspaces(OGYL_x_OGsat_rec);
  OGYL_x_rec0.print_subspace_details(OGYL_Tbndl_x.rec,"OGYL_Tbndl_x.rec");

  LD_Theta_bundle OGYL_Tbndl_xu(meta0,Sobs.ncrvs_tot,fspace0.ndof_full,LYbndl,-1);
  LD_svd_bundle::project_Theta_bundle(OGYL_Tbndl_xu,OGcode); // LD_svd_bundle OGYL_xu_svd(OGcode,OGYL_Tbndl_xu,true);
  LD_vspace_record  OGYL_xu_rec0(OGYL_Tbndl_xu.rec),
                    OGYL_xu_Osat_rec(OGYL_xu_rec0),
                    OGYL_xu_Gsat_rec(OGYL_xu_rec0);
  Rk_vspace_eval::evaluate_kth_ratio_condition<orthopolynomial_basis>(OGYL_xu_Osat_rec,OGYL_xu_rec0,
    Sobs,bases0,atol_use_R,rtol_use_R);
  G_vspace_eval::evaluate_infinitesimal_criterion<orthopolynomial_basis>(OGYL_xu_Gsat_rec,OGYL_xu_rec0,
    Sobs,bases0,tol_use_G);
  LD_vspace_record  OGYL_xu_OGsat_rec(OGYL_xu_Osat_rec,OGYL_xu_Gsat_rec); OGYL_Tbndl_xu.set_Vspaces(OGYL_xu_OGsat_rec);
  OGYL_xu_rec0.print_subspace_details(OGYL_Tbndl_xu.rec,"OGYL_Tbndl_xu.rec");

  OG_OGsat_rec.compare_subspaces(OGYL_x_OGsat_rec,"OGYL_x_OGsat_rec",OGYL_xu_OGsat_rec,"OGYL_xu_OGsat_rec");

  // LD_vector_bundle  Rnbndl(Rncode_svd.rec);
  //   Rn_vspace_eval::evaluate_nth_ratio_condition<LD_vector_bundle,orthopolynomial_basis>(Rnbndl,Rnbndl.rec,Sobs,bases0,atol_use_R,rtol_use_R);
  // r_xu_infgen Rn_rxu_ign(fspace0,Rnbndl.Vspaces); DoP853 intgr_rec_rxu(Rn_rxu_ign,intgr_rec_settings);
  // generated_ode_observations gen_Rn_rxu(Rn_rxu_ign,Sobs.ncrvs_tot,Sobs.min_npts_curve());
  // gen_Rn_rxu.set_solcurve_ICs(Sobs.curves);
  // gen_Rn_rxu.parallel_generate_solution_curves<r_xu_infgen,DoP853>(Rn_rxu_ign,intgr_rec_rxu,Sobs.get_default_IC_indep_range());

  rn_infgen OG_rn_ign(basis0,OGbndl0.Vspaces),
            OG_rn_OGsat_ign(basis0,OGbndl.Vspaces),
            OGYL_x_rn_OGsat_ign(basis0,OGYL_Tbndl_x.Vspaces),
            OGYL_xu_OGsat_rn_ign(basis0,OGYL_Tbndl_xu.Vspaces);
  DoP853_settings intgr_rec_settings; DoP853 intgr_rec_rn(OG_rn_ign,intgr_rec_settings);

  generated_ode_observations gen_OG_rn(OG_rn_ign,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  gen_OG_rn.set_solcurve_ICs(Sobs.curves);
  gen_OG_rn.parallel_generate_solution_curves<orthopolynomial_basis,rn_infgen,DoP853>(bases0,OG_rn_ign,intgr_rec_rn);

  generated_ode_observations gen_OG_OGsat_rn(OG_rn_OGsat_ign,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  gen_OG_OGsat_rn.set_solcurve_ICs(Sobs.curves);
  gen_OG_OGsat_rn.parallel_generate_solution_curves<orthopolynomial_basis,rn_infgen,DoP853>(bases0,OG_rn_OGsat_ign,intgr_rec_rn);

  generated_ode_observations gen_OGYL_x_OGsat_rn(OGYL_x_rn_OGsat_ign,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  gen_OGYL_x_OGsat_rn.set_solcurve_ICs(Sobs.curves);
  gen_OGYL_x_OGsat_rn.parallel_generate_solution_curves<orthopolynomial_basis,rn_infgen,DoP853>(bases0,OGYL_x_rn_OGsat_ign,intgr_rec_rn);

  generated_ode_observations gen_OGYL_xu_OGsat_rn(OGYL_xu_OGsat_rn_ign,Sobs.ncrvs_tot,Sobs.npts_per_crv);
  gen_OGYL_xu_OGsat_rn.set_solcurve_ICs(Sobs.curves);
  gen_OGYL_xu_OGsat_rn.parallel_generate_solution_curves<orthopolynomial_basis,rn_infgen,DoP853>(bases0,OGYL_xu_OGsat_rn_ign,intgr_rec_rn);

  nbuf2.load_name(obs_name.name,".",fam_name.name);
  gen_OG_rn.write_observed_solutions(nbuf0.name_file(nbuf1.load_name(nbuf2.name,".OG_null",".DoP853rnrec")));

  gen_OG_OGsat_rn.write_observed_solutions(nbuf0.name_file(nbuf1.load_name(nbuf2.name,".OG_OGsat",rec_suf)));
  gen_OGYL_x_OGsat_rn.write_observed_solutions(nbuf0.name_file(nbuf1.load_name(nbuf2.name,".OG_YLx_OGsat",rec_suf)));
  gen_OGYL_xu_OGsat_rn.write_observed_solutions(nbuf0.name_file(nbuf1.load_name(nbuf2.name,".OG_YLxu_OGsat",rec_suf)));

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;
}
