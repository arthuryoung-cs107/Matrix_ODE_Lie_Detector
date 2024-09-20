#include "LD_infinitesimal_generators.hh"
#include "LD_integrators.hh"
#include "LD_encodings.hh"
#include "LD_vspace_evaluators.hh"
// #include "solution_space_cokernals.hh"
#include "LD_cokernals.hh"

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

  LD_encoding_bundle OGcode(Ocode,Gcode);
  LD_svd_bundle OGcode_svd(OGcode); // OGcode_svd.print_details("OGcode_svd");

  // double  tol_use_L = 1e-10,
  //         tol_use_G = 1e-8,
  //         atol_use_R = 1e-12,
  //         rtol_use_R = 1e-08;
  double  tol_use_L = 1e-10,
          tol_use_G = 1e-06,
          // atol_use_R = 1e-06,
          atol_use_R = 1e-08,
          // atol_use_R = 1e-09,
          rtol_use_R = 0e-00;
          // rtol_use_R = 1e-08;
  const char rec_suf[] = ".DoP853rnrec";
  // const char rec_suf[] = "_1.DoP853rnrec";

  LD_vspace_record  L_rec0(Lcode_svd.rec);
  LD_vector_bundle LYbndl(L_rec0);
  L_vspace_eval::evaluate_lambda_signal_strength<LD_vector_bundle,orthopolynomial_basis>(LYbndl,L_rec0,Sobs,bases0,tol_use_L);
    L_rec0.print_subspace_details(LYbndl.rec,"LYbndl.rec");

  LD_vector_bundle  OGbndl0(OGcode_svd.rec);
    OGbndl0.rec.set_record_nullspc(OGcode_svd.rank_vec);
      OGbndl0.set_Vspaces();
        OGcode_svd.rec.print_subspace_details(OGbndl0.rec,"OGbndl0");

    LD_vector_bundle  OGbndl(OGcode_svd.rec);
    LD_vspace_record  OG_rec0(OGbndl.rec);
      OG_vspace_eval::evaluate_nthcond_infcrit(OGbndl,OG_rec0,Sobs,bases0,atol_use_R,rtol_use_R,tol_use_G);
    OG_rec0.print_subspace_details(OGbndl.rec,"OGbndl");

    LD_Theta_bundle OGYL_Tbndl_x(meta0,Sobs.ncrvs_tot,fspace0.ndof_full,LYbndl,0);
      LD_svd_bundle::project_Theta_bundle(OGYL_Tbndl_x,OGcode);
    LD_vspace_record  OGYL_x_rec0(OGYL_Tbndl_x.rec);
      OG_vspace_eval::evaluate_nthcond_infcrit(OGYL_Tbndl_x,OGYL_x_rec0,Sobs,bases0,atol_use_R,rtol_use_R,tol_use_G);
    OGYL_x_rec0.print_subspace_details(OGYL_Tbndl_x.rec,"OGYL_x");

    LD_Theta_bundle OGYL_Tbndl_xu(meta0,Sobs.ncrvs_tot,fspace0.ndof_full,LYbndl,-1);
      LD_svd_bundle::project_Theta_bundle(OGYL_Tbndl_xu,OGcode);
    LD_vspace_record  OGYL_xu_rec0(OGYL_Tbndl_xu.rec);
      OG_vspace_eval::evaluate_nthcond_infcrit(OGYL_Tbndl_xu,OGYL_xu_rec0,Sobs,bases0,atol_use_R,rtol_use_R,tol_use_G);
    OGYL_xu_rec0.print_subspace_details(OGYL_Tbndl_xu.rec,"OGYL_xu");

  OGbndl.rec.compare_subspaces(OGYL_Tbndl_x.rec,"OGYL_x",OGYL_Tbndl_xu.rec,"OGYL_xu");

  printf("\ntesting medoids\n\n");

  Frobenius_vspace_measure OGmsr_OGsat(OGbndl); OGmsr_OGsat.init_Frobenius_distances(OGbndl);
  // indepcomp_vspace_measure OGmsr_OGsat(OGbndl); OGmsr_OGsat.init_indepcomp_distances(OGbndl);
    k_medoids_package k_med_OG_OGsat(OGmsr_OGsat.dsym,OGmsr_OGsat.nset);
    int k_SC = k_med_OG_OGsat.comp_kSC_medoids();
    // const int kmin=2,kmax=10,klen=kmax-kmin+1;
    // int k_SC = k_med_OG_OGsat.comp_kSC_krange_medoids(kmin,kmax);
    // double SCs_OG_OGsat[klen]; k_med_OG_OGsat.comp_SC_krange_medoids(SCs_OG_OGsat,kmin,kmax);
    // for (size_t k = kmin, ik=0; k <= kmax; k++, ik++) printf("(%d, %e) ", k,SCs_OG_OGsat[ik]); printf("\n");

  // // debugging dataset
  // LD_encoding_bundle Ocode_check(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  // LD_R_encoder::encode_O_bundle<orthopolynomial_basis>(Ocode_check,Sobs,bases0,false);
  // LD_svd_bundle Ocode_check_svd(Ocode_check); Ocode_check_svd.print_details("Ocode_check_svd");
  // LD_vector_bundle Obndl_check_maxnull(Ocode_check_svd.rec);
  //   Obndl_check_maxnull.rec.set_record_uninull(Ocode_check_svd.max_nulldim()); Obndl_check_maxnull.set_Vspaces();
  //   Frobenius_vspace_measure Omsr_maxnull(Obndl_check_maxnull); Omsr_maxnull.init_Frobenius_distances(Obndl_check_maxnull);
  //     k_medoids_package k_med_O_maxnull(Omsr_maxnull.dsym,Omsr_maxnull.nset);
  //     // double sc3 = k_med_O_maxnull.comp_k_medoids(3); printf("sc3 = %.4f \n", sc3);
  //     // int k_SC = k_med_O_maxnull.find_k_SC_krange(10); printf("k_SC = %d, SC = %.4f \n", k_SC, k_med_O_maxnull.silh_coeff);
  //     // int k_SC = k_med_O_maxnull.find_kSC(); printf("k_SC = %d, SC = %.4f \n", k_SC, k_med_O_maxnull.silh_coeff);
  //     // int k_SC = k_med_O_maxnull.comp_kSC_medoids();
  //     int k_SC = k_med_O_maxnull.comp_kSC_krange_medoids(2,10);

  printf("\ntesting cokernals\n\n");

  const bool wdistance = true;
  // const bool wdistance = false;

  // double  atol_R_cok = atol_use_R,
  double  atol_R_cok = 1e-12,
          // rtol_R_cok = rtol_use_R,
          rtol_R_cok = 1e-6,
          tol_G_cok = tol_use_G;
          // tol_G_cok = 1e-3;

  LD_OG_encoder OGenc(meta0);

  // Jet_function_vector_space jfvs(Sobs,fspace0,OGenc,msr);
   // jfvs.encode_decompose_bundle<orthopolynomial_basis>(bases0,normalize_flag);
  // Jet_function_vector_space jfvs(Sobs,fspace0,OGenc,msr,bases0,normalize_flag);

  Jet_function_vector_space jfvs(Sobs,fspace0,OGenc,bases0,normalize_flag);

  // Frobenius_vspace_measure msr(fspace0.ndof_full,Sobs.ncrvs_tot);
  indepcomp_vspace_measure msr(fspace0.ndof_full,Sobs.ncrvs_tot);

  nullspace_clst_policy pol(msr,wdistance);
  // nullspace_near_policy pol(msr,wdistance);

  cokernal_sub_bundle ckrn(jfvs,pol,true);

  // OG_vspace_eval OG_evl(Sobs,fspace0.ndof_full,atol_R_cok,rtol_R_cok,tol_G_cok,false);
    // jfvs.evaluate_Vbndle0<OG_vspace_eval,orthopolynomial_basis>(OG_evl,bases0,true);
    // jfvs.init_encoded_clustering<OG_vspace_eval,orthopolynomial_basis>(OG_vspace_eval,);

  // LD_vector_bundle  OGbndl_cok(OGcode_svd.rec);
  // LD_vector_bundle  OGbndl_cok(OGYL_x_rec0);
  // LD_vector_bundle  OGbndl_cok(OGYL_xu_rec0);
  // OG_vspace_eval OG_evl(Sobs,OGbndl_cok.vlen_full,atol_R_cok,rtol_R_cok,tol_G_cok);
  // Frobenius_vspace_measure OGmsr_cok(OGbndl_cok);

  // solution_space_cokernals solspc_cok(Sobs.ncrvs_tot,Sobs.nobs);
  // solspc_cok.compute_solspc_cokernals<OG_vspace_eval,orthopolynomial_basis>(OGbndl_cok,Sobs,OG_evl,OGmsr_cok,bases0);

  bool reintegrate_data = false;
  if (reintegrate_data)
  {
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
  }

  printf("ending main\n"); free_evaluation_bases<orthopolynomial_basis>(bases0); return 0;
}
