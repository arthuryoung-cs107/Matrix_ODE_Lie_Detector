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

  LD_encoding_bundle Lcode(Sobs.ncrvs_tot,fspace0.perm_len,Sobs.npts_per_crv,1);
  LD_L_encoder::encode_L_bundle<orthopolynomial_basis>(Lcode,Sobs,bases0);
  LD_svd_bundle Lcode_svd(Lcode); Lcode_svd.print_details("Lcode_svd");
  // LD_matrix_svd_result Lmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L"))); Lmat_svd_check.print_details("Lmat_svd_check");

  LD_encoding_bundle Gcode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep);
  LD_G_encoder::encode_G_bundle<orthopolynomial_basis>(Gcode,Sobs,bases0);
  LD_svd_bundle Gcode_svd(Gcode); Gcode_svd.print_details("Gcode_svd");
  // LD_matrix_svd_result Gmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"G"))); Gmat_svd_check.print_details("Gmat_svd_check");

  LD_encoding_bundle Rncode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep*meta0.eor);
  LD_R_encoder::encode_Rn_bundle<orthopolynomial_basis>(Rncode,Sobs,bases0,meta0.eor);
  LD_svd_bundle Rncode_svd(Rncode); Rncode_svd.print_details("Rncode_svd");
  // LD_matrix_svd_result Rmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R"))); Rmat_svd_check.print_details("Rmat_svd_check");

  LD_encoding_bundle Qcode(Sobs.ncrvs_tot,fspace0.ndof_full,Sobs.npts_per_crv,meta0.ndep*(meta0.eor+1));
  LD_R_encoder::encode_Q_bundle<orthopolynomial_basis>(Qcode,Sobs,bases0);
  LD_svd_bundle Qcode_svd(Qcode); Qcode_svd.print_details("Qcode_svd");
  // LD_matrix_svd_result Qmat_svd_check(LD_svd_file(name.name_svd_file(obs_name,fam_name,"Q"))); Qmat_svd_check.print_details("Qmat_svd_check");

  LD_matrix Lmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"L")));
  LD_matrix_svd_result Lmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L"))); Lmat_svd.print_details("Lmat_svd");
    LD_vspace_record  L_Vrec(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns),
                      L_Vrec_svd(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns),
                      L_Vrec_rank(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns),
                      L_Vrec_minrank(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns),
                      L_Vrec_ss(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns);
                        L_Vrec_rank.set_record_rankvec(Lmat_svd.rank_vec);
                        L_Vrec_minrank.set_record_unirank(Lmat_svd.min_rank());
                        LD_L_matrix::eval_L_signal_strength(L_Vrec_ss,L_Vrec_svd,fspace0,Sobs,1e-10);
    // L_Vrec.copy_record(L_Vrec_rank); const char L_Vrec_name[] = "L_Vrec_rank";
    // L_Vrec.copy_record(L_Vrec_minrank); const char L_Vrec_name[] = "L_Vrec_minrank";
    L_Vrec.copy_record(L_Vrec_ss); const char L_Vrec_name[] = "L_Vrec_ss";
      L_Vrec.print_selected_details("L",false);
      L_Vrec_svd.compare_subspaces(L_Vrec_rank,"L_Vrec_rank",L_Vrec,L_Vrec_name);
    LD_vector_bundle L_Ybndl(L_Vrec);

  LD_matrix Rmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"R"))); double tol_use_R = 1e-5;
    LD_matrix_svd_result Rmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R"))); Rmat_svd.print_details("Rmat_svd");
      LD_vspace_record  R_Vrec(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.ncols,Rmat_svd.VTtns),
                        R_Vrec_svd(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.ncols,Rmat_svd.VTtns),
                        R_Vrec_null(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.ncols,Rmat_svd.VTtns),
                        R_Vrec_minnull(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.ncols,Rmat_svd.VTtns),
                        R_Vrec_nsat(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.ncols,Rmat_svd.VTtns);
                          R_Vrec_null.set_record_nullspc(Rmat_svd.rank_vec);
                          R_Vrec_minnull.set_record_uninull(Rmat_svd.min_nulldim());
                          LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(R_Vrec_nsat,R_Vrec_nsat,Sobs,bases0,tol_use_R);
      // R_Vrec.copy_record(R_Vrec_null); const char R_Vrec_name[] = "R_Vrec_null";
      // R_Vrec.copy_record(R_Vrec_minnull); const char R_Vrec_name[] = "R_Vrec_minnull";
      R_Vrec.copy_record(R_Vrec_nsat); const char R_Vrec_name[] = "R_Vrec_nsat";
        R_Vrec.print_selected_details("R",false); R_Vrec_svd.compare_subspaces(R_Vrec,R_Vrec_name,R_Vrec_null,"R_Vrec_null");
      LD_vector_bundle R_Kbndl(R_Vrec);

    LD_Theta_bundle RYL_Tbndl_x(meta0,Sobs.ncrvs_tot,Rmat.net_cols); RYL_Tbndl_x.set_Yspaces(L_Ybndl,0);
      LD_svd_bundle RYLmat_x_svd(Rmat,RYL_Tbndl_x); RYLmat_x_svd.print_details("RYLmat_x_svd");
        LD_vspace_record  &RYL_Trec_x = RYL_Tbndl_x.Trec,
                          &RYL_Vrec_x_svd = RYLmat_x_svd.rec,
                          RYL_Vrec_x(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_x),
                          RYL_Vrec_x_null(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_x),
                          RYL_Vrec_x_minnul(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_x),
                          RYL_Vrec_x_nsat(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_x);
                            RYL_Vrec_x_null.set_record_nullspc(RYLmat_x_svd.rank_vec);
                            RYL_Vrec_x_minnul.set_record_uninull(RYLmat_x_svd.min_nulldim());
                            LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(RYL_Vrec_x_nsat,RYL_Vrec_x_nsat,Sobs,bases0,tol_use_R);
        // RYL_Vrec_x.copy_record(RYL_Vrec_x_null); const char RYL_Vrec_x_name[] = "RYL_Vrec_x_null";
        // RYL_Vrec_x.copy_record(RYL_Vrec_x_minnul); const char RYL_Vrec_x_name[] = "RYL_Vrec_x_minnul";
        RYL_Vrec_x.copy_record(RYL_Vrec_x_nsat); const char RYL_Vrec_x_name[] = "RYL_Vrec_x_nsat";
          RYL_Vrec_x.print_selected_details("RYL_x",false);
          RYL_Vrec_x_svd.compare_subspaces(RYL_Vrec_x,RYL_Vrec_x_name,RYL_Vrec_x_null,"RYL_Vrec_x_null");
        RYL_Tbndl_x.set_Tspaces(RYL_Vrec_x);

    LD_Theta_bundle RYL_Tbndl_xu(meta0,Sobs.ncrvs_tot,Rmat.net_cols); RYL_Tbndl_xu.set_Yspaces(L_Ybndl,-1);
      LD_svd_bundle RYLmat_xu_svd(Rmat,RYL_Tbndl_xu); RYLmat_xu_svd.print_details("RYLmat_xu_svd");
        LD_vspace_record  &RYL_Trec_xu = RYL_Tbndl_xu.Trec,
                          &RYL_Vrec_xu_svd = RYLmat_xu_svd.rec,
                          RYL_Vrec_xu(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_xu),
                          RYL_Vrec_xu_null(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_xu),
                          RYL_Vrec_xu_minnul(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_xu),
                          RYL_Vrec_xu_nsat(Rmat.ncrvs_tot,Rmat.net_cols,Rmat.net_cols,RYL_Trec_xu);
                            RYL_Vrec_xu_null.set_record_nullspc(RYLmat_xu_svd.rank_vec);
                            RYL_Vrec_xu_minnul.set_record_uninull(RYLmat_xu_svd.min_nulldim());
                            LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(RYL_Vrec_xu_nsat,RYL_Vrec_xu_nsat,Sobs,bases0,tol_use_R);
        // RYL_Vrec_xu.copy_record(RYL_Vrec_xu_null); const char RYL_Vrec_xu_name[] = "RYL_Vrec_xu_null";
        // RYL_Vrec_xu.copy_record(RYL_Vrec_xu_minnul); const char RYL_Vrec_xu_name[] = "RYL_Vrec_xu_minnul";
        RYL_Vrec_xu.copy_record(RYL_Vrec_xu_nsat); const char RYL_Vrec_xu_name[] = "RYL_Vrec_xu_nsat";
          RYL_Vrec_xu.print_selected_details("RYL_xu",false);
          RYL_Vrec_xu_svd.compare_subspaces(RYL_Vrec_xu,RYL_Vrec_xu_name,RYL_Vrec_xu_null,"RYL_Vrec_xu_null");
        RYL_Tbndl_xu.set_Tspaces(RYL_Vrec_xu);

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
