#include "LD_matrices.hh"
#include "matrix_Lie_detector.hh"
#include "LD_infinitesimal_generators.hh"

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
  const int bor = 9;
  // const int bor = 10;

  LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",bor);
  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(name.name_domain_config_file(obs_name,fam_name)));
  fspace0.debugging_description();

  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);

  LD_matrix Lmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"L")));
  LD_matrix_svd_result Lmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L"))); // Lmat_svd.print_details("Lmat_svd");
  LD_vspace_record L_Vrec(Lmat_svd.ncrvs,Lmat_svd.ncol_use,Lmat_svd.ncols,Lmat_svd.VTtns);
  // double tol_use_L = 1e-12; LD_L_matrix::eval_L_signal_strength(Sobs,L_Vrec,fspace0,tol_use_L); // L_Vrec.print_selected_details("L");
  double tol_use_L = 0.0; LD_L_matrix::eval_L_signal_strength(Sobs,L_Vrec,fspace0,tol_use_L); L_Vrec.print_selected_details("L",false);
  int nnsat_s_L = matrix_Lie_detector::compare_relaxed_subspaces(Lmat_svd,L_Vrec,"L",tol_use_L);
  LD_vector_bundle L_Ybndl(L_Vrec);

  LD_matrix Rmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"R")));

  double tol_use_R = 1e-1;

  LD_matrix_svd_result Rmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R"))); Rmat_svd.print_details("Rmat_svd");
  LD_vector_bundle R_Vbndl(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.VTtns); LD_vspace_record &Rrec = R_Vbndl.rec;
  LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(R_Vbndl.rec,R_Vbndl,Sobs,bases0,tol_use_R); R_Vbndl.rec.print_selected_details("R",false);

  LD_matrix_svd_result RYLmat_svd0(Sobs.ncrvs_tot,fspace0.ndof_full);
  LD_Theta_bundle RYL_Tbndl0(meta0,Sobs.ncrvs_tot,fspace0.ndof_full); LD_vspace_record &RYLrec = RYL_Tbndl0.Vbndle.rec;
  RYL_Tbndl0.set_Yspaces(L_Ybndl,-1); // constrain just x parameter space
  matrix_Lie_detector::compute_AYmat_curve_svds(RYLmat_svd0,Rmat,RYL_Tbndl0); // RYLmat_svd0.print_details("RYLmat_svd");
  RYL_Tbndl0.init_Vbndle_premult(RYLmat_svd0.VTtns); // RYL_Tbndl0.Vbndle.rec.print_selected_details("RYL");
  LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(RYLrec,RYL_Tbndl0.Vbndle,Sobs,bases0,tol_use_R); RYLrec.print_selected_details("RYL",false);

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;

  // const int M_mat1 = 25,
  //           N_mat1 = 10,
  //           isub = 1,
  //           jsub = 1,
  //           M_mat2 = LD_linalg::min_T<int>(M_mat1 - isub,5),
  //           N_mat2 = LD_linalg::min_T<int>(N_mat1 - jsub,4);
  // double ** const mat_1 = Tmatrix<double>(M_mat1,N_mat1);
  // gsl_matrix * const mat_gsl1 = gsl_matrix_alloc(M_mat1,N_mat1);
  //
  // LD_linalg::fill_vec_012(mat_1[0],M_mat1*N_mat1);
  // LD_gsl::load_gsl_mat(mat_gsl1,mat_1);
  //
  // LD_linalg::print_A("mat_1",mat_1,M_mat1,N_mat1);
  // LD_gsl::print_A_gsl("mat_gsl1",mat_gsl1);
  //
  // {
  //   gsl_matrix_view mat_gsl2 = gsl_matrix_view_array(mat_1[0] + ((isub*N_mat1) + jsub),M_mat2,N_mat2),
  //                   mat_gsl2_gsl = gsl_matrix_submatrix(mat_gsl1,isub,jsub,M_mat2,N_mat2);
  //   gsl_matrix  submat_gsl_1 = mat_gsl2.matrix,
  //               submat_gsl_2 = mat_gsl2_gsl.matrix;
  //   LD_gsl::print_A_gsl("submat_gsl_1",&submat_gsl_1);
  //   LD_gsl::print_A_gsl("submat_gsl_2",&submat_gsl_2);
  // }
  //
  // free_Tmatrix<double>(mat_1);
  // gsl_matrix_free(mat_gsl1);

  // LD_matrix_svd_result Rmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"R")));
  // Rmat_svd.print_details("Rmat_svd");
  // LD_vector_bundle R_Vbndl(Rmat_svd.ncrvs,Rmat_svd.ncol_use,Rmat_svd.VTtns);
  //
  // if (false)
  // {
  //   const char Amat_letter[] = "R"; const char Amat_svd_name[] = "Rmat_svd"; double tol_use_A = 1e-3;
  //   // const char Amat_letter[] = "G"; const char Amat_svd_name[] = "Gmat_svd"; double tol_use_A = 1e-6;
  //
  //   LD_matrix Amat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,Amat_letter)));
  //   LD_matrix_svd_result Amat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,Amat_letter)));
  //   Amat_svd.print_details(Amat_svd_name);
  //
  //   LD_vspace_record A_Vrec(Amat_svd.ncrvs,Amat_svd.ncol_use,Amat_svd.ncols,Amat_svd.VTtns);
  //
  //   LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(Sobs,A_Vrec,bases0,tol_use_A);
  //   // LD_G_matrix::eval_inf_criterion<orthopolynomial_basis>(Sobs,A_Vrec,bases0,tol_use_A);
  //
  //   A_Vrec.print_selected_details(Amat_letter);
  //   bool nsat_bf_A[Sobs.ncrvs_tot];
  //   int nnsat_b_A = 0,
  //       nnsat_diff_A[Sobs.ncrvs_tot];
  //   printf("  %s nulldim\n", Amat_letter);
  //   for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++)
  //   {
  //     int nulldim_i = Amat_svd.ncol_use-Amat_svd.rank_vec[icrv];
  //     printf("%d ", nulldim_i);
  //     nnsat_diff_A[icrv] = A_Vrec.nV_spcvec[icrv]-nulldim_i;
  //     nnsat_b_A+=(int)(nsat_bf_A[icrv]=(nnsat_diff_A[icrv] >= 0));
  //   }
  //   printf("\n");
  //   printf("nsat_bigger_A = %d (of %d, tol = %.1e) \n", nnsat_b_A,Sobs.ncrvs_tot,tol_use_A);
  //   for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++) printf("%d ", (int) nsat_bf_A[icrv]);
  //   printf("\n");
  //   LD_linalg::print_xT("nsat_diff_A",nnsat_diff_A,Sobs.ncrvs_tot);
  // }

  // free_evaluation_bases<orthopolynomial_basis>(bases0);
  // return 0;
}
