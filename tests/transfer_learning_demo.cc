#include "LD_matrices.hh"
#include "matrix_Lie_detector.hh"

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
  LD_observations_set Sobs(meta0,ode_curve_observations(name.name_file(obs_name)));
  Sobs.load_additional_inputs(ode_curve_observations( name.name_file(obs_name),
                                                      name_dnp1xu.name_file(obs_name,"_dnp1xu"),
                                                      name_JFs.name_file(obs_name,"_JFs")));

  const int bor = 8;
  // const int bor = 9;
  // const int bor = 10;

  LD_name_buffer fam_name; fam_name.name_function_space("Chebyshev1",bor);
  orthopolynomial_space fspace0(meta0,orthopolynomial_config_file(name.name_domain_config_file(obs_name,fam_name)));
  fspace0.debugging_description();

  orthopolynomial_basis **bases0 = make_evaluation_bases<orthopolynomial_basis, orthopolynomial_space>(fspace0);

  LD_matrix Lmat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,"L")));
  LD_matrix_svd_result Lmat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,"L")));
  Lmat_svd.print_details("Lmat_svd");

  {
    double tol_use_L = 1e-12;
    LD_Theta_stack Theta_stack_L(Sobs.ncrvs_tot,fspace0.perm_len);
    Theta_stack_L.load_Theta_data_tns(Lmat_svd.VTtns);
    LD_L_matrix::eval_L_signal_strength(Sobs,Theta_stack_L,fspace0,tol_use_L);
    Theta_stack_L.print_satisfactory_details("L");
    bool nsat_sf_L[Sobs.ncrvs_tot];
    int nnsat_s_L = 0,
        nnsat_diff_L[Sobs.ncrvs_tot];
    printf("  %s rank\n", "L");
    for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++)
    {
      int rank_i = Lmat_svd.rank_vec[icrv];
      printf("%d ", rank_i);
      nnsat_diff_L[icrv] = rank_i-Theta_stack_L.ntheta_spcvec[icrv];
      nnsat_s_L+=(int)(nsat_sf_L[icrv]=(nnsat_diff_L[icrv] >= 0));
    }
    printf("\n");
    printf("nsat_smaller_L = %d (of %d, tol = %.1e) \n", nnsat_s_L,Sobs.ncrvs_tot,tol_use_L);
    for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++) printf("%d ", (int) nsat_sf_L[icrv]);
    printf("\n");
    LD_linalg::print_xT("nsat_diff_L",nnsat_diff_L,Sobs.ncrvs_tot);
    Theta_stack_L.set_Thetas(fspace0.perm_len);
  }

  {
    // const char Amat_letter[] = "R"; const char Amat_svd_name[] = "Rmat_svd"; double tol_use_A = 1e-3;
    const char Amat_letter[] = "G"; const char Amat_svd_name[] = "Gmat_svd"; double tol_use_A = 1e-6;

    LD_matrix Amat(fspace0,Sobs,LD_matrix_file(name.name_matrix_file(obs_name,fam_name,Amat_letter)));
    LD_matrix_svd_result Amat_svd(LD_svd_file(name.name_svd_file(obs_name,fam_name,Amat_letter)));
    Amat_svd.print_details("Amat_svd_name");

    LD_Theta_stack Theta_stack_A(Sobs.ncrvs_tot,fspace0.ndof_full);
    Theta_stack_A.load_Theta_data_tns(Amat_svd.VTtns);

    // LD_R_matrix::eval_Rn_condition<orthopolynomial_basis>(Sobs,Theta_stack_A,bases0,tol_use_A);
    LD_G_matrix::eval_inf_criterion<orthopolynomial_basis>(Sobs,Theta_stack_A,bases0,tol_use_A);

    Theta_stack_A.print_satisfactory_details(Amat_letter);
    bool nsat_bf_A[Sobs.ncrvs_tot];
    int nnsat_b_A = 0,
        nnsat_diff_A[Sobs.ncrvs_tot];
    printf("  %s nulldim\n", Amat_letter);
    for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++)
    {
      int nulldim_i = Amat_svd.ncol_use-Amat_svd.rank_vec[icrv];
      printf("%d ", nulldim_i);
      nnsat_diff_A[icrv] = Theta_stack_A.ntheta_spcvec[icrv]-nulldim_i;
      nnsat_b_A+=(int)(nsat_bf_A[icrv]=(nnsat_diff_A[icrv] >= 0));
    }
    printf("\n");
    printf("nsat_bigger_A = %d (of %d, tol = %.1e) \n", nnsat_b_A,Sobs.ncrvs_tot,tol_use_A);
    for (size_t icrv = 0; icrv < Sobs.ncrvs_tot; icrv++) printf("%d ", (int) nsat_bf_A[icrv]);
    printf("\n");
    LD_linalg::print_xT("nsat_diff_A",nnsat_diff_A,Sobs.ncrvs_tot);
    Theta_stack_A.set_Thetas(fspace0.ndof_full);
  }

  free_evaluation_bases<orthopolynomial_basis>(bases0);
  return 0;
}
