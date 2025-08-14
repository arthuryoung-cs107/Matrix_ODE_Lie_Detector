// #include <sys/stat.h>
// #include "named_odes.hh"
// #include "LD_trivial_flows.hh"

#include "LD_noise_aux.hh"
#include "Lie_detector.hh"
#include "LD_integrators.hh"

/*
  This file reads in noisy observed trajectories of the Riccati equation:

    dx u = (a*u)/x + b*(u^2)*(x^2),

    where x is the independent variable and u is the real valued dependent variable.
    Both a and b are real numbers, a*u denotes the scalar multiplication , and u^2 denotes the square of u.

  This is a quadratically nonlinear ordinary differential equation that is singular at x = 0.
  We observe trajectories and build a uniquely defined n=1 order model of the system.

  Assuming that the observations of the dependent variables, u, and their derivatives are corrupted by noise,
    we proceed with "denoising" the input data by searching for the underlying smooth vector field.
*/
// const char dir_name[] = "./denoise_data_directory/Gaussian_IC_perturbation/"; // data directory
const char dir_name[] = "./denoise_data_directory/Gaussian_IC_perturbation/rendering_data/"; // data directory
const char obs_name[] = "Riccati_xrange0_noise1_DoP853gen"; // name of observations file
// const char obs_name[] = "Riccati_xrange0_noise1_DoP853gen.jsol_R1"; // name of observations file
const char dat_suff[] = ".lddat"; // data file suffix (*.lddat)
char data_name[strlen(dir_name) + strlen(obs_name) + strlen(dat_suff) + 1];
const int data_name_len = sprintf(data_name,"%s%s%s", dir_name,obs_name,dat_suff);
// name of data file
// const char data_name[] =
// "./denoise_data_directory/Gaussian_IC_perturbation/"
// "Riccati_xrange0_noise1_DoP853gen" // Riccati_xrange0_noise0_DoP853gen Riccati_xrange0_true_DoP853gen
// ".lddat";

// ode_solspc_meta meta0(1,1); // Riccati equation : n = 1, q = 1

ode_curve_observations observations(data_name);
Lie_detector detector(observations);
  ode_solcurve_chunk observations_twin(detector.meta0,detector.ncrv,detector.npts_per_crv);
orthopolynomial_space function_space(detector.meta0,3); // initializes as unmapped multinomials

const int ndns_max = 10;

// basic experiment, suitable for any first order system
struct global_R1mat_experiment : public ode_solspc_meta
{
  Lie_detector &det;
  orthopolynomial_space &fspace0;

  const int nobs,
            nobs_h,
            nthread_wkspc,
            nsnap;

  LD_R_encoder R1enc;
  LD_svd  R1svd_global;

  // theta (parameter) space variables
  double * const svec_global,
         ** const VTmat_global,
         ** const theta_mat;

  // S (solution space) variables
  ode_solution ** const sols;
  ode_solcurve ** const curves,
               ** const curves_h;

  // Lie detectors on a single trivial flow
  curve_Lie_detector ** const cdets;

  // thread local workspaces
  double *** const xu_snap_wkspc;
  orthopolynomial_basis ** const fbases0;
  LD_spectral_tvfield ** const tvfields0;

  void denoise_data(ode_solcurve_chunk &Sobs_alt_)
  {
    ode_solution ** const sols_alt = Sobs_alt_.sols;
    ode_solcurve ** const curves_alt = Sobs_alt_.curves;

    int nsmooth=1;

    combine_trivial_flows(curves_alt,curves);
      write_curve_observations(Sobs_alt_,nsmooth);

    double t0 = LD_threads::tic();
    do
    {
     det.set_trivial_jets(curves_alt);

     encode_decompose_R1_matrix_global(sols_alt);
       R1svd_global.print_result("R1svd_global");

     smooth_trivial_flows(curves_alt);

     combine_trivial_flows(curves_alt,curves_alt);

     write_curve_observations(Sobs_alt_,++nsmooth);

     if (nsmooth >= 3) break;

    } while (true);
    double work_time = LD_threads::toc(t0);
    printf(
      "(global_R1mat_experiment::global_R1mat_experiment) "
      "encoded and SV-decomposed %d by %d (global) R1 matrix. "
      "Applied R1 pseudoinverse to Hermite jets, "
      "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
      "%d TIMES "
      "in %.3f seconds (%d threads)\n",
      R1svd_global.Muse, R1svd_global.Nuse,
      nobs_h, 2*nobs_h, ndep+1,
      nsmooth,
      work_time, LD_threads::numthreads()
    );
  }

  ~global_R1mat_experiment()
  {
    for (int i = 0; i < nthread_wkspc; i++) {delete tvfields0[i]; delete fbases0[i]; }
    delete [] tvfields0; delete [] fbases0;
    free_T3tensor<double>(xu_snap_wkspc);
    free_Tmatrix<double>(VTmat_global);
    free_Tmatrix<double>(theta_mat);
  }
  global_R1mat_experiment(Lie_detector &det_,orthopolynomial_space &fspace0_,int nsnap_=3) :
    ode_solspc_meta(det_.eor,det_.ndep),
    det(det_), fspace0(fspace0_),
      nobs(det_.nobs), nobs_h(det_.comp_nobs_h()),
        nthread_wkspc(LD_threads::numthreads()),
        nsnap(nsnap_),
    R1enc(det_.meta0,1),
    R1svd_global((det.ndep)*(det.nobs),fspace0.ndof_full),
      svec_global(R1svd_global.svec),
      VTmat_global(Tmatrix<double>(fspace0.ndof_full,fspace0.ndof_full)),
      theta_mat(Tmatrix<double>(det.nobs,fspace0.ndof_full)),
    sols(det.sols),
      curves(det.curves),
      curves_h(det.curves_h),
      cdets(det.cdets),
    xu_snap_wkspc(T3tensor<double>(nthread_wkspc,nsnap,fspace0.nvar)),
    fbases0(new orthopolynomial_basis*[nthread_wkspc]),
    tvfields0(new LD_spectral_tvfield*[nthread_wkspc])
    {
      for (int i = 0; i < nthread_wkspc; i++)
      {
        fbases0[i] = new orthopolynomial_basis(fspace0);
        tvfields0[i] = new LD_spectral_tvfield(*(fbases0[i]),svec_global,VTmat_global);
      }
      double t0 = LD_threads::tic();
      encode_decompose_R1_matrix_global(sols);
        R1svd_global.print_result("R1svd_global");
      double twork1 = LD_threads::toc(t0);
      smooth_trivial_flows(curves);
      double work_time = LD_threads::toc(t0);
      printf(
        "(global_R1mat_experiment::global_R1mat_experiment) "
        "encoded and SV-decomposed %d by %d (global) R1 matrix. "
        "Applied R1 pseudoinverse to Hermite jets, "
        "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
        "in %.3f seconds (%d threads)\n",
        R1svd_global.Muse, R1svd_global.Nuse,
        nobs_h, 2*nobs_h, ndep+1,
        work_time, LD_threads::numthreads()
      );
    }

  inline void combine_trivial_flows(ode_solcurve **crvs_alt_,ode_solcurve **crvs_)
  {
    #pragma omp parallel for
    for (int icrv = 0; icrv < det.ncrv; icrv++)
    {
      cdets[icrv]->recombine_flows(*(crvs_alt_[icrv]),*(crvs_[icrv]));
    }
  }
  void encode_decompose_R1_matrix_global(ode_solution **sols_)
  {
    #pragma omp parallel
    {
      orthopolynomial_basis &fbse_t = *( fbases0[LD_threads::thread_id()] );
      #pragma omp for
      for (int i = 0; i < nobs; i++)
      {
        R1enc.encode_normalize_rows(
            R1enc.submat_i(R1svd_global.Umat,i), // submatrix of R1 associated with observation i
            fbse_t, // thread local prolongation workspace
            *(sols_[i]), // jet space data associated with observation i
            false // normalization flag (off by default)
          );
      }
    }
    R1svd_global.decompose_U();
    R1svd_global.unpack_VTmat(VTmat_global);
  }
  void smooth_trivial_flows(ode_solcurve **crvs_)
  {
    #pragma omp parallel
    {
      const int t_id = LD_threads::thread_id();
      LD_spectral_tvfield &tvf_t = *(tvfields0[t_id]);
      LD_trivial_generator tgen_t(tvf_t);
      DoP853_settings set_t; DoP853 intgr_t(tgen_t,set_t);

      double ** const integr_wkspc_t = xu_snap_wkspc[t_id],
             * const s_state_t_f = integr_wkspc_t[nsnap-1],
             * const s_state_t = intgr_t.get_u_state();
      int icrv_t,
          isol_t;
      #pragma omp for
      for (int iobs = 0; iobs < nobs_h; iobs++)
      {
        det.get_icrvisol_h_given_iobs_h(icrv_t,isol_t,iobs);
        ode_solution &sol_h_i = *(curves_h[icrv_t]->sols[isol_t]),
                     &sol_0_i = *(crvs_[icrv_t]->sols[isol_t]),
                     &sol_1_i = *(crvs_[icrv_t]->sols[isol_t+1]);
        tjet_chart &tjchart_i = *(cdets[icrv_t]->tjcharts[isol_t]);

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of R1 matrix
        tvf_t.comp_spectral_theta_local(theta_mat[iobs],sol_h_i.x,sol_h_i.u);
        // the result is an improved estimate of the derivatives at interior colocation point
        tvf_t.eval_pr1_theta_image(tjchart_i.solh_alt.dxu,theta_mat[iobs]);

        // reevaluate forward flow
        tgen_t.flow_forward = true; // flow varies directly with x
        for (int i = 0; i <= ndep; i++) s_state_t[i] = sol_h_i.pts[i];
        intgr_t.init_curve_integration(0.5*(sol_1_i.x-sol_h_i.x),0); // curve id irrelevant here
        intgr_t.set_and_solve_time(0.0,sol_1_i.x-sol_h_i.x,nsnap,integr_wkspc_t);
        for (int i = 0; i <= ndep; i++) tjchart_i.sol1_alt.pts[i] = s_state_t_f[i];
        tgen_t.dudx_eval(tjchart_i.sol1_alt.pts);

        // reevaluate backward flow
        tgen_t.flow_forward = false; // flow varies oppositely with x
        for (int i = 0; i <= ndep; i++) s_state_t[i] = sol_h_i.pts[i];
        intgr_t.init_curve_integration(0.5*(sol_h_i.x-sol_0_i.x),0); // curve id irrelevant here
        intgr_t.set_and_solve_time(0.0,sol_h_i.x-sol_0_i.x,3,integr_wkspc_t);
        for (int i = 0; i <= ndep; i++) tjchart_i.sol0_alt.pts[i] = s_state_t_f[i];
        tgen_t.dudx_eval(tjchart_i.sol0_alt.pts);
      }
    }
  }
  // void write_sol_h_data(const char dir_name_[], const char obs_name_[], const char dat_suff_[])
  void write_sol_h_data()
  {
    const int len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_base = len_dir_name+len_obs_name+len_dat_suff;

    const char name_theta_mat[] = ".theta_mat";
      char fname_theta_mat[len_base+strlen(name_theta_mat)+1];
      sprintf(fname_theta_mat,"%s%s%s%s",dir_name,obs_name, name_theta_mat ,dat_suff);
      LD_io::write_Tmat<double>(fname_theta_mat,theta_mat,nobs_h,fspace0.ndof_full);

    int ode_meta[2],
        &eor_meta = ode_meta[0] = det.eor,
        &ndep_meta = ode_meta[1] = det.ndep,
        obs_meta[2],
        &ncrv_meta = obs_meta[0] = det.ncrv,
        &nobs_meta = obs_meta[1] = nobs_h,
        npts_per_crv_h[ncrv_meta];
    for (int i = 0; i < det.ncrv; i++) npts_per_crv_h[i] = det.npts_per_crv[i]-1;

    // write raw colocation point resulting from deterministic Hermite jet
    const char name_jsol_h[] = ".jsol_h";
      char fname_jsol_h[len_base+strlen(name_jsol_h)+1];
      sprintf(fname_jsol_h,"%s%s%s%s",dir_name,obs_name, name_jsol_h ,dat_suff);
      FILE * file_jsol_h = LD_io::fopen_SAFE(fname_jsol_h,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_h);
      fwrite(obs_meta,sizeof(int),2,file_jsol_h);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_h);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(curves_h[icrv]->sols[isol]->pts,sizeof(double),det.ndim,file_jsol_h);
      LD_io::fclose_SAFE(file_jsol_h);
      printf("(global_R1mat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h);

    // write smooth colocation point resulting from deterministic Hermite jet and R1 matrix pseudoinverse
    const char name_jsol_h_R1[] = ".jsol_h_R1";
      char fname_jsol_h_R1[len_base+strlen(name_jsol_h_R1)+1];
      sprintf(fname_jsol_h_R1,"%s%s%s%s",dir_name,obs_name, name_jsol_h_R1 ,dat_suff);
      FILE * file_jsol_h_R1 = LD_io::fopen_SAFE(fname_jsol_h_R1,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_h_R1);
      fwrite(obs_meta,sizeof(int),2,file_jsol_h_R1);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_h_R1);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->solh_alt.pts,sizeof(double),det.ndim,file_jsol_h_R1);
      LD_io::fclose_SAFE(file_jsol_h_R1);
      printf("(global_R1mat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h_R1);

    // write new forward observations generated by exponentiating from sol_h
    const char name_jsol_1_R1[] = ".jsol_1_R1";
      char fname_jsol_1_R1[len_base+strlen(name_jsol_1_R1)+1];
      sprintf(fname_jsol_1_R1,"%s%s%s%s",dir_name,obs_name, name_jsol_1_R1 ,dat_suff);
      FILE * file_jsol_1_R1 = LD_io::fopen_SAFE(fname_jsol_1_R1,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_1_R1);
      fwrite(obs_meta,sizeof(int),2,file_jsol_1_R1);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_1_R1);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->sol1_alt.pts,sizeof(double),det.ndim,file_jsol_1_R1);
      LD_io::fclose_SAFE(file_jsol_1_R1);
      printf("(global_R1mat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_1_R1);

    // write new backwards observations generated by exponentiating from sol_h
    const char name_jsol_0_R1[] = ".jsol_0_R1";
      char fname_jsol_0_R1[len_base+strlen(name_jsol_0_R1)+1];
      sprintf(fname_jsol_0_R1,"%s%s%s%s",dir_name,obs_name, name_jsol_0_R1 ,dat_suff);
      FILE * file_jsol_0_R1 = LD_io::fopen_SAFE(fname_jsol_0_R1,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_0_R1);
      fwrite(obs_meta,sizeof(int),2,file_jsol_0_R1);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_0_R1);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->sol0_alt.pts, sizeof(double),det.ndim,file_jsol_0_R1);
      LD_io::fclose_SAFE(file_jsol_0_R1);
      printf("(global_R1mat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_0_R1);
  }
  // void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_, const char dir_name_[], const char obs_name_[], const char dat_suff_[])
  void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_)
  {

    const int len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_base = len_dir_name+len_obs_name+len_dat_suff;
    int ode_meta[2],
        &eor_meta = ode_meta[0] = Sobs_.eor,
        &ndep_meta = ode_meta[1] = Sobs_.ndep,
        obs_meta[2],
        &ncrv_meta = obs_meta[0] = Sobs_.ncrv,
        &nobs_meta = obs_meta[1] = Sobs_.nobs;

    const char name_jsol_R1[] = ".jsol_R1";
      char fname_jsol_R1[len_base+strlen(name_jsol_R1)+1];
      sprintf(fname_jsol_R1,"%s%s%s_%d%s",dir_name,obs_name, name_jsol_R1,nsmooth_, dat_suff);
      FILE * file_jsol_R1 = LD_io::fopen_SAFE(fname_jsol_R1,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_R1);
      fwrite(obs_meta,sizeof(int),2,file_jsol_R1);
      fwrite(Sobs_.npts_per_crv,sizeof(int),ncrv_meta,file_jsol_R1);
      for (int iobs = 0; iobs < nobs_meta; iobs++)
        fwrite((Sobs_.sols[iobs])->pts, sizeof(double),Sobs_.ndim,file_jsol_R1);
      LD_io::fclose_SAFE(file_jsol_R1);
      printf("(global_R1mat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_R1);
  }

};

global_R1mat_experiment R1_experiment(detector,function_space,ndns_max);

int main()
{
  // observations.print_details();
  // detector.print_details();

  R1_experiment.write_sol_h_data();
  R1_experiment.denoise_data(observations_twin);

  return 0;
}
