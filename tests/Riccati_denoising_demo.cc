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
// const char dir_name[] = "./denoise_data_directory/Gaussian_IC_perturbation/rendering_data/"; // data directory
const char dir_name[] = "./denoise_data_directory/Uniform_IC_perturbation/rendering_data/"; // data directory
const char obs_name[] = "Riccati_xrange0_noise0_DoP853gen"; // name of observations file
// const char obs_name[] = "Riccati_xrange0_noise1_DoP853gen"; // name of observations file
// const char obs_name[] = "Riccati_xrange0_noise1_DoP853gen.jsol_R1"; // name of observations file

const char dat_suff[] = ".lddat"; // data file suffix (*.lddat)
char data_name[strlen(dir_name) + strlen(obs_name) + strlen(dat_suff) + 1];
const int data_name_len = sprintf(data_name,"%s%s%s", dir_name,obs_name,dat_suff);
// name of data file
// const char data_name[] =
// "./denoise_data_directory/Gaussian_IC_perturbation/"
// "Riccati_xrange0_noise1_DoP853gen" // Riccati_xrange0_noise0_DoP853gen Riccati_xrange0_true_DoP853gen
// ".lddat";

const char dnp1xu_suffix[] = "_dnp1xu";
char data_dnp1xu_name[strlen(dir_name) + strlen(obs_name) + strlen(dnp1xu_suffix) + strlen(dat_suff) + 1];
const int data_dnp1xu_name_len = sprintf(data_dnp1xu_name,"%s%s%s%s", dir_name,obs_name,dnp1xu_suffix,dat_suff);

// ode_solspc_meta meta0(1,1); // Riccati equation : n = 1, q = 1

const int curve_Lie_detector::combine_flag = 0; // 0, flag for updating smoothened coordinates, (0 lazy, 1 aggressive)
const int Hermite_exp = 1; // 0, flag for Hermite exponentiation technique (0 uses num. quadrature instead)

const bool v_verbose = false;
const bool Rmat_h_exp = false;
const bool stop_blowup = true;
// const bool stop_blowup = false;

const double res_ratio_tol = 1e-10;
const int write_sched_early = 5;

// const int ndns_max = 3; // max permitted denoising steps
// const int ndns_max = 5; // max permitted denoising steps
// const int ndns_max = 10;
// const int ndns_max = 20;
// const int ndns_max = 30;
// const int ndns_max = 40;
// const int ndns_max = 50;
// const int ndns_max = 100;
// const int write_sched = 2;

// const int ndns_max = 100;
// const int ndns_max = 200;
// const int ndns_max = 300;
const int ndns_max = 400;
const int write_sched = 10;

// const int ndns_max = 400;
// const int ndns_max = 500;
// const int write_sched = 20;

// const int ndns_max = 600;
// const int write_sched = 50;

// const int ndns_max = 999;
// const int write_sched = 100;

// ode_curve_observations observations(data_name);
  ode_curve_observations observations(data_name,data_dnp1xu_name);

Lie_detector detector(observations);
  ode_solcurve_chunk observations_twin(detector.meta0,detector.ncrv,detector.npts_per_crv);
orthopolynomial_space function_space(detector.meta0,3); // initializes as unmapped multinomials
  // double axlims[] = {0.0, 2.0, 0.0, 6.0}; // eyeballed the data set
  // const int set_fspace_flag = function_space.set_Legendre_coeffs(axlims,-0.95,0.95);
                                            // set_Legendre_coeffs(axlims,-0.95,0.95);
                                            // set_Chebyshev1_coeffs(axlims,-0.95,0.95);
                                            // set_Chebyshev2_coeffs(axlims,-0.95,0.95);

// basic experiment, suitable for any first order system
struct global_Rmat_experiment : public ode_solspc_meta
{
  Lie_detector &det;
  orthopolynomial_space &fspace0;

  const int nobs,
            nobs_h,
            nthread_wkspc,
            nsnap;

  LD_R_encoder Renc;
  LD_svd  Rsvd_global,
          Rsvd_h_global;

  // theta (parameter) space variables
  double * const svec_global, // singular values of primary R1 matrix
         * const svec_h_global, // singular values of half step R1 matrix
         ** const VTmat_global, // transposed row space of primary R1 matrix
         ** const VTmat_h_global, // transposed row space of half step R1 matrix
         ** const theta_mat; // parameter space image of every observation

  // S (solution space) members, points and sets (curves)
  ode_solcurve ** const curves,
               ** const curves_h;
  ode_solution ** const sols,
               ** const sols_h,
               ** const sols_h_alt;

  // trivial jet space jets, curves, and charts
  jet_chart ** const tjcharts;
  ode_trivial_curvejet ** const tcjets;
  ode_trivial_soljet ** const tsjets;

  // Lie detectors on a single trivial flow
  curve_Lie_detector ** const cdets;

  // thread local workspaces
  double *** const xu_snap_wkspc;
  orthopolynomial_basis ** const fbases0;
  LD_spectral_tvfield ** const tvfields0;

  /*
    [CONSTRUCTOR] on instantiation, this structure receives a Lie detector, which is a structure defined
    in the Lie_Detector.hh header file, that contains raw observational data. The following experimental
    structure proceeds to use the Lie detector data and the given function space to build a n'th order
    trivial flow model of the observed jet space. The evaluation of this constructor is a well defined
    linear program for order n >= 1 ordinary differential equations given any multinomial family of
    order k >= 1.
  */
  global_Rmat_experiment(Lie_detector &det_,orthopolynomial_space &fspace0_,int nsnap_=3) :
    ode_solspc_meta(det_.eor,det_.ndep),
    det(det_), fspace0(fspace0_),
      nobs(det_.nobs), nobs_h(det_.comp_nobs_h()),
        nthread_wkspc(LD_threads::numthreads()),
        nsnap(nsnap_),
    Renc(det_.meta0,det.kor), // Renc(det_.meta0,1),
    Rsvd_global((det.ndep*det.kor)*(det.nobs),fspace0.ndof_full),
    Rsvd_h_global((det.ndep*det.kor)*(det.nobs-det.ncrv),fspace0.ndof_full),
      svec_global(Rsvd_global.svec),
      svec_h_global(Rsvd_h_global.svec),
      VTmat_global(Tmatrix<double>(fspace0.ndof_full,fspace0.ndof_full)),
      VTmat_h_global(Tmatrix<double>(fspace0.ndof_full,fspace0.ndof_full)),
    theta_mat(Tmatrix<double>(det.nobs,fspace0.ndof_full)),
    curves(det.curves),
    curves_h(det.curves_h),
      sols(det.sols),
      sols_h(det.sols_h),
      sols_h_alt(det.sols_h_alt),
    tjcharts(det.tjcharts),
      tcjets(det.tcrvjets),
        tsjets(det.tsoljets),
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

      encode_decompose_R_matrix_global(VTmat_global,Rsvd_global,Renc,sols,nobs);
        Rsvd_global.print_result("Rsvd_global");

      // encode_decompose_R_matrix_global(VTmat_h_global,Rsvd_h_global,Renc,sols_h,nobs_h);
      smooth_enc_decomp_trivial_flows(VTmat_h_global,Rsvd_h_global,sols_h_alt, // output args
                                        nobs_h,Renc,sols_h,svec_global,VTmat_global); // input args
        Rsvd_h_global.print_result("Rsvd_h_global");

      // smooth_and_exp_trivial_flows(curves,curves_h,svec_global,VTmat_global);

      // if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_h_global,VTmat_h_global);
      // else exp_trivial_flows(tjcharts,nobs_h,svec_global,VTmat_global);

      exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_global,VTmat_global);

      double twork1 = LD_threads::toc(t0);

      // double res_1 = det.combine_trivial_flows(curves_alt,curves);

      double work_time = LD_threads::toc(t0);
      printf(
        "(global_Rmat_experiment::global_Rmat_experiment) "
        "encoded and SV-decomposed %d by %d (global) R%d matrix. "
        "Applied R1 pseudoinverse to Hermite jets, "
        "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
        "in %.3f seconds (%d threads)\n",
        Rsvd_global.Muse, Rsvd_global.Nuse, det.kor,
        nobs_h, 2*nobs_h, ndep+1,
        work_time, LD_threads::numthreads()
      );
      write_sol_h_data();
    }
    ~global_Rmat_experiment()
    {
      for (int i = 0; i < nthread_wkspc; i++) {delete tvfields0[i]; delete fbases0[i]; }
      delete [] tvfields0; delete [] fbases0;
      free_T3tensor<double>(xu_snap_wkspc);
      free_Tmatrix<double>(VTmat_global);
      free_Tmatrix<double>(VTmat_h_global);
      free_Tmatrix<double>(theta_mat);
    }

  void denoise_data(ode_solcurve_chunk &Sobs_alt_)
  {
    ode_solution ** const sols_alt = Sobs_alt_.sols;
    ode_solcurve ** const curves_alt = Sobs_alt_.curves;

    /*
     From the top : given solutions from det, we build Hermite jets over the colocation points.
     We smooth the derivatives of the colocation points using Rmat_obs pseudo inverse, then recompute
     an R matrix, this time over the colocation points.
    */
    det.set_trivial_jets(sols_h,curves); // this sets curves_h using curves as input

    encode_decompose_R_matrix_global(VTmat_global,Rsvd_global,Renc,sols,nobs);
      Rsvd_global.print_result("  Rsvd_global (0)");

    smooth_enc_decomp_trivial_flows(VTmat_h_global,Rsvd_h_global,sols_h_alt, // output args
                                      nobs_h,Renc,sols_h,svec_global,VTmat_global); // input args
      Rsvd_h_global.print_result("  Rsvd_h_global (0)");

    if (Hermite_exp)
    {
      if (Rmat_h_exp) exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_h_global,VTmat_h_global);
      else exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_global,VTmat_global);
    }
    else
    {
      if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_h_global,VTmat_h_global);
      else exp_trivial_flows(tjcharts,nobs_h,svec_global,VTmat_global);
    }

    double res_1 = det.combine_trivial_flows(curves_alt,curves);
    int nsmooth=1;
      write_curve_observations(Sobs_alt_,nsmooth);

    printf("   (global_Rmat_experiment::denoise_data)"
           " i=%d, res_i = %.8f "
           "\n", nsmooth, res_1);

    double  res_old = res_1,
            t0 = LD_threads::tic();
    do
    {
     det.set_trivial_jets(sols_h,curves_alt);
     encode_decompose_R_matrix_global(VTmat_global,Rsvd_global,Renc,sols_alt,nobs);
       if (v_verbose) Rsvd_global.print_result("    Rsvd_global");

     if (Rmat_h_exp)
     {
       smooth_enc_decomp_trivial_flows(VTmat_h_global,Rsvd_h_global,sols_h_alt, // output args
                                         nobs_h,Renc,sols_h,svec_global,VTmat_global); // input args
         if (v_verbose) Rsvd_h_global.print_result("    Rsvd_h_global");
     }

     if (Hermite_exp)
     {
       if (Rmat_h_exp) exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_h_global,VTmat_h_global);
       else exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_global,VTmat_global);
     }
     else
     {
       if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_h_global,VTmat_h_global);
       else exp_trivial_flows(tjcharts,nobs_h,svec_global,VTmat_global);
     }

     double res_i = det.combine_trivial_flows(curves_alt,curves_alt);

     if ( ((++nsmooth)<=write_sched_early) || ((nsmooth%write_sched) == 0) )
      write_curve_observations(Sobs_alt_,nsmooth,v_verbose);

     printf("   (global_Rmat_experiment::denoise_data)"
            " i=%d, res_i = %.8f ( res_i/res_0 = %.1e )"
            "\n", nsmooth, res_i, res_i/res_1);

     if ( nsmooth >= ndns_max ) break;
     if ( (stop_blowup)&&(res_i > res_old) ) break;
     if ( (res_i/res_1) < res_ratio_tol ) break;

     res_old = res_i;
    } while (true);
    double work_time = LD_threads::toc(t0);
    printf(
      "(global_Rmat_experiment::global_Rmat_experiment) "
      "encoded and SV-decomposed %d by %d (global, %s) R1 matrices. "
      "Applied R1 pseudoinverse to Hermite jets, "
      "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
      "\n%d TIMES\n "
      "in %.3f seconds (%d threads)\n",
      Rsvd_global.Muse, Rsvd_global.Nuse, (Rmat_h_exp)?("TWICE"):("once"),
      nobs_h, 2*nobs_h, ndep+1,
      nsmooth,
      work_time, LD_threads::numthreads()
    );
    write_curve_observations(Sobs_alt_,nsmooth,true);
  }

  void encode_decompose_R_matrix_global(double **VTmg_,LD_svd &Rsvdg_,LD_R_encoder &Renc_,ode_solution **sols_,int nobs_)
  {
    // #pragma omp parallel num_threads(1)
    #pragma omp parallel
    {
      orthopolynomial_basis &fbse_t = *( fbases0[LD_threads::thread_id()] );
      #pragma omp for
      for (int i = 0; i < nobs_; i++)
      {
        Renc_.encode_normalize_rows(
            Renc_.submat_i(Rsvdg_.Umat,i), // submatrix of R associated with observation i
            fbse_t, // thread local prolongation workspace
            *(sols_[i]), // jet space data associated with observation i
            false // normalization flag (off by default)
          );
      }
    }
    Rsvdg_.decompose_U();
    Rsvdg_.unpack_VTmat(VTmg_);
  }
  void smooth_trivial_flows(ode_solution **sols_out_,int nobs_,int kor_,ode_solution **sols_in_,double *sv_,double **VTm_)
  {
    #pragma omp parallel
    {
      LD_spectral_tvfield &tvf_t = *(tvfields0[LD_threads::thread_id()]);
      tvf_t.set_SVD_space(sv_,VTm_);
      #pragma omp for
      for (int iobs = 0; iobs < nobs_; iobs++)
      {
        if ( sols_out_[iobs] != sols_in_[iobs] ) sols_out_[iobs]->copy_xu(*(sols_in_[iobs]));

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of R1 matrix
        tvf_t.comp_spectral_theta_local(sols_out_[iobs]->x,sols_out_[iobs]->u);
        // the result is an improved estimate of the derivatives at interior colocation point
        tvf_t.eval_prn_theta_image(*(sols_out_[iobs]),kor_);
      }
    }
  }
  void smooth_enc_decomp_trivial_flows(double **VTm_out_,LD_svd &Rsvdg_,ode_solution **sout_,
    int nobs_,LD_R_encoder &Renc_,ode_solution **sin_,double *sv_,double **VTm_)
  {
    #pragma omp parallel
    {
      const int t_id = LD_threads::thread_id();
      orthopolynomial_basis &fbse_t = *(fbases0[t_id]);
      LD_spectral_tvfield &tvf_t = *(tvfields0[t_id]);
      tvf_t.set_SVD_space(sv_,VTm_);
      #pragma omp for
      for (int iobs = 0; iobs < nobs_; iobs++)
      {
        if ( sout_[iobs] != sin_[iobs] ) sout_[iobs]->copy_xu(*(sin_[iobs])); // make sure synced

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of R1 matrix
        tvf_t.comp_spectral_theta_local(sout_[iobs]->x,sout_[iobs]->u);
        // the result is an improved estimate of the derivatives at interior colocation point
        tvf_t.eval_prn_theta_image(*(sout_[iobs]),det.kor);

        Renc_.encode_normalize_rows(
            Renc_.submat_i(Rsvdg_.Umat,iobs), // submatrix of R associated with observation i
            fbse_t, // thread local prolongation workspace
            *(sout_[iobs]), // jet space data associated with observation i
            false // normalization flag (off by default)
          );
      }
    }
    Rsvdg_.decompose_U();
    Rsvdg_.unpack_VTmat(VTm_out_);
  }
  void exp_trivial_Rmat_Hermites(jet_chart **tjc_, ode_trivial_soljet **tsj_, int nobs_, double *sv_, double **VTm_)
  {
    #pragma omp parallel
    {
      const int t_id = LD_threads::thread_id();
      LD_spectral_tvfield &tvf_t = *(tvfields0[t_id]);
        tvf_t.set_SVD_space(sv_,VTm_);

      #pragma omp for
      for (int iobs = 0; iobs < nobs_; iobs++)
      {
       // reevaluate flows by R pseudoinverse corrected Hermite jets.
       jet_chart &tjchart_i = *(tjc_[iobs]);
       ode_trivial_soljet &tsjet_i = *(tsj_[iobs]);
       tvf_t.set_SVD_space(sv_,VTm_);

       tjchart_i.solh_alt.copy_xu(tjchart_i.solh); // make sure synced

       // use 0'th order Hermite jet info to compute local theta via pseudoinverse of R1 matrix
       tvf_t.comp_spectral_theta_local(tjchart_i.solh.x,tjchart_i.solh.u);
        tvf_t.eval_prn_theta_image(tjchart_i.solh_alt,det.kor);
       // the result is an improved estimate of the derivatives at interior colocation point

       // now modify the Hermite jet coefficients to reflect R pseudoinverse estimates
       tsjet_i.set_jet_given_solh(tjchart_i.solh_alt);

       // use corrected Hermite jet to predict new values for u0, u1
       tjchart_i.sol0_alt.x = tjchart_i.sol0.x; tjchart_i.sol1_alt.x = tjchart_i.sol1.x;
       tsjet_i.compute_u_01(tjchart_i.sol0_alt.u,tjchart_i.sol1_alt.u);

       // correct derivative estimates again
         tvf_t.comp_spectral_theta_local(tjchart_i.sol0_alt.x,tjchart_i.sol0_alt.u);
          tvf_t.eval_prn_theta_image(tjchart_i.sol0_alt,det.kor);
         tvf_t.comp_spectral_theta_local(tjchart_i.sol1_alt.x,tjchart_i.sol1_alt.u);
          tvf_t.eval_prn_theta_image(tjchart_i.sol1_alt,det.kor);

      }
    }
  }
  void exp_trivial_flows(jet_chart **tjc_, int nobs_, double *sv_, double **VTm_)
  {
    #pragma omp parallel
    {
      const int t_id = LD_threads::thread_id();
      LD_spectral_tvfield &tvf_t = *(tvfields0[t_id]);
        tvf_t.set_SVD_space(sv_,VTm_);
      LD_trivial_generator tgen_t(tvf_t);
      DoP853_settings set_t; DoP853 intgr_t(tgen_t,set_t);
      double ** const integr_wkspc_t = xu_snap_wkspc[t_id],
             * const s_state_t_f = integr_wkspc_t[nsnap-1],
             * const s_state_t = intgr_t.get_u_state();
      #pragma omp for
      for (int iobs = 0; iobs < nobs_; iobs++)
      {
        jet_chart &tjchart_i = *(tjc_[iobs]);

        // reevaluate forward flow by numerical quadrature
          tgen_t.flow_forward = true; // flow varies directly with x
          for (int i = 0; i <= ndep; i++) s_state_t[i] = tjchart_i.solh.pts[i];
          // specify del_x for dense output. Note that curve id irrelevant here
          intgr_t.init_curve_integration((tjchart_i.sol1.x-tjchart_i.solh.x)/((double) (nsnap-1)),0);
          intgr_t.set_and_solve_time(0.0,tjchart_i.sol1.x-tjchart_i.solh.x,nsnap,integr_wkspc_t);
          for (int i = 0; i <= ndep; i++) tjchart_i.sol1_alt.pts[i] = s_state_t_f[i];
          tvf_t.comp_spectral_theta_local(tjchart_i.sol1_alt.x,tjchart_i.sol1_alt.u);
          tvf_t.eval_prn_theta_image(tjchart_i.sol1_alt,det.kor);

        ////////////////

        // reevaluate backward flow by numerical quadrature
          tgen_t.flow_forward = false; // flow varies oppositely with x
          for (int i = 0; i <= ndep; i++) s_state_t[i] = tjchart_i.solh.pts[i];
          // specify del_x for dense output. Note that curve id irrelevant here
          intgr_t.init_curve_integration((tjchart_i.solh.x-tjchart_i.sol0.x)/((double) (nsnap-1)),0);
          intgr_t.set_and_solve_time(0.0,tjchart_i.solh.x-tjchart_i.sol0.x,nsnap,integr_wkspc_t);
          for (int i = 0; i <= ndep; i++) tjchart_i.sol0_alt.pts[i] = s_state_t_f[i];
          tvf_t.comp_spectral_theta_local(tjchart_i.sol0_alt.x,tjchart_i.sol0_alt.u);
          tvf_t.eval_prn_theta_image(tjchart_i.sol0_alt,det.kor);

      }
    }
  }
  void smooth_and_exp_trivial_flows(ode_solcurve **crvs_, ode_solcurve **crvs_h_, double *sv_, double **VTm_)
  {
    #pragma omp parallel
    {
      const int t_id = LD_threads::thread_id();
      LD_spectral_tvfield &tvf_t = *(tvfields0[t_id]);
      tvf_t.set_SVD_space(sv_,VTm_);
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
        ode_solution &sol_h_i = *(crvs_h_[icrv_t]->sols[isol_t]),
                     &sol_0_i = *(crvs_[icrv_t]->sols[isol_t]),
                     &sol_1_i = *(crvs_[icrv_t]->sols[isol_t+1]);
        jet_chart &tjchart_i = *(cdets[icrv_t]->tjcharts[isol_t]);

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of R1 matrix
        tvf_t.comp_spectral_theta_local(theta_mat[iobs],sol_h_i.x,sol_h_i.u);
        // the result is an improved estimate of the derivatives at interior colocation point
        tvf_t.eval_pr1_theta_image(tjchart_i.solh_alt.dxu,theta_mat[iobs]);

        // reevaluate forward flow
        tgen_t.flow_forward = true; // flow varies directly with x
        for (int i = 0; i <= ndep; i++) s_state_t[i] = sol_h_i.pts[i];
        // specify del_x for dense output. Note that curve id irrelevant here
        intgr_t.init_curve_integration((sol_1_i.x-sol_h_i.x)/((double) (nsnap-1)),0); // curve id irrelevant here
        intgr_t.set_and_solve_time(0.0,sol_1_i.x-sol_h_i.x,nsnap,integr_wkspc_t);
        for (int i = 0; i <= ndep; i++) tjchart_i.sol1_alt.pts[i] = s_state_t_f[i];
        // tgen_t.dudx_eval(tjchart_i.sol1_alt.pts);
        tvf_t.comp_spectral_theta_local(tjchart_i.sol1_alt.x,tjchart_i.sol1_alt.u);
        tvf_t.eval_prn_theta_image(tjchart_i.sol1_alt,det.kor);

        // reevaluate backward flow
        tgen_t.flow_forward = false; // flow varies oppositely with x
        for (int i = 0; i <= ndep; i++) s_state_t[i] = sol_h_i.pts[i];
        // specify del_x for dense output. Note that curve id irrelevant here
        intgr_t.init_curve_integration((sol_h_i.x-sol_0_i.x)/((double) (nsnap-1)),0); // curve id irrelevant here
        intgr_t.set_and_solve_time(0.0,sol_h_i.x-sol_0_i.x,nsnap,integr_wkspc_t);
        for (int i = 0; i <= ndep; i++) tjchart_i.sol0_alt.pts[i] = s_state_t_f[i];
        // tgen_t.dudx_eval(tjchart_i.sol0_alt.pts);
        tvf_t.comp_spectral_theta_local(tjchart_i.sol0_alt.x,tjchart_i.sol0_alt.u);
        tvf_t.eval_prn_theta_image(tjchart_i.sol0_alt,det.kor);
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

    const char name_Rsvd[] = ".Rsvd_g";
      char fname_Rsvd[len_base+strlen(name_Rsvd)+1];
      sprintf(fname_Rsvd,"%s%s%s%s",dir_name,obs_name, name_Rsvd ,dat_suff);
      Rsvd_global.write_LD_svd(fname_Rsvd);

    const char name_Rsvd_h[] = ".Rsvd_h_g";
      char fname_Rsvd_h[len_base+strlen(name_Rsvd_h)+1];
      sprintf(fname_Rsvd_h,"%s%s%s%s",dir_name,obs_name, name_Rsvd_h ,dat_suff);
      Rsvd_h_global.write_LD_svd(fname_Rsvd_h);

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
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h);

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
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h_R1);

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
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_1_R1);

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
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_0_R1);

  }
  // void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_, const char dir_name_[], const char obs_name_[], const char dat_suff_[])
  void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_, bool verbose_=true)
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
      char fname_jsol_R1[len_base+strlen(name_jsol_R1)+5];
      sprintf(fname_jsol_R1,"%s%s%s_%d%s",dir_name,obs_name, name_jsol_R1,nsmooth_, dat_suff);
      FILE * file_jsol_R1 = LD_io::fopen_SAFE(fname_jsol_R1,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_R1);
      fwrite(obs_meta,sizeof(int),2,file_jsol_R1);
      fwrite(Sobs_.npts_per_crv,sizeof(int),ncrv_meta,file_jsol_R1);
      for (int iobs = 0; iobs < nobs_meta; iobs++)
        fwrite((Sobs_.sols[iobs])->pts, sizeof(double),Sobs_.ndim,file_jsol_R1);
      LD_io::fclose_SAFE(file_jsol_R1);
      if (verbose_) printf("(global_Rmat_experiment::write_curve_observations) wrote %s\n",fname_jsol_R1);
  }

};

global_Rmat_experiment R1_experiment(detector,function_space);

int main()
{
  // observations.print_details();
  // detector.print_details();

  // R1_experiment.write_sol_h_data();
  R1_experiment.denoise_data(observations_twin);

  return 0;
}
