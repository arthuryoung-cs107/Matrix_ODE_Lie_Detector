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
// const char dir_name[] = "./denoise_data_directory/Uniform_IC_perturbation/rendering_data/"; // data directory
const char dir_name[] = "./denoise_data_directory/Uniform_IC_perturbation/"; // data directory
// const char obs_name[] = "Riccati_xrange0_true_DoP853gen"; // name of observations file
const char obs_name[] = "Riccati_xrange0_noise0_DoP853gen"; // name of observations file
// const char obs_name[] = "Riccati_xrange0_noise1_DoP853gen"; // name of observations file
// const char obs_name[] = "Riccati_xrange0_noise2_DoP853gen"; // name of observations file

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

const int Hermite_exp = 1; // 1, flag for Hermite exponentiation technique (0 uses num. quadrature instead)
const bool exp_staggered_Hermites = true; // use staggered Hermite jets for denoising exponentiation
const bool trunc_Hermite_exp = true; // flag for truncating Rmatrix corrected Hermite jet
const bool Rmat_telescoping_decomposition = true; // split up global SVD into curve based parallel SVD
const bool stepladder_update = true;

const int curve_Lie_detector::combine_flag = 1; // 0, flag for updating smoothened coordinates, (0 lazy, 1 aggressive)
const bool curve_Lie_detector::truncate_Hermite_exp = trunc_Hermite_exp;
// int curve_Lie_detector::kor_upd = 0;

const bool v_verbose = false;
const bool Rmat_h_exp = false;
const bool Rmat_Hermite_jets = false;
const bool stop_blowup = false;
const int nullity_stop = 2; // nullity dimension required to terminate (if 0, does not stop)

const double res_ratio_tol = 1e-12;
const double stepladder_ratio_tol = 1e-6;
const int write_sched_early = 5;

// const int ndns_max = 3; // max permitted denoising steps
// const int ndns_max = 5; // max permitted denoising steps
// const int ndns_max = 10;
// const int ndns_max = 50;
// const int write_sched = 1;

// const int ndns_max = 10;
// const int ndns_max = 20;
// const int ndns_max = 30;
// const int ndns_max = 40;
// const int ndns_max = 50;
// const int ndns_max = 100;
// const int write_sched = 2;

// const int ndns_max = 50;
const int ndns_max = 1000;
const int write_sched = 5;

// const int ndns_max = 2000;
// const int write_sched = 100;

// const int ndns_max = 100;
// const int ndns_max = 200;
// const int ndns_max = 300;
// const int ndns_max = 400;
// const int ndns_max = 500;
// const int write_sched = 10;

// const int ndns_max = 400;
// const int ndns_max = 500;
// const int write_sched = 20;

// const int ndns_max = 600;
// const int write_sched = 50;

// const int ndns_max = 1000;
// const int write_sched = 100;

// ode_curve_observations observations(data_name);
  ode_curve_observations observations(data_name,data_dnp1xu_name);

const int kor_obs = observations.kor_obs();

Lie_detector detector(observations);
  ode_solcurve_chunk observations_twin(detector.meta0,detector.ncrv,detector.npts_per_crv); // equal size data set for workspace
orthopolynomial_space function_space(detector.meta0,3); // initializes as unmapped multinomials
  // double axlims[] = {0.0, 2.0, 0.0, 6.0}; // eyeballed the data set
  // const int set_fspace_flag = function_space.set_Legendre_coeffs(axlims,-0.95,0.95);
                                            // set_Legendre_coeffs(axlims,-0.95,0.95);
                                            // set_Chebyshev1_coeffs(axlims,-0.95,0.95);
                                            // set_Chebyshev2_coeffs(axlims,-0.95,0.95);

// basic experiment, suitable for any first order system
struct global_Rmat_experiment : public global_multinomial_experiment
{
  const int nthread_wkspc,
            nsnap;

  LD_R_encoder Rkenc;
  LD_svd  &Rsvd_global,
          Rsvd_h_global;

  // theta (parameter) space variables
  double * const svec_Rk_global, // singular values of primary Rn matrix
         * const svec_Rk_h_global, // singular values of half step Rn matrix
         ** const VTmat_Rk_global, // transposed row space of primary Rn matrix
         ** const VTmat_Rk_h_global; // transposed row space of half step Rn matrix

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
    global_multinomial_experiment(det_,fspace0_,det_.kor),
      nthread_wkspc(LD_threads::numthreads()),
      nsnap(nsnap_),
    Rkenc(det_.meta0,det_.kor),
    Rsvd_global(Asvd_global),
    Rsvd_h_global((det.ndep*det.kor)*(det.nobs-det.ncrv),fspace0.ndof_full),
      svec_Rk_global(Rsvd_global.svec),
      svec_Rk_h_global(Rsvd_h_global.svec),
      VTmat_Rk_global(VTmat_global),
      VTmat_Rk_h_global(Tmatrix<double>(fspace0.ndof_full,fspace0.ndof_full)),
    xu_snap_wkspc(T3tensor<double>(nthread_wkspc,nsnap,fspace0.nvar)),
    fbases0(new orthopolynomial_basis*[nthread_wkspc]),
    tvfields0(new LD_spectral_tvfield*[nthread_wkspc])
    {
      for (int i = 0; i < nthread_wkspc; i++)
      {
        fbases0[i] = new orthopolynomial_basis(fspace0);
        tvfields0[i] = new LD_spectral_tvfield(*(fbases0[i]),svec_Rk_global,VTmat_Rk_global);
      }
      double t0 = LD_threads::tic();

      int rank0 = encode_decompose_R_matrix_global(VTmat_Rk_global,Rsvd_global,Rkenc,sols,nobs);
        Rsvd_global.print_result("Rsvd_global");

      if (Rmat_Hermite_jets) // use R matrix pseudo inverse for Hermite knots
        set_trivial_Rmat_Hermite_jets(sols_h,curves,svec_Rk_global,VTmat_Rk_global);

      smooth_enc_decomp_trivial_flows(VTmat_Rk_h_global,Rsvd_h_global,sols_h_alt, // output args
                                        nobs_h,Rkenc,sols_h,svec_Rk_global,VTmat_Rk_global); // input args
        Rsvd_h_global.print_result("Rsvd_h_global");

      double twork1 = LD_threads::toc(t0);

      if (Hermite_exp)
      {
        if (Rmat_h_exp) exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
        else exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_global,VTmat_Rk_global);
      }
      else
      {
        if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
        else exp_trivial_flows(tjcharts,nobs_h,svec_Rk_global,VTmat_Rk_global);
      }

      double work_time = LD_threads::toc(t0);
      printf(
        "(global_Rmat_experiment::global_Rmat_experiment) "
        "encoded and SV-decomposed %d by %d (global) R%d matrix. "
        "Applied Rk pseudoinverse to Hermite jets, "
        "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
        "in %.3f seconds (%d threads)\n",
        Rsvd_global.Muse, Rsvd_global.Nuse, det.kor,
        nobs_h, 2*nobs_h, ndep+1,
        work_time, LD_threads::numthreads()
      );
      // write_sol_h_data();
      write_initial_diagnostics();
    }
    ~global_Rmat_experiment()
    {
      for (int i = 0; i < nthread_wkspc; i++) {delete tvfields0[i]; delete fbases0[i]; }
      delete [] tvfields0; delete [] fbases0;
      free_T3tensor<double>(xu_snap_wkspc);
      free_Tmatrix<double>(VTmat_Rk_h_global);
    }

  void denoise_data(ode_solcurve_chunk &Sobs_alt_)
  {
    int kor_update = (stepladder_update)?(0):(kor_obs),
        rnk_vec[ndns_max],
        iwrite_vec[ndns_max],
        nwrite = 0;
    double res_vec[ndns_max];

    ode_solution ** const sols_alt = Sobs_alt_.sols;
    ode_solcurve ** const curves_alt = Sobs_alt_.curves;

    /*
     From the top : given solutions from det, we build Hermite jets over the colocation points.
     We smooth the derivatives of the colocation points using Rmat_obs pseudo inverse, then recompute
     an R matrix, this time over the colocation points.
    */

    int rank0 = rnk_vec[0] = encode_decompose_R_matrix_global(VTmat_Rk_global,Rsvd_global,Rkenc,sols,nobs);
      Rsvd_global.print_result("  Rsvd_global (0)");

    // double &res_1 = res_vec[0];
    double res_1;

    if (exp_staggered_Hermites)
      res_1 = exp_staggered_trivial_Rmat_Hermites(curves_alt,
        svec_Rk_global,VTmat_Rk_global,curves,kor_update);
    else
    {
      if (Rmat_Hermite_jets) // use R matrix pseudo inverse for Hermite knots
        set_trivial_Rmat_Hermite_jets(sols_h,curves,svec_Rk_global,VTmat_Rk_global);
      else det.set_trivial_jets(sols_h,curves); // this sets curves_h using curves as input

      smooth_enc_decomp_trivial_flows(VTmat_Rk_h_global,Rsvd_h_global,sols_h_alt, // output args
                                        nobs_h,Rkenc,sols_h,svec_Rk_global,VTmat_Rk_global); // input args
        Rsvd_h_global.print_result("  Rsvd_h_global (0)");

      if (Hermite_exp)
      {
       if (Rmat_h_exp) exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
       else exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_global,VTmat_Rk_global);
      }
      else
      {
       if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
       else exp_trivial_flows(tjcharts,nobs_h,svec_Rk_global,VTmat_Rk_global);
      }
      res_1 = det.combine_trivial_flows(curves_alt,curves_alt);
    }
      res_vec[0] = res_1;
      write_sol_h_data("_s");

    int nsmooth=1;
      write_curve_observations(Sobs_alt_,iwrite_vec[nwrite++] = nsmooth,true);


    printf("   (global_Rmat_experiment::denoise_data)"
           " i=1, r_1 = %.8f "
           "\n", res_1);

    // if (exp_staggered_Hermites)
    //   rnk_vec[nsmooth] = encode_decompose_R_matrix_global(VTmat_Rk_global,Rsvd_global,Rkenc,sols_alt,nobs);
    bool update_res_1 = false;
    double  res_old = res_1,
            t0 = LD_threads::tic();
    do
    {
      double ti = LD_threads::tic();
      int &rank_i =
        rnk_vec[nsmooth] = encode_decompose_R_matrix_global(VTmat_Rk_global,Rsvd_global,Rkenc,sols_alt,nobs);
       if (v_verbose) Rsvd_global.print_result("    Rsvd_global");
      double &res_i = res_vec[nsmooth];

      if (update_res_1)
      {
        res_1 = res_i;
        update_res_1 = false;
      }

      double ti_svd = LD_threads::tic();

      // employ staggered Hermite exponentiation technique
      if (exp_staggered_Hermites) res_i =
        exp_staggered_trivial_Rmat_Hermites(curves_alt,svec_Rk_global,VTmat_Rk_global,curves_alt,kor_update);
      else
      {
        if (Rmat_Hermite_jets) // use R matrix pseudo inverse for Hermite knots
          set_trivial_Rmat_Hermite_jets(sols_h,curves_alt,svec_Rk_global,VTmat_Rk_global);
        else det.set_trivial_jets(sols_h,curves_alt); // this sets curves_h using curves as input

        if (Rmat_h_exp)
        {
         smooth_enc_decomp_trivial_flows(VTmat_Rk_h_global,Rsvd_h_global,sols_h_alt, // output args
                                           nobs_h,Rkenc,sols_h,svec_Rk_global,VTmat_Rk_global); // input args
           if (v_verbose) Rsvd_h_global.print_result("    Rsvd_h_global");
        }
        if (Hermite_exp)
        {
         if (Rmat_h_exp) exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
         else exp_trivial_Rmat_Hermites(tjcharts,tsjets,nobs_h,svec_Rk_global,VTmat_Rk_global);
        }
        else
        {
         if (Rmat_h_exp) exp_trivial_flows(tjcharts,nobs_h,svec_Rk_h_global,VTmat_Rk_h_global);
         else exp_trivial_flows(tjcharts,nobs_h,svec_Rk_global,VTmat_Rk_global);
        }
        res_i = det.combine_trivial_flows(curves_alt,curves_alt);
      }

      double ti_exp = LD_threads::tic();

      if ( ((++nsmooth)<=write_sched_early) || ((nsmooth%write_sched) == 0) )
        write_curve_observations(Sobs_alt_,iwrite_vec[nwrite++] = nsmooth,v_verbose);

      printf("#" // "   (global_Rmat_experiment::denoise_data)"
            " i=%d, ri = %.3e ( ri/r0 = %.1e )."
            " Rk : rank = %d, kdim = %d, sN = %.1e, sN/s1 = %.1e."
            " dti=%.1e, dti_svd=%.1e (%.2f)" // , dti_exp=%.1e"
            "\n",
              nsmooth, res_i, res_i/res_vec[0],
              rank_i, Rsvd_global.Nuse-rank_i,
              Rsvd_global.sigma_min(), Rsvd_global.cond(),
              ti_exp-ti, ti_svd-ti, (ti_svd-ti)/(ti_exp-ti)
              // , ti_exp-ti_svd
            );

      if ( nsmooth >= ndns_max ) break;
      if ( (stop_blowup)&&(res_i > res_old) ) break;
      if ( (res_i/res_vec[0]) < res_ratio_tol ) break;
      if ( nullity_stop&&(rank_i<=( Rsvd_global.Nuse-nullity_stop )) ) break;

      if ( stepladder_update&&( (res_i/res_1) < stepladder_ratio_tol ) )
      {
        // if (kor_update==det.kor) break;
        // else update_kor_constraints(curves_alt,svec_Rk_global,VTmat_Rk_global,++kor_update);
        update_kor_constraints(curves_alt,svec_Rk_global,VTmat_Rk_global,det.kor);
        // kor_update++;
        // update_res_1 = true;
      }
      if (kor_update>det.kor) break;
      res_old = res_i;
    } while (true);
    double work_time = LD_threads::toc(t0);
    printf(
      "(global_Rmat_experiment::global_Rmat_experiment) "
      "encoded and SV-decomposed %d by %d (global, %s) Rk matrices. "
      "Applied Rk pseudoinverse to Hermite jets, "
      "then exponentiated 2*(%d)=%d flowouts over %d dimensions "
      "\n%d TIMES\n "
      "in %.3f seconds (%d threads)\n",
      Rsvd_global.Muse, Rsvd_global.Nuse, (Rmat_h_exp)?("TWICE"):("once"),
      nobs_h, 2*nobs_h, ndep+1,
      nsmooth,
      work_time, LD_threads::numthreads()
    );

    write_curve_observations(Sobs_alt_,nsmooth,true);
    encode_decompose_R_matrix_global(VTmat_Rk_global,Rsvd_global,Rkenc,sols_alt,nobs);
    write_sol_h_data("_f");

    const int hlen_sum = 4;
    int sum_meta[hlen_sum+1];
        sum_meta[0] = hlen_sum,
        sum_meta[1] = nwrite,
        sum_meta[2] = nsmooth,
        sum_meta[3] = (1*nwrite)+(1*nsmooth),
        sum_meta[4] = 1*nsmooth;

    const int len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_base = len_dir_name+len_obs_name+len_dat_suff;
    const char name_summary[] = ".denoise_summary";
      char fname_summary[len_base+strlen(name_summary)+1];
      sprintf(fname_summary,"%s%s%s%s",dir_name,obs_name, name_summary ,dat_suff);
      FILE * file_summary = LD_io::fopen_SAFE(fname_summary,"wb");
    fwrite(sum_meta,sizeof(int),hlen_sum+1,file_summary);

    fwrite(iwrite_vec,sizeof(int),nwrite,file_summary);
    fwrite(rnk_vec,sizeof(int),nsmooth,file_summary);

    fwrite(res_vec,sizeof(double),nsmooth,file_summary);

    LD_io::fclose_SAFE(file_summary);
    printf("(global_Rmat_experiment::denoise_data) wrote %s\n",fname_summary);

  }
  int encode_decompose_R_matrix_global(double **VTmg_,LD_svd &Rsvdg_,LD_R_encoder &Rkenc_,ode_solution **sols_,int nobs_)
  {
    #pragma omp parallel
    {
      orthopolynomial_basis &fbse_t = *( fbases0[LD_threads::thread_id()] );
      #pragma omp for
      for (int i = 0; i < nobs_; i++)
      {
        Rkenc_.encode_normalize_rows(
            Rkenc_.submat_i(Rsvdg_.Umat,i), // submatrix of R associated with observation i
            fbse_t, // thread local prolongation workspace
            *(sols_[i]), // jet space data associated with observation i
            false // normalization flag (off by default)
          );
      }
    }
    if (Rmat_telescoping_decomposition)
      return telescope_decompose_global_matrix(VTmg_,Rsvdg_,Rkenc_,nobs_);
    else
    {
      Rsvdg_.decompose_U();
      Rsvdg_.unpack_VTmat(VTmg_);
      return Rsvdg_.rank();
    }
  }
  double update_kor_constraints(ode_solcurve **crvs_o_,double *sv_, double **VTm_,int kor_upd_)
  {
    double res_out = 0.0;
    #pragma omp parallel reduction(+:res_out)
    {
      LD_spectral_tvfield &tvf_t = *(tvfields0[LD_threads::thread_id()]);
        tvf_t.set_SVD_space(sv_,VTm_);
      #pragma omp for
      for (int icrv = 0; icrv < ncrv; icrv++)
      {
        res_out += det.cdets[icrv]->combine_two_sols_dxu(
          *(crvs_o_[icrv]->sols[0]),
          *(crvs_o_[icrv]->sols[0]),
          (det.cdets[icrv]->tjchart_start())->sol1_alt,kor_upd_);
        for (int iobs = 1; iobs < det.cdets[icrv]->nobs_h; iobs++)
        {
          // apply vfield to resultant solution for improved derivatives
          // tvf_t.set_sol_dkxu(*(crvs_o_[icrv]->sols[iobs]), kor_upd_);
          res_out += det.cdets[icrv]->combine_two_sols_dxu(
            *(crvs_o_[icrv]->sols[iobs]),
            *(crvs_o_[icrv]->sols[iobs]),
            (det.cdets[icrv]->tjcharts[iobs])->sol1_alt,kor_upd_);
        }
        res_out += det.cdets[icrv]->combine_two_sols_dxu(
          *(crvs_o_[icrv]->sols[det.cdets[icrv]->nobs_h]),
          *(crvs_o_[icrv]->sols[det.cdets[icrv]->nobs_h]),
          (det.cdets[icrv]->tjchart_final())->sol1_alt,kor_upd_);
      }
    }
    return res_out;
  }
  double exp_staggered_trivial_Rmat_Hermites(ode_solcurve **crvs_o_,double *sv_, double **VTm_,ode_solcurve **crvs_i_,int kor_upd_=-1)
  {
    const int kor_upd = (kor_upd_>=0)?(kor_upd_):(det.kor);
    double res_out = 0.0;
    #pragma omp parallel reduction(+:res_out)
    {

      LD_spectral_tvfield &tvf_t = *(tvfields0[LD_threads::thread_id()]);
        tvf_t.set_SVD_space(sv_,VTm_);

      LD_lu lu_t(tjmeta.jor+1,tjmeta.jor+1);

      #pragma omp for
      for (int icrv = 0; icrv < ncrv; icrv++)
      {
        res_out += cdets[icrv]->compute_staggered_Hermite_flow_curve(crvs_o_[icrv]->sols, lu_t,tvf_t , crvs_i_[icrv]->sols, kor_upd);
      }
    }
    return res_out;
  }
  void set_trivial_Rmat_Hermite_jets(ode_solution ** sols_h_,ode_solcurve ** crvs_,double *sv_,double **VTm_)
  {
    // const int nobs_h = nobs-ncrv;
    #pragma omp parallel
    {
      LD_spectral_tvfield &tvf_t = *(tvfields0[LD_threads::thread_id()]);
        tvf_t.set_SVD_space(sv_,VTm_);

      ode_solchunk scw0(eor,ndep),
                   scw1(eor,ndep);
      ode_solution sw0(meta,scw0),
                   sw1(meta,scw1);

      LD_lu lu_t(tjmeta.jor+1,tjmeta.jor+1);
      double  ** const LUmat_t = lu_t.LUmat;
      int icrv_t,
          isol_t;

      #pragma omp for
      for (int iobs = 0; iobs < nobs_h; iobs++)
      {
        det.get_icrvisol_h_given_iobs_h(icrv_t,isol_t,iobs);
        ode_solution &obs0_i = *(crvs_[icrv_t]->sols[isol_t]),
                     &obs1_i = *(crvs_[icrv_t]->sols[isol_t+1]);
        ode_trivial_soljet &tjet_i = *(tcjets[icrv_t]->tjets[isol_t]);

        // build Hermites consistent with R matrix model
        sw0.copy_xu_dxun(obs0_i);
          tvf_t.set_sol_dkxu(sw0,det.kor);
        sw1.copy_xu_dxun(obs1_i);
          tvf_t.set_sol_dkxu(sw1,det.kor);

        tjet_i.set_and_solve_Hermite_jets( *(sols_h_[iobs]), lu_t, sw0,sw1 );

      }
    }
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

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of Rk matrix
        tvf_t.comp_spectral_theta_local(sols_out_[iobs]->x,sols_out_[iobs]->u);
        // the result is an improved estimate of the derivatives at interior colocation point
        tvf_t.eval_prn_theta_image(*(sols_out_[iobs]),kor_);
      }
    }
  }
  void smooth_enc_decomp_trivial_flows(double **VTm_out_,LD_svd &Rsvdg_,ode_solution **sout_,
    int nobs_,LD_R_encoder &Rkenc_,ode_solution **sin_,double *sv_,double **VTm_)
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

        // use 0'th order Hermite jet info to compute local theta via pseudoinverse of Rk matrix
        tvf_t.set_sol_dkxu( *(sout_[iobs]), det.kor );
        // the result is an improved estimate of the derivatives at interior colocation point

        // // use 0'th order Hermite jet info to compute local theta via pseudoinverse of Rk matrix
        // tvf_t.comp_spectral_theta_local(sout_[iobs]->x,sout_[iobs]->u);
        // // the result is an improved estimate of the derivatives at interior colocation point
        // tvf_t.eval_prn_theta_image(*(sout_[iobs]),det.kor);

        Rkenc_.encode_normalize_rows(
            Rkenc_.submat_i(Rsvdg_.Umat,iobs), // submatrix of R associated with observation i
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

       tvf_t.set_sol_dkxu(tjchart_i.solh_alt, det.kor, tjchart_i.solh);
        tsjet_i.set_jet_given_solh(tjchart_i.solh_alt);

       tjchart_i.sol0_alt.x = tjchart_i.sol0.x; tjchart_i.sol1_alt.x = tjchart_i.sol1.x;
       if (trunc_Hermite_exp) tsjet_i.compute_u_01(tjchart_i.sol0_alt.u,tjchart_i.sol1_alt.u,det.kor);
       else tsjet_i.compute_u_01(tjchart_i.sol0_alt.u,tjchart_i.sol1_alt.u);

       tvf_t.set_sol_dkxu(tjchart_i.sol0_alt, det.kor);
       tvf_t.set_sol_dkxu(tjchart_i.sol1_alt, det.kor);
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

  // void write_sol_h_data(const char dir_name_[], const char obs_name_[], const char dat_suff_[])
  void write_sol_h_data(const char ps_[] = "", bool write_loads_=true)
  {
    const int len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_ps = strlen(ps_),
              len_base = len_dir_name+len_obs_name+len_dat_suff+len_ps;

    const char name_Rsvd[] = ".Rsvd_g";
      char fname_Rsvd[len_base+strlen(name_Rsvd)+1];
      sprintf(fname_Rsvd,"%s%s%s%s%s",dir_name,obs_name, name_Rsvd, ps_ ,dat_suff);
      Rsvd_global.write_LD_svd(fname_Rsvd);

    // const char name_Rsvd_h[] = ".Rsvd_h_g";
    //   char fname_Rsvd_h[len_base+strlen(name_Rsvd_h)+1];
    //   sprintf(fname_Rsvd_h,"%s%s%s%s%s",dir_name,obs_name, name_Rsvd_h, ps_ ,dat_suff);
    //   Rsvd_h_global.write_LD_svd(fname_Rsvd_h);

    // const char name_theta_mat[] = ".theta_mat";
    //   char fname_theta_mat[len_base+strlen(name_theta_mat)+1];
    //   sprintf(fname_theta_mat,"%s%s%s%s%s",dir_name,obs_name, name_theta_mat, ps_ ,dat_suff);
    //   LD_io::write_Tmat<double>(fname_theta_mat,theta_mat,nobs_h,fspace0.ndof_full);

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
      sprintf(fname_jsol_h,"%s%s%s%s%s",dir_name,obs_name, name_jsol_h, ps_ ,dat_suff);
      FILE * file_jsol_h = LD_io::fopen_SAFE(fname_jsol_h,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_h);
      fwrite(obs_meta,sizeof(int),2,file_jsol_h);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_h);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(curves_h[icrv]->sols[isol]->pts,sizeof(double),det.ndim,file_jsol_h);
      LD_io::fclose_SAFE(file_jsol_h);
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h);

    // write smooth colocation point resulting from deterministic Hermite jet and Rk matrix pseudoinverse
    const char name_jsol_h_Rk[] = ".jsol_h_Rk";
      char fname_jsol_h_Rk[len_base+strlen(name_jsol_h_Rk)+1];
      sprintf(fname_jsol_h_Rk,"%s%s%s%s%s",dir_name,obs_name, name_jsol_h_Rk, ps_ ,dat_suff);
      FILE * file_jsol_h_Rk = LD_io::fopen_SAFE(fname_jsol_h_Rk,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_h_Rk);
      fwrite(obs_meta,sizeof(int),2,file_jsol_h_Rk);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_h_Rk);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->solh_alt.pts,sizeof(double),det.ndim,file_jsol_h_Rk);
      LD_io::fclose_SAFE(file_jsol_h_Rk);
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_h_Rk);

    // write new forward observations generated by exponentiating from sol_h
    const char name_jsol_1_Rk[] = ".jsol_1_Rk";
      char fname_jsol_1_Rk[len_base+strlen(name_jsol_1_Rk)+1];
      sprintf(fname_jsol_1_Rk,"%s%s%s%s%s",dir_name,obs_name, name_jsol_1_Rk, ps_,dat_suff);
      FILE * file_jsol_1_Rk = LD_io::fopen_SAFE(fname_jsol_1_Rk,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_1_Rk);
      fwrite(obs_meta,sizeof(int),2,file_jsol_1_Rk);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_1_Rk);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->sol1_alt.pts,sizeof(double),det.ndim,file_jsol_1_Rk);
      LD_io::fclose_SAFE(file_jsol_1_Rk);
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_1_Rk);

    // write new backwards observations generated by exponentiating from sol_h
    const char name_jsol_0_Rk[] = ".jsol_0_Rk";
      char fname_jsol_0_Rk[len_base+strlen(name_jsol_0_Rk)+1];
      sprintf(fname_jsol_0_Rk,"%s%s%s%s%s",dir_name,obs_name, name_jsol_0_Rk, ps_ ,dat_suff);
      FILE * file_jsol_0_Rk = LD_io::fopen_SAFE(fname_jsol_0_Rk,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_0_Rk);
      fwrite(obs_meta,sizeof(int),2,file_jsol_0_Rk);
      fwrite(npts_per_crv_h,sizeof(int),ncrv_meta,file_jsol_0_Rk);
      for (int icrv = 0; icrv < ncrv_meta; icrv++)
        for (int isol = 0; isol < npts_per_crv_h[icrv]; isol++)
          fwrite(cdets[icrv]->tjcharts[isol]->sol0_alt.pts, sizeof(double),det.ndim,file_jsol_0_Rk);
      LD_io::fclose_SAFE(file_jsol_0_Rk);
      printf("(global_Rmat_experiment::write_sol_h_data) wrote %s\n",fname_jsol_0_Rk);

  }
  // void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_, const char dir_name_[], const char obs_name_[], const char dat_suff_[])
  void write_curve_observations(ode_solcurve_chunk &Sobs_, int nsmooth_, bool verbose_=true)
  {

    const int len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_base = len_dir_name+len_obs_name+len_dat_suff;

    const size_t buf_pad = 6;

    const char name_Rsvd[] = ".Rsvd_g";
      char fname_Rsvd[len_base+strlen(name_Rsvd)+buf_pad];
      sprintf(fname_Rsvd,"%s%s%s_%d%s",dir_name,obs_name, name_Rsvd, nsmooth_ ,dat_suff);
      Rsvd_global.write_LD_svd(fname_Rsvd,false,verbose_);


    int ode_meta[2],
        &eor_meta = ode_meta[0] = Sobs_.eor,
        &ndep_meta = ode_meta[1] = Sobs_.ndep,
        obs_meta[2],
        &ncrv_meta = obs_meta[0] = Sobs_.ncrv,
        &nobs_meta = obs_meta[1] = Sobs_.nobs;

    const char name_jsol_Rk[] = ".jsol_Rk";
      char fname_jsol_Rk[len_base+strlen(name_jsol_Rk)+buf_pad];
      sprintf(fname_jsol_Rk,"%s%s%s_%d%s",dir_name,obs_name, name_jsol_Rk,nsmooth_, dat_suff);
      FILE * file_jsol_Rk = LD_io::fopen_SAFE(fname_jsol_Rk,"wb");
      fwrite(ode_meta,sizeof(int),2,file_jsol_Rk);
      fwrite(obs_meta,sizeof(int),2,file_jsol_Rk);
      fwrite(Sobs_.npts_per_crv,sizeof(int),ncrv_meta,file_jsol_Rk);
      for (int iobs = 0; iobs < nobs_meta; iobs++)
        fwrite((Sobs_.sols[iobs])->pts, sizeof(double),Sobs_.ndim,file_jsol_Rk);
      LD_io::fclose_SAFE(file_jsol_Rk);
      if (verbose_) printf("(global_Rmat_experiment::write_curve_observations) wrote %s\n",fname_jsol_Rk);
  }

  void write_initial_diagnostics()
  {
    write_sol_h_data();

    const int n_dxuk = ndep*(det.kor),
              len_dir_name = strlen(dir_name),
              len_obs_name = strlen(obs_name),
              len_dat_suff = strlen(dat_suff),
              len_base = len_dir_name+len_obs_name+len_dat_suff;
    int ode_meta[2],
        &kor_meta = ode_meta[0] = det.kor,
        &ndep_meta = ode_meta[1] = det.ndep,
        obs_meta[2],
        &ncrv_meta = obs_meta[0] = det.ncrv,
        &nobs_meta = obs_meta[1] = det.nobs;

    double * const svec_use = (Rmat_h_exp)?(svec_Rk_h_global):(svec_Rk_global),
           ** const VTmat_use = (Rmat_h_exp)?(VTmat_Rk_h_global):(VTmat_Rk_global);
    LD_spectral_tvfield &tvf = *(tvfields0[0]);
      tvf.set_SVD_space(svec_use,VTmat_use);

    ode_solchunk solchunk_work(eor,ndep);
    ode_solution sol_work(meta,solchunk_work);

    const size_t buf_pad = 6;
    const char name_dxuk_Rk[] = ".dxuk_Rk";
      char fname_dxuk_Rk[len_base+strlen(name_dxuk_Rk)+buf_pad];
      sprintf(fname_dxuk_Rk,"%s%s%s%s",dir_name,obs_name, name_dxuk_Rk, dat_suff);
      FILE * file_dxuk_Rk = LD_io::fopen_SAFE(fname_dxuk_Rk,"wb");
      fwrite(ode_meta,sizeof(int),2,file_dxuk_Rk);
      fwrite(obs_meta,sizeof(int),2,file_dxuk_Rk);
      fwrite(det.npts_per_crv,sizeof(int),ncrv_meta,file_dxuk_Rk);
      for (int iobs = 0; iobs < nobs_meta; iobs++)
      {
        sol_work.copy_xu(*(sols[iobs]));

        tvf.comp_spectral_theta_local(sol_work.x,sol_work.u);
          tvf.eval_prn_theta_image(sol_work,det.kor);
        // writing out the Rk matrix pseudoinverse prediction of dx u^(k) at noisy u
        fwrite(sol_work.dxu,sizeof(double),n_dxuk,file_dxuk_Rk);
      }
      LD_io::fclose_SAFE(file_dxuk_Rk);
      printf("(global_Rmat_experiment::write_initial_diagnostics) wrote %s\n",fname_dxuk_Rk);
  }
  // void write_denoising_summary(ode_solcurve_chunk &Sobs_, int nsmooth_, double *res_vec_)
  // {
  //   write_curve_observations(Sobs_,nsmooth_,true);
  //   write_sol_h_data("_f");
  //   // const int len_dir_name = strlen(dir_name),
  //   //           len_obs_name = strlen(obs_name),
  //   //           len_dat_suff = strlen(dat_suff),
  //   //           len_base = len_dir_name+len_obs_name+len_dat_suff;
  //
  // }

};

global_Rmat_experiment Rk_experiment(detector,function_space);

int main()
{

  // observations.print_details();
  // detector.print_details();

  // Rk_experiment.write_sol_h_data();
  Rk_experiment.denoise_data(observations_twin);

  return 0;
}
