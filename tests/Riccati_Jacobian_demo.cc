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

const char dir_name[] = "./denoise_data_directory/Uniform_IC_perturbation/"; // data directory
const char obs_name[] = "Riccati_xrange0_true_DoP853gen"; // name of observations file
const char dat_suff[] = ".lddat"; // data file suffix (*.lddat)

char data_name[strlen(dir_name) + strlen(obs_name) + strlen(dat_suff) + 1];
const int data_name_len = sprintf(data_name,"%s%s%s", dir_name,obs_name,dat_suff);

const char dnp1xu_suffix[] = "_dnp1xu";
char data_dnp1xu_name[strlen(dir_name) + strlen(obs_name) + strlen(dnp1xu_suffix) + strlen(dat_suff) + 1];
const int data_dnp1xu_name_len = sprintf(data_dnp1xu_name,"%s%s%s%s", dir_name,obs_name,dnp1xu_suffix,dat_suff);

const char JFs_suffix[] = "_JFs";
char data_JFs_name[strlen(dir_name) + strlen(obs_name) + strlen(JFs_suffix) + strlen(dat_suff) + 1];
const int data_JFs_name_len = sprintf(data_JFs_name,"%s%s%s%s", dir_name,obs_name,JFs_suffix,dat_suff);

const int curve_Lie_detector::combine_flag = 1; // 0, flag for updating smoothened coordinates, (0 lazy, 1 aggressive)
const bool curve_Lie_detector::truncate_Hermite_exp = true;

ode_curve_observations observations_full(data_name,data_dnp1xu_name,data_JFs_name);
ode_curve_observations observations(data_name);

const int kor_obs = 1;

Lie_detector detector(observations);
  ode_solcurve_chunk observations_twin(detector.meta0,detector.ncrv,detector.npts_per_crv); // equal size data set for workspace
orthopolynomial_space function_space(detector.meta0,3); // initializes as unmapped multinomials

Lie_detector_digital_twin LDtwin(detector,observations_twin,function_space);

// global_Rmat_experiment Rk_experiment(detector,function_space);

struct cubic_model : public global_multinomial_experiment
{
  const int nthread_wkspc;

  LD_R_encoder R1enc;
  LD_svd  &R1svd_global;

  // theta (parameter) space variables
  double * const svec_R1_global, // singular values of primary Rn matrix
         ** const VTmat_R1_global; // transposed row space of primary Rn matrix

  // thread local workspaces
  // double *** const xu_snap_wkspc;
  orthopolynomial_basis ** const fbases0;
  LD_spectral_tvfield ** const tvfields0;

  cubic_model(Lie_detector &det_,orthopolynomial_space &fspace0_) :
    global_multinomial_experiment(det_,fspace0_,det_.kor),
    nthread_wkspc(LD_threads::numthreads()),
    R1enc(det_.meta0,det_.kor),
    R1svd_global(Asvd_global),
    svec_R1_global(R1svd_global.svec),
    VTmat_R1_global(VTmat_global),
    fbases0(new orthopolynomial_basis*[nthread_wkspc]),
    tvfields0(new LD_spectral_tvfield*[nthread_wkspc])
  {
    for (int i = 0; i < nthread_wkspc; i++)
    {
      fbases0[i] = new orthopolynomial_basis(fspace0);
      tvfields0[i] = new LD_spectral_tvfield(*(fbases0[i]),svec_R1_global,VTmat_R1_global);
    }

    // generate cubic digital twin of observations
    double t0 = LD_threads::tic();

    encode_matrix(R1svd_global.Umat,R1enc,fbases0,sols,nobs,false );
    int rank_out = telescope_decompose_global_matrix( VTmat_R1_global,
                                                      R1svd_global,
                                                      R1enc.ncod, false);
    R1svd_global.print_result("R1svd_global");

    /*
      computing this singular value decomposition induces a trivial vector field model in the parameter space
    */

  }
  ~cubic_model()
  {
    for (int i = 0; i < nthread_wkspc; i++) {delete tvfields0[i]; delete fbases0[i]; }
    delete [] tvfields0; delete [] fbases0;
  }

  void evaluate_vartheta_set()
  {

    /*
      evaluate the image of the data set over its own parameter map, vartheta
      use these local parameters to estimate the n+1=2nd derivative
    */
    #pragma omp parallel num_threads(1)
    {
      LD_spectral_tvfield &tvf_t = *(tvfields0[LD_threads::thread_id()]);
        tvf_t.set_SVD_space(R1svd_global.svec,VTmat_R1_global);
      #pragma omp for
      for (int i = 0; i < nobs; i++)
      {
        ode_solution &sol_i = *(sols[i]),
                     &sol_i_twin = *(LDtwin.soltwin_i(i));
        double * const theta_i = LDtwin.theta_space[i];

        printf("sol_i_twin: "); sol_i_twin.print_sol();
        printf("sol_i_obs: "); sol_i.print_sol();
        printf("sol_i_obs_full: "); observations_full.print_sol_i(i);

        LDtwin.set_soltwin_i(i,sol_i);
        printf("sol_i_twin: "); sol_i_twin.print_sol();

        // compute vartheta at this value of x and u
        tvf_t.comp_spectral_theta_local(theta_i,sol_i_twin.x,sol_i_twin.u);

        // use vartheta alone to estimate dxu, then d2xu. Write to twin workspace
        tvf_t.eval_prn_theta_image(sol_i_twin,2,theta_i);
        printf("sol_i_twin: "); sol_i_twin.print_sol();

        // use vartheta, x, u, and dxu to estimate d2xu. (over)Write to minimal observations workspace
        tvf_t.eval_vdnxu_theta_image(sol_i,theta_i);
        printf("sol_i_obs: "); sol_i.print_sol();

        getchar();

      }
    }

  }

};

cubic_model cube(detector,function_space);

int main()
{
  cube.evaluate_vartheta_set();
  printf("done.\n");
  return 0;
}
