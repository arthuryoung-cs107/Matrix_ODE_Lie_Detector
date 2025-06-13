#ifndef LD_NOISE_AUX_HH
#define LD_NOISE_AUX_HH

#include "LD_framework.hh"


class LD_noise_aux
{

  public:

    LD_noise_aux() {}
    ~LD_noise_aux() {}

    static void perturb_pts(ode_curve_observations &obs_, int noise_level_, int seed_=1234, double nse_scl_=1.0)
    {
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim;
      double  sigma_noise_pts[ndim],
              sigma_noise_pts_scaling[ndim];
      LD_noise_aux::specify_sigma_noise_pts_scaling_median_mag(sigma_noise_pts_scaling,obs_,nse_scl_);

      // make corrections to noise scales
      sigma_noise_pts[0] = (nse_scl_*0.1)*fabs(obs_.pts_in[ndim] - obs_.pts_in[0]);

      LD_noise_aux::specify_sigma_noise_pts(sigma_noise_pts,obs_,sigma_noise_pts_scaling,noise_level_);
      LD_noise_aux::apply_noise_pts(obs_,noise_level_,sigma_noise_pts,seed_);
    }
    static void specify_sigma_noise_pts_scaling_median_mag(double *sigma_noise_pts_scaling_,ode_curve_observations &obs_,double scl_=1.0)
    {
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nobs = obs_.nobs;
      int i_med_1 = (int)(nobs/2),
          i_med_2 = i_med_1 + (nobs%2);
      double  sort_vec[nobs],
              &med_mag_1 = sort_vec[i_med_1],
              &med_mag_2 = sort_vec[i_med_2],
              * const pts_chunk = obs_.pts_in;
      for (size_t idim = 0; idim < ndim; idim++)
      {
        for (size_t i_obs = 0; i_obs < nobs; i_obs++) sort_vec[i_obs] = fabs(pts_chunk[(i_obs*ndim) + idim]);
        LD_linalg::sort_vec_inc<double>(sort_vec,nobs);
        sigma_noise_pts_scaling_[idim] = scl_*(0.5*(med_mag_1 + med_mag_2));
      }
    }
    static void specify_sigma_noise_pts(double *sigma_noise_pts_,ode_curve_observations &obs_,double *dim_scale_,int noise_level_)
    {
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim;
      if (noise_level_>=0)
      {
        double  sigma_base = pow(10.0,(double)(-noise_level_));
        for (size_t i = 0; i < ndim; i++) sigma_noise_pts_[i] = sigma_base*dim_scale_[i];
      }
      else for (int i = 0; i < ndim; i++) sigma_noise_pts_[i] = 0.0;
    }
    static void apply_noise_pts(ode_curve_observations &obs_,int noise_level_,double *sigma_noise_pts_,int seed_=13)
    {
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nobs = obs_.nobs;
      if (noise_level_>=0)
      {
        LD_rng ran(seed_);
        double *pts_j = obs_.pts_in;
        // noise points
        for (size_t j_obs = 0; j_obs < nobs; j_obs++, pts_j+=ndim)
          for (size_t i_dim = 0; i_dim < ndim; i_dim++)
            pts_j[i_dim] += ran.rand_gau(0.0,sigma_noise_pts_[i_dim]);
      }
    }

};

#endif
