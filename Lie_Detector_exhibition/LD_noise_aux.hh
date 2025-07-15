#ifndef LD_NOISE_AUX_HH
#define LD_NOISE_AUX_HH

#include "LD_framework.hh"


class LD_noise_aux
{

  ode_solspc_meta &meta;

  public:

    LD_noise_aux(ode_solspc_meta &meta_): meta(meta_),
      sigma_s(new double[meta_.ndim]),
      sigma_dnp1xu(new double[meta_.ndep]),
      sigma_JFs(new double[meta_.ndim*meta_.ndep])
      {}
    ~LD_noise_aux()
    {
      delete [] sigma_s;
      delete [] sigma_dnp1xu;
      delete [] sigma_JFs;
    }
    double  * const sigma_s,
            * const sigma_dnp1xu,
            * const sigma_JFs;
    void write_sigmas(const char name_[])
    {
      int hlen = 2,
          header[] = {hlen,meta.eor,meta.ndep};
      FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
      fwrite(header,sizeof(int),hlen+1,file_out);
      fwrite(sigma_s,sizeof(double),meta.ndim,file_out);
      fwrite(sigma_dnp1xu,sizeof(double),meta.ndep,file_out);
      fwrite(sigma_JFs,sizeof(double),meta.ndim*meta.ndep,file_out);
      LD_io::fclose_SAFE(file_out);
      printf("(LD_noise_aux::write_sigmas) wrote %s\n",name_);
    }


    static const char * dir_name;

    static void perturb_JFs(ode_curve_observations &obs_,double *sigma_noise_JFs_,int noise_level_,int seed_=1234, double nse_scl_=1.0)
    {
      const int eor = obs_.eor,
                ndep = obs_.ndep,
                nobs = obs_.nobs;
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nJFe = ndep*ndim;
      double  sigma_noise_JFs_scaling[nJFe];
      LD_noise_aux::specify_sigma_noise_JFs_scaling_median_mag(sigma_noise_JFs_scaling,obs_,nse_scl_);
      LD_noise_aux::specify_sigma_noise_JFs(sigma_noise_JFs_,obs_,sigma_noise_JFs_scaling,noise_level_);
      LD_noise_aux::apply_noise_JFs(obs_,noise_level_,sigma_noise_JFs_,seed_);
    }
    static void perturb_dnp1xu(ode_curve_observations &obs_,double *sigma_noise_dnp1xu_,int noise_level_, int seed_=1234, double nse_scl_=1.0)
    {
      const int ndep = obs_.ndep;
      double  sigma_noise_dnp1xu_scaling[ndep];
      LD_noise_aux::specify_sigma_noise_dnp1xu_scaling_median_mag(sigma_noise_dnp1xu_scaling,obs_,nse_scl_);
      LD_noise_aux::specify_sigma_noise_dnp1xu(sigma_noise_dnp1xu_,obs_,sigma_noise_dnp1xu_scaling,noise_level_);
      LD_noise_aux::apply_noise_dnp1xu(obs_,noise_level_,sigma_noise_dnp1xu_,seed_);
    }
    static void perturb_pts(ode_curve_observations &obs_, double sigma_noise_pts_[],
        int noise_level_, int seed_=1234, double nse_scl_=1.0, double xnse_scl_=1.0)
    {
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim;
      double  sigma_noise_pts_scaling[ndim];
      LD_noise_aux::specify_sigma_noise_pts_scaling_median_mag(sigma_noise_pts_scaling,obs_,nse_scl_);

      // make corrections to noise scales
      sigma_noise_pts_scaling[0] = (nse_scl_*xnse_scl_)*fabs(obs_.pts_in[ndim] - obs_.pts_in[0]);

      LD_noise_aux::specify_sigma_noise_pts(sigma_noise_pts_,obs_,sigma_noise_pts_scaling,noise_level_);

      // make corrections to final variances
      // sigma_noise_pts_[0] = 0.0;

      LD_noise_aux::apply_noise_pts(obs_,noise_level_,sigma_noise_pts_,seed_);
    }

    static void specify_sigma_noise_JFs_scaling_median_mag(double *sigma_noise_JFs_scaling_,ode_curve_observations &obs_,double scl_=1.0)
    {
      const int eor = obs_.eor,
                ndep = obs_.ndep,
                nobs = obs_.nobs;
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nJFe = ndep*ndim;
      int i_med_1 = (int)(nobs/2),
          i_med_2 = i_med_1 + (nobs%2);
      double  sort_vec[nobs],
              &med_mag_1 = sort_vec[i_med_1],
              &med_mag_2 = sort_vec[i_med_2],
              * const JFs_chunk = obs_.JFs_in;
      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t idim = 0; idim < ndim; idim++)
        {
          for (size_t i_obs = 0; i_obs < nobs; i_obs++)
            sort_vec[i_obs] = fabs(JFs_chunk[(i_obs*nJFe) + (idep*ndim) + idim]);
          LD_linalg::sort_vec_inc<double>(sort_vec,nobs);
          sigma_noise_JFs_scaling_[idep] = scl_*(0.5*(med_mag_1 + med_mag_2));
        }

    }
    static void specify_sigma_noise_dnp1xu_scaling_median_mag(double *sigma_noise_dnp1xu_scaling_,ode_curve_observations &obs_,double scl_=1.0)
    {
      const int ndep = obs_.ndep,
                nobs = obs_.nobs;
      int i_med_1 = (int)(nobs/2),
          i_med_2 = i_med_1 + (nobs%2);
      double  sort_vec[nobs],
              &med_mag_1 = sort_vec[i_med_1],
              &med_mag_2 = sort_vec[i_med_2],
              * const dnp1xu_chunk = obs_.dnp1xu_in;
      for (size_t idep = 0; idep < ndep; idep++)
      {
        for (size_t i_obs = 0; i_obs < nobs; i_obs++)
          sort_vec[i_obs] = fabs(dnp1xu_chunk[(i_obs*ndep) + idep]);
        LD_linalg::sort_vec_inc<double>(sort_vec,nobs);
        sigma_noise_dnp1xu_scaling_[idep] = scl_*(0.5*(med_mag_1 + med_mag_2));
      }
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

    static void specify_sigma_noise_JFs(double *sigma_noise_JFs_,ode_curve_observations &obs_,double *dim_scale_,int noise_level_)
    {
      const int eor = obs_.eor,
                ndep = obs_.ndep,
                nobs = obs_.nobs;
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nJFe = ndep*ndim;
      if (noise_level_>=0)
      {
        double  sigma_base = pow(10.0,(double)(-noise_level_));
        for (size_t i = 0; i < nJFe; i++) sigma_noise_JFs_[i] = sigma_base*dim_scale_[i];
      }
      else for (int i = 0; i < nJFe; i++) sigma_noise_JFs_[i] = 0.0;
    }
    static void specify_sigma_noise_dnp1xu(double *sigma_noise_dnp1xu_,ode_curve_observations &obs_,double *dim_scale_,int noise_level_)
    {
      const int ndep = obs_.ndep;
      if (noise_level_>=0)
      {
        double  sigma_base = pow(10.0,(double)(-noise_level_));
        for (size_t i = 0; i < ndep; i++) sigma_noise_dnp1xu_[i] = sigma_base*dim_scale_[i];
      }
      else for (int i = 0; i < ndep; i++) sigma_noise_dnp1xu_[i] = 0.0;
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

    static void apply_noise_JFs(ode_curve_observations &obs_,int noise_level_,double *sigma_noise_JFs_,int seed_=13)
    {
      const int eor = obs_.eor,
                ndep = obs_.ndep,
                nobs = obs_.nobs;
      ode_solspc_meta meta_local(obs_.eor,obs_.ndep);
      const int ndim = meta_local.ndim,
                nJFe = ndep*ndim;
      if (noise_level_>=0)
      {
        LD_rng ran(seed_);
        double *JFs_j = obs_.JFs_in;
        // noise points
        for (size_t j_obs = 0; j_obs < nobs; j_obs++, JFs_j+=nJFe)
          for (size_t iJFs = 0; iJFs < nJFe; iJFs++)
            JFs_j[iJFs] += ran.rand_gau(0.0,sigma_noise_JFs_[iJFs]);
      }
    }
    static void apply_noise_dnp1xu(ode_curve_observations &obs_,int noise_level_,double *sigma_noise_dnp1xu_,int seed_=13)
    {
      const int ndep = obs_.ndep,
                nobs = obs_.nobs;
      if (noise_level_>=0)
      {
        LD_rng ran(seed_);
        double *dnp1xu_j = obs_.dnp1xu_in;
        // noise points
        for (size_t j_obs = 0; j_obs < nobs; j_obs++, dnp1xu_j+=ndep)
          for (size_t i_dep = 0; i_dep < ndep; i_dep++)
            dnp1xu_j[i_dep] += ran.rand_gau(0.0,sigma_noise_dnp1xu_[i_dep]);
      }
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
