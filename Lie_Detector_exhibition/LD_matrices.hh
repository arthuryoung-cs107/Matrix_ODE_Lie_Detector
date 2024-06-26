#ifndef LD_MATS_HH
#define LD_MATS_HH

#include "LD_framework.hh"

struct LD_L_matrix: public LD_matrix
{
  LD_L_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,1,fspc_.perm_len) {}
  ~LD_L_matrix() {}

  static void eval_L_signal_strength(LD_observations_set &Sobs_,LD_Theta_stack &stack_,function_space &fspc_,double tol_=1e-13,bool verbose_=true)
  {
    const int ncrv = stack_.nspc,
              npts_max = Sobs_.max_npts_curve(),
              ncol_use = stack_.nrow_use;
    bool satisfy_crvmat[ncrv][ncol_use];
    int * const nsat_theta = stack_.ntheta_spcvec,
        ** const isat_theta = stack_.itheta_spcmat;
    double *** const Theta_tns = stack_.Theta_tns;

    ode_solcurve ** const crvs = Sobs_.curves;
    int npts_evl = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:npts_evl)
    {
      lamda_map_eval_workspace wkspc_t(ncol_use,fspc_,npts_max);
      bool  ** const Sat_t = wkspc_t.satisfy_flags_mat;
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrv; icrv++)
      {
        ode_solcurve &crv_i = *(crvs[icrv]);
        ode_solution ** const sols_i = crv_i.sols;

        const int nobs_crvi = crv_i.nobs;
        bool * const sat_crvi = satisfy_crvmat[icrv];

        npts_evl += nobs_crvi;
        wkspc_t.Theta = Theta_tns[icrv];
        LD_linalg::fill_vec<bool>(sat_crvi,ncol_use,false);

        int nsat_crvi = 0;
        for (size_t jsol = 0; jsol < nobs_crvi; jsol++)
        {
          if (Lie_detector::eval_lamfunc_signal_strength(Sat_t[jsol],*(sols_i[jsol]),wkspc_t,fspc_,tol_) != ncol_use)
          {
            nsat_crvi = 0;
            for (size_t ith = 0; ith < ncol_use; ith++)
              nsat_crvi += (int)(sat_crvi[ith] = (sat_crvi[ith]) || (Sat_t[jsol][ith]));
          }
          else LD_linalg::fill_vec<bool>(sat_crvi,nsat_crvi=ncol_use,true);
          if (nsat_crvi==ncol_use) break;
        }
        if (nsat_theta[icrv] = nsat_crvi)
        {
          for (size_t ith = 0, isat = 0; ith < ncol_use; ith++)
            if (sat_crvi[ith]) isat_theta[icrv][isat++] = ith;
        }
        else LD_linalg::fill_vec<int>(isat_theta[icrv],ncol_use,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_L_matrix::eval_L_signal_strength) evaluated %d conds. (sols) over %d bases (%d x %d) in %.4f seconds (%d threads)\n",
        npts_evl,
        ncrv,ncol_use,stack_.ndof_use,
        work_time,LD_threads::numthreads());
    }
  }

  template <class BSIS> void populate_L_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_L_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Lmat_i_)
  {
    const double L_row_i_mag = LD_linalg::norm_l2(Lmat_i_[0],perm_len);
    for (size_t i_L = 0; i_L < perm_len; i_L++) Lmat_i_[0][i_L] /= L_row_i_mag;
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_L_rows(chunk_,sol_,Amat_i_);}

    inline void fill_L_rows(partial_chunk &chunk_, ode_solution &sol_, double **Lmat_i_)
    {
      double  * const Lrow_i = Lmat_i_[0],
              * const lambda_x_vec = chunk_.Jac_mat[0];
      for (size_t i_L = 0; i_L < perm_len; i_L++) Lrow_i[i_L] = (dof_tun_flags[i_L])?(lambda_x_vec[i_L]):(0.0);
      if (normalization_flag) normalize_submat_rows(sol_,Lmat_i_);
    }
};

struct LD_R_matrix: public LD_matrix
{
  LD_R_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep*fspc_.eor) {}
  ~LD_R_matrix() {}

  template <class BSIS> static void eval_Rn_condition(LD_observations_set &Sobs_,LD_Theta_stack &stack_,BSIS **bases_,double tol_=1e-3,bool verbose_=true)
  {
    const int eor = Sobs_.eor,
              ncrv = stack_.nspc,
              ncol_use = stack_.nrow_use;
    bool satisfy_crvmat[ncrv][ncol_use];
    int * const nsat_theta = stack_.ntheta_spcvec,
        ** const isat_theta = stack_.itheta_spcmat;
    double *** const Theta_tns = stack_.Theta_tns;

    ode_solcurve ** const crvs = Sobs_.curves;

    int npts_evl = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:npts_evl)
    {
      BSIS &bsis_t = *(bases_[LD_threads::thread_id()]);
      Rn_cond_eval_workspace wkspc_t(ncol_use,Sobs_.meta,1);
      bool * const sat_t = wkspc_t.satisfy_flags;
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrv; icrv++)
      {
        ode_solcurve &crv_i = *(crvs[icrv]);
        ode_solution ** const sols_i = crv_i.sols;

        const int nobs_crvi = crv_i.nobs;
        bool * const sat_crvi = satisfy_crvmat[icrv];

        npts_evl += nobs_crvi;
        wkspc_t.Theta = Theta_tns[icrv];
        LD_linalg::fill_vec<bool>(sat_crvi,ncol_use,true);

        int nsat_crvi = 0;
        for (size_t jsol = 0; jsol < nobs_crvi; jsol++)
        {
          if (Lie_detector::eval_Theta_Rn_condition(sat_t,*(sols_i[jsol]),wkspc_t,bsis_t,eor,tol_))
          {
            nsat_crvi = 0;
            for (size_t ith = 0; ith < ncol_use; ith++)
              nsat_crvi += (int)(sat_crvi[ith] = (sat_crvi[ith])&&(sat_t[ith]));
          }
          else LD_linalg::fill_vec<bool>(sat_crvi,ncol_use,nsat_crvi = 0);
          if (!nsat_crvi) break;
        }
        if (nsat_theta[icrv] = nsat_crvi)
        {
          for (size_t ith = 0, isat = 0; ith < ncol_use; ith++)
            if (sat_crvi[ith]) isat_theta[icrv][isat++] = ith;
        }
        else LD_linalg::fill_vec<int>(isat_theta[icrv],ncol_use,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_R_matrix::eval_Rn_condition) evaluated %d (%d sols, n=%d, q=%d) conds. over %d bases (%d x %d) in %.4f seconds (%d threads)\n",
        npts_evl*eor*Sobs_.ndep,
        npts_evl,eor,Sobs_.ndep,
        ncrv,ncol_use,stack_.ndof_use,
        work_time,LD_threads::numthreads());
    }
  }

  template <class BSIS> void populate_R_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_Rn_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Rmat_i_)
  {
    for (size_t k = 0, idxu = 0; k < eor; k++)
      for (size_t idep = 0; idep < ndep; idep++, idxu++)
        for (size_t icol = 0; icol < ndof_tun; icol++)
          Rmat_i_[idxu][icol] /= sol_.dxu[idxu];
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_Rn_rows(chunk_,sol_,Amat_i_);}

    inline void fill_Rn_rows(partial_chunk &chunk_, ode_solution &sol_, double **Rmat_i_)
    {
      int i_dof = 0;
      double  * const dxu = sol_.dxu,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
          fill_x_Rn_columns(i_dof++,Rmat_i_,dxu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L]);

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_u_Rn_columns(i_dof++,idep,Rmat_i_,Lambda_xu_mat[idep+1][i_L],Jac_utheta_vdxu_mat[i_L]);

      if (normalization_flag) normalize_submat_rows(sol_,Rmat_i_);
    }
    inline void fill_x_Rn_columns(int i_dof, double **Rmat_, double *dxu_, double &lambda_, double **Jac_xtheta_i_vdxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++) Rmat_[idep][i_dof] = dxu_[idep]*lambda_; // 1st order contribution
      for (size_t k = 1, i_R = ndep; k < eor; k++)
        for (size_t idep = 0; idep < ndep; idep++, i_R++)
          Rmat_[i_R][i_dof] = dxu_[i_R]*lambda_ - Jac_xtheta_i_vdxu_[k-1][idep]; // subtract contribution to v1n
    }
    inline void fill_u_Rn_columns(int i_dof, int idep, double **Rmat_, double &lambda_, double *del_utheta_i)
    {
      Rmat_[idep][i_dof] = -(lambda_);
      for (size_t k = 1, i_R = idep + ndep; k < eor; k++, i_R += ndep)
        Rmat_[i_R][i_dof] = -(del_utheta_i[k-1]);
    }
};

struct LD_O_matrix: public LD_matrix
{
  LD_O_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep) {}
  ~LD_O_matrix() {}

  template <class BSIS> void populate_O_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_O_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Omat_i_)
  {
    for (size_t idep = 0; idep < ndep; idep++)
      for (size_t icol = 0; icol < ndof_tun; icol++)
        Omat_i_[idep][icol] /= sol_.dnxu[idep];
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_O_rows(chunk_,sol_,Amat_i_);}

    inline void fill_O_rows(partial_chunk &chunk_, ode_solution &sol_, double **Omat_i_)
    {
      int i_dof = 0;
      double  * const lambda_x_vec = chunk_.Jac_mat[0];

      if (eor==1)
      {
        double  ** const Lambda_xu_mat = chunk_.Jac_mat,
                * const dxu = sol_.dxu;

        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_x_O1_columns(i_dof++,Omat_i_,dxu,lambda_x_vec[i_L]);

        for (size_t idep = 0; idep < ndep; idep++)
          for (size_t i_L = 0; i_L < perm_len; i_L++)
            if (dof_tun_flags[i_L])
              fill_u_O1_columns(Omat_i_[idep][i_dof++],Lambda_xu_mat[idep+1][i_L]);
      }
      else
      {
        double  *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
                ** const Jac_utheta_vdxu_mat = chunk_.C_u,
                * const dnxu = sol_.dnxu;

        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_x_O_columns(i_dof++,Omat_i_,dnxu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L][eor-2]);

        for (size_t idep = 0; idep < ndep; idep++)
          for (size_t i_L = 0; i_L < perm_len; i_L++)
            if (dof_tun_flags[i_L])
              fill_u_O_columns(Omat_i_[idep][i_dof++],Jac_utheta_vdxu_mat[i_L][eor-2]);
      }

      if (normalization_flag) normalize_submat_rows(sol_,Omat_i_);
    }

    inline void fill_x_O1_columns(int i_dof, double **Omat_, double *dxu_, double &lambda_)
      {for (size_t idep = 0; idep < ndep; idep++) Omat_[idep][i_dof] = dxu_[idep]*lambda_;}
    inline void fill_u_O1_columns(double &Omat_idep_idof_, double &lambda_)
      {Omat_idep_idof_ = -(lambda_);}

    inline void fill_x_O_columns(int i_dof,double **Omat_,double *dnxu_,double &lambda_,double *del_xtheta_i_vdnxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++)
        Omat_[idep][i_dof] = dnxu_[idep]*lambda_ - del_xtheta_i_vdnxu_[idep];
    }
    inline void fill_u_O_columns(double &Omat_idep_idof_,double &del_untheta_i)
      {Omat_idep_idof_ = -(del_untheta_i);}
};

struct LD_P_matrix: public LD_matrix
{
  LD_P_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep) {}
  ~LD_P_matrix() {}

  template <class BSIS> void populate_P_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_P_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Pmat_i_)
  {
    for (size_t idep = 0; idep < ndep; idep++)
      for (size_t icol = 0; icol < ndof_tun; icol++)
        Pmat_i_[idep][icol] /= sol_.dnp1xu[idep];
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_P_rows(chunk_,sol_,Amat_i_);}

    inline void fill_P_rows(partial_chunk &chunk_, ode_solution &sol_, double **Pmat_i_)
    {
      int i_dof = 0;
      double  * const dnp1xu = sol_.dnp1xu,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
          fill_x_P_columns(i_dof++,Pmat_i_,dnp1xu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L][eor-1]);

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_u_P_columns(Pmat_i_[idep][i_dof++],Jac_utheta_vdxu_mat[i_L][eor-1]);

      if (normalization_flag) normalize_submat_rows(sol_,Pmat_i_);
    }
    inline void fill_x_P_columns(int i_dof,double **Pmat_,double *dnp1xu_,double &lambda_,double *Jac_xtheta_i_vdnxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++)
        Pmat_[idep][i_dof] = dnp1xu_[idep]*lambda_ - Jac_xtheta_i_vdnxu_[idep];
    }
    inline void fill_u_P_columns(double &Pmat_ij_,double par_utheta_i) {Pmat_ij_ = -(par_utheta_i);}
};

struct LD_Q_matrix: public LD_matrix
{
  LD_Q_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep*(fspc_.eor+1)) {}
  ~LD_Q_matrix() {}

  template <class BSIS> void populate_Q_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_Q_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Qmat_i_)
  {
    size_t idxu = 0;
    for (size_t k = 0; k < eor; k++)
      for (size_t idep = 0; idep < ndep; idep++, idxu++)
        for (size_t icol = 0; icol < ndof_tun; icol++)
          Qmat_i_[idxu][icol] /= sol_.dxu[idxu];
    for (size_t idep = 0; idep < ndep; idep++, idxu++)
      for (size_t icol = 0; icol < ndof_tun; icol++)
        Qmat_i_[idxu][icol] /= sol_.dnp1xu[idep];
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_Q_rows(chunk_,sol_,Amat_i_);}

    inline void fill_Q_rows(partial_chunk &chunk_, ode_solution &sol_, double **Qmat_i_)
    {
      int i_dof = 0;
      double  * const dxu = sol_.dxu,
              * const dnp1xu = sol_.dnp1xu,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
          fill_x_Q_columns(i_dof++,Qmat_i_,dxu,dnp1xu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L]);

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_u_Q_columns(i_dof++,idep,Qmat_i_,Lambda_xu_mat[idep+1][i_L],Jac_utheta_vdxu_mat[i_L]);

      if (normalization_flag) normalize_submat_rows(sol_,Qmat_i_);
    }
    inline void fill_x_Q_columns(int i_dof,double **Qmat_,double *dxu_,double *dnp1xu_,double &lambda_,double **Jac_xtheta_i_vdxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++) Qmat_[idep][i_dof] = dxu_[idep]*lambda_; // 1st order contribution
      for (size_t k = 1, i_Q = ndep; k < eor; k++)
        for (size_t idep = 0; idep < ndep; idep++, i_Q++)
          Qmat_[i_Q][i_dof] = dxu_[i_Q]*lambda_ - Jac_xtheta_i_vdxu_[k-1][idep]; // subtract contribution to v1n
      for (size_t idep = 0, i_Q = ndep*eor; idep < ndep; idep++, i_Q++)
        Qmat_[i_Q][i_dof] = dnp1xu_[idep]*lambda_ - Jac_xtheta_i_vdxu_[eor-1][idep];
    }
    inline void fill_u_Q_columns(int i_dof,int idep,double **Qmat_,double &lambda_,double *del_utheta_i)
    {
      Qmat_[idep][i_dof] = -(lambda_);
      for (size_t k = 1, i_Q = idep + ndep; k <= eor; k++, i_Q += ndep)
        Qmat_[i_Q][i_dof] = -(del_utheta_i[k-1]);
    }
};

struct LD_G_matrix: public LD_matrix
{
  LD_G_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep) {}
  ~LD_G_matrix() {}

  template <class BSIS> static void eval_inf_criterion(LD_observations_set &Sobs_,LD_Theta_stack &stack_,BSIS **bases_,double tol_=1e-3,bool verbose_=true)
  {
    const int ncrv = stack_.nspc,
              ncol_use = stack_.nrow_use;
    bool satisfy_crvmat[ncrv][ncol_use];
    int * const nsat_theta = stack_.ntheta_spcvec,
        ** const isat_theta = stack_.itheta_spcmat;
    double *** const Theta_tns = stack_.Theta_tns;

    ode_solcurve ** const crvs = Sobs_.curves;

    int npts_evl = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:npts_evl)
    {
      BSIS &bsis_t = *(bases_[LD_threads::thread_id()]);
      inf_crit_eval_workspace wkspc_t(ncol_use,Sobs_.meta,1);
      bool * const sat_t = wkspc_t.satisfy_flags;
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrv; icrv++)
      {
        ode_solcurve &crv_i = *(crvs[icrv]);
        ode_solution ** const sols_i = crv_i.sols;

        const int nobs_crvi = crv_i.nobs;
        bool * const sat_crvi = satisfy_crvmat[icrv];

        npts_evl += nobs_crvi;
        wkspc_t.Theta = Theta_tns[icrv];
        LD_linalg::fill_vec<bool>(sat_crvi,ncol_use,true);

        int nsat_crvi = 0;
        for (size_t jsol = 0; jsol < nobs_crvi; jsol++)
        {
          if (Lie_detector::eval_Theta_inf_criterion(sat_t,*(sols_i[jsol]),wkspc_t,bsis_t,tol_))
          {
            nsat_crvi = 0;
            for (size_t ith = 0; ith < ncol_use; ith++)
              nsat_crvi += (int)(sat_crvi[ith] = (sat_crvi[ith])&&(sat_t[ith]));
          }
          else LD_linalg::fill_vec<bool>(sat_crvi,ncol_use,nsat_crvi = 0);
          if (!nsat_crvi) break;
        }
        if (nsat_theta[icrv] = nsat_crvi)
        {
          for (size_t ith = 0, isat = 0; ith < ncol_use; ith++)
            if (sat_crvi[ith]) isat_theta[icrv][isat++] = ith;
        }
        else LD_linalg::fill_vec<int>(isat_theta[icrv],ncol_use,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_G_matrix::eval_inf_criterion) evaluated %d (%d sols, q=%d) conds. over %d bases (%d x %d) in %.4f seconds (%d threads)\n",
        npts_evl*Sobs_.ndep,
        npts_evl,Sobs_.ndep,
        ncrv,ncol_use,stack_.ndof_use,
        work_time,LD_threads::numthreads());
    }
  }

  template <class BSIS> void populate_G_matrix(BSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_G_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline void normalize_submat_rows(ode_solution &sol_,double **Gmat_i_)
  {
    for (size_t idep = 0; idep < ndep; idep++)
    {
      // const double normalizer_i = LD_linalg::norm_l2(sol_.JFs[idep],ndim);
      const double normalizer_i = LD_linalg::norm_l2(Gmat_i_[idep],ndof_tun);
      for (size_t icol = 0; icol < ndof_tun; icol++)
        Gmat_i_[idep][icol] /= normalizer_i;
    }
  }

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_G_rows(chunk_,sol_,Amat_i_);}

    inline void fill_G_rows(partial_chunk &chunk_, ode_solution &sol_, double **Gmat_i_)
    {
      int i_dof = 0;
      double  ** const JFs = sol_.JFs,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
          fill_x_G_columns(i_dof++,Gmat_i_,JFs,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L]);

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
            fill_u_G_columns(i_dof++,idep,Gmat_i_,JFs,Lambda_xu_mat[idep+1][i_L],Jac_utheta_vdxu_mat[i_L]);

      if (normalization_flag) normalize_submat_rows(sol_,Gmat_i_);
    }
    inline void fill_x_G_columns(int i_dof,double **Gmat_,double **JFs_,double &lambda_,double **Jac_xtheta_i_vdxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++) // tangent space constraint
      {
        Gmat_[idep][i_dof] = JFs_[idep][0]*lambda_;
        for (size_t k = 1, idim = nvar; k <= eor; k++)
          for (size_t jdep = 0; jdep < ndep; jdep++, idim++)
            Gmat_[idep][i_dof] += JFs_[idep][idim]*Jac_xtheta_i_vdxu_[k-1][jdep];
      }
    }
    inline void fill_u_G_columns(int i_dof,int idep,double **Gmat_,double **JFs_,double &lambda_,double *del_utheta_i)
    {
      for (size_t jdep = 0; jdep < ndep; jdep++) // tangent space contraint
      {
        Gmat_[jdep][i_dof] = JFs_[jdep][1+idep]*lambda_;
        for (size_t k = 1, idim = nvar+idep; k <= eor; k++, idim+=ndep)
          Gmat_[jdep][i_dof] += JFs_[jdep][idim]*del_utheta_i[k-1];
      }
    }
};

#endif
