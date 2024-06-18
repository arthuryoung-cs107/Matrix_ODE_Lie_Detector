#ifndef LD_MATS_HH
#define LD_MATS_HH

#include "LD_framework.hh"

struct LD_L_matrix: public LD_matrix
{
  LD_L_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,1,fspc_.perm_len) {}
  ~LD_L_matrix() {}

  template <class TBSIS> void populate_L_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

  protected:

    void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_)
      {fill_L_rows(chunk_,sol_,Amat_i_);}

    inline void fill_L_rows(partial_chunk &chunk_, ode_solution &sol_, double **Lmat_i_)
    {
      double  * const Lrow_i = Lmat_i_[0],
              * const lambda_x_vec = chunk_.Jac_mat[0];
      for (size_t i_L = 0; i_L < perm_len; i_L++) Lrow_i[i_L] = (dof_tun_flags[i_L])?(lambda_x_vec[i_L]):(0.0);

      if (normalization_flag)
      {
        const double L_row_i_mag = LD_linalg::norm_l2(Lrow_i,perm_len);
        for (size_t i_L = 0; i_L < perm_len; i_L++) Lrow_i[i_L] /= L_row_i_mag;
      }
    }
};

struct LD_R_matrix: public LD_matrix
{
  LD_R_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep*fspc_.eor) {}
  ~LD_R_matrix() {}

  template <class TBSIS> void populate_R_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

      if (normalization_flag)
      {
        for (size_t k = 0, idxu = 0; k < eor; k++)
          for (size_t idep = 0; idep < ndep; idep++, idxu++)
            for (size_t icol = 0; icol < ndof_tun; icol++)
              Rmat_i_[idxu][icol] /= sol_.dxu[idxu];
      }
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

  template <class TBSIS> void populate_O_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

      if (normalization_flag)
      {
        for (size_t idep = 0; idep < ndep; idep++)
          for (size_t icol = 0; icol < ndof_tun; icol++)
            Omat_i_[idep][icol] /= sol_.dnxu[idep];
      }
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

  template <class TBSIS> void populate_P_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

      if (normalization_flag)
      {
        for (size_t idep = 0; idep < ndep; idep++)
          for (size_t icol = 0; icol < ndof_tun; icol++)
            Pmat_i_[idep][icol] /= sol_.dnp1xu[idep];
      }
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

  template <class TBSIS> void populate_Q_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

      if (normalization_flag)
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

  template <class TBSIS> void populate_G_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
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

      if (normalization_flag)
      {
        for (size_t idep = 0; idep < ndep; idep++)
        {
          const double JFs_i_mag = sqrt(sol_.JFs_i_mag2(idep));
          for (size_t icol = 0; icol < ndof_tun; icol++)
            Gmat_i_[idep][icol] /= JFs_i_mag;
        }
      }
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
