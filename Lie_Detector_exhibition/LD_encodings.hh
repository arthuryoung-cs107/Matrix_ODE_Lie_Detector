#ifndef LD_ENCODE_HH
#define LD_ENCODE_HH

#include "LD_parameter_space.hh"

struct LD_L_encoder: public LD_encoder
{
  LD_L_encoder(ode_solspc_meta &meta_): LD_encoder(1,meta_) {}
  ~LD_L_encoder() {}

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_)
  {
    fbse_.fill_partial_chunk(sol_i_.pts,0);
    LD_L_encoder::encode_L_row(rows_i_[0],fbse_.Jac_mat[0],fbse_.dof_tun_flags,fbse_.perm_len);
    if (normalize_) LD_L_encoder::normalize_L_row(rows_i_[0],fbse_.perm_len);
  }

  /*
    encoding routines
  */
  template <class BSIS> static void encode_L_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      const int perm_len = base_t.perm_len;
      bool * const dof_tun_flags = base_t.dof_tun_flags;
      double * const lambda_x_vec = base_t.Jac_mat[0];
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        base_t.fill_partial_chunk(sols[i]->pts,0);
        double * const L_row_i = *(bndle_.get_iobs_encoding(i));
        LD_L_encoder::encode_L_row(L_row_i,lambda_x_vec,dof_tun_flags,perm_len);
        if (normalize_) LD_L_encoder::normalize_L_row(L_row_i,perm_len);
      }
    }
  }

  /*
    encoding techniques
  */
  static void encode_L_row(double *Lrow_,double *lambda_x_vec_,bool *dof_tun_flags_,int perm_len_)
    {for (size_t i_L = 0; i_L < perm_len_; i_L++) Lrow_[i_L] = (dof_tun_flags_[i_L])?(lambda_x_vec_[i_L]):(0.0);}

  /*
    normalization techniques
  */
  static void normalize_L_row(double *Lrow_,int perm_len_)
  {
    const double L_row_i_mag = LD_linalg::norm_l2(Lrow_,perm_len_);
    for (size_t i_L = 0; i_L < perm_len_; i_L++) Lrow_[i_L] /= L_row_i_mag;
  }
};

struct LD_G_encoder: public LD_encoder
{
  LD_G_encoder(ode_solspc_meta &meta_): LD_encoder(meta_.ndep,meta_) {}
  ~LD_G_encoder() {}

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_)
  {
    fbse_.fill_partial_chunk(sol_i_.pts,sol_i_.eor);
    LD_G_encoder::encode_G_rows(rows_i_,fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat);
    if (normalize_) LD_G_encoder::normalize_G_rows(rows_i_,sol_i_,fbse_.ndof_tun);
  }

  /*
    encoding routines
  */
  template <class BSIS> static void encode_G_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        base_t.fill_partial_chunk(sols[i]->pts);
        double ** const G_rows_i = bndle_.get_iobs_encoding(i);
        LD_G_encoder::encode_G_rows(G_rows_i,chunk_t,*(sols[i]),dof_tun_flags_mat);
        if (normalize_) LD_G_encoder::normalize_G_rows(G_rows_i,*(sols[i]),base_t.ndof_tun);
      }
    }
  }

  /*
    encoding techniques
  */
  static void encode_G_rows(double **Grows_,partial_chunk &chunk_,ode_solution &sol_,bool **dof_tun_flags_mat_)
  {
    const int eor = sol_.eor,
              ndep = sol_.ndep,
              nvar = sol_.nvar,
              perm_len = chunk_.perm_len;
    bool * const dof_tun_flags_x = *(dof_tun_flags_mat_);
    int i_dof = 0;
    double  ** const JFs = sol_.JFs,
            * const lambda_x_vec = chunk_.Jac_mat[0],
            ** const Lambda_xu_mat = chunk_.Jac_mat,
            *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
            ** const Jac_utheta_vdxu_mat = chunk_.C_u;

    for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tun_flags_mat_[0][i_L++]))
      if (dof_tun_flags_mat_[0][i_L])
        for (size_t idep = 0; idep < ndep; idep++)
        {
          Grows_[idep][i_dof] = JFs[idep][0]*lambda_x_vec[i_L];
          for (size_t k = 0, idim = nvar; k < eor; k++)
            for (size_t jdep = 0; jdep < ndep; jdep++, idim++)
              Grows_[idep][i_dof] += JFs[idep][idim]*Jac_xtheta_vdxu_tns[i_L][k][jdep];
        }

    for (size_t idep = 1; idep <= ndep; idep++)
      for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tun_flags_mat_[idep][i_L++]))
        if (dof_tun_flags_mat_[idep][i_L])
          for (size_t jdep = 0; jdep < ndep; jdep++)
          {
            Grows_[jdep][i_dof] = JFs[jdep][idep]*Lambda_xu_mat[idep][i_L];
            for (size_t k = 0, idim = ndep+idep; k < eor; k++, idim+=ndep)
              Grows_[jdep][i_dof] += JFs[jdep][idim]*Jac_utheta_vdxu_mat[i_L][k];
          }
  }

  /*
    normalization techniques
  */
  static void normalize_G_rows(double **Grows_,ode_solution &sol_,int ndof_)
  {
    for (size_t idep = 0; idep < sol_.ndep; idep++)
    {
      const double normalizer_i = LD_linalg::norm_l2(sol_.JFs[idep],sol_.ndim);
      for (size_t icol = 0; icol < ndof_; icol++) Grows_[idep][icol] /= normalizer_i;
      // LD_linalg::normalize_vec_l2(Grows_[idep],ndof_);
    }
  }
};

struct LD_R_encoder: public LD_encoder
{
  LD_R_encoder(ode_solspc_meta &meta_,int nor_=0): LD_encoder(meta_.ndep*((nor_)?(nor_):(meta_.eor)),meta_) {}
  ~LD_R_encoder() {}

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_)
  {
    const int nor = ncod/sol_i_.ndep;
    fbse_.fill_partial_chunk(sol_i_.pts,nor-1);
    LD_R_encoder::encode_R1_rows(rows_i_,fbse_.Jac_mat,sol_i_,fbse_.dof_tun_flags_mat,fbse_.perm_len);
    if (nor>1)
    {
      const int eor = sol_i_.eor;
      LD_R_encoder::encode_R2n_rows(rows_i_+sol_i_.ndep,fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat,((nor>eor)?(eor):(nor))-1);
      if (nor==(eor+1)) LD_R_encoder::encode_P_rows(rows_i_+(eor*sol_i_.ndep),fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat);
    }
    if (normalize_)
    {
      LD_R_encoder::normalize_Rn_rows(rows_i_,sol_i_,nor-1,fbse_.ndof_tun);
      if (nor==(sol_i_.eor+1)) LD_R_encoder::normalize_P_rows(rows_i_+(sol_i_.eor*sol_i_.ndep),sol_i_,fbse_.ndof_tun);
    }
  }

  /*
    encoding routines
  */
  template <class BSIS> static void encode_R1_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int ndep = set_.ndep,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      const int perm_len = base_t.perm_len;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      double ** const Lambda_xu_mat_t = base_t.Jac_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,0);
        double ** const R1_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_R1_rows(R1_rows_i,Lambda_xu_mat_t,sol_i,dof_tun_flags_mat,perm_len);
        if (normalize_) LD_R_encoder::normalize_Rk_rows(R1_rows_i,sol_i,1,base_t.ndof_tun);
      }
    }
  }
  template <class BSIS> static void encode_P_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int eor = set_.eor,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,eor);
        double ** const P_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_P_rows(P_rows_i,chunk_t,sol_i,dof_tun_flags_mat);
        if (normalize_) LD_R_encoder::normalize_P_rows(P_rows_i,sol_i,base_t.ndof_tun);
      }
    }
  }
  template <class BSIS> static void encode_Rk_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,int kor_,bool normalize_=false)
  {
    if (kor_==1) LD_R_encoder::encode_R1_bundle<BSIS>(bndle_,set_,bases_,normalize_);
    else if (kor_==(set_.eor+1)) LD_R_encoder::encode_P_bundle<BSIS>(bndle_,set_,bases_,normalize_);
    else
    {
      const int korm1 = kor_-1,
                nobs_max = bndle_.nobs();
      ode_solution ** const sols = set_.sols;
      #pragma omp parallel
      {
        BSIS &base_t = *(bases_[LD_threads::thread_id()]);
        partial_chunk &chunk_t = base_t.partials;
        bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
        #pragma omp for
        for (size_t i = 0; i < nobs_max; i++)
        {
          ode_solution &sol_i = *(sols[i]);
          base_t.fill_partial_chunk(sol_i.pts,korm1);
          double ** const R_rows_i = bndle_.get_iobs_encoding(i);
          LD_R_encoder::encode_Rk_rows(R_rows_i,chunk_t,sol_i,dof_tun_flags_mat,kor_);
          if (normalize_) LD_R_encoder::normalize_Rk_rows(R_rows_i,sol_i,kor_,base_t.ndof_tun);
        }
      }
    }
  }
  template <class BSIS> static void encode_O_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    LD_R_encoder::encode_Rk_bundle<BSIS>(bndle_,set_,bases_,set_.eor,normalize_);
  }
  template <class BSIS> static void encode_Q_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int eor = set_.eor,
              eorm1 = eor-1,
              ndep = set_.ndep,
              i_eorp1 = ndep*eor,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      const int perm_len = base_t.perm_len;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      double ** const Lambda_xu_mat_t = base_t.Jac_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,eor);
        double ** const Q_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_R1_rows(Q_rows_i,Lambda_xu_mat_t,sol_i,dof_tun_flags_mat,perm_len);
        if (eorm1>0) LD_R_encoder::encode_R2n_rows(Q_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat,eorm1);
        LD_R_encoder::encode_P_rows(Q_rows_i+i_eorp1,chunk_t,sol_i,dof_tun_flags_mat);
        if (normalize_)
        {
          LD_R_encoder::normalize_Rn_rows(Q_rows_i,sol_i,eorm1,base_t.ndof_tun);
          LD_R_encoder::normalize_P_rows(Q_rows_i+i_eorp1,sol_i,base_t.ndof_tun);
        }
      }
    }
  }
  template <class BSIS> static void encode_Rn_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,int eor_=0,bool normalize_=false)
  {
    if (eor_==1) LD_R_encoder::encode_R1_bundle<BSIS>(bndle_,set_,bases_,normalize_);
    else if (eor_==(set_.eor+1)) LD_R_encoder::encode_Q_bundle(bndle_,set_,bases_,normalize_);
    else
    {
      const int eorm1 = ((eor_)&&(eor_<=(set_.eor)))?(eor_-1):(set_.eor-1),
                ndep = set_.ndep,
                nobs_max = bndle_.nobs();
      ode_solution ** const sols = set_.sols;
      #pragma omp parallel
      {
        BSIS &base_t = *(bases_[LD_threads::thread_id()]);
        partial_chunk &chunk_t = base_t.partials;
        const int perm_len = base_t.perm_len;
        bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
        double ** const Lambda_xu_mat_t = base_t.Jac_mat;
        #pragma omp for
        for (size_t i = 0; i < nobs_max; i++)
        {
          ode_solution &sol_i = *(sols[i]);
          base_t.fill_partial_chunk(sol_i.pts,eorm1);
          double ** const Rn_rows_i = bndle_.get_iobs_encoding(i);
          LD_R_encoder::encode_R1_rows(Rn_rows_i,Lambda_xu_mat_t,sol_i,dof_tun_flags_mat,perm_len);
          if (eorm1>0) LD_R_encoder::encode_R2n_rows(Rn_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat,eorm1);
          if (normalize_) LD_R_encoder::normalize_Rn_rows(Rn_rows_i,sol_i,eorm1,base_t.ndof_tun);
        }
      }
    }
  }

  /*
    encoding techniques
  */
  static void encode_R1_rows(double **R1rows_,double **L_xu_mat_,ode_solution &sol_,bool **dof_tflags_mat_,int perm_len_)
  {
    const int ndep = sol_.ndep;
    int i_dof = 0;
    double  * const d1xu = sol_.dxu;
    for (size_t i_L = 0; i_L < perm_len_; i_dof += (int)(dof_tflags_mat_[0][i_L++]))
      if (dof_tflags_mat_[0][i_L])
        for (size_t idep = 0; idep < ndep; idep++) R1rows_[idep][i_dof] = d1xu[idep]*L_xu_mat_[0][i_L];
    for (size_t idep = 0; idep < ndep; idep++)
      for (size_t i_L = 0; i_L < perm_len_; i_dof += (int)(dof_tflags_mat_[idep+1][i_L++]))
        if (dof_tflags_mat_[idep+1][i_L]) R1rows_[idep][i_dof] = -(L_xu_mat_[idep+1][i_L]);
  }
  static void encode_P_rows(double **Prows_,partial_chunk &chunk_,ode_solution &sol_,bool **dof_tflags_mat_)
  {
    const int eorm1 = sol_.eor-1,
              ndep = sol_.ndep,
              perm_len = chunk_.perm_len;
    int i_dof = 0;
    double  * const dnp1xu = sol_.dnp1xu,
            * const lambda_x_vec = chunk_.Jac_mat[0],
            *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
            ** const Jac_utheta_vdxu_mat = chunk_.C_u;

    for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[0][i_L++]))
      if (dof_tflags_mat_[0][i_L])
        for (size_t idep = 0; idep < ndep; idep++)
          Prows_[idep][i_dof] = dnp1xu[idep]*lambda_x_vec[i_L] - Jac_xtheta_vdxu_tns[i_L][eorm1][idep];

    for (size_t idep = 0; idep < ndep; idep++)
      for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[idep+1][i_L++]))
        if (dof_tflags_mat_[idep+1][i_L])
          Prows_[idep][i_dof] = -(Jac_utheta_vdxu_mat[i_L][eorm1]);
  }
  static void encode_Rk_rows(double **Rrows_,partial_chunk &chunk_,ode_solution &sol_,bool **dof_tflags_mat_,int kor_)
  {
    if (kor_==1) LD_R_encoder::encode_R1_rows(Rrows_,chunk_.Jac_mat,sol_,dof_tflags_mat_,chunk_.perm_len);
    else if (kor_==(sol_.eor+1)) LD_R_encoder::encode_P_rows(Rrows_,chunk_,sol_,dof_tflags_mat_);
    else
    {
      const int korm2 = kor_-2,
                ndep = sol_.ndep,
                perm_len = chunk_.perm_len;
      int i_dof = 0;
      double  * const dkxu = sol_.u + (kor_*ndep),
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[0][i_L++]))
        if (dof_tflags_mat_[0][i_L])
          for (size_t idep = 0; idep < ndep; idep++)
            Rrows_[idep][i_dof] = dkxu[idep]*Lambda_xu_mat[0][i_L] - Jac_xtheta_vdxu_tns[i_L][korm2][idep];
      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[idep+1][i_L++]))
          if (dof_tflags_mat_[idep+1][i_L])
            Rrows_[idep][i_dof] = -(Jac_utheta_vdxu_mat[i_L][korm2]);
    }
  }
  static void encode_R2n_rows(double **R2rows_,partial_chunk &chunk_,ode_solution &sol_,bool **dof_tflags_mat_,int eorm1_)
  {
    if (eorm1_>0)
    {
      const int ndep = sol_.ndep,
                perm_len = chunk_.perm_len;
      int i_dof = 0;
      double  * const d2xu = sol_.dxu+ndep,
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[0][i_L++]))
        if (dof_tflags_mat_[0][i_L])
          for (size_t k = 0, idim = 0; k < eorm1_; k++)
            for (size_t idep = 0; idep < ndep; idep++, idim++)
              R2rows_[idim][i_dof] = d2xu[idim]*Lambda_xu_mat[0][i_L] - Jac_xtheta_vdxu_tns[i_L][k][idep];

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_dof += (int)(dof_tflags_mat_[idep+1][i_L++]))
          if (dof_tflags_mat_[idep+1][i_L])
            for (size_t k = 0, i_R = idep; k < eorm1_; k++, i_R+=ndep)
              R2rows_[i_R][i_dof] = -(Jac_utheta_vdxu_mat[i_L][k]);
    }
  }

  /*
    normalization techniques
  */
  static void normalize_P_rows(double **Prows_,ode_solution &sol_,int ndof_)
  {
    for (size_t idep = 0; idep < sol_.ndep; idep++)
    {
      for (size_t icol = 0; icol < ndof_; icol++) Prows_[idep][icol] /= sol_.dnp1xu[idep];
      // LD_linalg::normalize_vec_l2(Prows_[idep],ndof_);
    }
  }
  static void normalize_Rk_rows(double **Rkrows_,ode_solution &sol_,int kor_,int ndof_)
  {
    if (kor_==(sol_.eor+1)) LD_R_encoder::normalize_P_rows(Rkrows_,sol_,ndof_);
    else
      for (size_t idep = 0, idxu = kor_*sol_.ndep; idep < sol_.ndep; idep++, idxu++)
        for (size_t icol = 0; icol < ndof_; icol++) Rkrows_[idep][icol] /= sol_.u[idxu];
        // LD_linalg::normalize_vec_l2(Rkrows_[idep],ndof_);
  }
  static void normalize_Rn_rows(double **Rnrows_,ode_solution &sol_,int eorm1_,int ndof_)
  {
    for (size_t k = 0, idxu = 0; k <= eorm1_; k++)
      for (size_t idep = 0; idep < sol_.ndep; idep++, idxu++)
        for (size_t icol = 0; icol < ndof_; icol++) Rnrows_[idxu][icol] /= sol_.dxu[idxu];
        // LD_linalg::normalize_vec_l2(Rnrows_[idxu],ndof_);
  }
};

struct LD_O_encoder: public LD_encoder
{
  LD_O_encoder(ode_solspc_meta &meta_): LD_encoder(meta_.ndep,meta_) {}
  ~LD_O_encoder() {}

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_)
  {
    fbse_.fill_partial_chunk(sol_i_.pts,sol_i_.eor-1);
    LD_R_encoder::encode_Rk_rows(rows_i_,fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat,sol_i_.eor);
    if (normalize_) LD_R_encoder::normalize_Rk_rows(rows_i_,sol_i_,sol_i_.eor,fbse_.ndof_tun);
  }
};

struct LD_OG_encoder: public LD_encoder
{
  LD_OG_encoder(ode_solspc_meta &meta_): LD_encoder(2*meta_.ndep,meta_) {}
  ~LD_OG_encoder() {}

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_)
  {
    fbse_.fill_partial_chunk(sol_i_.pts,sol_i_.eor);
    LD_R_encoder::encode_Rk_rows(rows_i_,fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat,sol_i_.eor);
    LD_G_encoder::encode_G_rows(rows_i_+sol_i_.ndep,fbse_.partials,sol_i_,fbse_.dof_tun_flags_mat);
    if (normalize_)
    {
      LD_R_encoder::normalize_Rk_rows(rows_i_,sol_i_,sol_i_.eor,fbse_.ndof_tun);
      LD_G_encoder::normalize_G_rows(rows_i_+sol_i_.ndep,sol_i_,fbse_.ndof_tun);
    }
  }

  template <class BSIS> static void encode_OG_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,bool normalize_=false)
  {
    const int eor = set_.eor,
              ndep = set_.ndep,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,eor);
        double ** const OG_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_Rk_rows(OG_rows_i,chunk_t,sol_i,dof_tun_flags_mat,eor);
        LD_G_encoder::encode_G_rows(OG_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat);
        if (normalize_)
        {
          LD_R_encoder::normalize_Rk_rows(OG_rows_i,sol_i,eor,base_t.ndof_tun);
          LD_G_encoder::normalize_G_rows(OG_rows_i+ndep,sol_i,base_t.ndof_tun);
        }
      }
    }
  }
};

#endif
