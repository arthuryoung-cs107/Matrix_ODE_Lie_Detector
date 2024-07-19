#ifndef LD_ENCODE_HH
#define LD_ENCODE_HH

#include "LD_aux.hh"
#include "LD_function_space.hh"

struct LD_encoder
{
  LD_encoder(int ncod_): ncod(ncod_) {}
  ~LD_encoder() {}

  const int ncod;

};

struct LD_matrix_record
{
  LD_matrix_record(int ncol_, int nrow_, double **Amat_data_, bool init_=true):
    ncol(ncol_), nrow(nrow_), flgs_A(new bool[nrow_]), inds_A(new int[nrow_]), Amat_data(Amat_data_)
    { if (init_) init_default(); }
  ~LD_matrix_record() {delete [] flgs_A; delete [] inds_A;}

  const int ncol,
            nrow;
  bool * const flgs_A;
  int * const inds_A;
  double ** const Amat_data;

  inline void init_default() {for (size_t i = 0; i < nrow; i++) flgs_A[ inds_A[i] = (int) i ] = true;}

  inline void copy_record(LD_matrix_record &mrec_)
  {
    const int nrow_min = (nrow<mrec_.nrow)?(nrow):(mrec_.nrow);
    for (int i = 0, ii = 0; i < nrow_min; ii += (int)(mrec_.flgs_A[i++]))
      if ( flgs_A[i] = mrec_.flgs_A[i] ) inds_A[ii] = mrec_.flgs_A[ii];
  }
};

class LD_encoded_matrix
{
  const bool data_owner;

  protected:

    const int ncol_full,
              nobs_full,
              ncod_full,
              nrow_full;

    double *** const Atns_data,
            ** const Amat_data = Atns_data[0];

    LD_matrix_record * const mrec_ptr;

  public:

    LD_encoded_matrix(int ncol_,int nobs_,int ncod_=1): data_owner(true),
      ncol_full(ncol_), nobs_full(nobs_), ncod_full(ncod_), nrow_full(nobs_*ncod_),
      Atns_data(T3tensor<double>(nobs_,ncod_,ncol_)),
      mrec_ptr(new LD_matrix_record(ncol_,nrow_full,Amat_data)),
      ncol(ncol_full), nobs(nobs_full), ncod(ncod_full), nrow(nrow_full),
      Amat(new double*[nrow_full])
      {for (size_t i = 0; i < nrow_full; i++) Amat[i] = Amat_data[i];}
    ~LD_encoded_matrix()
    {
      if (data_owner)
      {
        free_T3tensor<double>(Atns_data);
        delete [] Amat;
        delete mrec_ptr;
      }
    }

    LD_matrix_record &mrec = *(mrec_ptr);
    int ncol,
        nobs,
        ncod,
        nrow,
        * const inds_A = mrec.inds_A;
    double  ** const Amat;

    inline double ** get_submat_i(int i_) {return Amat + i_*ncod;}

    inline int verify_nrow(int ncod_=0,int nobs_=0) {return nrow = ( (ncod_>0)?(ncod=ncod_):(ncod) )*( (nobs_>0)?(nobs=nobs_):(nobs) );}
};

class LD_encoding_bundle
{
  const bool data_owner;

  public:

    LD_encoding_bundle(int nset_,int ncol_,int *nobs_per_set_,int ncod_=1): data_owner(true),
      nset(nset_), ncol_full(ncol_), Amats(new double**[nset_]), Acodes(new LD_encoded_matrix*[nset_])
      {for (size_t i = 0; i < nset; i++) Amats[i] = ( Acodes[i] = new LD_encoded_matrix(ncol_,nobs_per_set_[i],ncod_) )->Amat;}
    LD_encoding_bundle(int nset_,int ncol_,int nobs_per_set_,int ncod_=1): data_owner(true),
      nset(nset_), ncol_full(ncol_), Amats(new double**[nset_]), Acodes(new LD_encoded_matrix*[nset_])
      {for (size_t i = 0; i < nset; i++) Amats[i] = ( Acodes[i] = new LD_encoded_matrix(ncol_,nobs_per_set_,ncod_) )->Amat;}
    ~LD_encoding_bundle()
    {
      if (data_owner)
      {
        for (size_t i = 0; i < nset; i++) delete Acodes[i];
        delete [] Acodes;
        delete [] Amats;
      }
    }

    const int nset,
              ncol_full;
    double *** const Amats;
    LD_encoded_matrix ** const Acodes;

    inline int nobs()
    {
      int acc = 0;
      for (size_t i = 0; i < nset; i++) acc += Acodes[i]->nobs;
      return acc;
    }
    inline int verify_nrow(int ncod_=0,int nobs_=0)
    {
      int acc = 0;
      for (size_t i = 0; i < nset; i++) acc += Acodes[i]->verify_nrow(ncod_,nobs_);
      return acc;
    }
    inline double ** get_iobs_encoding(int i_)
    {
      int iiAcode = 0,
          iiobs = i_;
      for (int ii = 0, idec = i_; ii < nset; ii++)
        if ((idec -= Acodes[iiAcode = ii]->nobs)>=0) iiobs = idec;
        else break;
      return Acodes[iiAcode]->get_submat_i(iiobs);
    }
    inline int max_nrows()
    {
      int mrows = 0;
      for (size_t i = 0; i < nset; i++) if (mrows < Acodes[i]->nrow) mrows = Acodes[i]->nrow;
      return mrows;
    }
    inline int min_nrows()
    {
      int mrows = Acodes[0]->nrow;
      for (size_t i = 1; i < nset; i++) if (mrows > Acodes[i]->nrow) mrows = Acodes[i]->nrow;
      return mrows;
    }
};

/*
  encoders
*/

struct LD_L_encoder: public LD_encoder
{
  LD_L_encoder(): LD_encoder(1) {}
  ~LD_L_encoder() {}

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
  static void encode_L_row(double *Lrow_,double *lambda_x_vec_,bool *dof_tun_flags_,int perm_len_)
    {for (size_t i_L = 0; i_L < perm_len_; i_L++) Lrow_[i_L] = (dof_tun_flags_[i_L])?(lambda_x_vec_[i_L]):(0.0);}
  static void normalize_L_row(double *Lrow_,int perm_len_)
  {
    const double L_row_i_mag = LD_linalg::norm_l2(Lrow_,perm_len_);
    for (size_t i_L = 0; i_L < perm_len_; i_L++) Lrow_[i_L] /= L_row_i_mag;
  }
};

struct LD_G_encoder: public LD_encoder
{
  LD_G_encoder(int ndep_): LD_encoder(ndep_) {}
  LD_G_encoder(ode_solspc_meta &meta_): LD_G_encoder(meta_.ndep) {}
  ~LD_G_encoder() {}

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

  static void normalize_G_rows(double **Grows_,ode_solution &sol_,int ndof_)
  {
    for (size_t idep = 0; idep < sol_.ndep; idep++)
    {
      // const double normalizer_i = LD_linalg::norm_l2(sol_.JFs[idep],sol_.ndim);
      const double normalizer_i = LD_linalg::norm_l2(Grows_[idep],ndof_);
      for (size_t icol = 0; icol < ndof_; icol++)
        Grows_[idep][icol] /= normalizer_i;
    }
  }
};

struct LD_R_encoder: public LD_encoder
{
  LD_R_encoder(int ndep_,int kor_): LD_encoder(ndep_*kor_) {}
  LD_R_encoder(ode_solspc_meta &meta_,int kor_=0): LD_R_encoder(meta_.ndep,(kor_)?(kor_):(meta_.eor)) {}
  ~LD_R_encoder() {}

  template <class BSIS> static void encode_Rn_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,int eor_=0,bool normalize_=false)
  {
    const int eorm1 = ((eor_)&&(eor_<=(set_.eor)))?(eor_-1):(set_.eor-1),
              ndep = set_.ndep,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel num_threads(1)
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      const int perm_len = base_t.perm_len;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      double ** const Lambda_xu_mat_t = chunk_t.Jac_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,eorm1);
        double ** const Rn_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_R1_rows(Rn_rows_i,Lambda_xu_mat_t,sol_i,dof_tun_flags_mat,perm_len);
        LD_R_encoder::encode_R2n_rows(Rn_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat,eorm1);
        if (normalize_) LD_R_encoder::normalize_Rn_rows(Rn_rows_i,sol_i,eorm1,base_t.ndof_tun);
      }
    }
  }

  template <class BSIS> static void encode_Q_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,int eor_=0,bool normalize_=false)
  {
    const int eor = set_.eor,
              eorm1 = eor-1,
              ndep = set_.ndep,
              i_eorp1 = ndep*eor,
              nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    #pragma omp parallel num_threads(1)
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      partial_chunk &chunk_t = base_t.partials;
      const int perm_len = base_t.perm_len;
      bool ** const dof_tun_flags_mat = base_t.dof_tun_flags_mat;
      double ** const Lambda_xu_mat_t = chunk_t.Jac_mat;
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        ode_solution &sol_i = *(sols[i]);
        base_t.fill_partial_chunk(sol_i.pts,eor);
        double ** const Q_rows_i = bndle_.get_iobs_encoding(i);
        LD_R_encoder::encode_R1_rows(Q_rows_i,Lambda_xu_mat_t,sol_i,dof_tun_flags_mat,perm_len);
        LD_R_encoder::encode_R2n_rows(Q_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat,eorm1);
        LD_R_encoder::encode_P_rows(Q_rows_i+i_eorp1,chunk_t,sol_i,dof_tun_flags_mat);
        if (normalize_)
        {
          LD_R_encoder::normalize_Rn_rows(Q_rows_i,sol_i,eorm1,base_t.ndof_tun);
          LD_R_encoder::normalize_P_rows(Q_rows_i+i_eorp1,sol_i,base_t.ndof_tun);
        }
      }
    }
  }

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
  static void encode_R2n_rows(double **R2rows_,partial_chunk &chunk_,ode_solution &sol_,bool **dof_tflags_mat_,int eorm1_)
  {
    const int ndep = sol_.ndep,
              nvar = sol_.nvar,
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

  static void normalize_Rn_rows(double **Rnrows_,ode_solution &sol_,int eorm1_,int ndof_)
  {
    for (size_t k = 0, idxu = 0; k <= eorm1_; k++)
      for (size_t idep = 0; idep < sol_.ndep; idep++, idxu++)
      // for (size_t icol = 0; icol < ndof_; icol++)
      //   Rmat_i_[idxu][icol] /= sol_.dxu[idxu];
      {
        for (size_t icol = 0; icol < ndof_; icol++) Rnrows_[idxu][icol] /= sol_.dxu[idxu];
        LD_linalg::normalize_vec_l2(Rnrows_[idxu],ndof_);
      }
  }
  static void normalize_P_rows(double **Prows_,ode_solution &sol_,int ndof_)
  {
    for (size_t idep = 0; idep < sol_.ndep; idep++)
      // for (size_t icol = 0; icol < ndof_; icol++)
        // Pmat_i_[idep][icol] /= sol_.dnp1xu[idep];
    {
      for (size_t icol = 0; icol < ndof_; icol++) Prows_[idep][icol] /= sol_.dnp1xu[idep];
      LD_linalg::normalize_vec_l2(Prows_[idep],ndof_);
    }
  }

};

#endif
