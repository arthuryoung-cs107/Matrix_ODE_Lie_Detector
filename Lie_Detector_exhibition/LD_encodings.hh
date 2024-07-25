#ifndef LD_ENCODE_HH
#define LD_ENCODE_HH

// #include "LD_aux.hh"
// #include "LD_function_space.hh"
#include "LD_parameter_space.hh"

struct vspace_evaluation_package
{
  vspace_evaluation_package(ode_solspc_subset &Sset_,int ncon_,int nvec_,int nsol_,double tol_): Sset(Sset_),
    ncon(ncon_), nvec_max(nvec_), nsol_max(nsol_),
    sat_flags_mat(Tmatrix<bool>(nsol_,nvec_)),
    tol(tol_), nvec_evl(nvec_), nsol_evl(nsol_) {}
  vspace_evaluation_package(vspace_evaluation_package &evl_,int nsol_): vspace_evaluation_package(evl_.Sset,evl_.ncon,evl_.nvec_evl,nsol_,evl_.tol) {}
  ~vspace_evaluation_package() {free_Tmatrix<bool>(sat_flags_mat);}

  ode_solspc_subset &Sset;
  const int ncon,
            nvec_max,
            nsol_max;
  bool ** const sat_flags_mat,
        * const sat_flags = sat_flags_mat[0];
  int nvec_evl,
      nsol_evl;
  double  tol,
          ** Vmat;

  inline ode_solution ** sols_iset(int i_) {return Sset.get_sol_subset_i(i_);}
  inline int nsol_iset(int i_) {return Sset.nobs_subset_i(i_);}
  inline int max_nsol_subset() {return Sset.max_nobs_subset();}

  protected:

    static double R_err(double vdkm1xu_, double vx_, double dkxu_)
    {
      return 1.0 - (vdkm1xu_/(vx_*dkxu_)); // relative error
      // return dkxu_ - (vdkm1xu_/vx_); // absolute error
    }

};

struct LD_encoder
{
  LD_encoder(int ncod_,ode_solspc_meta &meta_): ncod(ncod_), meta(meta_) {}
  ~LD_encoder() {}

  ode_solspc_meta &meta;
  const int ncod;

  template <class EVL, class BSE> static void leniently_evaluate_vspace(LD_vspace_record &reco_,LD_vspace_record &reci_,EVL evl_,BSE **bases_,bool verbose_)
  {
    const int nset = LD_linalg::min_T<int>(reco_.nspc,reci_.nspc),
              nvec = LD_linalg::min_T<int>(reco_.nvec,reci_.nvec);
    bool sat_setmat[nset][nvec];
    int nsol_acc = 0,
        nvec_acc = 0,
        * const nVevl = reci_.nV_spcvec,
        * const nVsat = reco_.nV_spcvec,
        ** const iVsat = reco_.iV_spcmat;
    double  *** const Vtns = reci_.Vtns_data,
            t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:nsol_acc,nvec_acc)
    {
      BSE &bse_t = *(bases_[LD_threads::thread_id()]);
      EVL evl_t(evl_,evl_.max_nsol_subset());
      bool ** const Sat_t = evl_t.sat_flags_mat;
      int &nvec_evl = evl_t.nvec_evl,
          &nsol_evl = evl_t.nsol_evl;
      #pragma omp for
      for (size_t iset = 0; iset < nset; iset++)
      {
        nvec_acc += nvec_evl = nVevl[iset];
        nsol_acc += nsol_evl = evl_t.nsol_iset(iset);
        evl_t.Vmat = Vtns[iset];

        ode_solution ** const sols_i = evl_t.sols_iset(iset);
        LD_linalg::fill_vec<bool>(sat_setmat[iset],nvec_evl,false);

        int nsat_seti = 0;
        for (size_t jsol = 0; jsol < nsol_evl; jsol++)
        {
          if (evl_t.nsat_eval_condition(Sat_t[jsol],*(sols_i[jsol]),bse_t) != nvec_evl)
          {
            nsat_seti = 0;
            for (size_t iV = 0; iV < nvec_evl; iV++)
              nsat_seti += (int)(sat_setmat[iset][iV] = (sat_setmat[iset][iV]) || (Sat_t[jsol][iV]));
          }
          else LD_linalg::fill_vec<bool>(sat_setmat[iset],nsat_seti=nvec_evl,true);
          if (nsat_seti==nvec_evl) break; // if all pass
        }
        if (nVsat[iset] = nsat_seti)
        {
          for (size_t iV = 0, isat = 0; iV < nvec_evl; iV++)
            if (sat_setmat[iset][iV]) iVsat[iset][isat++] = iV;
        }
        else LD_linalg::fill_vec<int>(iVsat[iset],nvec_evl,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_encoder::leniently_evaluate_vspace) evaluated %d conds. (%d sols) over %d bases (%.1f x %d, on avg.) in %.4f seconds (%d threads)\n",
        nsol_acc*evl_.ncon,nsol_acc,nset,((double)nvec_acc)/((double)nset),reci_.vlen,
        work_time,LD_threads::numthreads());
    }
  }

  template <class EVL, class BSE> static void strictly_evaluate_vspace(LD_vspace_record &reco_,LD_vspace_record &reci_,EVL evl_,BSE **bases_,bool verbose_)
  {
    const int nset = LD_linalg::min_T<int>(reco_.nspc,reci_.nspc),
              nvec = LD_linalg::min_T<int>(reco_.nvec,reci_.nvec);
    bool sat_setmat[nset][nvec];
    int nsol_acc = 0,
        nvec_acc = 0,
        * const nVevl = reci_.nV_spcvec,
        * const nVsat = reco_.nV_spcvec,
        ** const iVsat = reco_.iV_spcmat;
    double  *** const Vtns = reci_.Vtns_data,
            t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:nsol_acc,nvec_acc)
    {
      BSE &bse_t = *(bases_[LD_threads::thread_id()]);
      EVL evl_t(evl_,1);
      bool  * const sat_t = evl_t.sat_flags;
      int &nvec_evl = evl_t.nvec_evl,
          &nsol_evl = evl_t.nsol_evl;
      #pragma omp for
      for (size_t iset = 0; iset < nset; iset++)
      {
        nvec_acc += nvec_evl = nVevl[iset];
        nsol_acc += nsol_evl = evl_t.nsol_iset(iset);
        evl_t.Vmat = Vtns[iset];

        ode_solution ** const sols_i = evl_t.sols_iset(iset);
        LD_linalg::fill_vec<bool>(sat_setmat[iset],nvec_evl,true);

        int nsat_seti = 0;
        for (size_t jsol = 0; jsol < nsol_evl; jsol++)
        {
          if (evl_t.nsat_eval_condition(sat_t,*(sols_i[jsol]),bse_t))
          {
            nsat_seti = 0;
            for (size_t iV = 0; iV < nvec_evl; iV++)
              nsat_seti += (int)(sat_setmat[iset][iV] = (sat_setmat[iset][iV]) && (sat_t[iV]));
          }
          else LD_linalg::fill_vec<bool>(sat_setmat[iset],nvec_evl,nsat_seti = 0);
          if (!nsat_seti) break; // if none pass
        }
        if (nVsat[iset] = nsat_seti)
        {
          for (size_t iV = 0, isat = 0; iV < nvec_evl; iV++)
            if (sat_setmat[iset][iV]) iVsat[iset][isat++] = iV;
        }
        else LD_linalg::fill_vec<int>(iVsat[iset],nvec_evl,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_encoder::strictly_evaluate_vspace) evaluated %d conds. (%d sols) over %d bases (%.1f x %d, on avg.) in %.4f seconds (%d threads)\n",
        nsol_acc*evl_.ncon,nsol_acc,nset,((double)nvec_acc)/((double)nset),reci_.vlen,
        work_time,LD_threads::numthreads());
    }
  }
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
    LD_encoded_matrix(LD_encoded_matrix &Acode_,LD_encoded_matrix &Bcode_): LD_encoded_matrix(Acode_.ncol,Acode_.nobs,Acode_.ncod + Bcode_.ncod)
    {
      for (size_t iobs = 0, irowA = 0, irowB = 0; iobs < nobs_full; iobs++)
      {
        for (size_t icod = 0; icod < Acode_.ncod; icod++, irowA++)
          for (size_t icol = 0; icol < Acode_.ncol; icol++)
            Atns_data[iobs][icod][icol] = Acode_.Amat[irowA][icol];
        for (size_t icod = 0, iicod = Acode_.ncod; icod < Bcode_.ncod; icod++, iicod++, irowB++)
          for (size_t icol = 0; icol < Bcode_.ncol; icol++)
            Atns_data[iobs][iicod][icol] = Bcode_.Amat[irowB][icol];
      }
    }
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
    LD_encoding_bundle(LD_encoding_bundle &Abndle_, LD_encoding_bundle &Bbndle_): data_owner(true),
      nset(Abndle_.nset), ncol_full(Abndle_.ncol_full), Amats(new double**[nset]), Acodes(new LD_encoded_matrix*[nset])
      {for (size_t i = 0; i < nset; i++) Amats[i] = ( Acodes[i] = new LD_encoded_matrix(*(Abndle_.Acodes[i]),*(Bbndle_.Acodes[i])) )->Amat;}
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

    inline int nrows_mat_i(int i_) {return Acodes[i_]->nrow;}
};

/*
  encoders
*/

struct LD_L_encoder: public LD_encoder
{
  LD_L_encoder(ode_solspc_meta &meta_): LD_encoder(1,meta_) {}
  ~LD_L_encoder() {}

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
      // const double normalizer_i = LD_linalg::norm_l2(sol_.JFs[idep],sol_.ndim);
      const double normalizer_i = LD_linalg::norm_l2(Grows_[idep],ndof_);
      for (size_t icol = 0; icol < ndof_; icol++)
        Grows_[idep][icol] /= normalizer_i;
    }
  }
};

struct LD_R_encoder: public LD_encoder
{
  LD_R_encoder(ode_solspc_meta &meta_,int nor_=0): LD_encoder(meta_.ndep*((nor_)?(nor_):(meta_.eor)),meta_) {}
  ~LD_R_encoder() {}

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
          LD_R_encoder::encode_R2n_rows(Rn_rows_i+ndep,chunk_t,sol_i,dof_tun_flags_mat,eorm1);
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

  /*
    normalization techniques
  */
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
  static void normalize_Rk_rows(double **Rkrows_,ode_solution &sol_,int kor_,int ndof_)
  {
    if (kor_==(sol_.eor+1)) LD_R_encoder::normalize_P_rows(Rkrows_,sol_,ndof_);
    else
      for (size_t idep = 0, idxu = kor_*sol_.ndep; idep < sol_.ndep; idep++, idxu++)
        // for (size_t icol = 0; icol < ndof_; icol++)
        //   Rmat_i_[idxu][icol] /= sol_.dxu[idxu];
      {
        for (size_t icol = 0; icol < ndof_; icol++) Rkrows_[idep][icol] /= sol_.u[idxu];
        LD_linalg::normalize_vec_l2(Rkrows_[idep],ndof_);
      }
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

};

/*
  evaluation packages
*/

struct L_vspace_eval: public vspace_evaluation_package
{
  L_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-10): vspace_evaluation_package(Sset_,1,nvec_,1,tol_) {}
  L_vspace_eval(L_vspace_eval &evl_,int nsol_): vspace_evaluation_package(evl_,nsol_) {}
  ~L_vspace_eval() {}

  template <class BSE> static void evaluate_lambda_signal_strength(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_=1e-10,bool verbose_=true)
  {
    LD_encoder::leniently_evaluate_vspace<L_vspace_eval,BSE>(reco_,reci_,
      L_vspace_eval(Sset_,reco_.nvec,tol_),
      bases_,verbose_);
  }

  virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbasis_) // evaluate signal strength
  {
    const int vlen = fbasis_.perm_len;
    double * const lamvec = fbasis_.Jac_mat[0];

    fbasis_.fill_partial_chunk(sol_.pts,0); // LD_linalg::normalize_vec_l2(lamvec,vlen);

    int n_success = 0;
    for (size_t iV = 0; iV < nvec_evl; iV++)
    {
      double lam_i = 0.0;
      for (size_t iL = 0; iL < vlen; iL++) lam_i += Vmat[iV][iL]*lamvec[iL];
      n_success += (int)(sat_flags_[iV] = (fabs(lam_i) > tol));
    }
    return n_success;
  }
};

struct G_vspace_eval: public vspace_evaluation_package
{
  G_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-10): vspace_evaluation_package(Sset_,Sset_.ndep,nvec_,1,tol_), mags_JFs(new double[Sset.ndep]) {}
  G_vspace_eval(G_vspace_eval &evl_,int nsol_):
    vspace_evaluation_package(evl_,nsol_), mags_JFs(new double[Sset.ndep]) {}
  ~G_vspace_eval() {delete [] mags_JFs;}

  // infinitesimal criterion
  virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbasis_)
  {
    const int ndep = sol_.ndep,
              ndim = sol_.ndim;
    double ** const JFs = sol_.JFs;

    for (size_t idep = 0; idep < ndep; idep++)
    {
      double acc = 0.0;
      for (size_t idim = 0; idim < ndim; idim++) acc += JFs[idep][idim]*JFs[idep][idim];
      mags_JFs[idep] = sqrt(acc);
    }

    fbasis_.fill_partial_chunk(sol_.pts);

    int n_success = 0;
    for (size_t iV = 0; iV < nvec_evl; iV++)
    {
      double * const  viV = fbasis_.v_eval(Vmat[iV]),
                      acc = 0.0;
      for (size_t idim = 0; idim < ndim; idim++) acc += viV[idim]*viV[idim];
      double  mag_v = sqrt(acc);

      for (size_t idep = 0; idep < ndep; idep++)
      {
        acc = 0.0;
        for (size_t idim = 0; idim < ndim; idim++) acc += JFs[idep][idim]*viV[idim];
        if (!( sat_flags_[iV] = (fabs(acc / (mags_JFs[idep]*mag_v)) < tol) )) break;
      }
      if (sat_flags_[iV]) n_success++;
    }
    return n_success;
  }

  protected:

    double * const mags_JFs;
};

struct Rk_vspace_eval: public vspace_evaluation_package
{

  Rk_vspace_eval(ode_solspc_subset &Sset_,int nvec_,int kor_,double tol_=1e-10): vspace_evaluation_package(Sset_,Sset_.ndep,nvec_,1,tol_), kor(kor_) {}
  Rk_vspace_eval(Rk_vspace_eval &evl_,int nsol_):
    vspace_evaluation_package(evl_,nsol_), kor(evl_.kor) {}
  ~Rk_vspace_eval() {}

  virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbasis_) // k'th ratio condition
  {
    const int ndep = sol_.ndep;
    double * const dkxu = (kor <= sol_.eor)?(sol_.pts + idkxu):(sol_.dnp1xu);

    fbasis_.fill_partial_chunk(sol_.pts,korm1);

    int n_success = 0;
    for (size_t iV = 0; iV < nvec_evl; iV++)
    {
      double  * const viV = fbasis_.v_eval(Vmat[iV],korm1),
              * const vdkm1xu = viV + idkm1xu,
              &vx = viV[0];
      for (size_t idep = 0; idep < ndep; idep++)
        if (!( sat_flags_[iV] = (fabs(vspace_evaluation_package::R_err(vdkm1xu[idep],vx,dkxu[idep])) < tol) )) break;
      if (sat_flags_[iV]) n_success++;
    }
    return n_success;
  }

  protected:

    const int kor,
              korm1 = kor-1,
              idkxu = 1 + Sset.ndep*kor,
              idkm1xu = 1 + Sset.ndep*korm1;
};

struct Rn_vspace_eval: public vspace_evaluation_package
{
  Rn_vspace_eval(ode_solspc_subset &Sset_,int nvec_,int nor_,double tol_=1e-10):
    vspace_evaluation_package(Sset_,Sset_.ndep*nor_,nvec_,1,tol_), nor(nor_) {}
  Rn_vspace_eval(Rn_vspace_eval &evl_, int nsol_):
    vspace_evaluation_package(evl_,nsol_), nor(evl_.nor) {}
  ~Rn_vspace_eval() {}

  template <class BSE> static void evaluate_nth_ratio_condition(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_=1e-6,int nor_=0,bool verbose_=true)
  {
    LD_encoder::strictly_evaluate_vspace<Rn_vspace_eval,BSE>(reco_,reci_,
      Rn_vspace_eval(Sset_,reco_.nvec,(nor_)?(nor_):(Sset_.eor),tol_),
      bases_,verbose_);
  }

  // 1 - n'th ratio condition
  virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbasis_) // k'th ratio condition
  {
    const int eor = sol_.eor,
              ndep = sol_.ndep,
              ndxu = ndep*((nor<eor)?(nor):(eor));
    double  * const u = sol_.u,
            * const dxu = sol_.dxu;

    fbasis_.fill_partial_chunk(sol_.pts,norm1);

    int n_success = 0;
    for (size_t iV = 0; iV < nvec_evl; iV++)
    {
      double  * const viV = fbasis_.v_eval(Vmat[iV],norm1),
              * const vu = viV + 1,
              &vx = viV[0];
      for (size_t idxu = 0; idxu < ndxu; idxu++)
        if (!( sat_flags_[iV] = (fabs(vspace_evaluation_package::R_err(vu[idxu],vx,dxu[idxu])) < tol) )) break;
      if (sat_flags_[iV])
      {
        if (nor==(eor+1))
        {
          double * const dnp1xu = sol_.dnp1xu;
          for (size_t idep = 0; idep < ndep; idep++)
            if (!( sat_flags_[iV] = (fabs(vspace_evaluation_package::R_err(vu[idep+ndxu],vx,dnp1xu[idep])) < tol) )) break;
          if (sat_flags_[iV]) n_success++;
        }
        else n_success++;
      }
    }
    return n_success;
  }

  protected:

    const int nor,
              norm1 = nor-1;
};

#endif
