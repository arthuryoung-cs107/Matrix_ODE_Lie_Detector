#ifndef LD_PARSPC_HH
#define LD_PARSPC_HH

#include "LD_aux.hh"
#include "LD_function_space.hh"

struct LD_vspace_record
{
  LD_vspace_record(int nspc_,int nvec_, int vlen_, double ***Vtns_data_,bool init_=true):
    nspc(nspc_), nvec(nvec_), vlen(vlen_),
    nV_spcvec(new int[nspc_]), iV_spcmat(Tmatrix<int>(nspc_,nvec_)), Vtns_data(Vtns_data_)
    {if (init_) init_default();}
  LD_vspace_record(int nspc_,int nvec_, int vlen_, LD_vspace_record &rec_):
    nspc(nspc_), nvec(nvec_), vlen(vlen_),
    nV_spcvec(new int[nspc_]), iV_spcmat(Tmatrix<int>(nspc_,nvec_)), Vtns_data(rec_.Vtns_data)
    {copy_record(rec_);}
  LD_vspace_record(LD_vspace_record &rec_):
    nspc(rec_.nspc), nvec(rec_.nvec), vlen(rec_.vlen),
    nV_spcvec(new int[rec_.nspc]), iV_spcmat(Tmatrix<int>(rec_.nspc,rec_.nvec)), Vtns_data(rec_.Vtns_data)
    {copy_record(rec_);}
  LD_vspace_record(LD_vspace_record &rec1_,LD_vspace_record &rec2_): // combine records
    LD_vspace_record(LD_linalg::max_T<int>(rec1_.nspc,rec2_.nspc),
      LD_linalg::max_T<int>(rec1_.nvec,rec2_.nvec),LD_linalg::max_T<int>(rec1_.vlen,rec2_.vlen),
      rec1_.Vtns_data,false) {init_combined_records(rec1_,rec2_);}
  ~LD_vspace_record() { free_Tmatrix<int>(iV_spcmat); delete [] nV_spcvec; }

  const int nspc,
            nvec,
            vlen;
  int * const nV_spcvec,
      ** const iV_spcmat;
  double *** const Vtns_data;

  void init_combined_records(LD_vspace_record &rec1_,LD_vspace_record &rec2_);
  inline void init_default()
    {for (size_t i = 0; i < nspc; i++) LD_linalg::fill_vec_012<int>(iV_spcmat[i],nV_spcvec[i]=nvec);}
  inline void copy_record(LD_vspace_record &rec_)
  {
    LD_linalg::copy_vec<int>(nV_spcvec,rec_.nV_spcvec,nspc);
    for (size_t i = 0; i < nspc; i++)
      for (size_t j = 0; j < nV_spcvec[i]; j++) iV_spcmat[i][j] = rec_.iV_spcmat[i][j];
  }
  inline void set_record_nullspc(int *rvec_)
    {for (size_t i = 0; i < nspc; i++) LD_linalg::fill_vec_012<int>(iV_spcmat[i],nV_spcvec[i] -= rvec_[i],rvec_[i]);}
  inline void set_record_rankvec(int *rvec_)
    {for (size_t i = 0; i < nspc; i++) LD_linalg::fill_vec_012<int>(iV_spcmat[i],nV_spcvec[i] = rvec_[i]);}
  inline void set_record_unirank(int nvec_)
    {for (size_t i = 0; i < nspc; i++) LD_linalg::fill_vec_012<int>(iV_spcmat[i],nV_spcvec[i] = nvec_);}
  inline void set_record_uninull(int nvec_)
    {for (size_t i = 0; i < nspc; nV_spcvec[i++] = nvec_) LD_linalg::fill_vec_012<int>(iV_spcmat[i],nvec_,nV_spcvec[i]-nvec_);}

  void print_subspace_details(const char name_[], bool longwinded_=false, bool print_comp_=false);
  void print_subspace_details(LD_vspace_record &rec1_,const char name_[],bool longwinded_=false,bool print_comp_=false);
  void compare_subspaces(LD_vspace_record &rec1_,const char name1_[],LD_vspace_record &rec2_,const char name2_[]);
  void write_vspace_record(const char name_[],bool write_Vtns_data_=false);

  inline int comp_iV_spcmat_len()
  {
    int len_out = 0;
    for (size_t i = 0; i < nspc; i++) len_out += nV_spcvec[i];
    return len_out;
  }
};

class LD_vector_space
{

  const bool data_owner;
  double  ** const Vmat_data;

  public:

    LD_vector_space(int vlen_);
    LD_vector_space(int vlen_,int *inds_V_,double **Vmat_data_,double **Vmat_);
    LD_vector_space(LD_vector_space &Vspc_, bool deepcopy_);
    ~LD_vector_space();

    const int vlen_full;

    int nvec_use,
        vlen_use,
        * const inds_V;

    double ** const Vmat;

    // inline int set_Vspace(LD_vector_space &Vspc_,bool set_complements_=true)
    // {
    //   nvec_use = Vspc_.nvec_use; vlen_use = Vspc_.vlen_use;
    //
    //   double ** const Vdata_in = Vspc_.get_Vmat_data();
    //   for (size_t i = 0; i < nvec_use; i++)
    //     for (size_t j = 0; j < vlen_use; j++) Vmat_data[i][j] = Vdata_in[i][j];
    //   for (size_t i = 0; i < nvec_use; i++) Vmat[i] = Vmat_data[inds_V[i] = Vspc_.inds_V[i]];
    //   return nvec_use;
    // }
    inline void reset_Vspace()
    {
      nvec_use = vlen_use = vlen_full;
      for (size_t i = 0; i < vlen_full; i++) Vmat[i] = Vmat_data[inds_V[i] = i];
    }
    // inline int load_Vmat_data(double **Vmat_, bool fsat_[], int nvec_)
    // {
    //   int iv = 0;
    //   for (size_t i = 0, inv = 0; i < nvec_; i++)
    //     if (fsat_[i]) LD_linalg::copy_vec<double>(Vmat[iv++] = Vmat_data[inds_V[iv] = iv],Vmat_[i],vlen_use);
    //   return nvec_use = iv;
    // }
    inline void set_Vmat(int nvec_use_)
    {
      bool fevl[vlen_full]; for (size_t i = 0; i < vlen_full; i++) fevl[i] = false;
      for (size_t i = 0; i < nvec_use_; i++) fevl[inds_V[i]] = true;
      for (size_t i = 0, iV = 0, inV = nvec_use_; i < vlen_full; i++)
        if (fevl[i]) Vmat[iV++] = Vmat_data[inds_V[iV] = i]; // reorder if necessary
        else Vmat[inV++] = Vmat_data[inds_V[inV] = i]; // place extraneous indices outside
      nvec_use = nvec_use_;
    }
    inline bool check_default_configuration()
    {
      if (nvec_use==vlen_use)
      {
        for (size_t i = 0; i < nvec_use; i++) if (Vmat[i]!=Vmat_data[i]) return false;
        return true;
      }
      else return false;
    }
    inline void post_multiply_rowvec(double *AVrow_,double *Arow_)
    {
      for (size_t jvec = 0; jvec < nvec_use; jvec++)
      {
        double AVj_acc = 0.0;
        for (size_t iele = 0; iele < vlen_use; iele++) AVj_acc += Arow_[iele]*Vmat[jvec][iele];
        AVrow_[jvec] = AVj_acc;
      }
    }
    inline void pre_multiply_colvec(double *AVcol_,double *Acol_)
    {
      for (size_t iele = 0; iele < vlen_use; iele++) AVcol_[iele] = Acol_[0]*Vmat[0][iele];
      for (size_t jvec = 1; jvec < nvec_use; jvec++)
        for (size_t iele = 0; iele < vlen_use; iele++)
          AVcol_[iele] += Acol_[jvec]*Vmat[jvec][iele];
    }

    inline double * Vrowi_data(int i_) {return Vmat_data[i_];}
    inline double ** get_Vmat_data() {return Vmat_data;}
};

class LD_vector_bundle
{

  const bool data_owner;

  protected:

    double  *** const Vtns_data,
            ** const Vrows;

    LD_vspace_record * const rec_ptr;

  public:

    LD_vector_bundle(int nspc_,int vlen_);
    LD_vector_bundle(LD_vector_bundle &bndle_);
    LD_vector_bundle(int nspc_,int vlen_, double ***Vtns_in_): LD_vector_bundle(nspc_,vlen_) {load_Vtns_data(Vtns_in_);}
    LD_vector_bundle(LD_vspace_record &rec_): LD_vector_bundle(rec_.nspc,rec_.vlen,rec_.Vtns_data) {set_Vspaces(rec_);}
    ~LD_vector_bundle();

    inline void load_Vtns_data(double ***Vtns_,bool load_complements_=true)
    {
      if (load_complements_)
        for (size_t i = 0; i < nspc; i++)
          for (size_t j = 0; j < vlen_full; j++)
            for (size_t k = 0; k < vlen_full; k++)
              Vtns_data[i][j][k] = Vtns_[i][j][k];
      else
        for (size_t i = 0; i < nspc; i++)
          for (size_t j = 0; j < nV_spcvec[i]; j++)
            for (size_t k = 0; k < Vspaces[i]->vlen_use; k++)
              Vtns_data[i][j][k] = Vtns_[i][j][k];
    }

    LD_vspace_record &rec = *(rec_ptr);

    const int nspc,
              vlen_full;
    int *  const nV_spcvec = rec.nV_spcvec,
        ** const iV_spcmat = rec.iV_spcmat;
    double *** const Vtns;

    LD_vector_space ** const Vspaces;

    inline void reset_Vspaces()
    {
      for (size_t i = 0; i < nspc; i++) Vspaces[i]->reset_Vspace();
      LD_linalg::fill_vec<int>(nV_spcvec,nspc,vlen_full);
    }
    inline void set_Vspaces() {for (size_t i = 0; i < nspc; i++) Vspaces[i]->set_Vmat(nV_spcvec[i]);}
    inline void set_Vspaces(LD_vspace_record &rec_)
      {rec.copy_record(rec_); for (size_t i = 0; i < nspc; i++) Vspaces[i]->set_Vmat(nV_spcvec[i]);}
    inline bool check_default_configuration()
    {
      for (size_t i = 0; i < nspc; i++) if (!(Vspaces[i]->check_default_configuration())) return false;
      return true;
    }

    void write_LD_vector_bundle(const char name_[],bool write_Vtns_=false);

    inline int comp_Vtns_len()
    {
      int len_out = 0;
      for (size_t i = 0; i < nspc; i++) len_out += nV_spcvec[i]*(Vspaces[i]->vlen_use);
      return len_out;
    }
};

class LD_Theta_space: public ode_solspc_element
{
  const bool  vspc_owner,
              data_owner;

  const int perm_len;

  LD_vector_space * const spc_ptr,
                  ** const bse_ptrs;

  public:

    LD_Theta_space(ode_solspc_meta &meta_, int ndof_);
    LD_Theta_space(ode_solspc_meta &meta_, int ndof_, LD_vector_space *spc_ptr_, LD_vector_space **bse_ptrs_, int *pbse_);
    ~LD_Theta_space();

    LD_vector_space &spc = *(spc_ptr);

    int ndof_use,
        * const pbse_nvar_use;
    double  ** const Tmat = spc.Vmat;

    int init_Vspce_premult(double **Wmat_);

    inline int set_ivar_Yspace(LD_vector_space *Yspace_, int ivar_=0)
    {
      if ((0<=ivar_)&&(ivar_<nvar)) {bse_ptrs[ivar_] = Yspace_; return spc.nvec_use = verify_Yspaces(true);}
      else return spc.nvec_use = set_evry_Yspace(Yspace_);
    }
    inline int set_evry_Yspace(LD_vector_space *Yspace_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspace_;
      return spc.nvec_use = verify_Yspaces(true);
    }
    inline int set_evry_Yspace(LD_vector_space **Yspaces_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspaces_[ivar];
      return spc.nvec_use = verify_Yspaces(true);
    }

    inline void post_multiply_bases(double **AYmat_, double **Amat_, int mrows_, bool force_check_=false)
    {
      const int ivar_bse_case = check_bse_case(force_check_);
      switch (ivar_bse_case)
      {
        case -2: // multiple restricted dimensions, proceed by evaluating all
          post_multiply_evry_basis(AYmat_,Amat_,mrows_);
          break;
        case -1: // no restricted dimensions, treat as identity map
          for (size_t i = 0; i < mrows_; i++)
            for (size_t j = 0; j < ndof_use; j++)
              AYmat_[i][j] = Amat_[i][j];
          break;
        default: // single restricted dimnsion, proceed by evaluating one dimension
          post_multiply_ivar_basis(AYmat_,Amat_,mrows_,ivar_bse_case);
      }
    }

  protected:

    inline int verify_Yspaces(bool fset_complements_=false)
    {
      int ndof_acc = 0;
      for (size_t ivar = 0; ivar < nvar; ivar++)
        ndof_acc += (pbse_nvar_use[ivar] = (bse_ptrs[ivar]==NULL)?(perm_len):(bse_ptrs[ivar]->nvec_use));
      if (fset_complements_) set_Yspace_complements(ndof_acc);
      return ndof_use = ndof_acc;
    }
    inline void set_Yspace_complements(int ndof_)
    {
      const int ndof_full = perm_len*nvar;
      for (size_t ivar = 0, icmp = ndof_, iKiv = 0; ivar < nvar; ivar++, iKiv+=perm_len)
        if (bse_ptrs[ivar] != NULL)
          for (size_t ivc = bse_ptrs[ivar]->nvec_use; ivc < perm_len; ivc++, icmp++)
          {
            double  * const Vrow_cmp_i = spc.Vrowi_data(icmp);
            for (size_t jcmp = 0; jcmp < iKiv; jcmp++) Vrow_cmp_i[jcmp] = 0.0;
            LD_linalg::copy_vec<double>(Vrow_cmp_i+iKiv,spc.Vmat[ivc],perm_len);
            for (size_t jcmp = iKiv+perm_len; jcmp < ndof_full; jcmp++) Vrow_cmp_i[jcmp] = 0.0;
          }
    }
    inline int check_bse_case(bool force_check_=false)
    {
      if (force_check_) return -2; // check dimension to be safe
      bool ivar_flag = false;
      int ivar_out = -1; // if stays -1, then there are no restricted dimensons
      for (size_t ivar = 0; ivar < nvar; ivar++)
        if (bse_ptrs[ivar] != NULL)
        {
          ivar_out = ivar;
          if (ivar_flag) return -2; // if already one dimension is restricted, then multiple are restricted
          else ivar_flag = true; // flag that at least one dimension is restricted
        }
      return ivar_out;
    }

    // safe but disrupted
    void post_multiply_evry_basis(double **AYmat_, double **Amat_, int mrows_);
    // fast special case
    void post_multiply_ivar_basis(double **AYmat_, double **Amat_, int mrows_, int ivar_);
};

class LD_Theta_bundle: public ode_solspc_element
{

  const bool  Vbndl_owner,
              data_owner;

  protected:

    LD_vector_bundle * const Vbndle_ptr;
    LD_vector_space ** const Yspaces_nvar;

    const int perm_len;

  public:

    LD_Theta_bundle(ode_solspc_meta &meta_,int nspc_,int ndof_);
    LD_Theta_bundle(LD_Theta_bundle &Tbndle_);
    LD_Theta_bundle(ode_solspc_meta &meta_,int nspc_,int ndof_,LD_vector_bundle &Ybndle_,int ivar_=-1):
      LD_Theta_bundle(meta_,nspc_,ndof_) {set_Yspaces(Ybndle_,ivar_);}
    ~LD_Theta_bundle();

    LD_vector_bundle &Vbndle = *(Vbndle_ptr);
    LD_vspace_record &rec = Vbndle.rec;

    const int nspc = Vbndle.nspc;

    int ** const pbse_nvar_spcmat;
    double  *** const Ttns = Vbndle.Vtns;

    LD_vector_space ** const Vspaces = Vbndle.Vspaces;
    LD_Theta_space ** const Tspaces;

    inline void set_Yspaces(LD_vector_bundle &Ybndle_,int ivar_=-1)
    {
      LD_vector_space ** const Yspaces = Ybndle_.Vspaces;
      if (ivar_<0) // provide same basis for all dimensions
        for (size_t ispc = 0; ispc < nspc; ispc++)
          Vbndle.nV_spcvec[ispc] = Tspaces[ispc]->set_evry_Yspace(Yspaces[ispc]);
      else // provide basis for just one dimension
        for (size_t ispc = 0; ispc < nspc; ispc++)
          Vbndle.nV_spcvec[ispc] = Tspaces[ispc]->set_ivar_Yspace(Yspaces[ispc],ivar_);
    }
    inline void init_Vbndle_premult(double ***Wtns_)
    {
      for (size_t ispc = 0; ispc < nspc; ispc++)
        Vbndle.nV_spcvec[ispc] = Tspaces[ispc]->init_Vspce_premult(Wtns_[ispc]);
    }
    inline void set_Vspaces() {Vbndle.set_Vspaces();}
    inline void set_Vspaces(LD_vspace_record &rec_) {Vbndle.set_Vspaces(rec_);}
    // inline void set_Vspace_data_i(int i_) {Vbndle.set_Vspace_data_i(i_);}
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

    inline int verify_nrow(int ncod_=0,int nobs_=0)
      {return nrow = ( (ncod_>0)?(ncod=ncod_):(ncod) )*( (nobs_>0)?(nobs=nobs_):(nobs) );}
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
      {int acc = 0; for (size_t i = 0; i < nset; i++) acc += Acodes[i]->nobs; return acc;}
    inline int verify_nrows(int ncod_=0,int nobs_=0)
      {int acc = 0; for (size_t i = 0; i < nset; i++) acc += Acodes[i]->verify_nrow(ncod_,nobs_); return acc;}
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
      {int mrows = 0; for (size_t i = 0; i < nset; i++) if (mrows < Acodes[i]->nrow) mrows = Acodes[i]->nrow; return mrows;}
    inline int min_nrows()
      {int mrows = Acodes[0]->nrow; for (size_t i = 1; i < nset; i++) if (mrows > Acodes[i]->nrow) mrows = Acodes[i]->nrow; return mrows;}
    inline int nrows_mat_i(int i_) {return Acodes[i_]->nrow;}
};

struct LD_encoder
{
  LD_encoder(int ncod_,ode_solspc_meta &meta_): ncod(ncod_), meta(meta_) {}
  ~LD_encoder() {}

  ode_solspc_meta &meta;
  const int ncod;

  virtual void encode_normalize_rows(double **rows_i_,function_space_basis &fbse_,ode_solution &sol_i_,bool normalize_) = 0;

  template <class BSIS> static void encode_bundle(LD_encoding_bundle &bndle_,ode_solspc_subset &set_,BSIS **bases_,LD_encoder &enc_,bool normalize_,bool verbose_=true)
  {
    const int nobs_max = bndle_.nobs();
    ode_solution ** const sols = set_.sols;
    double t0 = LD_threads::tic();
    #pragma omp parallel
    {
      BSIS &base_t = *(bases_[LD_threads::thread_id()]);
      #pragma omp for
      for (size_t i = 0; i < nobs_max; i++)
      {
        enc_.encode_normalize_rows(bndle_.get_iobs_encoding(i),base_t,*(sols[i]),normalize_);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
      printf("(LD_encoder::encode_bundle) encoded %d row constraints (nobs = %d, ncod = %d) in %.2f seconds (%d threads).\n",
        nobs_max*(enc_.ncod), nobs_max, enc_.ncod, work_time, LD_threads::numthreads());
  }
};

struct LD_vspace_evaluator
{
  LD_vspace_evaluator(ode_solspc_subset &Sset_,int ncon_,int nvec_,int nsol_,double tol_);
  LD_vspace_evaluator(LD_vspace_evaluator &evl_,int nsol_): LD_vspace_evaluator(evl_.Sset,evl_.ncon,evl_.nvec_evl,nsol_,evl_.tol) {}
  ~LD_vspace_evaluator() {free_Tmatrix<bool>(sat_flags_mat);}

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
          for (size_t iV = 0, isat = 0, insat = nsat_seti; iV < nvec_evl; iV++)
            if (sat_setmat[iset][iV]) iVsat[iset][isat++] = iV;
            else iVsat[iset][insat++] = iV;
        }
        // else LD_linalg::fill_vec<int>(iVsat[iset],nvec_evl,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_vspace_evaluator::leniently_evaluate_vspace) evaluated %d conds. (%d sols) over %d bases (%.1f x %d, on avg.) in %.4f seconds (%d threads)\n",
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
          for (size_t iV = 0, isat = 0, insat = nsat_seti; iV < nvec_evl; iV++)
            if (sat_setmat[iset][iV]) iVsat[iset][isat++] = iV;
            else iVsat[iset][insat++] = iV;
        }
        // else LD_linalg::fill_vec<int>(iVsat[iset],nvec_evl,-1);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
    {
      printf("(LD_vspace_evaluator::strictly_evaluate_vspace) evaluated %d conds. (%d sols) over %d bases (%.1f x %d, on avg.) in %.4f seconds (%d threads)\n",
        nsol_acc*evl_.ncon,nsol_acc,nset,((double)nvec_acc)/((double)nset),reci_.vlen,
        work_time,LD_threads::numthreads());
    }
  }
};

struct LD_vspace_measure
{
  LD_vspace_measure(int vlen_, int nset0_):
    vlen(vlen_), nset0(nset0_), nset(nset0_),
    dsym0(Tsym<double>(nset0_-1)), dsym(new double*[nset0_-1]) {}
  ~LD_vspace_measure() {free_Tmatrix<double>(dsym0); delete [] dsym;}

  const int vlen;
  int nset;
  double ** const dsym;

  inline void init_distances(LD_vector_bundle &vbndle_, int nset_=0, bool verbose_=true)
    {init_distances(vbndle_.Vtns,vbndle_.nV_spcvec,nset_,verbose_);}
  virtual void init_distances(double ***Vtns_, int *nV_spcvec_, int nset_=0, bool verbose_=true) = 0;
  virtual double comp_distance(double **Va_,int na_,double **Vb_,int nb_) = 0;

  void update_distances(int nrdc_,int ij_prnts_[],LD_vector_space *vspcs_new_[],bool verbose_=true);

  inline double dij(int i_,int j_) {return (i_<j_)?(dsym[i_][j_-i_-1]):((i_>j_)?(dsym[j_][i_-j_-1]):(0.0));}

  protected:

    const int nset0,
              nset0m1 = nset0-1,
              dsym0_len = (nset0*nset0m1)/2;

    double  ** const dsym0,
             * const dsym0_vec = dsym0[0];

    static int get_isym_rowcol(int row_, int col_, int nsetm1_)
    {
      if (row_==col_) return -1;
      const int rw = (row_<col_)?(row_):(col_);
      int lrws = 0;
      for (size_t i = 0, lrw = nsetm1_; i < rw; i++, lrw--) lrws += lrw;
      return lrws + ((row_<col_)?(col_-row_):(row_-col_)) - 1;
    }
    static void get_rowcol_isym(int &row, int &col, int isym_, int nsetm1_)
    {
      row = 0;
      int lenr_dec = nsetm1_,
          isym_dec = isym_;
      while ((isym_dec -= (lenr_dec--)) >= 0) row++;
      col = isym_dec+lenr_dec+row+2;
    }

};

struct indepcomp_vspace_measure: public LD_vspace_measure
{
  indepcomp_vspace_measure(int vlen_, int nset0_): LD_vspace_measure(vlen_,nset0_) {}
  indepcomp_vspace_measure(LD_vector_bundle &vbndle_): LD_vspace_measure(vbndle_.vlen_full,vbndle_.nspc) {}
  ~indepcomp_vspace_measure() {}

  virtual void init_distances(double ***Vtns_, int *nV_spcvec_, int nset_=0, bool verbose_=true)
    {init_indepcomp_distances(Vtns_,nV_spcvec_,nset_,verbose_);}
  virtual double comp_distance(double **Va_,int na_,double **Vb_,int nb_)
    {return indepcomp_vspace_measure::comp_indp_distance(Va_,na_,Vb_,nb_,vlen);}

  inline void init_indepcomp_distances(LD_vector_bundle &vbndle_, int nset_=0, bool verbose_=true)
    {init_indepcomp_distances(vbndle_.Vtns,vbndle_.nV_spcvec,nset_,verbose_);}
  inline void init_indepcomp_distances(double ***Vtns_, int *nV_spcvec_, int nset_, bool verbose_)
    {indepcomp_vspace_measure::comp_indp_distances(dsym,dsym0_vec,(nset_)?(nset = nset_):(nset),Vtns_,nV_spcvec_,vlen,verbose_);}

  static void comp_indp_distances(double **dsym_,double *dvec_,int nset_,double ***Vtns_,int *nV_,int vlen_,bool verbose_=true)
  {
    const int nsetm1 = nset_-1,
              len_sym = (nsetm1*nset_)/2;
    for (size_t i = 0, jcap = nsetm1, isym = 0; i < nsetm1; i++, isym += jcap--) dsym_[i] = dvec_+isym;
    int ncol_acc = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:ncol_acc)
    {
      int i,j;
      #pragma omp for
      for (size_t isym = 0; isym < len_sym; isym++)
      {
        LD_vspace_measure::get_rowcol_isym(i,j,isym,nsetm1);
        dvec_[isym] = indepcomp_vspace_measure::comp_indp_distance(Vtns_[i],nV_[i],Vtns_[j],nV_[j],vlen_);

        ncol_acc += nV_[i] + nV_[j];
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_) printf("(indepcomp_vspace_measure::comp_indp_distances) computed distance matrix (%d vspaces, %d matrix products of avg. dims. %d x %.1f) in %.3f seconds\n",
      nset_, len_sym, vlen_, ncol_acc/((double) (2*len_sym)),work_time);
  }
  static double comp_indp_distance(double **Va_,int na_,double **Vb_,int nb_,int vlen_)
  {
    const int ma = vlen_-na_,
              mb = vlen_-nb_;
    double  ** const  Ya = Va_ + na_, // assumes that complementary bases stored behind
            ** const  Yb = Vb_ + nb_,
                      acc_indep = 0.0;

    // compute | KbT Ya |^2_F
    for (size_t i = 0; i < nb_; i++)
      for (size_t j = 0; j < ma; j++)
      {
        double acc_ij = 0.0;
        for (size_t l = 0; l < vlen_; l++) acc_ij += Vb_[i][l]*Ya[j][l];
        acc_indep += acc_ij*acc_ij;
      }

    // compute | YbT Ka |^2_F
    for (size_t i = 0; i < mb; i++)
      for (size_t j = 0; j < na_; j++)
      {
        double acc_ij = 0.0;
        for (size_t l = 0; l < vlen_; l++) acc_ij += Yb[i][l]*Va_[j][l];
        acc_indep += acc_ij*acc_ij;
      }

    return acc_indep;
  }


};

struct Frobenius_vspace_measure: public LD_vspace_measure
{
  Frobenius_vspace_measure(int vlen_, int nset0_): LD_vspace_measure(vlen_,nset0_) {}
  Frobenius_vspace_measure(LD_vector_bundle &vbndle_): LD_vspace_measure(vbndle_.vlen_full,vbndle_.nspc) {}
  ~Frobenius_vspace_measure() {}

  virtual void init_distances(double ***Vtns_, int *nV_spcvec_, int nset_=0, bool verbose_=true)
    {init_Frobenius_distances(Vtns_,nV_spcvec_,nset_,verbose_);}
  virtual double comp_distance(double **Va_,int na_,double **Vb_,int nb_)
    {return Frobenius_vspace_measure::comp_Frob_distance(Va_,na_,Vb_,nb_,vlen);}

  inline void init_Frobenius_distances(LD_vector_bundle &vbndle_, int nset_=0, bool verbose_=true)
    {init_Frobenius_distances(vbndle_.Vtns,vbndle_.nV_spcvec,nset_,verbose_);}
  inline void init_Frobenius_distances(double ***Vtns_, int *nV_spcvec_, int nset_, bool verbose_)
    {Frobenius_vspace_measure::comp_Frob_distances(dsym,dsym0_vec,(nset_)?(nset = nset_):(nset),Vtns_,nV_spcvec_,vlen,verbose_);}

  static void comp_Frob_distances(double **dsym_,double *dvec_,int nset_,double ***Vtns_,int *nV_,int vlen_,bool verbose_=true)
  {
    const int nsetm1 = nset_-1,
              len_sym = (nsetm1*nset_)/2;
    for (size_t i = 0, jcap = nsetm1, isym = 0; i < nsetm1; i++, isym += jcap--) dsym_[i] = dvec_+isym;
    int ncol_acc = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:ncol_acc)
    {
      int i,j;
      #pragma omp for
      for (size_t isym = 0; isym < len_sym; isym++)
      {
        LD_vspace_measure::get_rowcol_isym(i,j,isym,nsetm1);
        dvec_[isym] = Frobenius_vspace_measure::comp_Frob_distance(Vtns_[i],nV_[i],Vtns_[j],nV_[j],vlen_);

        ncol_acc += nV_[i] + nV_[j];
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_) printf("(Frobenius_vspace_measure::comp_Frob_distances) computed distance matrix (%d vspaces, %d matrix products of avg. dims. %d x %.1f) in %.3f seconds\n",
      nset_, len_sym, vlen_, ncol_acc/((double) (2*len_sym)),work_time);
  }
  static double comp_Frob_distance(double **Va_,int na_,double **Vb_,int nb_,int vlen_)
    {return 1.0 - sqrt(Frobenius_vspace_measure::comp_Frob_closeness(Va_,na_,Vb_,nb_,vlen_)/((double) ((na_<nb_)?(na_):(nb_))));}
  static double comp_Frob_closeness(double **Va_,int na_,double **Vb_,int nb_,int vlen_)
  {
    double acc_fro = 0.0;
    for (size_t i = 0; i < na_; i++)
      for (size_t j = 0; j < nb_; j++)
      {
        double acc_ij = 0.0;
        for (size_t l = 0; l < vlen_; l++) acc_ij += Va_[i][l]*Vb_[j][l];
        acc_fro += acc_ij*acc_ij;
      }
    return acc_fro;
  }
};

struct LD_svd_bundle: public LD_vector_bundle
{
  LD_svd_bundle(int nspc_,int vlen_): LD_vector_bundle(nspc_,vlen_),
    rank_vec(new int[nspc_]), Smat(Tmatrix<double>(nspc_,vlen_)) {}
  LD_svd_bundle(LD_vector_bundle &Vbndl_): LD_vector_bundle(Vbndl_),
    rank_vec(new int[Vbndl_.nspc]), Smat(Tmatrix<double>(Vbndl_.nspc,Vbndl_.vlen_full)) {}
  LD_svd_bundle(LD_encoding_bundle &Abndl_,bool verbose_=true):
    LD_svd_bundle(Abndl_.nset,Abndl_.ncol_full)
    {compute_Acode_curve_svds(Abndl_,verbose_);}
  LD_svd_bundle(LD_encoding_bundle &Abndle_,LD_Theta_bundle &Tbndle_,bool init_=false,bool verbose_=true):
    LD_svd_bundle(Abndle_.nset,Abndle_.ncol_full)
  {
    compute_AYmat_curve_svds(Abndle_,Tbndle_.Tspaces,verbose_);
    if (init_) Tbndle_.init_Vbndle_premult(VTtns);
  }
  ~LD_svd_bundle()
  {
    free_Tmatrix<double>(Smat);
    delete [] rank_vec;
  }

  int * const rank_vec;
  double  ** const Smat,
          *** const VTtns = Vtns_data;

  static void project_Theta_bundle(LD_Theta_bundle &Tbndle_,LD_encoding_bundle &Abndle_,bool verbose_=true)
    {LD_svd_bundle svd_bndle(Abndle_,Tbndle_,true,verbose_);}

  void compute_Acode_curve_svds(LD_encoding_bundle &Abndle_,bool verbose_=true);
  void compute_AYmat_curve_svds(LD_encoding_bundle &Abndle_,LD_Theta_space ** const Tspaces_,bool verbose_=true);

  inline void evaluate_iscl_ranks(int iscl_)
    {for (size_t i = 0; i < nspc; i++) rank_vec[i] = LD_svd::rank(Smat[i],nV_spcvec[i],iscl_);}

  void print_details(const char name_[] =  "Amat_SVD");

  inline int nulldim_i(int i_) {return nV_spcvec[i_]-rank_vec[i_];}
  inline int min_rank() {return LD_linalg::min_val<int>(rank_vec,nspc);}
  inline int max_rank() {return LD_linalg::max_val<int>(rank_vec,nspc);}
  inline int min_nulldim()
  {
    int nd_out = vlen_full, nd_i;
    for (size_t i = 0; i < nspc; i++)
      if (nd_out > (nd_i=nulldim_i(i))) nd_out = nd_i;
    return nd_out;
  }
  inline int max_nulldim()
  {
    int nd_out = 0, nd_i;
    for (size_t i = 0; i < nspc; i++)
      if (nd_out < (nd_i=nulldim_i(i))) nd_out = nd_i;
    return nd_out;
  }

  inline void set_Vspaces_nullspc() {rec.set_record_nullspc(rank_vec); set_Vspaces();}
  inline void set_Vspaces_rankvec() {rec.set_record_rankvec(rank_vec); set_Vspaces();}

  void write_LD_svd_bundle(const char name_[]);
};

#endif
