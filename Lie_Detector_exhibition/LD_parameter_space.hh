#ifndef LD_PARSPC_HH
#define LD_PARSPC_HH

#include "LD_encodings.hh"

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
  ~LD_vspace_record() { free_Tmatrix<int>(iV_spcmat); delete [] nV_spcvec; }

  const int nspc,
            nvec,
            vlen;
  int * const nV_spcvec,
      ** const iV_spcmat;
  double *** const Vtns_data;

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
  {
    for (size_t i = 0; i < nspc; i++)
    {
      LD_linalg::fill_vec_012<int>(iV_spcmat[i],nvec_,nV_spcvec[i]-nvec_);
      nV_spcvec[i] = nvec_;
    }
  }

  void print_selected_details(const char name_[], bool longwinded_=true);
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
    ~LD_vector_space();

    const int vlen_full;

    int nvec_use,
        vlen_use,
        * const inds_V;

    double ** const Vmat;

    inline void set_Vmat(int nvec_use_)
      {nvec_use = nvec_use_; for (size_t i = 0; i < nvec_use; i++) Vmat[i] = Vmat_data[inds_V[i]];}
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

    inline void load_Vtns_data(double ***Vtns_)
    {
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
      if ((0<=ivar_)&&(ivar_<nvar)) {bse_ptrs[ivar_] = Yspace_; return spc.nvec_use = verify_Yspaces();}
      else return spc.nvec_use = set_evry_Yspace(Yspace_);
    }
    inline int set_evry_Yspace(LD_vector_space *Yspace_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspace_;
      return spc.nvec_use = verify_Yspaces();
    }
    inline int set_evry_Yspace(LD_vector_space **Yspaces_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspaces_[ivar];
      return spc.nvec_use = verify_Yspaces();
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

    inline int verify_Yspaces()
    {
      int ndof_acc = 0;
      for (size_t ivar = 0; ivar < nvar; ivar++)
        ndof_acc += (pbse_nvar_use[ivar] = (bse_ptrs[ivar]==NULL)?(perm_len):(bse_ptrs[ivar]->nvec_use));
      return ndof_use = ndof_acc;
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
    ~LD_Theta_bundle();

    LD_vector_bundle &Vbndle = *(Vbndle_ptr);
    LD_vspace_record &Trec = Vbndle.rec;

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
    inline void set_Tspaces() {Vbndle.set_Vspaces();}
    inline void set_Tspaces(LD_vspace_record &rec_) {Vbndle.set_Vspaces(rec_);}
};

struct Theta_eval_workspace
{
  Theta_eval_workspace(int ntheta_,int nsols_=1): ntheta(ntheta_), nsols(nsols_),
    ntheta_use(ntheta_), Theta(NULL),
    satisfy_flags_mat(Tmatrix<bool>(nsols_,ntheta_)) {}
  ~Theta_eval_workspace() {free_Tmatrix<bool>(satisfy_flags_mat);}

  const int ntheta,
            nsols;
  int ntheta_use;
  bool  ** const satisfy_flags_mat,
        * satisfy_flags = satisfy_flags_mat[0];
  double  ** Theta;
};

struct lambda_map_eval_workspace: public Theta_eval_workspace
{
  lambda_map_eval_workspace(int ntheta_,function_space &fspc_,int nsols_=1):
    Theta_eval_workspace(ntheta_,nsols_),
    lamvec(new double[fspc_.perm_len]), vxu_wkspc(vxu_workspace(fspc_.nvar,fspc_.comp_ord_len())) {}
  ~lambda_map_eval_workspace() {delete [] lamvec;}

  double * const lamvec;
  vxu_workspace vxu_wkspc;
};

struct Rn_cond_eval_workspace: public Theta_eval_workspace
{
  Rn_cond_eval_workspace(int ntheta_,ode_solspc_meta &meta_,int nsols_=1):
    Theta_eval_workspace(ntheta_,nsols_), v(new double[meta_.ndim]) {}
  ~Rn_cond_eval_workspace() {delete [] v;}

  double  * const v;
};

struct inf_crit_eval_workspace: public Rn_cond_eval_workspace
{
  inf_crit_eval_workspace(int ntheta_,ode_solspc_meta &meta_,int nsols_=1):
    Rn_cond_eval_workspace(ntheta_,meta_,nsols_), mags_JFs(new double[meta_.ndep]) {}
  ~inf_crit_eval_workspace() {delete [] mags_JFs;}

  double  * const mags_JFs;
};

#endif
