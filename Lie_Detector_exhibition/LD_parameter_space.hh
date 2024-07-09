#ifndef LD_PARSPC_HH
#define LD_PARSPC_HH

#include "LD_aux.hh"
#include "LD_function_space.hh"

struct LD_vspace_record
{
  LD_vspace_record(int nspc_,int nvec_, int vlen_, double ***Vtns_data_): nspc(nspc_), nvec(nvec_), vlen(vlen_),
    nV_spcvec(new int[nspc_]), iV_spcmat(Tmatrix<int>(nspc_,nvec_)), Vtns_data(Vtns_data_) {}
  ~LD_vspace_record() { free_Tmatrix<int>(iV_spcmat); delete [] nV_spcvec; }

  const int nspc,
            nvec,
            vlen;
  int * const nV_spcvec,
      ** const iV_spcmat;
  double *** const Vtns_data;

  inline void copy_record(LD_vspace_record &rec_)
  {
    LD_linalg::copy_vec<int>(nV_spcvec,rec_.nV_spcvec,nspc);
    for (size_t i = 0; i < nspc; i++)
      for (size_t j = 0; j < nV_spcvec[i]; j++) iV_spcmat[i][j] = rec_.iV_spcmat[i][j];
  }

  void print_selected_details(const char name_[], bool longwinded_=true);
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

    inline void set_Vmat(int nvec_use_,int vlen_use_)
    {
      nvec_use = nvec_use_; vlen_use=vlen_use_;
      for (size_t i = 0; i < nvec_use; i++) Vmat[i] = Vmat_data[inds_V[i]];
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
    LD_vector_bundle(int nspc_,int vlen_, double ***Vtns_in_): LD_vector_bundle(nspc_,vlen_) {load_Vtns_data(Vtns_in_);}
    LD_vector_bundle(LD_vspace_record &rec_): LD_vector_bundle(rec_.nspc,rec_.vlen,rec_.Vtns_data) {set_Vspaces(rec_,rec_.vlen);}
    LD_vector_bundle(LD_vector_bundle &bndle_);
    ~LD_vector_bundle();

    inline void load_Vtns_data(double ***Vtns_)
    {
      for (size_t i = 0; i < nspc; i++)
        for (size_t j = 0; j < nvec_use; j++)
          for (size_t k = 0; k < vlen_use; k++)
            Vtns_data[i][j][k] = Vtns_[i][j][k];
    }

    LD_vspace_record &rec = *(rec_ptr);

    const int nspc,
              vlen_full;
    int nvec_use,
        vlen_use,
        *  const nV_spcvec = rec.nV_spcvec,
        ** const iV_spcmat = rec.iV_spcmat;
    double *** const Vtns;

    LD_vector_space ** const Vspaces;

    inline void set_Vspaces(LD_vspace_record &rec_, int vlen_use_=0)
    {
      vlen_use = (vlen_use_)?(vlen_use_):(rec_.vlen); rec.copy_record(rec_);
      for (size_t i = 0; i < nspc; i++) Vspaces[i]->set_Vmat(nV_spcvec[i],vlen_use);
    }
    inline void set_Vspaces(int vlen_use_=0)
    {
      vlen_use = (vlen_use_)?(vlen_use_):(vlen_full);
      for (size_t i = 0; i < nspc; i++) Vspaces[i]->set_Vmat(nV_spcvec[i],vlen_use);
    }
    inline bool check_default_configuration()
    {
      for (size_t i = 0; i < nspc; i++) if (!(Vspaces[i]->check_default_configuration())) return false;
      return true;
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
      if ((0<=ivar_)&&(ivar_<nvar)) {bse_ptrs[ivar_] = Yspace_; return verify_Yspaces();}
      else return set_evry_Yspace(Yspace_);
    }
    inline int set_evry_Yspace(LD_vector_space *Yspace_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspace_;
      return verify_Yspaces();
    }
    inline int set_evry_Yspace(LD_vector_space **Yspaces_)
    {
      for (size_t ivar = 0; ivar < nvar; ivar++) bse_ptrs[ivar] = Yspaces_[ivar];
      return verify_Yspaces();
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
    LD_vector_space ** const Vspaces = Vbndle_ptr->Vspaces,
                    ** const Yspaces_nvar;

    const int perm_len;

  public:

    LD_Theta_bundle(ode_solspc_meta &meta_,int nspc_,int ndof_);
    LD_Theta_bundle(LD_Theta_bundle &Tbndle_);
    ~LD_Theta_bundle();

    LD_vector_bundle &Vbndle = *(Vbndle_ptr);

    const int nspc = Vbndle.nspc;

    int * const ndof_spcvec,
        ** const pbse_nvar_spcmat;
    double *** const Ttns = Vbndle.Vtns;

    LD_Theta_space ** const Tspaces;

    inline void set_Yspaces(LD_vector_bundle &Ybndle_,int ivar_=-1)
    {
      LD_vector_space ** const Yspaces = Ybndle_.Vspaces;
      if (ivar_<0) // provide same basis for all dimensions
        for (size_t ispc = 0; ispc < nspc; ispc++)
          ndof_spcvec[ispc] = Tspaces[ispc]->set_evry_Yspace(Yspaces[ispc]);
      else // provide basis for just one dimension
        for (size_t ispc = 0; ispc < nspc; ispc++)
          ndof_spcvec[ispc] = Tspaces[ispc]->set_ivar_Yspace(Yspaces[ispc],ivar_);
    }
    inline void init_Vbndle_premult(double ***Wtns_)
    {
      for (size_t ispc = 0; ispc < nspc; ispc++)
        Vbndle.nV_spcvec[ispc] = ndof_spcvec[ispc] = Tspaces[ispc]->init_Vspce_premult(Wtns_[ispc]);
    }
};

struct Theta_eval_workspace
{
  Theta_eval_workspace(int ntheta_,int nsols_=1): ntheta(ntheta_), nsols(nsols_), Theta(NULL),
    satisfy_flags_mat(Tmatrix<bool>(nsols_,ntheta_)) {}
  ~Theta_eval_workspace() {free_Tmatrix<bool>(satisfy_flags_mat);}

  const int ntheta,
            nsols;
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
    Theta_eval_workspace(ntheta_,nsols_), ntheta_use(ntheta_), v(new double[meta_.ndim]) {}
  ~Rn_cond_eval_workspace() {delete [] v;}

  int ntheta_use;
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