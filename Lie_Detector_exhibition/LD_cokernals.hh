#ifndef LD_COKERN_HH
#define LD_COKERN_HH

#include "LD_framework.hh"
// #include "LD_vspace_evaluators.hh"

struct ode_solspc_setbundle: public ode_solspc_subset
{
  ode_solspc_setbundle(ode_solspc_meta &meta_,int nobs_full_,double **pts_,ode_solution **sols_,int nset_,int *nsol_per_set_):
    ode_solspc_subset(meta_,nobs_full_,pts_,sols_), nset(nset_), nsol_per_set(nsol_per_set_) {}
  ~ode_solspc_setbundle() {}

  const int nset;
  int * const nsol_per_set;

  virtual int nobs_subset_i(int i_) {return nsol_per_set[i_];}
  virtual int max_nobs_subset() {return LD_linalg::max_val<int>(nsol_per_set,nset);}
  virtual ode_solution ** get_sol_subset_i(int i_) {return sols + cumnobs_subset_i(i_);}

  inline int cumnobs_subset_i(int i_)
    {int acc = 0; for (size_t i = 0; i < i_; i++) acc+=nsol_per_set[i]; return acc;}
};

class cokernal_workspace
{
  const bool data_owner;

  public:
    cokernal_workspace(int iwkspc_,bool &factive_,double *wvec_, LD_vector_space *vspc_): data_owner(false),
      iwkspc(iwkspc_), factive(factive_), wvec(wvec_), vspc(vspc_) {}
    cokernal_workspace(int iwkspc_,bool &factive_,int vlen_): data_owner(true),
      iwkspc(iwkspc_), factive(factive_), wvec(new double[vlen_]), vspc(new LD_vector_space(vlen_)) {}
    ~cokernal_workspace() {if (data_owner) {factive=false; delete [] wvec; delete vspc;}}

    const int iwkspc;
    bool &factive;
    double * const wvec;
    LD_vector_space * const vspc;
};

struct cokernal_space
{
  // initializing cokernal space to original curve kernals
  cokernal_space(cokernal_workspace *wkspc_):
    factive(wkspc_->factive),
    gen_spawn(0), n_0Vkrn(1),
    iwkspc(wkspc_->iwkspc), ickrn_parent(iwkspc), jckrn_parent(iwkspc),
    ickrn(iwkspc), nvec_use((wkspc_->vspc)->nvec_use), inds_V((wkspc_->vspc)->inds_V),
    ivec_0Vkrn(new int[n_0Vkrn]),
    wvec(wkspc_->wvec), Vmat((wkspc_->vspc)->Vmat), Vmat_data((wkspc_->vspc)->get_Vmat_data())
    {ivec_0Vkrn[0] = iwkspc;}
  // allocating new cokernal space by overwriting given cokernal workspace
  cokernal_space(int gen_spawn_,cokernal_workspace *wkspc_,cokernal_space * ckrni_,cokernal_space * &ckrnj_):
    factive(wkspc_->factive),
    gen_spawn(gen_spawn_), n_0Vkrn(ckrni_->n_0Vkrn + ckrnj_->n_0Vkrn),
    iwkspc(wkspc_->iwkspc), ickrn_parent(ckrni_->ickrn), jckrn_parent(ckrnj_->ickrn),
    ickrn(ckrni_->ickrn), nvec_use((wkspc_->vspc)->nvec_use), inds_V((wkspc_->vspc)->inds_V),
    ivec_0Vkrn(new int[n_0Vkrn]),
    wvec(wkspc_->wvec), Vmat((wkspc_->vspc)->Vmat), Vmat_data((wkspc_->vspc)->get_Vmat_data())
    {
      set_ivec_0Vkrn(ckrni_->ivec_0Vkrn,ckrni_->n_0Vkrn,ckrnj_->ivec_0Vkrn,ckrnj_->n_0Vkrn);
      factive = true;
      ckrni_->factive = ckrnj_->factive = false;
      delete ckrni_; delete ckrnj_; ckrnj_ = NULL;
    }
  // overwriting cokernal space with ckrni_ workspace
  cokernal_space(int gen_spawn_,cokernal_space * ckrni_,cokernal_space * &ckrnj_):
    factive(ckrni_->factive),
    gen_spawn(gen_spawn_), n_0Vkrn(ckrni_->n_0Vkrn + ckrnj_->n_0Vkrn),
    iwkspc(ckrni_->iwkspc), ickrn_parent(ckrni_->ickrn), jckrn_parent(ckrnj_->ickrn),
    ickrn(ckrni_->ickrn), nvec_use(ckrni_->nvec_use), inds_V(ckrni_->inds_V),
    ivec_0Vkrn(new int[n_0Vkrn]),
    wvec(ckrni_->wvec), Vmat(ckrni_->Vmat), Vmat_data(ckrni_->Vmat_data)
    {
      set_ivec_0Vkrn(ckrni_->ivec_0Vkrn,ckrni_->n_0Vkrn,ckrnj_->ivec_0Vkrn,ckrnj_->n_0Vkrn);
      factive = true;
      ckrnj_->factive = false;
      delete ckrni_; delete ckrnj_; ckrnj_ = NULL;
    }

  ~cokernal_space() {delete [] ivec_0Vkrn;}

  bool  &factive;
  const int gen_spawn,
            n_0Vkrn,
            iwkspc,
            ickrn_parent,
            jckrn_parent;
  int ickrn,
      &nvec_use,
      * const inds_V,
      * const ivec_0Vkrn;
  double  * const wvec,
          ** const Vmat,
          ** const Vmat_data;

  inline void init_cokernal_space(int vlen_,int nsat_,double *wvec_,double *Vm_d_,int *isat_)
  {
    const int vl2 = vlen_*vlen_;
    for (size_t i = 0; i < vlen_; i++) wvec[i] = wvec_[i];
    for (size_t i = 0; i < vl2; i++) Vmat_data[0][i] = Vm_d_[i];
    for (size_t i = 0; i < vlen_; i++) Vmat[i] = Vmat_data[inds_V[i] = isat_[i]];
    nvec_use = nsat_;
  }

  private:

    inline void set_ivec_0Vkrn(int *ivec_i_,int nveci_,int *ivec_j_,int nvecj_)
    {
      for (size_t i = 0; i < nveci_; i++) ivec_0Vkrn[i] = ivec_i_[i];
      for (size_t i = 0; i < nvecj_; i++) ivec_0Vkrn[i+nveci_] = ivec_j_[i];
      int imin, iimin;
      sort_ivec_0Vkrn(ivec_0Vkrn,nveci_+nvecj_,imin,iimin);
    }
    inline void sort_ivec_0Vkrn(int *ivec_,int ni_,int &im_,int &iim_)
    {
      im_ = ivec_[iim_=0];
      for (size_t i = 1; i < ni_; i++) if (im_>ivec_[i]) im_ = ivec_[iim_=i];
      ivec_[iim_] = ivec_[0]; ivec_[0] = im_;
      if (ni_>2) sort_ivec_0Vkrn(ivec_+1,ni_-1,im_,iim_);
    }

};

class cokernal_bundle
{

  const int nset0,
            vlen;

  bool  * const factive_ckrn_wkspc,
        * const factive_ckrn_wkspc_xtra = factive_ckrn_wkspc + nset0;

  int nxtra_wkspc,
      * const cokern_dimvec;

  cokernal_workspace  ** const ckrn_wkspcs,
                      ** const ckrn_wkspcs_xtra = ckrn_wkspcs + nset0;

  public:

    cokernal_bundle(LD_vector_bundle &Vb0_,int *cokern_dimvec_,double **Wmat0_):
      nset0(Vb0_.nspc), vlen(Vb0_.vlen_full),
      factive_ckrn_wkspc(new bool[2*nset0]),
      nxtra_wkspc(0),
      cokern_dimvec(cokern_dimvec_),
      ckrn_wkspcs(new cokernal_workspace*[2*nset0]),
      ckrn_spcs(new cokernal_space*[nset0])
      {
        LD_vector_space ** const Vspcs0 = Vb0_.Vspaces;
        for (size_t i = 0; i < nset0; i++)
        {
          ckrn_spcs[i] = new cokernal_space( ckrn_wkspcs[i] =
                          new cokernal_workspace(i, factive_ckrn_wkspc[i] = true, Wmat0_[i], Vspcs0[i]) );

          ckrn_wkspcs_xtra[i] = NULL;
          factive_ckrn_wkspc_xtra[i] = false;
        }
      }

    ~cokernal_bundle()
    {
      delete [] factive_ckrn_wkspc;

      for (size_t i = 0; i < nset0; i++) if (ckrn_spcs[i] != NULL) delete ckrn_spcs[i];
      delete [] ckrn_spcs;
      for (size_t i = 0; i < nset0; i++) {delete ckrn_wkspcs[i]; if (ckrn_wkspcs_xtra[i] != NULL) delete ckrn_wkspcs_xtra[i];}
      delete [] ckrn_wkspcs;
    }

    int &nset = cokern_dimvec[0],
        &kSC = cokern_dimvec[1],
        &nred_succ = cokern_dimvec[2];

    cokernal_space ** const ckrn_spcs;

    int collapse_cokernals(LD_vspace_measure &msr_,int gen_spawn_,bool *frdc_succ_,int *nsat_ckrn_,int *ij_prs_,
      double *w_ckrn_,double *vvec_ckrn_,int *isat_ckrn_,bool wdistance_,bool verbose_);
};

struct Jet_function_vector_space
{

  Jet_function_vector_space(LD_observations_set &Sobs_,function_space &fspc_,LD_encoder &enc_):
    nset0(Sobs_.ncrvs_tot), vlen_full(fspc_.ndof_full),
    Sobs(Sobs_), fspc(fspc_), enc(enc_),
    Acode(LD_encoding_bundle(nset0,vlen_full,Sobs_.npts_per_crv,enc_.ncod)),
    Vbndle0(LD_vector_bundle(nset0,vlen_full)),
    svd0(LD_svd_bundle(Vbndle0))
    {}

  template <class BSE> Jet_function_vector_space(LD_observations_set &Sobs_,function_space &fspc_,LD_encoder &enc_,
    BSE **bases_, bool normalize_, bool verbose_=true):
    Jet_function_vector_space(Sobs_,fspc_,enc_)
  {
    LD_encoder::encode_bundle<BSE>(Acode,Sobs,bases_,enc,normalize_,verbose_);
  }
  ~Jet_function_vector_space() {}

  const int nset0,
            vlen_full;

  LD_observations_set &Sobs;
  function_space &fspc;
  LD_encoder &enc;

  LD_encoding_bundle Acode;
  LD_vector_bundle Vbndle0;
  LD_svd_bundle svd0; // shared data with Vbndle0

  LD_vspace_record &rec0 = Vbndle0.rec;

};

class cokernal_policy
{

  public:

    cokernal_policy(LD_vspace_measure &msr_,bool wdistance_): wdistance(wdistance_), msr(msr_) {}
    ~cokernal_policy() {}

    const bool wdistance;
    LD_vspace_measure &msr;

    virtual double ** init_cokernal_collapse(Jet_function_vector_space &jfvs_,bool verbose_=true) = 0;
    virtual int iterate_cokernal_reduction(cokernal_bundle &cokern_,int gen_,bool verbose_=true) = 0;
};

class nullspace_ckrn_policy: public cokernal_policy
{
  public:

    nullspace_ckrn_policy(LD_vspace_measure &msr_,bool wdistance_): cokernal_policy(msr_,wdistance_) {}
    ~nullspace_ckrn_policy() {}

    virtual double ** init_cokernal_collapse(Jet_function_vector_space &jfvs_,bool verbose_=true)
    {
      jfvs_.svd0.compute_Acode_curve_svds(jfvs_.Acode,verbose_);
      jfvs_.svd0.evaluate_iscl_ranks(2*(jfvs_.vlen_full));
      if (verbose_) jfvs_.svd0.print_details("svd0");

      jfvs_.svd0.set_Vspaces_nullspc();
      if (verbose_) jfvs_.rec0.print_subspace_details("rec0",false,false);

      if (wdistance) msr.init_distances(jfvs_.Vbndle0.Vspaces,jfvs_.svd0.Smat,jfvs_.nset0,verbose_);
      else msr.init_distances(jfvs_.Vbndle0.Vspaces,jfvs_.nset0,verbose_);

      return jfvs_.svd0.Smat;
    }

  protected:

    int iterate_nullspace_cokernals(cokernal_bundle &cokern_,int gen_,k_medoids_package &kmed_,bool verbose_);

};

class nullspace_clst_policy: public nullspace_ckrn_policy
{
  public:

    nullspace_clst_policy(LD_vspace_measure &msr_,bool wdistance_): nullspace_ckrn_policy(msr_,wdistance_) {}
    ~nullspace_clst_policy() {}

    virtual int iterate_cokernal_reduction(cokernal_bundle &cokern_,int gen_,bool verbose_=true)
    {
      k_medoids_package kmed(msr.dsym,cokern_.nset);
      printf("\n");
      const int kSC0 = cokern_.kSC = kmed.comp_kSC_medoids(2,verbose_);
      if (iterate_nullspace_cokernals(cokern_,gen_,kmed,verbose_)) return cokern_.nred_succ;
      else
      {
        if (kSC0 > 2)
        {
          printf("\n(nullspace_clst_policy::iterate_cokernal_reduction) gen %d - IRREDUCABLE kSC0=%d clusters (nset=%d). Attempting to reduce kSC < %d clusters.\n",gen_,kSC0,cokern_.nset,kSC0);
          int kSCm1 = kSC0-1; cokern_.kSC = kmed.comp_kSC_krange_medoids(2,kSCm1,verbose_);
          while ( (kSCm1 >= 2) && !(iterate_nullspace_cokernals(cokern_,gen_,kmed,verbose_)) )
            if (kSCm1 == 2) break;
            else cokern_.kSC = kmed.comp_kSC_krange_medoids(2,--kSCm1,verbose_);
          if (cokern_.nred_succ) return cokern_.nred_succ;
        }
        printf("\n(nullspace_clst_policy::iterate_cokernal_reduction) gen %d - IRREDUCABLE CLUSTERS (kSC0=%d, nset=%d). Attempting to consolidate remaining cokernal spaces\n",gen_,kSC0,cokern_.nset);
        // if all else fails, just try consolidating one pair
        kmed.set_one_medoid(); cokern_.kSC = 1;
        return iterate_nullspace_cokernals(cokern_,gen_,kmed,verbose_);
      }
    }

};

class nullspace_near_policy: public nullspace_ckrn_policy
{
  public:

    nullspace_near_policy(LD_vspace_measure &msr_,bool wdistance_): nullspace_ckrn_policy(msr_,wdistance_) {}
    ~nullspace_near_policy() {}

    virtual int iterate_cokernal_reduction(cokernal_bundle &cokern_,int gen_,bool verbose_=true)
    {
      k_medoids_package kmed(msr.dsym,cokern_.nset);
      printf("\n");
      kmed.set_one_medoid(); cokern_.kSC = 1;
      return iterate_nullspace_cokernals(cokern_,gen_,kmed,verbose_);
    }

};

class cokernal_sub_bundle
{
  const int vlen_full;
  int nset,
      kSC,
      nred_succ,
      * const cokern_dimvec = &nset;

  cokernal_bundle cokern;
  cokernal_space ** const ckrn_spcs =  cokern.ckrn_spcs;

  public:

    cokernal_sub_bundle(Jet_function_vector_space &jfvs_,cokernal_policy &pol_,bool verbose_);
    ~cokernal_sub_bundle() {}

};


#endif
