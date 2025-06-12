#ifndef LD_COKERN_HH
#define LD_COKERN_HH

#include "LD_framework.hh"

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

class cokernal_family
{
  cokernal_workspace ** const wkspcs;

  public:

    cokernal_family(cokernal_space &ckrn_med_,int nckrn_,cokernal_workspace **wkspcs_):
      wkspcs(wkspcs_),
      ckrn_med(ckrn_med_),
      nckrn(nckrn_), nkrn0_ckrn(new int[nckrn_]),
      ckrn_spcs(new cokernal_space*[nckrn_]) {}
    cokernal_family(int nckrn_,int imed_,int ipts_med_[],cokernal_space **ckrn_spcs_,cokernal_workspace **wkspcs_):
      cokernal_family(*(ckrn_spcs_[imed_]),nckrn_,wkspcs_)
      {for (size_t i = 0; i < nckrn_; i++) nkrn0_ckrn[i] = (ckrn_spcs[i] = ckrn_spcs_[ipts_med_[i]])->n_0Vkrn;}

    ~cokernal_family()
    {
      delete [] nkrn0_ckrn;
      delete [] ckrn_spcs;
    }

    cokernal_space &ckrn_med;

    const int nckrn;
    int * const nkrn0_ckrn;
    cokernal_space ** const ckrn_spcs;

    void print_details(int k_)
    {
      printf("\n(cokernal_family::print_details) k=%d cokernal family: imed = %d, nkcrn = %d, nkrn0=%d, nvec0=%d -> %d=nvecf \n", k_,ckrn_med.ickrn,nckrn,nkrn0(),nvec0(),nvecf());
      for (int i = 0; i < nckrn; i++)
      {
        printf("  %s ckrn %d, ickrn=%d, nkrn0_i=%d, nvec0=%d -> %d=nvecf %s jkrn0: ",
          (ckrn_med.ickrn == ckrn_spcs[i]->ickrn)?("["):("("),
            i,ckrn_spcs[i]->ickrn,ckrn_spcs[i]->n_0Vkrn, nvec0_i(i), ckrn_spcs[i]->nvec_use,
          (ckrn_med.ickrn == ckrn_spcs[i]->ickrn)?("]"):(")"));
        for (int j = 0; j < ckrn_spcs[i]->n_0Vkrn; j++) printf("%d ", ckrn_spcs[i]->ivec_0Vkrn[j]);
        printf("\n");
      }
    }

    inline int nvecf() {int nvecf_out = 0; for (size_t i = 0; i < nckrn; i++) nvecf_out += ckrn_spcs[i]->nvec_use; return nvecf_out;}
    inline int nkrn0() {int nkrn0_out = 0; for (size_t i = 0; i < nckrn; i++) nkrn0_out += ckrn_spcs[i]->n_0Vkrn; return nkrn0_out;}

    inline int nvec0()
    {
      int nvec_acc = 0;
      for (size_t i = 0; i < nckrn; i++)
        for (size_t j = 0; j < nkrn0_ckrn[i]; j++)
          nvec_acc += wkspcs[ckrn_spcs[i]->ivec_0Vkrn[j]]->vspc->nvec_use;
      return nvec_acc;

    }

    inline int nvec0_i(int i_)
    {
      int nvec_acc = 0;
      for (size_t j = 0; j < nkrn0_ckrn[i_]; j++) nvec_acc += wkspcs[ckrn_spcs[i_]->ivec_0Vkrn[j]]->vspc->nvec_use;
      return nvec_acc;
    }
};

class cokernal_bundle
{

  const int nset0,
            vlen;

  bool  * const factive_ckrn_wkspc,
        * const factive_ckrn_wkspc_xtra = factive_ckrn_wkspc + nset0;

  int nxtra_wkspc;

  double ** const Wmat;

  cokernal_workspace  ** const ckrn_wkspcs,
                      ** const ckrn_wkspcs_xtra = ckrn_wkspcs + nset0;

  LD_vector_space ** const vspcs;

  public:

    cokernal_bundle(LD_vector_bundle &Vb0_,double **Wmat0_):
      nset0(Vb0_.nspc), vlen(Vb0_.vlen_full),
      factive_ckrn_wkspc(new bool[2*nset0]),
      nxtra_wkspc(0),
      Wmat(new double*[nset0]),
      ckrn_wkspcs(new cokernal_workspace*[2*nset0]),
      vspcs(new LD_vector_space*[nset0]),
      nset(nset0), kSC(2), nred_succ(0), generation(0),
      ckrn_spcs(new cokernal_space*[nset0])
      {
        LD_vector_space ** const Vspcs0 = Vb0_.Vspaces;
        for (size_t i = 0; i < nset0; i++)
        {
          ckrn_spcs[i] = new cokernal_space( ckrn_wkspcs[i] =
                          new cokernal_workspace(i, factive_ckrn_wkspc[i] = true, Wmat[i] = Wmat0_[i], vspcs[i] = Vspcs0[i]));
          ckrn_wkspcs_xtra[i] = NULL;
          factive_ckrn_wkspc_xtra[i] = false;
        }
      }

    ~cokernal_bundle()
    {
      delete [] factive_ckrn_wkspc;
      delete [] Wmat; delete [] vspcs;

      for (size_t i = 0; i < nset0; i++) if (ckrn_spcs[i] != NULL) delete ckrn_spcs[i];
      delete [] ckrn_spcs;
      for (size_t i = 0; i < nset0; i++) {delete ckrn_wkspcs[i]; if (ckrn_wkspcs_xtra[i] != NULL) delete ckrn_wkspcs_xtra[i];}
      delete [] ckrn_wkspcs;
    }

    int nset,
        kSC,
        nred_succ,
        generation;

    cokernal_space ** const ckrn_spcs;

    inline cokernal_workspace * kern0_i(int i_) {return ckrn_wkspcs[i_];}
    inline int nvK(int nK_, int i_vK_[])
    {
      int nvK_out = 0;
      for (size_t i = 0; i < nK_; i++) nvK_out += ckrn_spcs[i_vK_[i]]->nvec_use;
      return nvK_out;
    }
    inline int nvK()
    {
      int nvK_out = 0;
      for (size_t i = 0; i < nset; i++) nvK_out += ckrn_spcs[i]->nvec_use;
      return nvK_out;
    }

    cokernal_family ** spawn_cokernal_families(int **i_vspc0_fam_[],int *i_vspc0_set_[],k_medoids_results &Kf_res_,bool verbose_=true)
    {
      cokernal_family ** const ckfams = new cokernal_family*[Kf_res_.Kmed];
      for (size_t k = 0, kk = 0; k < Kf_res_.Kmed; k++)
      {
        ckfams[k] = new cokernal_family(Kf_res_.npts_med[k],Kf_res_.i_meds[k],Kf_res_.ipts_med[k],ckrn_spcs,ckrn_wkspcs);
        i_vspc0_fam_[k] = i_vspc0_set_+kk;
        for (size_t imem = 0; imem < ckfams[k]->nckrn; kk++, imem++)
          i_vspc0_set_[kk] = (ckfams[k]->ckrn_spcs[imem])->ivec_0Vkrn;
        if (verbose_) ckfams[k]->print_details(k);
      }
      return ckfams;
    }

    int collapse_cokernals(LD_vspace_measure &msr_,bool *frdc_succ_,int *nsat_ckrn_,int *ij_prs_,
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

class LD_cokernal_policy
{

  public:

    LD_cokernal_policy(LD_vspace_measure &msr_,bool wdistance_): wdistance(wdistance_), msr(msr_) {}
    ~LD_cokernal_policy() {}

    const bool wdistance;
    LD_vspace_measure &msr;

    virtual double ** init_cokernal_collapse(Jet_function_vector_space &jfvs_,bool verbose_=true) = 0;
    virtual int reduce_cokernal_bundle(cokernal_bundle &cokern_,bool verbose_=true) = 0;

};

struct cokernal_refinement
{
  cokernal_refinement(Jet_function_vector_space &jfvs_,LD_cokernal_policy &pol_,bool verbose_=true);
  ~cokernal_refinement()
  {
    delete [] ckrn_V0_assign;
    delete [] ckfm_V0_assign;

    for (size_t i = 0; i < nfam; i++) delete ckfams[i];
    delete [] ckfams;
  }

  const int vlen_full,
            nset0;

  int * const ckrn_V0_assign,
      * const ckfm_V0_assign;

  cokernal_bundle cokern;
  cokernal_space ** const ckrn_spcs = cokern.ckrn_spcs;

  const int nvKf,
            nset,
            nfam;

  int ** const i_vspc0_set,
      *** const i_vspc0_fam;

  k_medoids_results Kf_res;
  // LD_svd  Kf_svd;

  cokernal_family ** const ckfams;

  void print_refinement_diagnostics(Jet_function_vector_space &jfvs_,LD_cokernal_policy &pol_)
  {
    LD_vector_space ** const vspcs0 = jfvs_.Vbndle0.Vspaces;
    printf("(cokernal_refinement::print_refinement_diagnostics) k = %d kernal families of %d cokernals:\n",
      nfam, nset);
    for (int k = 0; k < nfam; k++)
    {
      cokernal_family &ckfk = *(ckfams[k]);
      const int imed_k = Kf_res.i_meds[k];
      int * const ickrn_fam_k = Kf_res.ipts_med[k];
      printf("\n(k=%d) imed = %d, nckrn = %d.\n  i_ckrns : ", k, imed_k, ckfk.nckrn);

      for (int i = 0; i < ckfams[k]->nckrn; i++)
        printf("%s%d%s ",(ickrn_fam_k[i]==imed_k)?("["):(""),ickrn_fam_k[i],(ickrn_fam_k[i]==imed_k)?("]"):(""));

      printf("--> iV0 : ");
      for (int i = 0; i < ckfams[k]->nckrn; i++)
        for (int j = 0; j < ckfams[k]->nkrn0_ckrn[i]; j++)
          printf("%d ", ckfams[k]->ckrn_spcs[i]->ivec_0Vkrn[j]);

      printf("\n  (tot. nV0 = %d ; nvec0=%d -> %d=nvecf )\n  [ ickrn : iV0 | (nV0, nvec0->nvecf) ] - \n",
              ckfams[k]->nkrn0(), ckfams[k]->nvec0(), ckfams[k]->nvecf()
            );
      for (int i = 0; i < ckfams[k]->nckrn; i++)
      {
        printf("  [ %d : ", ickrn_fam_k[i]);
        for (int j = 0; j < ckfams[k]->ckrn_spcs[i]->n_0Vkrn; j++)
          printf("%d ", ckfams[k]->ckrn_spcs[i]->ivec_0Vkrn[j]);
        printf(" | (%d, %d->%d) ] %s\n",
          ckfams[k]->ckrn_spcs[i]->n_0Vkrn, ckfams[k]->nvec0_i(i), ckfams[k]->ckrn_spcs[i]->nvec_use,
                                          (ickrn_fam_k[i]==imed_k)?(" (<-MED)"):(""));
      }
    }
    printf("\nFinal kernal assignments:\n");
    for (int ikrn = 0; ikrn < nset0; ikrn++)
      printf("krn %d -> ckrn %d (ckrn_fam %d) \n",ikrn,ckrn_V0_assign[ikrn],ckfm_V0_assign[ikrn]);
  }
};

#endif
