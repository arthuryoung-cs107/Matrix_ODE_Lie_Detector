#ifndef LD_SSCOK_HH
#define LD_SSCOK_HH

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
      * const cokern_dimvec,
      &nset = cokern_dimvec[0],
      &kSC = cokern_dimvec[1],
      &nred_succ = cokern_dimvec[2];

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

    cokernal_space ** const ckrn_spcs;

    int collapse_cokernals(LD_vspace_measure &msr_,int gen_spawn_,bool *frdc_succ_,int *nsat_ckrn_,int *ij_prs_,
      double *w_ckrn_,double *vvec_ckrn_,int *isat_ckrn_,bool verbose_=true)
    {
      const int vl2 = vlen*vlen;
      int i_i, i_j,
          ij_prnts[nred_succ][2],
          gen_iprnt[nred_succ],
          gen_jprnt[nred_succ],
          nvec_use_iprnt[nred_succ],
          nvec_use_jprnt[nred_succ],
          nV0krn_iprnt[nred_succ],
          nV0krn_jprnt[nred_succ];
      cokernal_space * ckrn_rdc[nred_succ];
      for (size_t i = 0, irdc = 0, i_ij = 0, iw = 0, iv = 0; irdc < nred_succ; irdc += (int)(frdc_succ_[i++]), i_ij+=2, iw+=vlen, iv+=vl2)
        if (frdc_succ_[i])
        {
          cokernal_space  * const ckrn_i_old = ckrn_spcs[i_i = ij_prnts[irdc][0] = ij_prs_[i_ij]],
                          * &ckrn_j_old = ckrn_spcs[i_j = ij_prnts[irdc][1] = ij_prs_[i_ij+1]];
          gen_iprnt[irdc] = ckrn_i_old->gen_spawn; gen_jprnt[irdc] = ckrn_j_old->gen_spawn;
          nvec_use_iprnt[irdc] = ckrn_i_old->nvec_use; nvec_use_jprnt[irdc] = ckrn_j_old->nvec_use;
          nV0krn_iprnt[irdc] = ckrn_i_old->n_0Vkrn; nV0krn_jprnt[irdc] = ckrn_j_old->n_0Vkrn;

          if ( ckrn_i_old->gen_spawn ) // if NOT an original kernal, can be overwritten
            ckrn_spcs[i_i] = new cokernal_space(gen_spawn_,ckrn_i_old,ckrn_j_old); // spawn new cokernal space by overwriting ckrn_i and its workspace
          else // if an original kernal, extra workspace must be overwritten or allocated
          {
            int i_check = 0;
            while ( (i_check<nxtra_wkspc) && (factive_ckrn_wkspc_xtra[i_check]) ) i_check++;
            if (i_check==nxtra_wkspc) // if all allocated workspaces are active, allocate a new one.
              ckrn_wkspcs_xtra[nxtra_wkspc++] = new cokernal_workspace(nset0+nxtra_wkspc,factive_ckrn_wkspc_xtra[nxtra_wkspc]=true,vlen);

            // spawn new cokernal space at ckrn_i using an allocated workspace
            ckrn_spcs[i_i] = new cokernal_space(gen_spawn_,ckrn_wkspcs_xtra[i_check],ckrn_i_old,ckrn_j_old);
          }
          // set the underlying vector space
          (ckrn_rdc[irdc] = ckrn_spcs[i_i])->init_cokernal_space(vlen,nsat_ckrn_[i],w_ckrn_+iw,vvec_ckrn_+iv,isat_ckrn_+iw);
        }

      const int nset_new = nset-nred_succ;
      int nV_setvec_new[nset_new];
      double ** Vtns_new[nset_new];
      LD_vector_space * vspcs_new[nset_new];
      for (size_t i = 0, ickrn_new = 0; ickrn_new < nset_new; i++)
        if (ckrn_spcs[i] != NULL)
        {
          ckrn_spcs[i]->ickrn = ickrn_new;
          nV_setvec_new[ickrn_new] = ckrn_spcs[i]->nvec_use;
          Vtns_new[ickrn_new] = ckrn_spcs[i]->Vmat;
          vspcs_new[ickrn_new] = ckrn_wkspcs[ckrn_spcs[i]->iwkspc]->vspc;
          ckrn_spcs[ickrn_new++] = ckrn_spcs[i];
        }
      for (size_t i = nset_new; i < nset0; i++) ckrn_spcs[i] = NULL;

      if (verbose_)
      {
        printf("  (cokernal_bundle::collapse_cokernals) reassigned %d cokernals:\n",nset_new);
        for (size_t i = 0; i < nset_new; i++)
        {
          cokernal_space &cki = *(ckrn_spcs[i]);
          printf("    (ickrn=%d, iwkspc=%d): nvec=%d, n_0Vkrn=%d, (gen,iprnt,jprnt)=(%d,%d,%d)",
            cki.ickrn, cki.iwkspc,
            cki.nvec_use, cki.n_0Vkrn,
            cki.gen_spawn,cki.ickrn_parent,cki.jckrn_parent);
          if (cki.n_0Vkrn > 1)
          {
            printf(", i_0Vkrn: ( ");
            for (size_t iV0 = 0; iV0 < cki.n_0Vkrn; iV0++) printf("%d ", cki.ivec_0Vkrn[iV0]);
            printf(")\n");
          }
          else printf(".\n");
        }

        printf("    New cokernals: \n");
        for (size_t irdc = 0; irdc < nred_succ; irdc++)
        {
          cokernal_space &ckrni = *(ckrn_rdc[irdc]);
          printf("      (irdc=%d) [iVprnt=%d, jVprnt=%d -> ickrn=%d] : gen [%d,%d -> %d] ; nV [%d,%d -> %d] ; nV0krn [%d,%d -> %d]. iV0krn = ( ",
            irdc,
            ckrni.ickrn_parent,ckrni.jckrn_parent,ckrni.ickrn,
            gen_iprnt[irdc],gen_jprnt[irdc],ckrni.gen_spawn,
            nvec_use_iprnt[irdc],nvec_use_jprnt[irdc],ckrni.nvec_use,
            nV0krn_iprnt[irdc],nV0krn_jprnt[irdc],ckrni.n_0Vkrn);
          for (size_t i = 0; i < ckrni.n_0Vkrn; i++) printf("%d ", ckrni.ivec_0Vkrn[i]);
          printf(").\n");
        }
      }
      // msr_.update_distances(nred_succ,ij_prnts[0],vspcs_new,verbose_);
      msr_.init_distances(Vtns_new,nV_setvec_new,nset_new,verbose_);

      return nset_new;
    }
};

class Jet_function_vector_space
{
    LD_observations_set &Sobs;
    function_space &fspc;
    LD_encoder &enc;
    LD_vspace_measure &msr;

    const int nset0,
              vlen_full;

  public:

    int nset,
        kSC,
        nred_succ,
        * const cokern_dimvec = &nset;

    LD_encoding_bundle Acode;
    LD_vector_bundle Vbndle0;
    LD_svd_bundle svd0; // shared data with Vbndle0

    LD_vspace_record &rec0 = Vbndle0.rec;

    Jet_function_vector_space(LD_observations_set &Sobs_,function_space &fspc_,LD_encoder &enc_,LD_vspace_measure &msr_):
      Sobs(Sobs_), fspc(fspc_), enc(enc_), msr(msr_), nset0(Sobs_.ncrvs_tot), vlen_full(fspc_.ndof_full),
      nset(nset0),
      Acode(LD_encoding_bundle(nset0,vlen_full,Sobs_.npts_per_crv,enc_.ncod)),
      Vbndle0(LD_vector_bundle(nset0,vlen_full)),
      svd0(LD_svd_bundle(Vbndle0))
      {}
    ~Jet_function_vector_space() {}

    template <class BSE> Jet_function_vector_space(LD_observations_set &Sobs_,function_space &fspc_,LD_encoder &enc_,LD_vspace_measure &msr_,
      BSE **bases_, bool normalize_, bool verbose_=true):
      Jet_function_vector_space(Sobs_,fspc_,enc_,msr_)
    {
      LD_encoder::encode_bundle<BSE>(Acode,Sobs,bases_,enc,normalize_,verbose_);
      svd0.compute_Acode_curve_svds(Acode,verbose_);
        svd0.evaluate_iscl_ranks(2*vlen_full);
        svd0.print_details("svd0");

      svd0.set_Vspaces_nullspc();
        rec0.print_subspace_details("rec0",false,false);

      cokernal_bundle cokern(Vbndle0,cokern_dimvec,svd0.Smat);
      msr.init_distances(Vbndle0,nset,verbose_);

      int generation = 0;
      double t0 = LD_threads::tic();
      while (iterate_nullspace_cokernals(cokern,msr,generation,verbose_)) generation++;
      double work_time = LD_threads::toc(t0);

      if (verbose_)
        printf("(Jet_function_vector_space::Jet_function_vector_space) collapsed %d kernal spaces (%d x %d) into %d in %.1f seconds (%d threads)\n",
          nset0, vlen_full, vlen_full, nset,
          work_time, LD_threads::numthreads() );

      LD_encoded_matrix ** const Aenc_mats = Acode.Acodes;
      LD_vector_space ** const krn0_spcs = Vbndle0.Vspaces;
      cokernal_space ** const ckrn_spcs = cokern.ckrn_spcs;

      int nKf = 0,
          ckrn_V0_assign[nset0];
      for (size_t ickrn = 0; ickrn < nset; nKf+=ckrn_spcs[ickrn++]->nvec_use)
        for (size_t iset = 0; iset < ckrn_spcs[ickrn]->n_0Vkrn; iset++)
          ckrn_V0_assign[ckrn_spcs[ickrn]->ivec_0Vkrn[iset]] = ickrn;

      int mrow_Amat[nset0],
          nvec_krn0[nset0],
          nvec_ckrn[nset0];
      double  ** A_mats[nset0],
              ** V_krn0[nset0],
              ** V_ckrn[nset0],
              krn0_stats[nset0][4],
              ckrn_stats[nset0][4];

      #pragma omp parallel for
      for (size_t iset = 0; iset < nset0; iset++)
      {
        const int iickrn = ckrn_V0_assign[iset],
                  mr_Amat_i = mrow_Amat[iset] = Aenc_mats[iset]->nrow,
                  nV_krn0_i = nvec_krn0[iset] = krn0_spcs[iset]->nvec_use,
                  nV_ckrn_i = nvec_ckrn[iset] = ckrn_spcs[iickrn]->nvec_use,
                  len_Akrn0_i = mr_Amat_i*nV_krn0_i,
                  len_Ackrn_i = mr_Amat_i*nV_ckrn_i;
        double  Ai_wkspc[len_Akrn0_i],
               *A_krn0_i[mr_Amat_i],
               *A_ckrn_i[mr_Amat_i];
        for (size_t i = 0, ii_krn0 = 0, ii_ckrn = 0; i < mr_Amat_i; i++, ii_krn0+=nV_krn0_i, ii_ckrn+=nV_ckrn_i)
        {
          A_krn0_i[i] = Ai_wkspc + ii_krn0;
          A_ckrn_i[i] = Ai_wkspc + ii_ckrn;
        }

        LD_linalg::A__B_CT( A_krn0_i,
                            A_mats[iset]=Aenc_mats[iset]->Amat,
                            V_krn0[iset]=krn0_spcs[iset]->Vmat,
                            mr_Amat_i,vlen_full,nV_krn0_i,true);
        LD_linalg::abs_vec<double>(A_krn0_i[0],len_Akrn0_i);
        LD_linalg::comp_Tvec_stats<double>(krn0_stats[iset],A_krn0_i[0],len_Akrn0_i);

        LD_linalg::A__B_CT( A_ckrn_i,
                            A_mats[iset],
                            V_ckrn[iset]=ckrn_spcs[iickrn]->Vmat,
                            mr_Amat_i,vlen_full,nV_ckrn_i,true);
        LD_linalg::abs_vec<double>(A_ckrn_i[0],len_Ackrn_i);
        LD_linalg::comp_Tvec_stats<double>(ckrn_stats[iset],A_ckrn_i[0],len_Ackrn_i);
      }

      printf("  original kernal vs. cokernal inner product performance\n");
      for (size_t iset = 0; iset < nset0; iset++)
        printf("    (crv %d, iickrn %d, krn0 vs. ckrn):  nV (%d, %d) ; min (%.2e, %.2e) ; med (%.2e, %.2e) ; avg (%.2e, %.2e) ; max (%.2e, %.2e) \n",
                    iset, ckrn_V0_assign[iset], krn0_spcs[iset]->nvec_use, ckrn_spcs[ckrn_V0_assign[iset]]->nvec_use,
                    krn0_stats[iset][0], ckrn_stats[iset][0],
                    krn0_stats[iset][1], ckrn_stats[iset][1],
                    krn0_stats[iset][2], ckrn_stats[iset][2],
                    krn0_stats[iset][3], ckrn_stats[iset][3]
                  );

      LD_svd  Kf_svd(vlen_full,nKf);
      for (size_t iset = 0, jK = 0; iset < nset; iset++)
        for (size_t lcol = 0; lcol < ckrn_spcs[iset]->nvec_use; lcol++, jK++)
          for (size_t i = 0; i < vlen_full; i++) Kf_svd.Umat[i][jK] = ckrn_spcs[iset]->Vmat[lcol][i];

      Kf_svd.decompose_U(vlen_full,nKf);

      Kf_svd.print_result("Kf");
      int rank_Kf = Kf_svd.rank();
      double  nuke_norm_Kf = Kf_svd.norm_nuke(),
              effrank_Kf = nuke_norm_Kf/Kf_svd.s0;
      printf("nKf = %d, rank_Kf = %d, nuke_norm_Kf = %.2f (%.2f), effective rank_Kf = %.2f (%.2f)\n",
        nKf,rank_Kf,
        nuke_norm_Kf, nuke_norm_Kf/((double)nKf),
        effrank_Kf,effrank_Kf/((double)nKf));
    }

    int iterate_nullspace_cokernals(cokernal_bundle &cokern_,LD_vspace_measure &msr_,int gen_,bool verbose_=true)
    {
      k_medoids_package kmed(msr_.dsym,nset);
      const int kSC0 = kSC = kmed.comp_kSC_medoids(2,verbose_);
      if (iterate_nullspace_cokernals(cokern_,msr_,gen_,kmed,verbose_)) return nred_succ;
      else
      {
        if (kSC0 > 2)
        {
          int kSCm1 = kSC0-1; kSC = kmed.comp_kSC_krange_medoids(2,kSCm1,verbose_);
          while ( (kSCm1 >= 2) && !(iterate_nullspace_cokernals(cokern_,msr_,gen_,kmed,verbose_)) )
            if (kSCm1 == 2) break;
            else kSC = kmed.comp_kSC_krange_medoids(2,--kSCm1,verbose_);
          if (nred_succ) return nred_succ;
        }
        printf("(Jet_function_vector_space::iterate_nullspace_cokernals) gen %d - irreducable clusters (nset=%d). Attempting to consolidate remaining cokernal spaces\n",gen_,nset);
        // if all else fails, just try consolidating one pair
        kmed.set_one_medoid(); kSC = 1;
        return iterate_nullspace_cokernals(cokern_,msr_,gen_,kmed,verbose_);
      }
    }

    int iterate_nullspace_cokernals(cokernal_bundle &cokern_,LD_vspace_measure &msr_,int gen_,k_medoids_package &kmed_,bool verbose_)
    {
      int nmem_clst_vec[kSC],
          imem_clst_vec[nset],
          *imem_clst_mat[kSC];

      kmed_.assign_clusters(nmem_clst_vec,imem_clst_vec,imem_clst_mat,kSC,verbose_);

      cokernal_space ** const ckrn_spcs = cokern_.ckrn_spcs;
      bool f_reduce_success[kSC];
      int nsvd_cmps = nred_succ = 0,
          nuld_ckrn[kSC],
          isat_ckrn[kSC][vlen_full],
          i_n_pairs_ckrn[kSC][2],
          ij_pairs_ckrn[kSC][2],
          ij_nn_pairs[kSC][2];
      double  dpairs_ckrn[kSC][2],
              svec_ckrns[kSC][vlen_full],
              Vmat_ckrns[kSC][vlen_full][vlen_full],
              t0 = LD_threads::tic();
      #pragma omp parallel reduction(+:nsvd_cmps,nred_succ)
      {
        LD_svd  Vsvd_t(2*vlen_full,vlen_full);
        int i_setk,j_setk;
        double  * const svec_V_t = Vsvd_t.svec,
                ** const Vmat_V_t = Vsvd_t.Vmat,
                ** const Umat_V_t = Vsvd_t.Umat;
        #pragma omp for
        for (size_t k = 0; k < kSC; k++)
        {
          if (f_reduce_success[k] = (nmem_clst_vec[k] > 1))
          {
            const int npairs_k = i_n_pairs_ckrn[k][1] = (nmem_clst_vec[k])*(nmem_clst_vec[k]-1)/2;
            int ipair = i_n_pairs_ckrn[k][0] = 0,
                ipairs_k[npairs_k],
                iipairs_k[npairs_k][2];
            double dpairs_k[npairs_k];
            kmed_.comp_sort_cluster_distances(dpairs_k,ipairs_k,iipairs_k[0],imem_clst_mat[k],nmem_clst_vec[k]);
            ij_nn_pairs[k][0] = iipairs_k[0][0]; ij_nn_pairs[k][1] = iipairs_k[0][1];
            dpairs_ckrn[k][1] = dpairs_k[0];
            do
            {
              i_setk = iipairs_k[ipair][0]; j_setk = iipairs_k[ipair][1];

              for (size_t i = 0; i < vlen_full; i++)
                for (size_t j = 0; j < vlen_full; j++)
                  Umat_V_t[i][j] = (ckrn_spcs[i_setk]->wvec[i])*(ckrn_spcs[i_setk]->Vmat_data[i][j]);

              for (size_t i = 0; i < vlen_full; i++)
                for (size_t j = 0; j < vlen_full; j++)
                  Umat_V_t[i+vlen_full][j] = (ckrn_spcs[j_setk]->wvec[i])*(ckrn_spcs[j_setk]->Vmat_data[i][j]);

              Vsvd_t.decompose_U(2*vlen_full,vlen_full); nsvd_cmps++;

              if (f_reduce_success[k] = (bool)( nuld_ckrn[k] = (vlen_full-Vsvd_t.rank()) ))
              {
                nred_succ++;
                for (size_t i = 0, isat = vlen_full-nuld_ckrn[k]; i < nuld_ckrn[k]; i++, isat++) isat_ckrn[k][i] = isat;
                for (size_t i = nuld_ckrn[k], insat = 0; i < vlen_full; i++, insat++) isat_ckrn[k][i] = insat;
                i_n_pairs_ckrn[k][0] = ipair;
                ij_pairs_ckrn[k][0] = i_setk; ij_pairs_ckrn[k][1] = j_setk;
                dpairs_ckrn[k][0] = dpairs_k[ipair];
                for (size_t i = 0; i < vlen_full; i++) svec_ckrns[k][i] = svec_V_t[i];
                for (size_t i = 0; i < vlen_full; i++)
                  for (size_t j = 0; j < vlen_full; j++)
                    Vmat_ckrns[k][i][j] = Vmat_V_t[j][i];
                break;
              }
              else if (ipair == (npairs_k-1)) break;
              else ipair++;
            } while (true);
          }
          else
          {
            ij_nn_pairs[k][0] = ij_nn_pairs[k][1] = kmed_.i_meds[k];
            dpairs_ckrn[k][0] = dpairs_ckrn[k][1] = 0.0;
            i_n_pairs_ckrn[k][0] = -1;
            i_n_pairs_ckrn[k][1] = 0;
          }
        }
      }
      double work_time = LD_threads::toc(t0);

      if (verbose_)
      {
        printf("  (Jet_function_vector_space::iterate_nullspace_cokernals) gen %d: reduced %d (of kSC = %d) cokernals via %d svds in %.2f seconds (%d threads).\n  Reduction details: \n",
          gen_,nred_succ, kSC,
          nsvd_cmps,work_time,LD_threads::numthreads() );
        for (size_t k = 0; k < kSC; k++)
        {
          int imd=kmed_.i_meds[k], inn=ij_nn_pairs[k][0], jnn=ij_nn_pairs[k][1];
          printf("    k=%d: imed=%d (nv=%d), nclst=%d. inn=%d, jnn=%d (nv = %d, %d ; dij=%.2e) of %d pairs. Intersection %s",
            k,imd,ckrn_spcs[imd]->nvec_use,nmem_clst_vec[k],
            inn,jnn,ckrn_spcs[inn]->nvec_use,ckrn_spcs[jnn]->nvec_use,dpairs_ckrn[k][1],
            i_n_pairs_ckrn[k][1],
            (f_reduce_success[k])?( (i_n_pairs_ckrn[k][0])?("FOUND: "):("FOUND (default).\n") ):("NOT found.\n"));
          if (f_reduce_success[k] && i_n_pairs_ckrn[k][0])
            printf("i_cok=%d, j_cok=%d (d=%.2e, %d).\n",
              ij_pairs_ckrn[k][0],ij_pairs_ckrn[k][1],dpairs_ckrn[k][0],i_n_pairs_ckrn[k][0]);
        }
      }

      if (nred_succ)
        nset = cokern_.collapse_cokernals(msr_,gen_+1,
          f_reduce_success,nuld_ckrn,ij_pairs_ckrn[0],svec_ckrns[0],Vmat_ckrns[0][0],isat_ckrn[0],verbose_);

      return nred_succ;
    }

    template <class BSE> void encode_decompose_bundle(BSE **bases_,bool normalize_,bool verbose_=true)
    {
      LD_encoder::encode_bundle<BSE>(Acode,Sobs,bases_,enc,normalize_,verbose_);
      svd0.compute_Acode_curve_svds(Acode,verbose_);
        if (verbose_) svd0.print_details("svd0");
      svd0.set_Vspaces_nullspc();
        if (verbose_) rec0.print_subspace_details("rec0",false,false);
    }
    template <class EVL,class BSE> void evaluate_Vbndle0(EVL &evl_, BSE **bases_, bool reset_=true, bool verbose_=true)
    {
      if (reset_) Vbndle0.reset_Vspaces();
      evl_.template evaluate_vspace_record<BSE>(rec0,evl_.Sset,bases_,verbose_); Vbndle0.set_Vspaces();
        rec0.print_subspace_details("rec0",false,false);
    }

    // inline int * inn_k(int k_) {return inn_space + (2*k_);}
    // inline void ijnn_k(int k_, int &i_,int &j_)
    //   {int * const ijnn_k_vec = inn_space + (2*k_); i_ = ijnn_k_vec[0]; j_ = ijnn_k_vec[1];}

  protected:

};

class cokernal_reduction
{

  const int nset0;
  LD_vector_bundle  &bndle0;
  LD_vector_space ** const Vckrns;

  public:

    cokernal_reduction(LD_vector_bundle &bndle0_,int *nppc_,ode_solution **s0_):
      nset0(bndle0_.nspc), bndle0(bndle0_), Vckrns(new LD_vector_space*[bndle0_.nspc]),
      nobs(LD_linalg::sum_vec<int>(nppc_,nset0)), nset(nset0),
      nvrs(new int[nset0]), nppc(new int[nset0]), inn_space(new int[nset0]), sf(new ode_solution*[nobs])
      {
        for (size_t i = 0; i < nset0; i++) Vckrns[i] = NULL;
        set_sols(nset,nppc_,s0_);
      }
    ~cokernal_reduction()
    {
      delete [] nvrs; delete [] nppc; delete [] inn_space; delete [] sf;
      for (size_t i = 0; i < nset0; i++) if (Vckrns[i] != NULL) delete Vckrns[i];
      delete [] Vckrns;
    }

    const int nobs;
    int nset;

    int init_nearest_neighbors(LD_vspace_measure &msr_,LD_vector_bundle &Vb_,int nset_=0,bool verbose_=true)
    {
      if (nset_) nset = nset_;

      msr_.init_distances(Vb_,nset,verbose_);
      k_medoids_package kmed(msr_.dsym,nset);

      const int kSC = kmed.comp_kSC_medoids(2,verbose_);
      int nmem_clst0_vec[kSC],
          imem_clst0_vec[nset],
          *imem_clst0_mat[kSC];

      kmed.assign_clusters(nmem_clst0_vec,imem_clst0_vec,imem_clst0_mat,kSC);

      double dnn_clst[kSC];
      for (size_t k = 0; k < kSC; k++)
        // kmed.get_nearest_neighbors(inn_k(k),imem_clst0_mat[k],nmem_clst0_vec[k]);
        dnn_clst[k] = kmed.get_nearest_neighbors(inn_k(k),imem_clst0_mat[k],nmem_clst0_vec[k]);

      printf("clusters\n");
      for (size_t k = 0; k < kSC; k++)
      {
        printf("(%d, %d): ", k, kmed.i_meds[k]);
        for (size_t i = 0; i < nmem_clst0_vec[k]; i++) printf("%d ", imem_clst0_mat[k][i]);
        printf("\n");
      }
      printf("nearest neigbors\n");
      for (size_t k = 0; k < kSC; k++) printf("(%d,%d): %d, %d, %e\n", k,kmed.i_meds[k],
        (inn_k(k))[0],(inn_k(k))[1],dnn_clst[k]);

      return kSC;
    }

    inline int * inn_k(int k_) {return inn_space + (2*k_);}
    inline void ijnn_k(int k_, int &i_,int &j_)
      {int * const ijnn_k_vec = inn_space + (2*k_); i_ = ijnn_k_vec[0]; j_ = ijnn_k_vec[1];}

    inline void set_sols(int nset_,int nsol_per_set_[],ode_solution ** sols_)
    {
      for (size_t ic = 0, is = 0; ic < nset_; ic++)
      {
        nppc[ic] = nsol_per_set_[ic];
        for (size_t ip = 0; ip < nppc[ic]; ip++, is++) sf[is]=sols_[is];
      }
    }

    template <class EVL,class BSE> void reduce_kernals(int k_,LD_vector_bundle &Vb_,ode_solspc_setbundle &Sb_,
      EVL &evl_,BSE **bases_,bool verbose_=true)
    {
      init_k_Vspaces(k_);
      const int vlen = Vb_.vlen_full;
      int nvec_acc = 0,
          nsol_acc = 0,
          nred_succ = 0;
      #pragma omp parallel num_threads(1) reduction(+:nvec_acc,nsol_acc,nred_succ)
      {
        LD_svd svd_t(vlen,vlen);
        BSE &bse_t = *(bases_[LD_threads::thread_id()]);
        EVL evl_t(evl_,1);
        bool  sat_k[vlen],
              * const sat_t = evl_t.sat_flags;
        int i_setk, j_setk,
            nvec_i, nvec_j,
            nsol_i, nsol_j,
            &nvec_evl = evl_t.nvec_evl,
            &nsol_evl = evl_t.nsol_evl;
        double  ** const Umat_t = svd_t.Umat,
                ** const Vmat_t = svd_t.Vmat;
        #pragma omp for
        for (size_t k = 0; k < 1; k++)
        // for (size_t k = 0; k < k_; k++)
        {
          ijnn_k(k,i_setk,j_setk);
          LD_vector_space &Vspc_i = *(Vb_.Vspaces[i_setk]),
                          &Vspc_j = *(Vb_.Vspaces[j_setk]);
          nvec_acc += (nvec_evl = (nvec_i = Vspc_i.nvec_use) + (nvec_j = Vspc_j.nvec_use));
          nsol_acc += (nsol_evl = (nsol_i = Sb_.nsol_per_set[i_setk]) + (nsol_j = Sb_.nsol_per_set[j_setk]));

          printf("k = %d: i = %d (%d sols), j = %d (%d sols)\n", k, i_setk, nsol_i, j_setk, nsol_j);

          LD_io::write_Tmat<int>("iV_Vmati.lddat",Vb_.iV_spcmat+i_setk,1,vlen);
          LD_io::write_Tmat<int>("iV_Vmatj.lddat",Vb_.iV_spcmat+j_setk,1,vlen);
          LD_io::write_Tmat<double>("Vmat_data_i.lddat",Vb_.rec.Vtns_data[i_setk],vlen,vlen);
          LD_io::write_Tmat<double>("Vmat_data_j.lddat",Vb_.rec.Vtns_data[j_setk],vlen,vlen);
          LD_io::write_Tmat<double>("Vmat_full_i.lddat",Vspc_i.Vmat,vlen,vlen);
          LD_io::write_Tmat<double>("Vmat_full_j.lddat",Vspc_j.Vmat,vlen,vlen);
          LD_io::write_Tmat<double>("Vmat_i.lddat",Vspc_i.Vmat,nvec_i,vlen);
          LD_io::write_Tmat<double>("Vmat_j.lddat",Vspc_j.Vmat,nvec_j,vlen);

          cokernal_reduction::load_Psi(Umat_t,Vspc_i.Vmat,nvec_i,Vspc_j.Vmat,nvec_j,vlen);
          LD_io::write_Tmat<double>("Umat_t.lddat",Umat_t,vlen,nvec_evl);
          svd_t.decompose_U(vlen,nvec_evl); svd_t.print_result("Psi");
          LD_io::write_Tmat<double>("Vmat_t.lddat",Vmat_t,nvec_evl,nvec_evl);
          // cokernal_reduction::load_W(evl_t.Vmat = Vckrns[k]->Vmat,Vmat_t,Vspc_i.Vmat,nvec_i,Vspc_j.Vmat,nvec_j,vlen);
          nvec_evl = cokernal_reduction::load_W(evl_t.Vmat = Vckrns[k]->Vmat,Vmat_t,Vspc_i.Vmat,nvec_i,Vspc_j.Vmat,nvec_j,vlen);
          // evl_t.Vmat = Vspc_i.Vmat; nvec_evl = nvec_i;

          ode_solution ** sols_i = Sb_.get_sol_subset_i(i_setk);
          int nsat_setk = 0;
          LD_linalg::fill_vec<bool>(sat_k,nvec_evl,true);
          for (size_t isol = 0; isol < nsol_i; isol++) // evaluate first set of solutions
          {
            if (evl_t.nsat_eval_condition(sat_t,*(sols_i[isol]),bse_t))
            {
              nsat_setk = 0;
              for (size_t iV = 0; iV < nvec_evl; iV++)
                nsat_setk += (int)(sat_k[iV] = (sat_k[iV]) && (sat_t[iV]));
              printf("i (%d of %d) nsat_setk = %d pass\n", isol, nsol_i, nsat_setk);
            }
            else LD_linalg::fill_vec<bool>(sat_k,nvec_evl,nsat_setk = 0);
            if (!nsat_setk) break; // if none pass
          }
          // if (nvrs[k] = nsat_setk) // then evaluate second set of solutions
          // {
          //   printf("nsat_setk = %d initially pass\n", nsat_setk);
          if (true) // then evaluate second set of solutions
          {
            LD_linalg::fill_vec<bool>(sat_k,nvec_evl,true);

            ode_solution ** sols_j = Sb_.get_sol_subset_i(j_setk);
            for (size_t jsol = 0; jsol < nsol_j; jsol++)
            {
              if (evl_t.nsat_eval_condition(sat_t,*(sols_j[jsol]),bse_t))
              {
                nsat_setk = 0;
                for (size_t iV = 0; iV < nvec_evl; iV++)
                  nsat_setk += (int)(sat_k[iV] = (sat_k[iV]) && (sat_t[iV]));
                printf("j (%d of %d) nsat_setk = %d pass\n", jsol, nsol_j, nsat_setk);
              }
              else LD_linalg::fill_vec<bool>(sat_k,nvec_evl,nsat_setk = 0);
              if (!nsat_setk) break; // if none pass
            }
            if (nvrs[k] = nsat_setk) nred_succ++;
          }
        }
      }
      nred_success = nred_succ;
      printf("nred_success = %d\n", nred_success);
      LD_linalg::print_x("nvrs", nvrs, k_);
    }

  protected:

    int * const nvrs,
        * const nppc,
        * const inn_space;

    int nred_success;

    ode_solution ** const sf;


    inline void init_k_Vspaces(int k_)
    {
      for (size_t k = 0; k < k_; k++)
        if (Vckrns[k] == NULL) Vckrns[k] = new LD_vector_space(bndle0.vlen_full);
        else Vckrns[k]->reset_Vspace();
    }
    // inline void init_k_Vspaces(int k_,LD_vector_bundle &Vb_)
    // {
    //   for (size_t k = 0; k < k_; k++)
    //     if (Vckrns[k] == NULL) Vckrns[k] = new LD_vector_space(*(Vb_.Vspaces[*(inn_k(k))]),true);
    //     else Vckrns[k]->set_Vspace(*(Vb_.Vspaces[*(inn_k(k))]));
    // }

    static void load_Psi(double **Psi_,double **Vi_,int nvi_,double **Vj_, int nvj_, int vlen_)
    {
      for (size_t ii = 0; ii < vlen_; ii++)
      {
        for (size_t jj = 0; jj < nvi_; jj++) Psi_[ii][jj] = Vi_[jj][ii];
        for (size_t jj = 0, ll = nvi_; jj < nvj_; jj++, ll++) Psi_[ii][ll] = -Vj_[jj][ii];
      }
    }

    // static void load_W(double **Wm_,double **VPsi_,double **Vi_,int nvi_,double **Vj_, int nvj_, int vlen_)
    // {
    //   const int N_Psi = nvi_+nvj_;
    //   for (size_t i = 0; i < N_Psi; i++)
    //     for (size_t j = 0; j < vlen_; j++)
    //       Wm_[i][j] = VPsi_[0][i]*Vi_[0][j];
    //
    //   for (size_t l = 1; l < nvi_; l++)
    //     for (size_t i = 0; i < N_Psi; i++)
    //       for (size_t j = 0; j < vlen_; j++)
    //         Wm_[i][j] += VPsi_[l][i]*Vi_[l][j];
    //
    //   for (size_t l = 0, ll = nvi_; l < nvj_; l++, ll++)
    //     for (size_t i = 0; i < N_Psi; i++)
    //       for (size_t j = 0; j < vlen_; j++)
    //         Wm_[i][j] += VPsi_[ll][i]*Vj_[l][j];
    // }

    static int load_W(double **Wm_,double **VPsi_,double **Vi_,int nvi_,double **Vj_,int nvj_,int vlen_)
    {
      const bool ichoice = nvi_<=nvj_;
      const int nvo = (ichoice)?(nvi_):(nvj_);
      double  ** const Vo = (ichoice)?(Vi_):(Vj_),
              * Km[nvo];

      if (ichoice) for (size_t i = 0; i < nvi_; i++) Km[i] = VPsi_[i] + nvj_;
      else for (size_t i = 0; i < nvj_; i++) Km[i] = VPsi_[i+nvi_] + nvi_;

      LD_io::write_Tmat<double>("Km.lddat",Km,nvo,nvo);

      cokernal_reduction::normalize_submat_cols(Km,nvo);

      LD_io::write_Tmat<double>("Km_nrm.lddat",Km,nvo,nvo);

      for (size_t i = 0; i < nvo; i++)
        for (size_t j = 0; j < vlen_; j++)
          Wm_[i][j] = Km[0][i]*Vo[0][j];
      for (size_t l = 1; l < nvo; l++)
        for (size_t i = 0; i < nvo; i++)
          for (size_t j = 0; j < vlen_; j++)
            Wm_[i][j] += Km[l][i]*Vo[l][j];

      LD_io::write_Tmat<double>("Wm.lddat",Wm_,nvo,vlen_);

      return nvo;
    }

    // static int load_W(double **Wm_,double **VPsi_,double **Vi_,int nvi_,double **Vj_,int nvj_,int vlen_)
    // {
    //   const bool ichoice = nvi_<=nvj_;
    //   const int nvo = (ichoice)?(nvi_):(nvj_);
    //   double  ** const Vo = (ichoice)?(Vi_):(Vj_),
    //           * Ym[nvo],
    //           * Km[nvo];
    //
    //   if (ichoice) for (size_t i = 0; i < nvi_; i++) Km[i] = (Ym[i] = VPsi_[i]) + nvj_;
    //   else for (size_t i = 0; i < nvj_; i++) Km[i] = (Ym[i] = VPsi_[i+nvi_]) + nvi_;
    //
    //   LD_io::write_Tmat<double>("Km.lddat",Km,nvo,nvo);
    //   LD_io::write_Tmat<double>("Ym.lddat",Ym,nvo,nvo);
    //
    //   cokernal_reduction::normalize_submat_cols(Ym,nvo);
    //   cokernal_reduction::normalize_submat_cols(Km,nvo);
    //
    //   LD_io::write_Tmat<double>("Km_nrm.lddat",Km,nvo,nvo);
    //   LD_io::write_Tmat<double>("Ym_nrm.lddat",Ym,nvo,nvo);
    //
    //   for (size_t i = 0; i < nvo; i++)
    //     for (size_t j = 0; j < vlen_; j++)
    //       Wm_[i][j] = Ym[0][i]*Vo[0][j];
    //   for (size_t l = 1; l < nvo; l++)
    //     for (size_t i = 0; i < nvo; i++)
    //       for (size_t j = 0; j < vlen_; j++)
    //         Wm_[i][j] += Ym[l][i]*Vo[l][j];
    //
    //   for (size_t i = 0; i < nvo; i++)
    //     for (size_t j = 0; j < vlen_; j++)
    //       Wm_[i+nvo][j] = Km[0][i]*Vo[0][j];
    //   for (size_t l = 1; l < nvo; l++)
    //     for (size_t i = 0; i < nvo; i++)
    //       for (size_t j = 0; j < vlen_; j++)
    //         Wm_[i+nvo][j] += Km[l][i]*Vo[l][j];
    //
    //   LD_io::write_Tmat<double>("Wm.lddat",Wm_,2*nvo,vlen_);
    //
    //   return 2*nvo;
    // }

    static void normalize_submat_cols(double *Km_[], int nvo_)
    {
      double accv[nvo_];
      for (size_t j = 0; j < nvo_; j++) accv[j] = Km_[0][j]*Km_[0][j];
      for (size_t i = 1; i < nvo_; i++)
        for (size_t j = 0; j < nvo_; j++)
          accv[j] += Km_[i][j]*Km_[i][j];
      for (size_t j = 0; j < nvo_; j++) accv[j] = 1.0/sqrt(accv[j]);
      for (size_t i = 0; i < nvo_; i++)
        for (size_t j = 0; j < nvo_; j++)
          Km_[i][j] *= accv[j];
    }

};

class solution_space_cokernals
{

  int * const nsol_per_set;
  double ** const pts_full;
  ode_solution ** const sols_full;

  public:

    solution_space_cokernals(int ncrv_,int nobs_):
      nset(ncrv_),
      nsol_per_set(new int[ncrv_]), pts_full(new double*[nobs_]), sols_full(new ode_solution*[nobs_])
      {}
    ~solution_space_cokernals() {delete [] nsol_per_set; delete [] pts_full; delete [] sols_full;}

    int nset;

    template <class EVL,class BSE> void compute_solspc_cokernals(LD_vector_bundle &bndl0_,LD_observations_set &Sobs_,
      EVL &evl_,LD_vspace_measure &msr_,BSE **bases_, bool verbose_=true)
    {
      ode_solspc_meta &meta0 = Sobs_.meta;
      cokernal_reduction red(bndl0_,Sobs_.npts_per_crv,Sobs_.sols);
      set_pts_sols(nset = red.nset,Sobs_.npts_per_crv,Sobs_.sols);
      const int nobs = red.nobs,
                kSC0 = init_cokernals<EVL,BSE>( red,bndl0_,
                                                ode_solspc_setbundle(meta0,nobs,pts_full,sols_full,nset,nsol_per_set),
                                                evl_,msr_,bases_,verbose_);
    }

    template <class EVL,class BSE> int init_cokernals(cokernal_reduction &red_,LD_vector_bundle &Vb_,
      ode_solspc_setbundle Sb_,
      EVL &evl_,LD_vspace_measure &msr_,BSE **bases_,bool verbose_=true)
    {
      evl_.template evaluate_vspace_record<BSE>(Vb_.rec,Sb_,bases_,verbose_); Vb_.set_Vspaces();
        // Vb_.rec.print_subspace_details("Vb_.rec");

      const int kSC = red_.init_nearest_neighbors(msr_,Vb_,nset,verbose_);

      red_.reduce_kernals<EVL,BSE>(kSC,Vb_,Sb_,evl_,bases_);

      return kSC;
    }

  protected:

    inline void set_pts_sols(int nset_,int nsol_per_set_[],ode_solution ** sols_)
    {
      for (size_t ic = 0, is = 0; ic < nset_; ic++)
      {
        nsol_per_set[ic] = nsol_per_set_[ic];
        for (size_t ip = 0; ip < nsol_per_set[ic]; ip++, is++) pts_full[is] = (sols_full[is]=sols_[is])->pts;
      }
    }
};

#endif
