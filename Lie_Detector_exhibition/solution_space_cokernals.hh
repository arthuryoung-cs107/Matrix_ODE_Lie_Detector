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
      int nmem_clst_vec[kSC],
          imem_clst_vec[nset],
          *imem_clst_mat[kSC];

      kmed.assign_clusters(nmem_clst_vec,imem_clst_vec,imem_clst_mat,kSC);

      double dnn_clst[kSC];
      for (size_t k = 0; k < kSC; k++)
        // kmed.get_nearest_neighbors(inn_k(k),imem_clst_mat[k],nmem_clst_vec[k]);
        dnn_clst[k] = kmed.get_nearest_neighbors(inn_k(k),imem_clst_mat[k],nmem_clst_vec[k]);

      printf("clusters\n");
      for (size_t k = 0; k < kSC; k++)
      {
        printf("(%d, %d): ", k, kmed.i_meds[k]);
        for (size_t i = 0; i < nmem_clst_vec[k]; i++) printf("%d ", imem_clst_mat[k][i]);
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
