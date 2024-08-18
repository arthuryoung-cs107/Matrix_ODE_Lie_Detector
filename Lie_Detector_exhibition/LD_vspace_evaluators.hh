#ifndef LD_VEVLS_HH
#define LD_VEVLS_HH

#include "LD_parameter_space.hh"

struct L_vspace_eval: public LD_vspace_evaluator
{
  L_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-10): LD_vspace_evaluator(Sset_,1,nvec_,1,tol_) {}
  L_vspace_eval(L_vspace_eval &evl_,int nsol_): LD_vspace_evaluator(evl_,nsol_) {}
  ~L_vspace_eval() {}

  template <class BND,class BSE> static void evaluate_lambda_signal_strength(BND &bndle_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    L_vspace_eval::evaluate_lambda_signal_strength<BSE>(bndle_.rec,reci_,Sset_,bases_,tol_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_lambda_signal_strength(LD_vspace_record &reco_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    LD_vspace_evaluator::leniently_evaluate_vspace<L_vspace_eval,BSE>(reco_,reci_,
      L_vspace_eval(Sset_,reco_.nvec,tol_),bases_,verbose_);
  }

  // lambda map signal strength
  virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbasis_)
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

struct G_vspace_eval: public LD_vspace_evaluator
{
  G_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-6): LD_vspace_evaluator(Sset_,Sset_.ndep,nvec_,1,tol_), mags_JFs(new double[Sset.ndep]) {}
  G_vspace_eval(G_vspace_eval &evl_,int nsol_):
    LD_vspace_evaluator(evl_,nsol_), mags_JFs(new double[Sset.ndep]) {}
  ~G_vspace_eval() {delete [] mags_JFs;}

  template <class BND,class BSE> static void evaluate_infinitesimal_criterion(BND &bndle_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    G_vspace_eval::evaluate_infinitesimal_criterion<BSE>(bndle_.rec,reci_,Sset_,bases_,tol_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_infinitesimal_criterion(LD_vspace_record &reco_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    LD_vspace_evaluator::strictly_evaluate_vspace<G_vspace_eval,BSE>(reco_,reci_,
      G_vspace_eval(Sset_,reco_.nvec,tol_),bases_,verbose_);
  }

  // infinitesimal criterion
  virtual int nsat_eval_condition(bool sat_flags_[],ode_solution &sol_,function_space_basis &fbse_)
  {
    fbse_.fill_partial_chunk(sol_.pts);
    return count_n_success(sat_flags_,fbse_,sol_.JFs,comp_JFs_mags(sol_.JFs,sol_.ndep,sol_.ndim),nvec_evl);
  }
  inline double * comp_JFs_mags(double **JFs_,int ndep_,int ndim_)
  {
    for (size_t idep = 0; idep < ndep_; idep++)
    {
      double acc = 0.0;
      for (size_t idim = 0; idim < ndim_; idim++) acc += JFs_[idep][idim]*JFs_[idep][idim];
      mags_JFs[idep] = sqrt(acc);
    }
    return mags_JFs;
  }

  inline int count_n_success(bool sat_flags_[],function_space_basis &fbse_,double **JFs_,double mags_[],int nvec_)
  {
    int n_success = 0;
    for (size_t iV = 0; iV < nvec_; iV++)
      if (sat_flags_[iV] = check_inf_criterion(JFs_,fbse_.v_eval(Vmat[iV]),mags_,fbse_.ndep,fbse_.ndim))
        n_success++;
    return n_success;
  }
  inline bool check_inf_criterion(double **JFs_,double viV_[],double mags_[],int ndep_,int ndim_)
  {
    const double mag_v = LD_linalg::norm_l2(viV_,ndim_);
    bool sat_out;
    for (size_t idep = 0; idep < ndep_; idep++)
    {
      double acc = 0.0;
      for (size_t idim = 0; idim < ndim_; idim++) acc += JFs_[idep][idim]*viV_[idim];
      if (!( sat_out = (fabs(acc / (mags_[idep]*mag_v)) < tol) )) break;
    }
    return sat_out;
  }

  protected:

    double * const mags_JFs;

};

struct R_vspace_eval: public LD_vspace_evaluator
{
  R_vspace_eval(ode_solspc_subset &Sset_,int ncon_,int nvec_,int ord_,double atol_=1e-8,double rtol_=0e-0): LD_vspace_evaluator(Sset_,ncon_,nvec_,1,atol_), ord(ord_), rtol(rtol_) {}
  R_vspace_eval(R_vspace_eval &evl_,int nsol_):
    LD_vspace_evaluator(evl_,nsol_), ord(evl_.ord), rtol(evl_.rtol) {}
  ~R_vspace_eval() {}

  const int ord;
  double  &atol = tol,
          rtol;

  inline bool check_ratio_condition(double viV_[],double dkxu_[],int len_,int offs_)
  {
    bool sat_out;
    for (size_t ievl = 0, idim = offs_; ievl < len_; ievl++, idim++)
      if (!( sat_out = R_err_sat(viV_[idim],viV_[0],dkxu_[ievl]) )) break;
    return sat_out;
  }

  protected:

    static double R_err(double vdkm1xu_, double vx_, double dkxu_)
    {
      return dkxu_ - (vdkm1xu_/vx_); // absolute error
    }
    inline double R_err_sat(double vdkm1xu_, double vx_, double dkxu_)
      {return fabs(dkxu_ - (vdkm1xu_/vx_)) < max_d(atol,fabs(rtol*dkxu_));}
    inline double max_d(double a_,double b_)
      {return (a_>b_)?(a_):(b_);}
};

struct Rk_vspace_eval: public R_vspace_eval
{
  Rk_vspace_eval(ode_solspc_subset &Sset_,int nvec_,int kor_,double atol_,double rtol_):
    R_vspace_eval(Sset_,Sset_.ndep,nvec_,kor_,atol_,rtol_) {}
  Rk_vspace_eval(Rk_vspace_eval &evl_,int nsol_):
    R_vspace_eval(evl_,nsol_) {}
  ~Rk_vspace_eval() {}

  template <class BND,class BSE> static void evaluate_kth_ratio_condition(BND &bndle_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int kor_=0,bool verbose_=true)
  {
    Rk_vspace_eval::evaluate_kth_ratio_condition<BSE>(bndle_.rec,reci_,Sset_,bases_,atol_,rtol_,kor_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_kth_ratio_condition(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int kor_=0,bool verbose_=true)
  {
    LD_vspace_evaluator::strictly_evaluate_vspace<Rk_vspace_eval,BSE>(reco_,reci_,
      Rk_vspace_eval(Sset_,reco_.nvec,(kor_)?(kor_):(Sset_.eor),atol_,rtol_),
      bases_,verbose_);
  }

  virtual int nsat_eval_condition(bool sat_flags_[],ode_solution &sol_,function_space_basis &fbse_) // k'th ratio condition
  {
    fbse_.fill_partial_chunk(sol_.pts,korm1);
    return count_n_success(sat_flags_,fbse_,(kor <= sol_.eor)?(sol_.pts + idkxu):(sol_.dnp1xu),nvec_evl);
  }
  inline int count_n_success(bool sat_flags_[],function_space_basis &fbse_,double dkxu_[],int nvec_)
  {
    int n_success = 0;
    for (size_t iV = 0; iV < nvec_; iV++)
      if (sat_flags_[iV] = check_ratio_condition(fbse_.v_eval(Vmat[iV],korm1),dkxu_,fbse_.ndep,idkm1xu))
        n_success++;
    return n_success;
  }

  protected:

    const int kor = ord,
              korm1 = kor-1,
              idkxu = 1 + Sset.ndep*kor,
              idkm1xu = 1 + Sset.ndep*korm1;
};

struct Rn_vspace_eval: public R_vspace_eval
{
  Rn_vspace_eval(ode_solspc_subset &Sset_,int nvec_,int nor_,double atol_,double rtol_):
    R_vspace_eval(Sset_,Sset_.ndep*nor_,nvec_,nor_,atol_,rtol_) {}
  Rn_vspace_eval(Rn_vspace_eval &evl_, int nsol_):
    R_vspace_eval(evl_,nsol_) {}
  ~Rn_vspace_eval() {}

  template <class BND,class BSE> static void evaluate_nth_ratio_condition(BND &bndle_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int nor_=0,bool verbose_=true)
  {
    Rn_vspace_eval::evaluate_nth_ratio_condition<BSE>(bndle_.rec,reci_,Sset_,bases_,atol_,rtol_,nor_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_nth_ratio_condition(LD_vspace_record &reco_,
    LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int nor_=0,bool verbose_=true)
  {
    LD_vspace_evaluator::strictly_evaluate_vspace<Rn_vspace_eval,BSE>(reco_,reci_,
      Rn_vspace_eval(Sset_,reco_.nvec,(nor_)?(nor_):(Sset_.eor),atol_,rtol_),
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
        if (!( sat_flags_[iV] = R_err_sat(vu[idxu],vx,dxu[idxu]) )) break;

      if (sat_flags_[iV])
      {
        if (nor==(eor+1))
        {
          double * const dnp1xu = sol_.dnp1xu;
          for (size_t idep = 0; idep < ndep; idep++)
            if (!( sat_flags_[iV] = R_err_sat(vu[idep+ndxu],vx,dnp1xu[idep]) )) break;
          if (sat_flags_[iV]) n_success++;
        }
        else n_success++;
      }
    }
    return n_success;
  }

  protected:

    const int nor = ord,
              norm1 = nor-1;
};

struct OG_vspace_eval: public LD_vspace_evaluator
{
  protected:

    Rk_vspace_eval O_evl;
    G_vspace_eval G_evl;

  public:

    OG_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double atol_,double rtol_,double gtol_):
      LD_vspace_evaluator(Sset_,2*(Sset_.ndep),nvec_,1,atol_),
      O_evl(Rk_vspace_eval(Sset_,nvec_,Sset_.eor,atol_,rtol_)), G_evl(G_vspace_eval(Sset_,nvec_,gtol_)) {}
    OG_vspace_eval(OG_vspace_eval &evl_,int nsol_): // sharing solution space and settings
      LD_vspace_evaluator(evl_,nsol_),
      O_evl(Rk_vspace_eval(evl_.O_evl,1)), G_evl(G_vspace_eval(evl_.G_evl,1)) {}
    ~OG_vspace_eval() {}

    double  &gtol = G_evl.tol,
            &atol = O_evl.atol,
            &rtol = O_evl.rtol;

    template <class BSE> void evaluate_vspace_record(LD_vspace_record &reco_,ode_solspc_subset &Sset_,BSE **bases_,
      bool verbose_=true)
    {
      LD_vspace_record reci_(reco_);
      OG_vspace_eval::evaluate_nthcond_infcrit<BSE>(reco_,reci_,Sset_,bases_,atol,rtol,gtol);
    }
    template <class BND,class BSE> static void evaluate_nthcond_infcrit(BND &bndle_,
      LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,double gtol_,bool verbose_=true)
    {
      OG_vspace_eval::evaluate_nthcond_infcrit<BSE>(bndle_.rec,reci_,Sset_,bases_,atol_,rtol_,gtol_,verbose_);
      bndle_.set_Vspaces();
    }
    template <class BSE> static void evaluate_nthcond_infcrit(LD_vspace_record &reco_,
      LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,double gtol_,bool verbose_=true)
    {
      LD_vspace_evaluator::strictly_evaluate_vspace<OG_vspace_eval,BSE>(reco_,reci_,
        OG_vspace_eval(Sset_,reco_.nvec,atol_,rtol_,gtol_),bases_,verbose_);
    }
    template <class BSE> static void evaluate_nthcond_infcrit(LD_vspace_record &recO_,LD_vspace_record &recG_,
      LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,double gtol_,bool verbose_=true)
    {
      LD_vspace_evaluator::strictly_evaluate_vspace<Rk_vspace_eval,BSE>(recO_,reci_,
        Rk_vspace_eval(Sset_,recO_.nvec,Sset_.eor,atol_,rtol_),bases_,verbose_);
      LD_vspace_evaluator::strictly_evaluate_vspace<G_vspace_eval,BSE>(recG_,reci_,
        G_vspace_eval(Sset_,recG_.nvec,gtol_),bases_,verbose_);
    }

    virtual int nsat_eval_condition(bool *sat_flags_,ode_solution &sol_,function_space_basis &fbse_)
    {
      fbse_.fill_partial_chunk(sol_.pts);
      double * const JFs_mags = G_evl.comp_JFs_mags(sol_.JFs,sol_.ndep,sol_.ndim);
      const int idkm1xu = 1 + sol_.ndep*(sol_.eor-1);
      int n_success = 0;
      for (size_t iV = 0; iV < nvec_evl; iV++)
      {
        double * const viV = fbse_.v_eval(Vmat[iV]);
        if (sat_flags[iV] = O_evl.check_ratio_condition(viV,sol_.dnxu,sol_.ndep,idkm1xu)) // ratio cond successful
          if (sat_flags[iV] = G_evl.check_inf_criterion(sol_.JFs,viV,JFs_mags,sol_.ndep,sol_.ndim)) // inf crit successful
            n_success++;
      }
      return n_success;
    }

};

#endif
