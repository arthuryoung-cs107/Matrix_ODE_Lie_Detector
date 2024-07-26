#ifndef LD_VEVLS_HH
#define LD_VEVLS_HH

#include "LD_parameter_space.hh"

struct L_vspace_eval: public vspace_evaluation_package
{
  L_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-10): vspace_evaluation_package(Sset_,1,nvec_,1,tol_) {}
  L_vspace_eval(L_vspace_eval &evl_,int nsol_): vspace_evaluation_package(evl_,nsol_) {}
  ~L_vspace_eval() {}

  template <class BND,class BSE> static void evaluate_lambda_signal_strength(BND &bndle_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    L_vspace_eval::evaluate_lambda_signal_strength<BSE>(bndle_.rec,reci_,Sset_,bases_,tol_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_lambda_signal_strength(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    vspace_evaluation_package::leniently_evaluate_vspace<L_vspace_eval,BSE>(reco_,reci_,
      L_vspace_eval(Sset_,reco_.nvec,tol_),
      bases_,verbose_);
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

struct G_vspace_eval: public vspace_evaluation_package
{
  G_vspace_eval(ode_solspc_subset &Sset_,int nvec_,double tol_=1e-6): vspace_evaluation_package(Sset_,Sset_.ndep,nvec_,1,tol_), mags_JFs(new double[Sset.ndep]) {}
  G_vspace_eval(G_vspace_eval &evl_,int nsol_):
    vspace_evaluation_package(evl_,nsol_), mags_JFs(new double[Sset.ndep]) {}
  ~G_vspace_eval() {delete [] mags_JFs;}

  template <class BND,class BSE> static void evaluate_infinitesimal_criterion(BND &bndle_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    G_vspace_eval::evaluate_infinitesimal_criterion<BSE>(bndle_.rec,reci_,Sset_,bases_,tol_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_infinitesimal_criterion(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double tol_,bool verbose_=true)
  {
    vspace_evaluation_package::strictly_evaluate_vspace<G_vspace_eval,BSE>(reco_,reci_,
      G_vspace_eval(Sset_,reco_.nvec,tol_),
      bases_,verbose_);
  }

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

struct R_vspace_eval: public vspace_evaluation_package
{
  R_vspace_eval(ode_solspc_subset &Sset_,int ncon_,int nvec_,int ord_,double atol_=1e-8,double rtol_=0e-0): vspace_evaluation_package(Sset_,ncon_,nvec_,1,atol_), ord(ord_), rtol(rtol_) {}
  R_vspace_eval(R_vspace_eval &evl_,int nsol_):
    vspace_evaluation_package(evl_,nsol_), ord(evl_.ord), rtol(evl_.rtol) {}
  ~R_vspace_eval() {}

  protected:

    const int ord;
    double  &atol = tol,
            rtol;

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

  template <class BND,class BSE> static void evaluate_kth_ratio_condition(BND &bndle_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int kor_=0,bool verbose_=true)
  {
    Rk_vspace_eval::evaluate_kth_ratio_condition<BSE>(bndle_.rec,reci_,Sset_,bases_,atol_,rtol_,kor_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_kth_ratio_condition(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int kor_=0,bool verbose_=true)
  {
    vspace_evaluation_package::strictly_evaluate_vspace<Rk_vspace_eval,BSE>(reco_,reci_,
      Rk_vspace_eval(Sset_,reco_.nvec,(kor_)?(kor_):(Sset_.eor),atol_,rtol_),
      bases_,verbose_);
  }

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
        if (!( sat_flags_[iV] = R_err_sat(vdkm1xu[idep],vx,dkxu[idep]) )) break;
      if (sat_flags_[iV]) n_success++;
    }
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

  template <class BND,class BSE> static void evaluate_nth_ratio_condition(BND &bndle_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int nor_=0,bool verbose_=true)
  {
    Rn_vspace_eval::evaluate_nth_ratio_condition<BSE>(bndle_.rec,reci_,Sset_,bases_,atol_,rtol_,nor_,verbose_);
    bndle_.set_Vspaces();
  }
  template <class BSE> static void evaluate_nth_ratio_condition(LD_vspace_record &reco_,LD_vspace_record &reci_,ode_solspc_subset &Sset_,BSE **bases_,double atol_,double rtol_,int nor_=0,bool verbose_=true)
  {
    vspace_evaluation_package::strictly_evaluate_vspace<Rn_vspace_eval,BSE>(reco_,reci_,
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

#endif
