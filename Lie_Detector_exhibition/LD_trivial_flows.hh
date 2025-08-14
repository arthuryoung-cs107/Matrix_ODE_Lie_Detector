#ifndef LD_TRV_FLW_HH
#define LD_TRV_FLW_HH

#include "LD_framework.hh"
#include "LD_encodings.hh"

struct ode_trivial_1jet : public ode_jetspc_element
{
  ode_trivial_1jet(ode_jetspc_meta &jmeta_,ode_solution &sol0_,ode_solution &sol1_) :
    ode_jetspc_element(jmeta_), sol0(sol0_), sol1(sol1_),
    hx(0.5*(sol1_.x - sol0_.x)),
    acube_coeffs(Tmatrix<double>(ndep,4)),
    schnk(eor,ndep), solh(meta,schnk)
    {
      const double  h_sqrd = hx*hx,
                    h_cubd = h_sqrd*hx;
      solh.x = 0.5*(sol1_.x + sol0_.x);
      for (int i = 0; i < ndep; i++)
      {
        const double  u00 = sol0_.u[i],
                      u10 = sol0_.dxu[i],
                      u01 = sol1_.u[i],
                      u11 = sol1_.dxu[i];
        acube_coeffs[i][0] = (0.5)*u00 + (0.25*hx)*u10 + (0.5)*u01 + (-0.25*hx)*u11;
        acube_coeffs[i][1] = (-0.75/hx)*u00 + (-0.25)*u10 + (0.75/hx)*u01 + (-0.25)*u11;
        acube_coeffs[i][2] = (-0.25/hx)*u10 + (0.25/hx)*u11; // + (0.0)*u00 + (0.0)*u01
        acube_coeffs[i][3] = (0.25/h_cubd)*u00 + (0.25/h_sqrd)*u10 + (-0.25/h_cubd)*u01 + (0.25/h_sqrd)*u11;

        solh.u[i] = acube_coeffs[i][0]; solh.dxu[i] = acube_coeffs[i][1];
        for (int k = 2; k <= eor+1; k++)
          solh.u[i+(k*ndep)] = dkxui_hat(k,i);
      }
    }
  ~ode_trivial_1jet()
  {
    free_Tmatrix<double>(acube_coeffs);
  }

  ode_solution  &sol0,
                &sol1;
  ode_solchunk  schnk;
  ode_solution   solh;
  const double hx;
  double ** const acube_coeffs;

  inline double dkxui_hat(int k_, int i_) {return (k_<=jor)?( acube_coeffs[i_][k_]*F_lk(k_,k_) ):(0.0);}
};

struct ode_trivial_sjet : public ode_sjet
{
  ode_trivial_sjet(ode_jetspc_meta &jmeta_,double *avec_,ode_solution &sol0_,ode_solution &sol1_) :
    ode_sjet(jmeta_,avec_,0.5*(sol0_.x + sol1_.x)),
    t1jet(jmeta_,sol0_,sol1_),
    sol0(sol0_), sol1(sol1_)
    {}
  ~ode_trivial_sjet() {}

  ode_trivial_1jet t1jet;

  double &xh = e0;
  ode_solution  &sol0,
                &sol1;

  void print_jet_details()
  {
    printf("(ode_trivial_sjet::print_jet_details) jor = %d, alpha coefficients: \n  ", jor);
    for (int i = 0; i < ndep; i++)
    {
      double *avec_i = a_i(i);
      for (int k = 0; k <= jor; k++) printf("%.3e ", avec_i[k]);
      printf("\n  ");
    }
    printf("s0, sh, s1:\n");
    for (int i = 0; i < ndim; i++) printf("   %.3e, %.3e, %.3e \n", sol0.s_idim(i),s_hat_idim(i),sol1.s_idim(i));
  }
  inline void stage_regularized_trivial_jet(double **bm_,double ***At_,double *sv_,double *lv_)
  {
    const int kor = comp_kor();
    const double hx_local = hx();
    double  ** const Amat0 = At_[0],
             * const Avec0 = Amat0[0];
  }
  inline void set_trivial_Amat(double **Amat_)
    {set_trivial_Amat(Amat_, hx() );}
  inline void set_trivial_Amat(double **Amat_,ode_solution &sol0_,ode_solution &sol1_)
    {set_trivial_Amat(Amat_, 0.5*(sol1_.x-sol0_.x) );}
  inline void set_trivial_Amat(double **Amat_,double hx_)
  {
    const int kor = comp_kor();
    double  ** const Amat0 = Amat_,
             * const Avec0 = Amat0[0],
            ** const Amat1 = Amat0+(kor+1),
             * const Avec1 = Amat1[0];
    Amat0[0][0] = Amat1[0][0] = 1.0;
    for (int k = 1; k <= jor; k++)
    {
      Amat0[0][k] = -hx_*Amat0[0][k-1];
      Amat1[0][k] =  hx_*Amat1[0][k-1];
    }
    for (int l = 1; l <= kor; l++)
    {
      for (int k = 0; k < l; k++)
        Amat0[l][k] = Amat1[l][k] = 0.0;
      for (int k = l, ih = 0; k <= jor; k++, ih++)
      {
        const double F_lk_local = (double)(F_lk(l,k));
        Amat0[l][k] = F_lk_local*Amat0[0][ih];
        Amat1[l][k] = F_lk_local*Amat1[0][ih];
      }
    }
  }
  inline double * set_trivial_bvec(int i_)
    {return set_trivial_bvec(i_,sol0,sol1);}
  inline double * set_trivial_bvec(int i_,ode_solution &sol0_,ode_solution &sol1_)
  {
    const int kor = comp_kor();
    double * const avec_i = a_i(i_);

    /*
      Set 2*(kor+1) entries of constraint vector in each dimension.
      Return address of assigned coefficients for solve in place
    */
    for (int k = 0; k <= kor; k++)
      avec_i[k] = sol0_.dkui(k,i_);
    for (int k = 0, ik = kor+1; k <= kor; k++, ik++)
      avec_i[ik] = sol1_.dkui(k,i_);
    return avec_i;
  }

  inline void set_sol_h(ode_solution &sol_h_)
  {
    // for (int i = 0; i < sol_h_.ndim; i++) sol_h_.pts[i] = s_hat_idim(i);
    sol_h_.x = xh;
    for (int k = 0, iu = 0; k <= sol_h_.eor; k++)
      for (int i = 0; i < ndep; i++, iu++)
        sol_h_.u[iu] = dkxui_hat(k,i);
  }

  inline double ui_hat(int i_) {return a_ij(i_,0);}
  inline double dkxui_hat(int k_, int i_) {return (k_<=jor)?( a_ij(i_,k_)*F_lk(k_,k_) ):(0.0);}
  inline double s_hat_idim(int i_) {return (i_==0)?(xh):( dkxui_hat((i_-1)/ndep,(i_-1)%ndep) );}
  // inline double hx() {return 0.5*(sol1.x-sol0.x);}
  inline double hx() {return t1jet.hx;}

  /*
    Follows from jor = 2*(kor+1)-1, as 2*(kor+1) constraints determine jor order polynomial.
    Remains accurate in the regularized case due to integer division.
  */
  inline int comp_kor() {return (((jor+1)/2)-1);}
};

struct ode_trivial_curvejet : public ode_solspc_element // assumes input points lie on common trivial flow
{
  ode_trivial_curvejet(ode_solcurve &crv_,ode_jetspc_meta &jmeta_trivial_) : ode_solspc_element(crv_.meta),
    sols(crv_.sols),
    jmeta_trivial(jmeta_trivial_),
    nxh(crv_.nobs-1),
    achunk_tjets(new double[nxh*ndep*(jmeta_trivial.jor+1)]),
    tjets(new ode_trivial_sjet*[nxh]),
    meta_1jet(1,ndep),
    fspc_1jet(meta_1jet,3),
      fbse_1jet(fspc_1jet),
    Lsvd_f1jet(crv_.nobs,fspc_1jet.perm_len),
      R1svd_f1jet(ndep*(crv_.nobs),fspc_1jet.ndof_full)
    {
      for (int i = 0, ia = 0, alen = ndep*(jmeta_trivial.jor+1); i < nxh; i++, ia+=alen)
        tjets[i] = new ode_trivial_sjet(jmeta_trivial,
                                        achunk_tjets+ia,
                                        *(crv_.sols[i]),
                                        *(crv_.sols[i+1]));
      LD_L_encoder Lenc(meta_1jet);
      LD_R_encoder R1enc(meta_1jet,1);
      for (int i = 0; i <= nxh; i++)
      {
        Lenc.encode_normalize_rows(Lenc.submat_i(Lsvd_f1jet.Umat,i),fbse_1jet,*(crv_.sols[i]),false);
          R1enc.encode_normalize_rows(R1enc.submat_i(R1svd_f1jet.Umat,i),fbse_1jet,*(crv_.sols[i]),false);
      }
      Lsvd_f1jet.decompose_U();
        R1svd_f1jet.decompose_U();
    }
  ~ode_trivial_curvejet()
  {
    for (size_t i = 0; i < nxh; i++) delete tjets[i];
    delete [] tjets;
    delete [] achunk_tjets;
  }

  // ode_solcurve &scrv;
  ode_solution ** const sols;
  ode_jetspc_meta &jmeta_trivial;

  const int nxh;
  double * const achunk_tjets;

  ode_trivial_sjet ** const tjets;

  ode_solspc_meta meta_1jet;
  orthopolynomial_space fspc_1jet; // unmapped multinomials by default
  orthopolynomial_basis fbse_1jet;
  LD_svd  Lsvd_f1jet,
          R1svd_f1jet;

  inline void print_details()
  {
    Lsvd_f1jet.print_result("Lsvd_f1jet");
    R1svd_f1jet.print_result("R1svd_f1jet");
  }

  inline void stage_regularized_trivial_jet(double **bm_,double ***At_,int j_,double *sv_,double *lv_)
    {tjets[j_]->stage_regularized_trivial_jet(bm_,At_,sv_,lv_);}
  inline void set_trivial_Amat(double **Amat_,int j_) {tjets[j_]->set_trivial_Amat(Amat_);}
  inline double * set_trivial_bvec(int j_,int i_) {return tjets[j_]->set_trivial_bvec(i_);}

  inline void print_jet_j_details(int j_) {tjets[j_]->print_jet_details();}

};

struct LD_trivial_jet_experiment : public LD_experiment
{
  LD_trivial_jet_experiment(LD_observations_set &Sset_) : LD_experiment(Sset_),
    kor(Sset.eor + ((Sset.dnp1xu_tns==NULL)?(0):(1)) ), jor(2*(kor+1) - 1),
    jmeta_trivial(Sset.meta, jor),
    tcrvjets( new ode_trivial_curvejet*[ncrvs_tot] ),
    sigma_chunk( new double[(kor+1)*(Sset_.ndep)] ),
    lambda_chunk( new double[Sset_.ndep] )
    {
      for (int i = 0, is = 0; i < Sset.ndep; i++)
      {
        lambda_chunk[i] = 1.0;
        for (int k = 0; k <= kor; k++, is++) sigma_chunk[is] = 1.0;
      }
      for (int i = 0; i < ncrvs_tot; i++)
        tcrvjets[i] = new ode_trivial_curvejet(*(curves[i]),jmeta_trivial);
    }
  LD_trivial_jet_experiment(LD_observations_set &Sset_,double *svec_,double *lvec_) :
    LD_experiment(Sset_),
    kor(Sset.eor + ((Sset.dnp1xu_tns==NULL)?(0):(1)) ), jor(2*(kor+1) ), // determined jor + 1
    jmeta_trivial(Sset.meta, jor),
    tcrvjets( new ode_trivial_curvejet*[ncrvs_tot] ),
    sigma_chunk( new double[(kor+1)*(Sset_.ndep)] ),
    lambda_chunk( new double[Sset_.ndep] )
    {
      for (int i = 0, is = 0; i < Sset.ndep; i++)
      {
        lambda_chunk[i] = lvec_[i];
        for (int k = 0; k <= kor; k++, is++) sigma_chunk[is] = svec_[is];
      }
      for (int i = 0; i < ncrvs_tot; i++)
      {
        tcrvjets[i] = new ode_trivial_curvejet(*(curves[i]),jmeta_trivial);
      }
    }
  ~LD_trivial_jet_experiment()
  {
    delete [] sigma_chunk;
    delete [] lambda_chunk;
    for (size_t i = 0; i < ncrvs_tot; i++) delete tcrvjets[i];
    delete [] tcrvjets;
  }

  const int kor,
            jor;

  ode_jetspc_meta jmeta_trivial;

  ode_trivial_curvejet ** const tcrvjets;

  // scratch debugging printout
  void print_details()
  {
    printf("(LD_trivial_jet_experiment::print_details) \n");
    int icrv = 0,
        ipts = 0;
    tcrvjets[icrv]->print_jet_j_details(ipts);
  }

  void determine_regularized_trivial_jets()
  {
    int net_clcp = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:net_clcp)
    {
      const int ndep = Sset.ndep,
                jorp1 = jor+1;
      LD_lu lu_t(jorp1,jorp1); // might be better to use cholesky
      double  ** const LUmat_t = lu_t.LUmat,
              ** const bmat_t = Tmatrix<double>(ndep,jorp1),
              *** const Atns_t = T3tensor<double>(ndep,jorp1,jorp1);
      #pragma omp for
      for (int icrv = 0; icrv < ncrvs_tot; icrv++)
      {
        ode_trivial_curvejet &tcj_icrv = *(tcrvjets[icrv]);
        for (int j = 0; j < tcj_icrv.nxh; j++)
        {
          tcj_icrv.stage_regularized_trivial_jet(bmat_t,Atns_t,j,sigma_chunk,lambda_chunk);
          // for (int i = 0; i < ndep; i++)
          //   lu_t.solve_system();
        }
        net_clcp += tcj_icrv.nxh;
      }
      free_Tmatrix<double>(bmat_t);
      free_T3tensor<double>(Atns_t);
    }
    double work_time = LD_threads::toc(t0);
    printf("(LD_trivial_jet_experiment::determine_regularized_trivial_jets) determined %d jets (%d dimensions, %d collocation points, order %d) over %d flows (%.1f points, on avg.) in %.4f seconds (%d threads)\n",
      (Sset.ndep)*net_clcp,
      Sset.ndep,net_clcp,jor,
      ncrvs_tot, ((double)net_clcp)/((double)ncrvs_tot),
      work_time,LD_threads::numthreads());
  }

  void determine_trivial_jets()
  {
    int net_clcp = 0;
    double t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:net_clcp)
    {
      const int ndep = Sset.ndep;
      LD_lu lu_t(jor+1,jor+1);
      double  ** const LUmat_t = lu_t.LUmat;
      #pragma omp for
      for (int icrv = 0; icrv < ncrvs_tot; icrv++)
      {
        ode_trivial_curvejet &tcj_icrv = *(tcrvjets[icrv]);
        for (int j = 0; j < tcj_icrv.nxh; j++)
        {
          tcj_icrv.set_trivial_Amat(LUmat_t,j);
          lu_t.decompose_A();
          for (int i = 0; i < ndep; i++)
            lu_t.solve_system(tcj_icrv.set_trivial_bvec(j,i));
        }
        net_clcp += tcj_icrv.nxh;
      }
    }
    double work_time = LD_threads::toc(t0);
    printf("(LD_trivial_jet_experiment::determine_trivial_jets) determined %d jets (%d dimensions, %d collocation points, order %d) over %d flows (%.1f points, on avg.) in %.4f seconds (%d threads)\n",
      (Sset.ndep)*net_clcp,
      Sset.ndep,net_clcp,jor,
      ncrvs_tot, ((double)net_clcp)/((double)ncrvs_tot),
      work_time,LD_threads::numthreads());
  }

  protected:

    double  * const sigma_chunk,
            * const lambda_chunk;

};

struct LD_global_tjet_experiment : public LD_trivial_jet_experiment
{
  LD_global_tjet_experiment(LD_observations_set &Sset_) :
    LD_trivial_jet_experiment(Sset_),
    meta_global_1jet(1,Sset_.ndep),
    fspc_global_1jet(meta_global_1jet,3),
    Lsvd_global_f1jet(Sset_.nobs,fspc_global_1jet.perm_len),
      R1svd_global_f1jet((Sset_.ndep)*(Sset_.nobs),fspc_global_1jet.ndof_full),
    Lsvd_global_f1jet_h(Sset_.nobs-Sset_.ncrvs_tot,fspc_global_1jet.perm_len),
      R1svd_global_f1jet_h((Sset_.ndep)*(Sset_.nobs-Sset_.ncrvs_tot),fspc_global_1jet.ndof_full)
    {
      // fspc_global_1jet.debugging_description();
      LD_L_encoder Lenc(meta_global_1jet);
      LD_R_encoder R1enc(meta_global_1jet,1);
      #pragma omp parallel for
      for (int icrv = 0; icrv < Sset_.ncrvs_tot; icrv++)
      {
        ode_trivial_curvejet &tcjet_i = *(tcrvjets[icrv]);
        orthopolynomial_basis &fbse_global_1jet_i = tcjet_i.fbse_1jet;
        const int iobs0 = Sset_.get_iobs(icrv,0);
        Lenc.encode_normalize_rows(Lenc.submat_i(Lsvd_global_f1jet.Umat,iobs0),
          fbse_global_1jet_i,*(sols_full[iobs0]),false);
        R1enc.encode_normalize_rows(R1enc.submat_i(R1svd_global_f1jet.Umat,iobs0),
          fbse_global_1jet_i,*(sols_full[iobs0]),false);
        for (int isol=1, iobs=iobs0+1, iobsh=iobs0-icrv; isol <= tcjet_i.nxh; isol++, iobs++, iobsh++)
        {
          Lenc.encode_normalize_rows(Lenc.submat_i(Lsvd_global_f1jet.Umat,iobs),
            fbse_global_1jet_i,*(sols_full[iobs]),false);
          R1enc.encode_normalize_rows(R1enc.submat_i(R1svd_global_f1jet.Umat,iobs),
            fbse_global_1jet_i,*(sols_full[iobs]),false);

          Lenc.encode_normalize_rows(Lenc.submat_i(Lsvd_global_f1jet_h.Umat,iobsh),
            fbse_global_1jet_i,(tcjet_i.tjets[isol-1])->t1jet.solh,false);
          R1enc.encode_normalize_rows(R1enc.submat_i(R1svd_global_f1jet_h.Umat,iobsh),
            fbse_global_1jet_i,(tcjet_i.tjets[isol-1])->t1jet.solh,false);
        }
      }
      Lsvd_global_f1jet.decompose_U();
      R1svd_global_f1jet.decompose_U();
      Lsvd_global_f1jet_h.decompose_U();
      R1svd_global_f1jet_h.decompose_U();

    }
  ~LD_global_tjet_experiment() {}

  ode_solspc_meta meta_global_1jet;
  orthopolynomial_space fspc_global_1jet; // unmapped multinomials by default

  LD_svd  Lsvd_global_f1jet,
          R1svd_global_f1jet,
          Lsvd_global_f1jet_h,
          R1svd_global_f1jet_h;

};

struct LD_svd_vector_field : public ode_system
{
  LD_svd_vector_field(function_space &fspc_,double *sv_,double **VTm_,int ndep_,int eor_) :
    ode_system(eor_,ndep_),
    fspc(fspc_),
    ndof_ODE(eor_*ndep_), rank(0),
    sv(sv_), VTm(VTm_),
    sigma0(sv_[0]), sigmaN(sv_[fspc_.ndof_full-1]),
    ispc(0)
    {}
  ~LD_svd_vector_field()  {}

  function_space &fspc;

  const int ndof_ODE;
  int rank;

  void init_dudx_eval(int ispc_) {ispc = ispc_;}
  void JacF_eval(double x_, double *u_, double **dls_out_) {} // do later
  void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) {} // do later

  inline void set_rank(int iscl_=0)
    {rank = LD_svd::rank(sv,fspc.ndof_full,(iscl_)?(iscl_):(fspc.ndof_full));}
  inline double condition_number() { return sigmaN/sigma0; }

  protected:

    double *  const sv,
           ** const VTm,
           &sigma0,
           &sigmaN;
    int ispc;

};

struct LD_trivial_vector_field : public LD_svd_vector_field
{
  LD_trivial_vector_field(function_space_basis &fbse_,double *sv_,double **VTm_,int ndep_=0) :
    LD_svd_vector_field(fbse_.fspc,sv_,VTm_,(ndep_)?(ndep_):(fbse_.ndep),fbse_.eor),
    fbse(fbse_),
    vxu_wkspc(fspc.nvar,fspc.comp_ord_len()),
    pts_local(new double[fspc.ndim]),
    vn_local(new double[fspc.ndim]),
    lamvec_local(new double[fspc.perm_len]),
    theta_local(new double[fspc.ndof_full])
    {}
  ~LD_trivial_vector_field()
  {
    delete [] pts_local;
    delete [] vn_local;
    delete [] lamvec_local;
    delete [] theta_local;
  }

  function_space_basis &fbse;

  virtual void dudx_eval(double x_,double *u_,double *dudx_)
  {
    // compute local parameter values
    comp_nullspace_theta_local(theta_local,x_,u_);
    eval_theta_image(u_,dudx_);
  }

  inline void comp_nullspace_theta_local(double *theta_,double x_,double *u_)
  {
    if (!rank) set_rank();
    pts_local[0] = x_;
    for (int i = 0; i < fspc.ndep; i++) pts_local[i+1] = u_[i];
    fspc.lamvec_eval(pts_local,lamvec_local,vxu_wkspc);
    const int len_lam = fspc.perm_len,
              len_theta = fspc.ndof_full,
              nvec_use = fspc.ndof_full;
    double Wnet = 0.0;
    for (int i = 0; i < len_theta; i++) theta_[i] = 0.0;
    for (int ith = rank; ith < nvec_use; ith++) // compute local parameter values via specified rank
    {
      double  vx_ith = 0.0;
      for (int i = 0; i < len_lam; i++) vx_ith += VTm[ith][i]*lamvec_local[i];
      for (int i = 0; i < len_theta; i++) theta_[i] += vx_ith*VTm[ith][i];
      Wnet += vx_ith*vx_ith;
    }
    for (int i = 0; i < len_theta; i++) theta_[i] /= Wnet; // normalize local parameters
  }
  inline void eval_pr1_theta_image(double *dudx_) {eval_pr1_theta_image(dudx_,theta_local);}
  inline void eval_pr1_theta_image(double *dudx_, double *theta_)
  {
    for (int i = 0, ideltheta = fspc.perm_len; i < ndep; i++, ideltheta+=fspc.perm_len)
    {
      double * const theta_ui = theta_ + ideltheta;
      dudx_[i] = 0.0;
      for (int j = 0; j < fspc.perm_len; j++)
        dudx_[i] += lamvec_local[j]*theta_ui[j];
    }
  }
  inline void eval_theta_image(double *u_,double *dudx_)
  {
     // load and set first thru n-1 derivative terms
     for (int i = ndep; i < ndof_ODE; i++) pts_local[i+1] = dudx_[i-ndep] = u_[i];

     // pad n'th order terms with zeroes (probably unnecessary with eorcap)
     for (int i = ndof_ODE+1; i < ndim; i++) pts_local[i] = 0.0;

     // evaluate n'th prolongation of vector field parameterized by theta_local
     fbse.v_eval(pts_local,vn_local,theta_local,eor-1);

     // unpack n+1 derivative terms via prolongation
     for (int i = ndof_ODE-ndep; i < ndof_ODE; i++) dudx_[i] = vn_local[i+1];
  }

  protected:

    vxu_workspace vxu_wkspc;
    double * const pts_local,
           * const vn_local,
           * const lamvec_local,
           * const theta_local;

};

struct LD_spectral_tvfield : public LD_trivial_vector_field
{
  LD_spectral_tvfield(function_space_basis &fbse_,double *sv_,double **VTm_,int ndep_=0) :
    LD_trivial_vector_field(fbse_,sv_,VTm_,ndep_)
    {}
  ~LD_spectral_tvfield() {}

  virtual void dudx_eval(double x_,double *u_,double *dudx_)
  {
    // compute local parameter values
    comp_spectral_theta_local(theta_local,x_,u_);
    eval_theta_image(u_,dudx_);
  }
  inline void comp_spectral_theta_local(double x_,double *u_)
    {comp_spectral_theta_local(theta_local,x_,u_);}
  inline void comp_spectral_theta_local(double *theta_,double x_,double *u_)
  {
    pts_local[0] = x_;
    for (int i = 0; i < fspc.ndep; i++) pts_local[i+1] = u_[i];
    fspc.lamvec_eval(pts_local,lamvec_local,vxu_wkspc);
    const int len_lam = fspc.perm_len,
              len_theta = fspc.ndof_full,
              nvec_use = fspc.ndof_full;
    double Wnet = 0.0;
    for (int i = 0; i < len_theta; i++) theta_[i] = 0.0;
    for (int ith = 0; ith < nvec_use; ith++) // compute local parameter values via spectral weighting scheme
    {
      double  vx_ith = 0.0;
      for (int i = 0; i < len_lam; i++) vx_ith += VTm[ith][i]*lamvec_local[i];

      const double  ww_i = sigmaN/sv[ith], // scaled Moore-Penrose pseudoinverse
                    w_i = ww_i*ww_i; // squared to reflect squared scaling of vx_ith
      // if (ww_i>1e-14)
      // {
        for (int i = 0; i < len_theta; i++) theta_[i] += w_i*vx_ith*VTm[ith][i];
        Wnet += w_i*vx_ith*vx_ith;
      // }
    }
    for (int i = 0; i < len_theta; i++) theta_[i] /= Wnet; // normalize local parameters
  }
  inline void comp_rowspace_image(double *theta_,double *vnu_,double *vnu_syn_,double x_,double *u_)
  {
    comp_spectral_theta_local(theta_local,x_,u_);
    for (int i = fspc.nvar; i < fspc.ndim; i++) pts_local[i] = u_[i-1];
    fbse.v_eval(pts_local,vn_local,theta_local,fspc.eor); // yields an estimate for eor+1 derivative
    for (int i = 1; i < fspc.ndim; i++) vnu_[i-1] = vn_local[i];

    // compute remaining n'th prolonged terms
    for (int k = 1, ipts_dkm1xu = fspc.nvar; k <= fspc.eor; k++)
    {
      /*
        Compute k'th prolongation, using k'th order derivatives already computed in k-1'th coordinates.
        Yields estimate for k+1 order derivative.
      */
      for (int i = 0; i < fspc.ndep; i++, ipts_dkm1xu++)
        pts_local[ipts_dkm1xu] = vn_local[ipts_dkm1xu-fspc.ndep];
      fbse.v_eval(pts_local,vn_local,theta_local,k); // yields an estimate for k+1 derivatives
    }
    for (int i = 0; i < fspc.ndof_full; i++) theta_[i] = theta_local[i];
    for (int i = 1; i < fspc.ndim; i++) vnu_syn_[i-1] = vn_local[i];
  }
};

class LD_trivial_generator : public ode_system
{
  LD_spectral_tvfield &vf;

  public:

    LD_trivial_generator(LD_spectral_tvfield &vf_) :
      ode_system(1,vf_.nvar),
      vf(vf_), flow_forward(1)
      {}
    ~LD_trivial_generator() {}

    bool flow_forward;

    inline void dudx_eval(double *pts_)
    {
      vf.comp_spectral_theta_local(pts_[0],pts_+1);
      vf.eval_pr1_theta_image(pts_+vf.nvar);
    }
    void dudx_eval(double eps_, double *s_, double *v_)
    {
      vf.comp_spectral_theta_local(s_[0],s_+1);
      vf.eval_pr1_theta_image(v_+1);
      if (flow_forward) v_[0] = 1.0;
      else // flip orientation of trivial vector field
      {
        v_[0] = -1.0;
        for (int i = 1; i < vf.nvar; i++) v_[i] = -v_[i];
      }
    }

    virtual void JacF_eval(double x_, double *u_, double **dls_out_) {} // do later
    virtual void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) {} // do later
};

struct LD_gtjet_Rmat_experiment : public LD_global_tjet_experiment
{
  LD_gtjet_Rmat_experiment(LD_observations_set &Sset_,orthopolynomial_space &fspc_) :
    LD_global_tjet_experiment(Sset_),
    fspc(fspc_),
    Rnsvd_global( (Sset_.eor*Sset_.ndep)*(Sset_.nobs), fspc.ndof_full ),
      Rnp1svd_global( ((Sset_.eor+1)*Sset_.ndep)*(Sset_.nobs), fspc.ndof_full ),
    Rnsvd_global_h( (Sset_.eor*Sset_.ndep)*(Sset_.nobs-Sset_.ncrvs_tot), fspc.ndof_full ),
      Rnp1svd_global_h( ((Sset_.eor+1)*Sset_.ndep)*(Sset_.nobs-Sset_.ncrvs_tot), fspc.ndof_full ),
    theta_chunk(new double[ (Sset_.nobs)*(fspc.ndof_full) ]),
    vnu_chunk(new double[ (Sset_.nobs)*(Sset_.ndep*(Sset_.eor+1)) ]),
    vnu_syn_chunk(new double[ (Sset_.nobs)*(Sset_.ndep*(Sset_.eor+1)) ])
    {
      const int nrm_flag = false;
      LD_R_encoder  Rnenc(fspc.meta,fspc.eor),
                    Rnp1enc(fspc.meta,fspc.eor+1);
      #pragma omp parallel
      {
        orthopolynomial_basis fbse_t(fspc_);
        #pragma omp for
        for (int icrv = 0; icrv < Sset_.ncrvs_tot; icrv++)
        {
          ode_trivial_curvejet &tcjet_i = *(tcrvjets[icrv]);
          const int iobs0 = Sset_.get_iobs(icrv,0);
          Rnenc.encode_normalize_rows(Rnenc.submat_i(Rnsvd_global.Umat,iobs0),
            fbse_t,*(sols_full[iobs0]),nrm_flag);
          Rnp1enc.encode_normalize_rows(Rnp1enc.submat_i(Rnp1svd_global.Umat,iobs0),
            fbse_t,*(sols_full[iobs0]),nrm_flag);
          for (int isol=1, iobs=iobs0+1, iobsh=iobs0-icrv; isol <= tcjet_i.nxh; isol++, iobs++, iobsh++)
          {
            Rnenc.encode_normalize_rows(Rnenc.submat_i(Rnsvd_global.Umat,iobs),
              fbse_t,*(sols_full[iobs]),nrm_flag);
            Rnp1enc.encode_normalize_rows(Rnp1enc.submat_i(Rnp1svd_global.Umat,iobs),
              fbse_t,*(sols_full[iobs]),nrm_flag);

            Rnenc.encode_normalize_rows(Rnenc.submat_i(Rnsvd_global_h.Umat,iobsh),
              fbse_t,(tcjet_i.tjets[isol-1])->t1jet.solh,nrm_flag);
            Rnp1enc.encode_normalize_rows(Rnp1enc.submat_i(Rnp1svd_global_h.Umat,iobsh),
              fbse_t,(tcjet_i.tjets[isol-1])->t1jet.solh,nrm_flag);
          }
        }
      }

      Rnsvd_global.decompose_U();
        Rnp1svd_global.decompose_U();
      Rnsvd_global_h.decompose_U();
        Rnp1svd_global_h.decompose_U();
    }
  ~LD_gtjet_Rmat_experiment()
  {
    delete [] theta_chunk;
    delete [] vnu_chunk;
    delete [] vnu_syn_chunk;
  }

  orthopolynomial_space &fspc;

  LD_svd  Rnsvd_global,
          Rnp1svd_global,
          Rnsvd_global_h,
          Rnp1svd_global_h;

  void compute_rowspace_image(LD_svd &Rsvd_,LD_observations_set &Sset_,orthopolynomial_basis **fbse_)
  {
    double ** const VTmat = Tmatrix<double>(Rsvd_.Nuse,Rsvd_.Nuse);
    Rsvd_.unpack_VTmat(VTmat);

    ode_solution ** const sols = Sset_.sols;

    double t0 = LD_threads::tic();
    #pragma omp parallel
    {
      // orthopolynomial_basis fbse_t(fspc);
      orthopolynomial_basis &fbse_t = *(fbse_[LD_threads::thread_id()]);
      LD_spectral_tvfield tvf_t(fbse_t,Rsvd_.svec,VTmat,fbse_t.ndep);
      const int vnu_len = (fbse_t.ndep)*(fbse_t.eor+1);
      #pragma omp for
      for (int iobs = 0; iobs < Sset_.nobs; iobs++)
      {
        tvf_t.comp_rowspace_image(theta_chunk + iobs*(fbse_t.ndof_full),
                                  vnu_chunk + iobs*vnu_len,
                                  vnu_syn_chunk + iobs*vnu_len,
                                  sols[iobs]->x, sols[iobs]->u);
      }
    }
    double twork = LD_threads::toc(t0);
    printf("(LD_gtjet_Rmat_experiment::compute_rowspace_image) evaluated %d Rmat row spaces (%d coordinates, N=%d parameter dimensions) in %.4f seconds (%d threads)\n",
    Sset_.nobs, (Sset_.ndep)*(Sset_.eor+1),fspc.ndof_full, twork, LD_threads::numthreads());

    free_Tmatrix<double>(VTmat);
  }

  void write_theta_chunk(const char name_[])
  {
    FILE * file = fopen(name_,"wb");
    if (file != NULL)
    {
      int hlen = 2,
          header[] = {hlen,Sset.nobs,fspc.ndof_full};
      fwrite(header, sizeof(int), hlen+1, file);
      fwrite(theta_chunk, sizeof(double), header[1]*header[2], file);
      fclose(file);
      printf("(LD_gtjet_Rmat_experiment::write_theta_chunk) wrote %s\n",name_);
    }
    else
    {
      printf("(LD_gtjet_Rmat_experiment::write_theta_chunk) fopen error with %s (wb attempted). Exiting\n", name_);
      exit(1);
    }
  }

  void write_vnu_chunk(const char name_[])
  {
    FILE * file = fopen(name_,"wb");
    if (file != NULL)
    {
      int hlen = 2,
          header[] = {hlen,Sset.nobs,fspc.ndep*(fspc.eor+1)};
      fwrite(header, sizeof(int), hlen+1, file);
      fwrite(vnu_chunk, sizeof(double), header[1]*header[2], file);
      fclose(file);
      printf("(LD_gtjet_Rmat_experiment::write_vnu_chunk) wrote %s\n",name_);
    }
    else
    {
      printf("(LD_gtjet_Rmat_experiment::write_vnu_chunk) fopen error with %s (wb attempted). Exiting\n", name_);
      exit(1);
    }
  }
  void write_vnu_syn_chunk(const char name_[])
  {
    FILE * file = fopen(name_,"wb");
    if (file != NULL)
    {
      int hlen = 2,
          header[] = {hlen,Sset.nobs,fspc.ndep*(fspc.eor+1)};
      fwrite(header, sizeof(int), hlen+1, file);
      fwrite(vnu_syn_chunk, sizeof(double), header[1]*header[2], file);
      fclose(file);
      printf("(LD_gtjet_Rmat_experiment::write_vnu_syn_chunk) wrote %s\n",name_);
    }
    else
    {
      printf("(LD_gtjet_Rmat_experiment::write_vnu_syn_chunk) fopen error with %s (wb attempted). Exiting\n", name_);
      exit(1);
    }
  }

  protected:

     double * const theta_chunk,
            * const vnu_chunk,
            * const vnu_syn_chunk;
};

struct LD_gtjet_comparison
{
  LD_gtjet_comparison(LD_observations_set &Sref_,LD_gtjet_Rmat_experiment &tjexp_nse_) :
    Sref(Sref_), Snse(tjexp_nse_.Sset),
    fspc(tjexp_nse_.fspc),
    tjexp_ref(Sref_,tjexp_nse_.fspc), tjexp_nse(tjexp_nse_),
    sols_ref(tjexp_ref.sols_full), sols_nse(tjexp_nse.sols_full)
    {
      tjexp_ref.determine_trivial_jets();
      // tjexp_ref.compute_rowspace_image(tjexp_ref.Rnsvd_global);
    }
  ~LD_gtjet_comparison() {}

  LD_observations_set &Sref,
                      &Snse;
  orthopolynomial_space &fspc;
  LD_gtjet_Rmat_experiment tjexp_ref,
                           &tjexp_nse;
  ode_solution ** const sols_ref,
               ** const sols_nse;
};

#endif
