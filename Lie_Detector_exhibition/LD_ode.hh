#ifndef LD_ODE_HH
#define LD_ODE_HH
#include "LD_util.hh"

#include <cstdio>

struct ode_solspc_meta
{
  ode_solspc_meta(int eor_,int ndep_): eor(eor_), ndep(ndep_) {}
  ~ode_solspc_meta() {}

  const int eor,
            ndep,
            ndim = 1 + ndep*(eor+1),
            nvar = 1 + ndep;
};

struct ode_solspc_element
{
  ode_solspc_element(ode_solspc_meta &meta_): meta(meta_) {}
  ~ode_solspc_element() {}

  ode_solspc_meta &meta;

  const int &eor = meta.eor,
            &ndep = meta.ndep,
            &ndim = meta.ndim,
            &nvar = meta.nvar;
};

struct ode_solution: public ode_solspc_element
{
  ode_solution(ode_solspc_meta &meta_, double *pts_): ode_solspc_element(meta_), pts(pts_) {}
  ~ode_solution() {}

  double  * const pts,
          &x = pts[0],
          * const u = pts + 1,
          * const dxu = u + ndep,
          * const dnxu = u + ndep*eor,
          * dnp1xu = NULL,
          ** JFs = NULL;

  void print_sol();

  inline double s_idim(int i_)
    {return (i_>=0)?( (i_<ndim)?(pts[i_]):( dnp1xu[ i_-ndim ] ) ):(0.0);}
  inline double dkui(int k_,int i_)
    {return ( k_<=eor )?( u[i_ + ndep*k_] ):( dnp1xu[i_+ndep*(k_-(eor+1))] );}
};

struct ode_solspc_subset: public ode_solspc_element // when the data is not necessarily contiguous
{
  ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,bool palloc_=true,bool Jalloc_=true);
  ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_=NULL,double ***JFs_tns_=NULL);
  ode_solspc_subset(ode_solspc_subset &s_):
    ode_solspc_subset(s_.meta,s_.nobs,s_.pts_mat,s_.sols,s_.dnp1xu_mat,s_.JFs_tns) {}
  ~ode_solspc_subset();

  const bool local_data;
  const int nobs;
  double  ** const pts_mat;
  ode_solution ** const sols;

  double  ** dnp1xu_mat,
          *** JFs_tns;

  inline  void initialize_solutions(bool force_alloc_=false)
  {
    for (size_t i = 0; i < nobs; i++)
      if (force_alloc_||(sols[i] == NULL)) sols[i] = new ode_solution(meta,pts_mat[i]);
    if (dnp1xu_mat != NULL)
      for (size_t i = 0; i < nobs; i++) sols[i]->dnp1xu = dnp1xu_mat[i];
    if (JFs_tns != NULL)
      for (size_t i = 0; i < nobs; i++) sols[i]->JFs = JFs_tns[i];
  }
  inline void print_jsol(int j_) {sols[j_]->print_sol();}
  inline void print_subset() {for (size_t i = 0; i < nobs; i++) print_jsol(i);}

  virtual int nobs_subset_i(int i_) {return 1;}
  virtual int max_nobs_subset() {return 1;}
  virtual ode_solution ** get_sol_subset_i(int i_) {return sols + i_;}
};

struct solspc_data_chunk: public ode_solspc_subset // when data is vectorized
{
  solspc_data_chunk(ode_solspc_meta &meta_,int nobs_,bool palloc_=true,bool Jalloc_=true);
  solspc_data_chunk(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_=NULL,double ***JFs_tns_=NULL);
  solspc_data_chunk(ode_solspc_subset &s_):
    solspc_data_chunk(s_.meta,s_.nobs,s_.pts_mat,s_.sols,s_.dnp1xu_mat,s_.JFs_tns) {}
  ~solspc_data_chunk();

  const bool data_owner;
  double  * const pts_chunk = pts_mat[0],
          * dnp1xu_chunk = NULL,
          * JFs_chunk = NULL,
          ** JFs_rows = NULL;

  inline void alloc_dnp1xu_safe(double *dnp1xu_in_=NULL)
  {
    if ((dnp1xu_mat==NULL)&&data_owner&&(!local_data)) dnp1xu_mat=Tmatrix<double>(nobs,ndep);
    initialize_solutions(false);
    initialize_additional_data();
    if (dnp1xu_in_!=NULL) for (size_t i = 0; i < nobs*ndep; i++) dnp1xu_chunk[i] = dnp1xu_in_[i];
  }
  inline void alloc_JFs_safe(double *JFs_in_=NULL)
  {
    if ((JFs_tns==NULL)&&data_owner&&(!local_data)) JFs_tns=T3tensor<double>(nobs,ndep,ndim);
    initialize_solutions(false);
    initialize_additional_data();
    if (JFs_in_!=NULL) for (size_t i = 0; i < nobs*ndep*ndim; i++) JFs_chunk[i] = JFs_in_[i];
  }
  inline void initialize_additional_data()
  {
    if (dnp1xu_mat != NULL) dnp1xu_chunk = dnp1xu_mat[0];
    if (JFs_tns != NULL) {JFs_rows = JFs_tns[0]; JFs_chunk = JFs_rows[0];}
  }
};

struct ode_solcurve: public solspc_data_chunk // data is vectorized, and on a one dimensional flow
{
  ode_solcurve(int icrv_,ode_solspc_meta &meta_,int nobs_,bool palloc_=true,bool Jalloc_=true);
  ode_solcurve(int icrv_,ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_=NULL,double ***JFs_tns_=NULL);
  ode_solcurve(int icrv_,ode_solspc_subset &s_):
    ode_solcurve(icrv_,s_.meta,s_.nobs,s_.pts_mat,s_.sols,s_.dnp1xu_mat,s_.JFs_tns) {}
  ode_solcurve(ode_solcurve &c_):
    ode_solcurve(c_.icrv,c_.meta,c_.nobs,c_.pts_mat,c_.sols,c_.dnp1xu_mat,c_.JFs_tns) {}
  ~ode_solcurve();

  const int icrv;
  double  * const pts0 = pts_chunk,
          * const eps_vec;

  inline ode_solution * sol_0() {return sols[0];}
  inline ode_solution * sol_f() {return sols[nobs-1];}

  inline void print_curve() {print_subset();}
};

struct ode_jetspc_meta : public ode_solspc_element
{
  ode_jetspc_meta(ode_solspc_meta &smeta_,int jor_) : ode_solspc_element(smeta_),
    jor(jor_), Fjor(Tsym<int>(jor+1))
  {
    for (int i = 0; i <= jor; i++) Fjor[0][i] = 1;
    for (int i = 1, jend = jor-1; i <= jor; i++, jend--)
      for (int j = 0; j <= jend; j++)
        Fjor[i][j] = (j+1)*Fjor[i-1][j+1];
  }
  ~ode_jetspc_meta()
    {free_Tmatrix<int>(Fjor);}

  const int jor;

  inline int get_F_lder_kcol(int l_,int k_) {return (l_<=k_)?( Fjor[l_][k_-l_] ):(0);}

  inline void print_details()
  {
    printf("(ode_jetspc_meta::print_details) jor = %d. F coefficients: \n", jor);
    for (int l = 0; l <= jor; l++)
    {
      for (int k = 0; k <= jor; k++) printf("%.1e ", (double) get_F_lder_kcol(l,k) );
      printf("\n");
    }
  }
  inline void print_F()
  {
    printf("(ode_jetspc_meta::print_F) F coefficients : \n");
    for (int l = 0; l <= jor; l++)
    {
      for (int k = 0; k <= jor; k++) printf("%.1e ", (double) get_F_lder_kcol(l,k) );
      printf("\n");
    }
  }

  private:

    int ** const Fjor;

};

struct ode_jetspc_element : public ode_solspc_element
{
  ode_jetspc_element(ode_jetspc_meta &jmeta_) : ode_solspc_element(jmeta_.meta),
    jmeta(jmeta_) {}
  ~ode_jetspc_element() {}

  ode_jetspc_meta &jmeta;
  const int &jor = jmeta.jor;

  inline int F_lk(int l_,int k_) {return jmeta.get_F_lder_kcol(l_,k_);}
};

struct ode_sjet : public ode_jetspc_element
{
  ode_sjet(ode_jetspc_meta &jmeta_, double *avec_, double e0_=0.0) : ode_jetspc_element(jmeta_),
    avec(avec_), e0(e0_) {}
  ~ode_sjet() {}
  double  * const avec,
                  e0;

  inline double * a_i(int i_) {return avec + i_*(jor+1);}
  inline double a_ij(int i_,int j_) {return avec[j_ + i_*(jor+1)];}
};

struct ode_trivial_sjet : public ode_sjet
{
  ode_trivial_sjet(ode_jetspc_meta &jmeta_,double *avec_,ode_solution &sol0_,ode_solution &sol1_) :
    ode_sjet(jmeta_,avec_,0.5*(sol0_.x + sol1_.x)),
    sol0(sol0_), sol1(sol1_) // jor = 2*(kor+1)-1, as 2*(kor+1) constraints determine jor order polynomial
    {}
  ~ode_trivial_sjet() {}

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
  inline void set_trivial_Amat(double **Amat_)
  {
    const int kor = comp_kor();
    const double hx_local = hx();
    double  ** const Amat0 = Amat_,
             * const Avec0 = Amat0[0],
            ** const Amat1 = Amat0+(kor+1),
             * const Avec1 = Amat1[0];

    Amat0[0][0] = Amat1[0][0] = 1.0;
    for (int k = 1; k <= jor; k++)
    {
      Amat0[0][k] = -hx_local*Amat0[0][k-1];
      Amat1[0][k] =  hx_local*Amat1[0][k-1];
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
  {
    const int kor = comp_kor();
    double * const avec_i = a_i(i_);

    /*
      Set 2*(kor+1) entries of constraint vector in each dimension.
      Return address of assigned coefficients for solve in place
    */
    for (int k = 0; k <= kor; k++)
      avec_i[k] = sol0.dkui(k,i_);
    for (int k = 0, ik = kor+1; k <= kor; k++, ik++)
      avec_i[ik] = sol1.dkui(k,i_);
    return avec_i;
  }

  inline double s_hat_idim(int i_) {return (i_>0)?( dkxui_hat((i_-1)/ndep,(i_-1)%ndep) ):(xh);}
  inline double dkxui_hat(int k_, int i_) {return (k_<=jor)?( a_ij(i_,k_)*F_lk(k_,k_) ):(0.0);}

  inline double hx() {return 0.5*(sol1.x-sol0.x);}
  inline int comp_kor() {return (((jor+1)/2)-1);}
};

struct ode_trivial_curvejet : public ode_solspc_element
{
  ode_trivial_curvejet(ode_solcurve &crv_,ode_jetspc_meta &jmeta_trivial_) : ode_solspc_element(crv_.meta),
    jmeta_trivial(jmeta_trivial_),
    nxh(crv_.nobs-1),
    achunk_tjets(new double[nxh*ndep*(jmeta_trivial.jor+1)]),
    tjets(new ode_trivial_sjet*[nxh])
    {
      for (int i = 0, ia = 0, alen = ndep*(jmeta_trivial.jor+1); i < nxh; i++, ia+=alen)
        tjets[i] = new ode_trivial_sjet(jmeta_trivial,achunk_tjets+ia,*(crv_.sols[i]),*(crv_.sols[i+1]));
    }
  ~ode_trivial_curvejet()
  {
    for (size_t i = 0; i < nxh; i++) delete tjets[i];
    delete [] tjets;
    delete [] achunk_tjets;
  }

  ode_jetspc_meta &jmeta_trivial;

  const int nxh;
  double * const achunk_tjets;

  ode_trivial_sjet ** const tjets;

  void print_jet_j_details(int j_)
  {
    tjets[j_]->print_jet_details();
  }

  inline void set_trivial_Amat(double **Amat_,int j_)
    {tjets[j_]->set_trivial_Amat(Amat_);}
  inline double * set_trivial_bvec(int j_,int i_)
    {return tjets[j_]->set_trivial_bvec(i_);}
};

struct ode_system: public ode_solspc_meta
{
  ode_system(int eor_, int ndep_): ode_solspc_meta(eor_,ndep_) {}
  ~ode_system() {}

  virtual void init_dudx_eval(int crv_count_) {}

  virtual void dudx_eval(double x_, double *u_, double *dudx_) = 0;
  virtual void JacF_eval(double x_, double *u_, double **dls_out_) = 0;
  virtual void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) = 0;
};

struct ode_integrator
{
  ode_integrator(ode_system &ode_): ode(ode_) {}
  ~ode_integrator() {}

  ode_system &ode;

  char name[50];

  const int ndof_ODE = ode.ndep*ode.eor;

  double del_x;

  inline void ff(double tt_, double *in_, double *out_) {ode.dudx_eval(tt_,in_,out_);}
  inline void init_curve_integration(double del_x_, int crv_count_)
  {
    del_x = del_x_;
    ode.init_dudx_eval(crv_count_);
  }

  inline void unpack_time_sol(double xstart_, int snaps_, double **wkspc_, double *pts_chunk_)
  {
    const int ndep = ode.ndep,
              nvar = 1 + ndep,
              ndim = ode.ndim;
    for (size_t j_obs = 0, jj = 0; j_obs < snaps_; j_obs++, jj+=ndim)
    {
      double  * const pts_j = pts_chunk_+jj,
              x_j = pts_j[0] = xstart_ + del_x*((double) j_obs);
      for (size_t l = 0; l < ndof_ODE; l++) pts_j[l+1] = wkspc_[j_obs][l];
      ode.dudx_eval(x_j,wkspc_[j_obs],pts_j+nvar);
    }
  }

  virtual double * get_u_state() = 0;
  virtual void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_) = 0;
};

#endif
