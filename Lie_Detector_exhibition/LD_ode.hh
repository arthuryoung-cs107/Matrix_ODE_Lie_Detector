#ifndef LD_ODE_HH
#define LD_ODE_HH
#include "LD_util.hh"

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
struct ode_solchunk
{
  ode_solchunk(int eor_,int ndep_) :
    data_owner(true),
    pts(new double[1+ndep_*(eor_+2)]), JFs(Tmatrix<double>(ndep_,1+ndep_*(eor_+1)))
    {}
  ode_solchunk(double *pts_,double **JFs_) :
    data_owner(false),
    pts(pts_), JFs(JFs_)
    {}
  ~ode_solchunk()
  {
    if (data_owner)
    {
      free_Tmatrix<double>(JFs);
      delete [] pts;
    }
  }
  const bool data_owner;
  double  * const pts,
          ** const JFs;
};

struct ode_solution: public ode_solspc_element
{
  ode_solution(ode_solspc_meta &meta_, double *pts_): ode_solspc_element(meta_), pts(pts_) {}
  ode_solution(ode_solspc_meta &meta_,ode_solchunk &sc_): ode_solspc_element(meta_),
    pts(sc_.pts), dnp1xu(sc_.pts+ndim), JFs(sc_.JFs)
    {}
  ~ode_solution() {}

  double  * const pts,
          &x = pts[0],
          * const u = pts + 1,
          * const dxu = u + ndep,
          * const dnxu = u + ndep*eor,
          * dnp1xu = NULL,
          ** JFs = NULL;

  inline void print_sol(int ndim_print_=0)
  {
    const int ndim_print = (ndim_print_)?(ndim_print_):(ndim);
    for (int idim = 0; idim < ndim_print; idim++)
      printf("%.2e ", pts[idim]);
    printf("\n");
  }
  inline double s_idim(int i_)
    {return (i_>=0)?( (i_<ndim)?(pts[i_]):( dnp1xu[ i_-ndim ] ) ):(0.0);}
  inline double dkui(int k_,int i_)
    {return ( k_<=eor )?( u[i_ + ndep*k_] ):( dnp1xu[i_+ndep*(k_-(eor+1))] );}
  inline void copy_pts(ode_solution &sol_, int len_=0)
  {
    const int len = (len_)?(len_):(ndim);
    x = sol_.x;
    for (int i = 1; i <= len; i++) pts[i] = sol_.pts[i];
  }

};

class tjet_chart
{
  double * const pts_h_alt,
         * const pts_0_alt,
         * const pts_1_alt;

  public:

    ode_solution &solh,
                 &sol0,
                 &sol1,
                 solh_alt,
                 sol0_alt,
                 sol1_alt;

  tjet_chart(ode_solution &solh_,ode_solution &sol0_,ode_solution &sol1_,bool init_=false) :
    pts_h_alt(new double[solh_.ndim]),
    pts_0_alt(new double[solh_.ndim]),
    pts_1_alt(new double[solh_.ndim]),
    solh(solh_), sol0(sol0_), sol1(sol1_),
    solh_alt(solh_.meta,pts_h_alt),
    sol0_alt(solh_.meta,pts_0_alt),
    sol1_alt(solh_.meta,pts_1_alt)
    {
      if (init_)
      {
        solh_alt.copy_pts(solh);
        sol0_alt.copy_pts(sol0);
        sol1_alt.copy_pts(sol1);
      }
    }
  ~tjet_chart()
    { delete [] pts_h_alt;
      delete [] pts_0_alt;
      delete [] pts_1_alt; }

  inline void reload_sol_h() {solh_alt.copy_pts(solh);}

};

struct ode_solspc_subset: public ode_solspc_element // when the data is not necessarily contiguous
{
  // allocates pointers to solutions and points, does not initialize.
  ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,bool palloc_=true,bool Jalloc_=true);
  // inherits solutions and points from elsewhere, does not initialize.
  ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_=NULL,double ***JFs_tns_=NULL);
  // inherits solutions and points from another set, does not initialize.
  ode_solspc_subset(ode_solspc_subset &s_):
    ode_solspc_subset(s_.meta,s_.nobs,s_.pts_mat,s_.sols,s_.dnp1xu_mat,s_.JFs_tns) {}

  // receives points as contiguous memory from elswhere, allocates, initializes.
  ode_solspc_subset(ode_solspc_meta &meta_,int nobs_,double *pts_,double *dnp1xu_=NULL,double **JFs_=NULL) :
    ode_solspc_subset(meta_,nobs_,(dnp1xu_==NULL)?(false):(true),(JFs_==NULL)?(false):(true))
    {
      for (int i = 0; i < nobs; i++) pts_mat[i] = pts_ + (i*ndim);
      if (dnp1xu_!=NULL) for (int i = 0; i < nobs; i++) dnp1xu_mat[i] = dnp1xu_ + (i*ndep);
      if (JFs_!=NULL) for (int i = 0; i < nobs; i++) JFs_tns[i] = JFs_ + (i*ndep);
      initialize_solutions(false);
    }

  ~ode_solspc_subset();

  const bool local_data;
  const int nobs;
  double  ** const pts_mat;
  ode_solution ** const sols;

  double  ** dnp1xu_mat,
          *** JFs_tns;

  inline  void initialize_solutions(bool force_alloc_=false)
  {
    for (int i = 0; i < nobs; i++)
      if (force_alloc_||(sols[i] == NULL)) sols[i] = new ode_solution(meta,pts_mat[i]);
    if (dnp1xu_mat != NULL)
      for (int i = 0; i < nobs; i++) sols[i]->dnp1xu = dnp1xu_mat[i];
    if (JFs_tns != NULL)
      for (int i = 0; i < nobs; i++) sols[i]->JFs = JFs_tns[i];
  }
  inline void print_jsol(int j_) {sols[j_]->print_sol();}
  inline void print_subset() {for (int i = 0; i < nobs; i++) print_jsol(i);}

  virtual int nobs_subset_i(int i_) {return 1;}
  virtual int max_nobs_subset() {return 1;}
  virtual ode_solution ** get_sol_subset_i(int i_) {return sols + i_;}
};

struct solspc_data_chunk: public ode_solspc_subset // when data is vectorized
{
  solspc_data_chunk(ode_solspc_meta &meta_,int nobs_,bool palloc_=true,bool Jalloc_=true);
  solspc_data_chunk(ode_solspc_meta &meta_,int nobs_,double **pts_mat_,ode_solution **sols_,double **dnp1xu_mat_=NULL,double ***JFs_tns_=NULL);
  solspc_data_chunk(ode_solspc_subset &s_) :
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

struct ode_solcurve_chunk : public solspc_data_chunk
{
  ode_solcurve_chunk(ode_solspc_meta &meta_,int ncrv_,int *npts_per_crv_) :
    solspc_data_chunk(meta_,sum_Tvec<int>(npts_per_crv_,ncrv_),true,true),
    ncrv(ncrv_),
    npts_per_crv(Tvec_copy<int>(npts_per_crv_,ncrv_)),
    curves(new ode_solcurve*[ncrv_])
    {
      for (int icrv = 0, ipts=0; icrv < ncrv; ipts+=npts_per_crv[icrv++])
        curves[icrv] = new ode_solcurve(icrv,meta,
                                        npts_per_crv[icrv],
                                        pts_mat+ipts,
                                        sols+ipts,
                                        dnp1xu_mat+ipts,
                                        JFs_tns+ipts
                                        );
    }
  ~ode_solcurve_chunk()
  {
    for (int i = 0; i < ncrv; i++) delete curves[i];
    delete [] curves;
    delete [] npts_per_crv;
  }

  const int ncrv;
  int * const npts_per_crv;
  ode_solcurve ** const curves;
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

  static int comp_kor(int jor_) {return (((jor_+1)/2)-1);}

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
