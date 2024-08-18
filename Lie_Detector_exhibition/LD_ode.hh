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
