#ifndef LD_ODE_HH
#define LD_ODE_HH
#include <cstddef>

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
          * const dxu = pts + nvar;
  void print_sol();
};

struct ode_solspc_subset: public ode_solspc_element
{
  ode_solspc_subset(ode_solspc_meta &meta_, int nobs_);
  ~ode_solspc_subset();

  const int nobs;
  double  ** const pts_mat;
  ode_solution ** const sols;

  inline void print_jsol(int j_) {sols[j_]->print_sol();}
  inline void print_subset() {for (size_t i = 0; i < nobs; i++) print_jsol(i);}
};

struct ode_solcurve: public ode_solspc_subset
{
  ode_solcurve(ode_solspc_meta &meta_, int nobs_, double *pts_chunk_, int icrv_);
  ~ode_solcurve();

  const int icrv;
  double  * const pts_chunk,
          * const pts0 = pts_chunk;

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

  double del_t;

  inline void ff(double tt_, double *in_, double *out_) {ode.dudx_eval(tt_,in_,out_);}
  inline void init_curve_integration(int crv_count_) {ode.init_dudx_eval(crv_count_);}

  inline void unpack_time_sol(double xstart_, int snaps_, double **wkspc_, double *pts_chunk_)
  {
    const int ndep = ode.ndep,
              nvar = 1 + ndep,
              ndim = ode.ndim;
    for (size_t j_obs = 0, jj = 0; j_obs < snaps_; j_obs++, jj+=ndim)
    {
      double  * const pts_j = pts_chunk_+jj,
              x_j = pts_j[0] = xstart_ + del_t*((double) j_obs);
      for (size_t l = 0; l < ndof_ODE; l++) pts_j[l+1] = wkspc_[j_obs][l];
      ode.dudx_eval(x_j,wkspc_[j_obs],pts_j+nvar);
    }
  }

  virtual double * get_u_state() = 0;
  virtual void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_) = 0;

};

#endif
