#ifndef LD_INTEGR_HH
#define LD_INTEGR_HH

#include "LD_ode.hh"

struct ode_exponentiator
{
  ode_exponentiator(ode_system &ode_): ode(ode_) {}
  ~ode_exponentiator() {}

  char name[50];

  ode_system &ode;

  const int ndof_ODE = ode.ndep*ode.eor;

  inline void ff(double tt_, double *in_, double *out_) {ode.dudx_eval(tt_,in_,out_);}

  virtual void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_) = 0;
};

// struct runge_kutta_integrator: public ode_exponentiator
struct runge_kutta_integrator: public ode_integrator
{
  runge_kutta_integrator(ode_system &ode_, const char name_[], int nk_);
  ~runge_kutta_integrator();

  double * get_u_state() {return u_state;};

  protected:

    int nff;

    double  t,
            * const u_state;

    double ** const k_wkspc;
};

struct rk_fixed: public runge_kutta_integrator
{
  rk_fixed(ode_system &ode_, const char name_[], int nk_);
  ~rk_fixed();

  int nstep_int = 1000;
  double del_t_fixed = 0.0;

  void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_);

  protected:

    inline int pos_round_up(double val_)
      {int val_trunc = (int)val_; return (val_>((double) val_trunc))?(val_trunc+1):(val_trunc);}
    inline int determine_nmin(int nstep_, double tdur_)
      {return (del_t_fixed)?(pos_round_up(tdur_/del_t_fixed)):(nstep_int);}

    virtual void step(double dt_) = 0;
};

struct rk_adaptive_settings
{
  rk_adaptive_settings( int nmax_,int nstiff_,
                        double atoli_,double rtoli_,double fac1_,double fac2_,double safe_,double beta_,double uround_):
    nmax(nmax_),nstiff(nstiff_),
    atoli(atoli_), rtoli(rtoli_), fac1(fac1_), fac2(fac2_), safe(safe_), beta(beta_), uround(uround_) {}
  rk_adaptive_settings(rk_adaptive_settings &s_):
    rk_adaptive_settings(s_.nmax,s_.nstiff,s_.atoli,s_.rtoli,s_.fac1,s_.fac2,s_.safe,s_.beta,s_.uround) {}
  ~rk_adaptive_settings() {}

  const size_t  ilen = 2,
                dlen = 7;
  const int nmax,
            nstiff,
            * const ivec = &nmax;
  const double  atoli,
                rtoli,
                fac1,
                fac2,
                safe,
                beta,
                uround,
                * const dvec = &atoli;
};

class rk_adaptive: public runge_kutta_integrator
{
  public:
    rk_adaptive(ode_system &ode_, const char name_[], int nk_, int nrc_);
    ~rk_adaptive();

    int nmax,
        nstiff,
        * const ivec = &nmax;
    double  atoli,
            rtoli,
            fac1,
            fac2,
            safe,
            beta,
            uround,
            * const dvec = &atoli;

    inline void init_settings(rk_adaptive_settings &s_)
    {
      for (size_t i = 0; i < s_.ilen; i++) ivec[i] = s_.ivec[i];
      for (size_t i = 0; i < s_.dlen; i++) dvec[i] = s_.dvec[i];
    }

    void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_);

  protected:

    int nstep,
        naccpt,
        nrejct,
        nff,
        nonsti,
        iasti,
        csn;
    double  h,
            told;

    double * u_state_new,
            ** rc_wkspc;

    virtual int solve(double tstart_, double tend_, int snaps_, double **wkspc_) = 0;

    inline void init_counters() {nstep=naccpt=nrejct=nff=nonsti=iasti=0;}
    inline double max_d(double a_, double b_) {return (a_>b_)?(a_):(b_);}
    inline double min_d(double a_, double b_) {return (a_<b_)?(a_):(b_);}
};

class Sun5: public rk_fixed
{
  public:

    Sun5(ode_system &ode_);
    ~Sun5();

  protected:

    void step(double dt_);

  private:

    const double  r6,
                  a11=(16-r6)/72.,a12=(328-167*r6)/1800,a13=(-2+3*r6)/450,
                  a21=(328+167*r6)/1800,a22=(16+r6)/72,a23=(-2-3*r6)/450,
                  a31=(85-10*r6)/180,a32=(85+10*r6)/180,a33=1/18.;

    double  * const q = u_state,
            * const dq,
            * k1 = k_wkspc[0],
            * k2 = k_wkspc[1],
            * k3 = k_wkspc[2],
            ** const kb_wkspc,
            * k1b = kb_wkspc[0],
            * k2b = kb_wkspc[1],
            * k3b = kb_wkspc[2];
};

struct DoPri5_settings: public rk_adaptive_settings
{
  // DoPri5_settings(  int nmax_=100000,int nstiff_=1000,
  //                   double atoli_=1e-7,double rtoli_=1e-7, // default tolerances by Ashby
  //                   double fac1_=0.2,double fac2_=10.0,
  //                   double safe_=0.9,double beta_=0.04,double uround_=1e-16):
  //     rk_adaptive_settings(nmax_,nstiff_,atoli_,rtoli_,fac1_,fac2_,safe_,beta_,uround_) {}
  DoPri5_settings(  int nmax_=100000,int nstiff_=1000,
                    double atoli_=1e-12,double rtoli_=1e-10,
                    double fac1_=0.2,double fac2_=10.0,
                    double safe_=0.9,double beta_=0.04,double uround_=1e-16):
      rk_adaptive_settings(nmax_,nstiff_,atoli_,rtoli_,fac1_,fac2_,safe_,beta_,uround_) {}
  DoPri5_settings(rk_adaptive_settings &s_): rk_adaptive_settings(s_) {}
  ~DoPri5_settings() {}
};

class DoPri5: public rk_adaptive
{
  public:
    DoPri5(ode_system &ode_);
    DoPri5(ode_system &ode_,rk_adaptive_settings &s_): DoPri5(ode_) {init_settings(s_);}
    ~DoPri5();

  protected:

    double  * k1 = k_wkspc[0],
            * k2 = k_wkspc[1],
            * k3 = k_wkspc[2],
            * k4 = k_wkspc[3],
            * k5 = k_wkspc[4],
            * k6 = k_wkspc[5];
    double  * rc1 = rc_wkspc[0],
            * rc2 = rc_wkspc[1],
            * rc3 = rc_wkspc[2],
            * rc4 = rc_wkspc[3],
            * rc5 = rc_wkspc[4];

    double  * y = u_state,
            * yy1 = u_state_new,
            * ysti;

  private:

    const double  c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0,
			            a21=0.2, a31=3.0/40.0, a32=9.0/40.0,
	                a41=44.0/45.0, a42=-56.0/15.0, a43=32.0/9.0,
	                a51=19372.0/6561.0, a52=-25360.0/2187.0,
			            a53=64448.0/6561.0, a54=-212.0/729.0,
		              a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0,
			            a64=49.0/176.0, a65=-5103.0/18656.0,
			            a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0,
			            a75=-2187.0/6784.0, a76=11.0/84.0,
			            e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0,
			            e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0,
			            d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0,
			            d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0,
			            d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;

    int solve(double tstart_,double tend_,int snaps_, double **wkspc_);
    double hinit(double hmax_, double posneg_);
};

struct dop853_settings
{
  dop853_settings(long nmax_=16777216,
                  long nstiff_=4,
                  double fac1_=1.0/3.0,
                  double fac2_=6.0,
                  double safe_=0.9,
                  double beta_=0.0,
                  double uround_=2.3e-16):
    nmax(nmax_), nstiff(nstiff_), fac1(fac1_), fac2(fac2_), safe(safe_), beta(beta_), uround(uround_) {}
  dop853_settings(dop853_settings &set_):
    nmax(set_.nmax), nstiff(set_.nstiff), fac1(set_.fac1), fac2(set_.fac2), safe(set_.safe), beta(set_.beta), uround(set_.uround) {}
  ~dop853_settings() {}

  const long  nmax, // The maximum number of integration steps.
              nstiff; // The number of timesteps between which to test for stiffness.
  const double  fac1, // A factor used in estimating the length of the next timestep.
                fac2, // A factor used in estimating the length of the next timestep.
                safe, // A safety factor used in estimating the length of the next timestep.
                beta, // An exponent used in estimating the length of the next timestep.
                uround; // An estimate of the smallest number that can be added to 1 to obtain a different number, used to determine when the timestep has become too small.

};

class dop853_integrator: public ode_integrator, public dop853_settings
{
  public:
    dop853_integrator(ode_system &ode_, dop853_settings &set_);
    ~dop853_integrator();

    // const int n = ndof_integration;
    const int n = ndof_ODE;

    /** The absolute local error tolerance. */
    double atoli = 1e-14;
    /** The relative local error tolerance. */
    double rtoli = 0e0;
    /** The total number of integration steps. */
    long nstep;
    /** The number of accepted integration steps. */
    long naccpt;
    /** The number of rejected integration steps. */
    long nrejct;
    /** The number of integration steps of fixed size. */
    long nfixed;
    /** The number of evaluations of the differential equation. */
    long nff;
    /** A counter for the number of non-stiff integration steps. */
    long nonsti;
    /** A counter for the number of stiff integration steps. */
    long iasti;
    /** A counter for the current snapshot number. */
    int csn;
    /** The current time. */
    double t;
    /** The previous time. */
    double told;
    /** The current time plus the integration step size. */
    double tph;
    /** The integration step size. */
    double h;
    /** An array holding the current state vector. */
    double *w;
    /** A temporary array for assembling a new state vector. */
    double *ww1;
    /** Storage for the Runge--Kutta intermediate steps. */
    double *k1,*k2,*k3,*k4,*k5,*k6,*k7,*k8,*k9,*k10;
    /** Storage for the dense output. */
    double *rc1,*rc2,*rc3,*rc4,*rc5,*rc6,*rc7,*rc8;

    int solve(double tstart,double tend,int snaps,double **ws);
    inline int solve(double tstart,double tend) {return solve(tstart,tend,0,NULL);}

    void fixed_step(double h_,bool last);
    bool detect_stiffness();
    void dense_output();
    void dense(double *ws,double ti);
    void init_counters();
    virtual void output() {};

    // double * get_wvec() {return w;}
    double * get_u_state() {return w;}
    void set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_);

  private:
      void step12(double *p1);
      double error_estimation();
      double hinit(double hmax,double posneg);
      /** Returns the minimum of two numbers.
       * \param[in] (a,b) the two numbers.
       * \return The minimum */
      inline double min_d(double a,double b) {return a<b?a:b;}
      /** Returns the maximum of two numbers.
       * \param[in] (a,b) the two numbers.
       * \return The maximum */
      inline double max_d(double a,double b) {return a>b?a:b;}
};


#endif
