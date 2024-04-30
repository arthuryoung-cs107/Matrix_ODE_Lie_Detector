#ifndef LD_INTEGR_HH
#define LD_INTEGR_HH

#include "LD_ode.hh"

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

    double * get_wvec() {return w;}
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
