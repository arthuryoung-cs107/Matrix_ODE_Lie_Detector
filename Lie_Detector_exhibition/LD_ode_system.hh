#ifndef LD_ODE_SYS_HH
#define LD_ODE_SYS_HH

#include "LD_ode.hh"
#include "LD_aux.hh"

class Duffing_ode: public ode_system
{
  public:
    // default parameters inducing chaotic trajectories
    Duffing_ode(double alpha_=-1.0, double beta_=1.0, double gamma_=0.5, double delta_=0.3, double omega_=1.2);
    ~Duffing_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_);
    void JacF_eval(double x_, double *u_, double **Jac_out_);

    inline const double * get_default_IC_indep_range() {return Duffing_x_range_def;}
    inline const double * get_default_IC_range() {return Duffing_IC_range_def[0];}

  private:

    const int dparam_len = 5;
    const double  alpha,
                  beta,
                  gamma,
                  delta,
                  omega,
                  * const dparams = &alpha;

    const double  pi_private = 3.14159265358979323846;
    const double Duffing_x_range_def[2] = {0.0, 5.0*(2.0*pi_private)/omega};
    const double Duffing_IC_range_def[2][2] = { {-1.25, 1.25}, {-1.25, 1.25} };
};

#endif
