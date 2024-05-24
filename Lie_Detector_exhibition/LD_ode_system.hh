#ifndef LD_ODE_SYS_HH
#define LD_ODE_SYS_HH

#include "LD_ode.hh"
#include "LD_aux.hh"
#include <cstring>

struct known_ode: public ode_system
{
  known_ode(int eor_, int ndep_, const char name_[]): ode_system(eor_,ndep_),
  name(new char[strlen(name_)+1])
  {strcpy(name,name_);}
  ~known_ode() {delete [] name;}

  char * const name;
};

class Duffing_ode: public known_ode
{
  public:
    // default parameters inducing chaotic trajectories
    Duffing_ode(double alpha_=-1.0, double beta_=1.0, double gamma_=0.5, double delta_=0.3, double omega_=1.2):
      known_ode(2,1,"Duffing"),
      alpha(alpha_), beta(beta_), gamma(gamma_), delta(delta_), omega(omega_) {}
    ~Duffing_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      double  u_in = u_[0], dudx_in = u_[1];
      dudx_[0] = dudx_in;
      dudx_[1] = (gamma*cos(omega*x_))-((delta*dudx_in) + u_in*(alpha + (beta*(u_in*u_in))));
    }
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      Jac_out_[0][0] = gamma*omega*sin(omega*x_);
      Jac_out_[0][1] = alpha + (3.0*beta*u_[0]*u_[0]);
      Jac_out_[0][2] = delta;
      Jac_out_[0][3] = 1.0;
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
      {dnp1xu_[0] = -1.0*(delta*u_[2] + u_[1]*(alpha + 3.0*beta*u_[0]*u_[0]) + gamma*omega*sin(omega*x_));}


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
