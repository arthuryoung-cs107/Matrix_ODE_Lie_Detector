#include "LD_ode_system.hh"
#include <cstring>

known_ode::known_ode(int eor_, int ndep_, const char name_[]): ode_system(eor_,ndep_),
name(new char[strlen(name_)+1])
{strcpy(name,name_);}
known_ode::~known_ode() {delete [] name;}

Duffing_ode::Duffing_ode(double alpha_,double beta_,double gamma_,double delta_,double omega_): known_ode(2,1,"Duffing"),
alpha(alpha_), beta(beta_), gamma(gamma_), delta(delta_), omega(omega_) {}
void Duffing_ode::dudx_eval(double x_, double *u_, double *dudx_)
{
  double  u_in = u_[0],
          dudx_in = u_[1];
  dudx_[0] = dudx_in;
  dudx_[1] = (gamma*cos(omega*x_))-((delta*dudx_in) + u_in*(alpha + (beta*(u_in*u_in))));
}
void Duffing_ode::JacF_eval(double x_, double *u_, double **Jac_out_)
{
  Jac_out_[0][0] = gamma*omega*sin(omega*x_);
  Jac_out_[0][1] = alpha + (3.0*beta*u_[0]*u_[0]);
  Jac_out_[0][2] = delta;
  Jac_out_[0][3] = 1.0;
}
