#ifndef LD_ODE_SYS_HH
#define LD_ODE_SYS_HH

#include "LD_ode.hh"
#include "LD_aux.hh"
#include <cstring>

struct named_ode: public ode_system
{
  named_ode(int eor_, int ndep_, const char name_[]): ode_system(eor_,ndep_),
  name(new char[strlen(name_)+1])
  {strcpy(name,name_);}
  ~named_ode() {delete [] name;}

  char * const name;

  virtual double nse_scl() {return 1.0;}
};

class Riccati_ode: public named_ode
{
  public:
    Riccati_ode(double coeff1_=2.0, double coeff2_=-1.0):
      named_ode(1,1,"Riccati"),
      coeff1(coeff1_), coeff2(coeff2_) {}
    ~Riccati_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
      {dudx_[0] = (coeff1*(u_[0]/x_)) + (coeff2*(x_*x_)*(u_[0]*u_[0]));}
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      double x_sqrd = x_*x_;
      Jac_out_[0][0] = ((coeff1*u_[0])/x_sqrd) - ((2.0*coeff2)*x_*(u_[0]*u_[0]));
      Jac_out_[0][1] = ((-1.0*coeff1)/x_) - ((2.0*coeff2)*(x_sqrd)*u_[0]);
      Jac_out_[0][2] = 1.0;
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
      {dnp1xu_[0] = (coeff1*((u_[1]*x_)-u_[0])/(x_*x_)) + (coeff2*(2.0*x_*x_*u_[0]*u_[0])*((u_[1]/u_[0])+(1.0/x_)));}

    inline const double * get_default_IC_indep_range(int xrange_=0) {return Riccati_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return Riccati_IC_range_def[icrange_];}

    virtual double nse_scl() {return 0.4;}

  private:

    const int dparam_len = 2;
    const double  coeff1,
                  coeff2,
                  * const dparams = &coeff1;

    const double Riccati_x_range_def[1][2] =
    {
      {0.2, 2.0}
    };
    const double Riccati_IC_range_def[1][2] =
    {
      // u
      {0.02, 1.0}
    };
};

class Pendulum_ode: public named_ode
{
  public:
    Pendulum_ode(double r_=2.0, double mass_=1.0, double gamma_=1.0, double Aforce_=0.0, double Fforce_=1.0, double gravity_=-9.81):
      named_ode(2,1,"Pendulum"),
      radius(r_), mass(mass_), gamma(gamma_), Aforce(Aforce_), Fforce(Fforce_), gravity(gravity_) {}
    ~Pendulum_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      const double  u_in = u_[0],
                    dudx_in = u_[1];
      dudx_[0] = dudx_in;
      dudx_[1] = (c1*sin(u_in)) - (c2*dudx_in) + c3*(cos(Fforce*x_));
    }
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      Jac_out_[0][0] = Fforce*c3*sin(Fforce*x_);
      Jac_out_[0][1] = -1.0*c1*cos(u_[0]);
      Jac_out_[0][2] = c2;
      Jac_out_[0][3] = 1.0;
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
      {dnp1xu_[0] = (c1*cos(u_[0]))*u_[1] - (c2*u_[2]) - c3*(sin(Fforce*x_)*(Fforce));}

    inline const double * get_default_IC_indep_range(int xrange_=0) {return Pendulum_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return Pendulum_IC_range_def[icrange_][0];}

  private:

    const int dparam_len = 6;
    const double  radius,
                  mass,
                  gamma,
                  Aforce,
                  Fforce,
                  gravity,
                  * const dparams = &radius;

    const double  pi_private = 3.14159265358979323846,
                  c1 = gravity/radius,
                  c2 = gamma/(mass*radius*radius),
                  c3 = Aforce/(mass*radius*radius);

    const double Pendulum_x_range_def[1][2] =
    {
      {0.0, 4.0*pi_private}
    };
    const double Pendulum_IC_range_def[1][2][2] =
    {
      // u,              dxu
      { {-1.25, 1.25}, {-1.0, 1.0} }
    };
};

class Brusselator_ode: public named_ode
{
  public:
    Brusselator_ode(double c1_a_=1.0, double c1_b_=-4.0, double c1_c_=1.0, double c2_a_=0.0, double c2_b_=3.0, double c2_c_=-1.0):
      named_ode(1,2,"Brusselator"),
      c1_a(c1_a_), c1_b(c1_b_), c1_c(c1_c_),
      c2_a(c2_a_), c2_b(c2_b_), c2_c(c2_c_) {}
    ~Brusselator_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      const double  u1_sqrd_u2 = u_[0]*u_[0]*u_[1];
      dudx_[0] = c1_a + c1_b*u_[0] + c1_c*u1_sqrd_u2;
      dudx_[1] = c2_a + c2_b*u_[0] + c2_c*u1_sqrd_u2;
    }
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      const double  u1_u2 = u_[0]*u_[1],
                    u1_sqrd = u_[0]*u_[0];
      // del u1'
      Jac_out_[0][0] = 0.0; // x
      Jac_out_[0][1] = -1.0*(c1_b + c1_c*(2.0*u1_u2)); // u1
      Jac_out_[0][2] = -1.0*(c1_c*(u1_sqrd)); // u2
      Jac_out_[0][3] = 1.0; // u1'
      Jac_out_[0][4] = 0.0; // u2'

      // del u2'
      Jac_out_[1][0] = 0.0; // x
      Jac_out_[1][1] = -1.0*(c2_b + c2_c*(2.0*u1_u2)); // u1
      Jac_out_[1][2] = -1.0*(c2_c*(u1_sqrd)); // u2
      Jac_out_[1][3] = 0.0; // u1'
      Jac_out_[1][4] = 1.0; // u2'
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
    {
      const double d2xu_fac = 2.0*u_[0]*u_[2]*u_[1] + (u_[0]*u_[0])*u_[3];
      dnp1xu_[0] = c1_b*u_[2] + c1_c*d2xu_fac; // u1'
      dnp1xu_[1] = c2_b*u_[2] + c2_c*d2xu_fac; // u2'
    }

    inline const double * get_default_IC_indep_range(int xrange_=0) {return Brusselator_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return Brusselator_IC_range_def[icrange_][0];}

  private:

    const int dparam_len = 6;
    const double  c1_a,
                  c1_b,
                  c1_c,
                  c2_a,
                  c2_b,
                  c2_c,
                  * const dparams = &c1_a;

    const double  Brusselator_u1ic = 1.5,
                  Brusselator_u2ic = 3.0,
                  Brusselator_u1delic = 0.25*Brusselator_u1ic,
                  Brusselator_u2delic = 0.25*Brusselator_u2ic;

    const double Brusselator_x_range_def[1][2] =
    {
      {0.0, 20.0}
    };
    const double Brusselator_IC_range_def[1][2][2] =
    {
      // u,
      {
        {Brusselator_u1ic-(Brusselator_u1delic), Brusselator_u1ic+(Brusselator_u1delic)}, // u1
        {Brusselator_u2ic-(Brusselator_u2delic), Brusselator_u2ic+(Brusselator_u2delic)} // u2
      }
    };
};

class Bessel_ode: public named_ode
{
  public:
    Bessel_ode(double alpha_=1.0):
      named_ode(2,1,"Bessel"),
      alpha(alpha_) {}
    ~Bessel_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      const double  u_in = u_[0],
                    dudx_in = u_[1];
      dudx_[0] = dudx_in;
      dudx_[1] = -((dudx_in/x_) + (1.0 - (alpha/x_)*(alpha/x_))*u_in);
    }
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      Jac_out_[0][0] = -((u_[1]/x_) + 2.0*(alpha*alpha*u_[0])/(x_*x_))/x_ ;
      Jac_out_[0][1] = 1.0 - (alpha/x_)*(alpha/x_);
      Jac_out_[0][2] = 1.0/x_;
      Jac_out_[0][3] = 1.0;
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
      {dnp1xu_[0] = -((u_[2] - (u_[1]/x_) + 2.0*u_[0]*(alpha/x_)*(alpha/x_))/x_ + (1.0 - (alpha/x_)*(alpha/x_))*u_[1]);}

    inline const double * get_default_IC_indep_range(int xrange_=0) {return Bessel_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return Bessel_IC_range_def[icrange_][0];}

  private:

    const int dparam_len = 1;
    const double  alpha,
                  * const dparams = &alpha;

    const double  j1_root = 3.83170597020751e+000,
                  jp1_1 = -402.759395702553e-003;

    const double Bessel_x_range_def[1][2] =
    {
      {j1_root, 15.0}
    };
    const double Bessel_IC_range_def[1][2][2] =
    {
      // u,              dxu
      { {-1.0,1.0}, {-1.1*jp1_1, 1.1*jp1_1} }
    };
};

class VanDerPol_ode: public named_ode
{
  public:
    VanDerPol_ode(double mu_=1.0):
      named_ode(2,1,"VanDerPol"),
      mu(mu_) {}
    ~VanDerPol_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      const double  u_in = u_[0],
                    dudx_in = u_[1];
      dudx_[0] = dudx_in;
      dudx_[1] = mu*(1.0 - (u_in*u_in))*dudx_in - u_in;
    }
    void JacF_eval(double x_, double *u_, double **Jac_out_)
    {
      Jac_out_[0][0] = 0.0;
      Jac_out_[0][1] = 2.0*mu*u_[0]*u_[1] + 1.0;
      Jac_out_[0][2] = -mu*(1.0-(u_[0]*u_[0]));
      Jac_out_[0][3] = 1.0;
    }
    void dnp1xu_eval(double x_, double *u_, double *dnp1xu_)
      {dnp1xu_[0] = mu*((1.0-(u_[0]*u_[0]))*u_[2] - 2.0*u_[0]*u_[1]*u_[1]) - u_[1];}

    inline const double * get_default_IC_indep_range(int xrange_=0) {return VanDerPol_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return VanDerPol_IC_range_def[icrange_][0];}

  private:

    const int dparam_len = 1;
    const double  mu,
                  * const dparams = &mu;

    const double VanDerPol_x_range_def[1][2] =
    {
      {0.0, 20.0}
    };
    const double VanDerPol_IC_range_def[1][2][2] =
    {
      // u,              dxu
      { {-1.0,1.0}, {-1.0, 1.0} }
    };
};

class Duffing_ode: public named_ode
{
  public:
    // default parameters inducing chaotic trajectories
    Duffing_ode(double alpha_=-1.0, double beta_=1.0, double gamma_=0.5, double delta_=0.3, double omega_=1.2):
      named_ode(2,1,"Duffing"),
      alpha(alpha_), beta(beta_), gamma(gamma_), delta(delta_), omega(omega_) {}
    ~Duffing_ode() {}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      const double  u_in = u_[0],
                    dudx_in = u_[1];
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

    inline const double * get_default_IC_indep_range(int xrange_=0) {return Duffing_x_range_def[xrange_];}
    inline const double * get_default_IC_range(int icrange_=0) {return Duffing_IC_range_def[icrange_][0];}

  private:

    const int dparam_len = 5;
    const double  alpha,
                  beta,
                  gamma,
                  delta,
                  omega,
                  * const dparams = &alpha;

    const double  pi_private = 3.14159265358979323846;
    const double Duffing_x_range_def[5][2] =
    {
      {0.0, 5.0*(2.0*pi_private)/omega},
      {0.0, 4.0*(2.0*pi_private)/omega},
      {0.0, 6.0*(2.0*pi_private)/omega},
      {0.0, 7.0*(2.0*pi_private)/omega},
      {0.0, 7.5*(2.0*pi_private)/omega}
    };
    const double Duffing_IC_range_def[1][2][2] =
    {
      // u,              dxu
      { {-1.25, 1.25}, {-1.25, 1.25} }
    };
};

#endif
