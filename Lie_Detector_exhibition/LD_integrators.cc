#include "LD_integrators.hh"
#include "LD_util.hh"
#include <cstdio>
#include <cstring>
#include <cmath>

runge_kutta_integrator::runge_kutta_integrator(ode_system &ode_, const char name_[], int nk_):
// ode_exponentiator(ode_),
ode_integrator(ode_),
u_state(new double[ndof_ODE]), k_wkspc(Tmatrix<double>(nk_,ndof_ODE))
{strcpy(name,name_);}
runge_kutta_integrator::~runge_kutta_integrator()
{delete [] u_state; free_Tmatrix<double>(k_wkspc);}

rk_fixed::rk_fixed(ode_system &ode_, const char name_[], int nk_): runge_kutta_integrator(ode_,name_,nk_)
{}
rk_fixed::~rk_fixed() {}

rk_adaptive::rk_adaptive(ode_system &ode_, const char name_[], int nk_, int nrc_):
runge_kutta_integrator(ode_,name_,nk_), u_state_new(new double[ndof_ODE]), rc_wkspc(Tmatrix<double>(nrc_,ndof_ODE))
{}
rk_adaptive::~rk_adaptive()
{delete [] u_state_new;free_Tmatrix<double>(rc_wkspc);}

void rk_fixed::set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_)
{
  if (snaps_>1)
  {
    memcpy(*wkspc_,u_state,ndof_ODE*sizeof(double));
    int nsteps = snaps_-1,
        nmin = determine_nmin(nsteps,tend_-tstart_);
    if  (nsteps<nmin)
    {
      bool step_after = (bool)( nmin%nsteps );
      int nstep_int_sub = nmin/(nsteps),
          snap_count = 0;
      double  t_dur = tend_-tstart_,
              del_t_I = t_dur/((double)(nmin)),
              del_t_O = t_dur/((double)(nsteps)),
              del_t_diff = del_t_O-(del_t_I*((double) nstep_int_sub));
      t = tstart_;
      for (size_t o_step = 0; o_step < nsteps; o_step++)
      {
        for (size_t i_step = 0; i_step < nstep_int_sub; i_step++) step(del_t_I);
        if (step_after)
        {
          step(del_t_diff);
          memcpy(wkspc_[++snap_count],u_state,ndof_ODE*sizeof(double));
        }
        else memcpy(wkspc_[++snap_count],u_state,ndof_ODE*sizeof(double));
      }
    }
    else
    {
      int snap_count = 0;
      double del_t_I = (tend_-tstart_)/((double)(nsteps));
      for (size_t o_step = 0; o_step < nsteps; o_step++)
      {
        step(del_t_I);
        memcpy(wkspc_[++snap_count],u_state,ndof_ODE*sizeof(double));
      }
    }
  }
  else if (snaps_ == 1) memcpy(*wkspc_,u_state,ndof_ODE*sizeof(double));
}

void rk_adaptive::set_and_solve_time(double tstart_, double tend_, int snaps_, double **wkspc_)
{
  int exit_code = solve(tstart_,tend_,snaps_,wkspc_);
  if (exit_code==1) printf("(rk_adaptive::set_and_solve_time) Warning - very stiff diffeq. Max stiffness increments reached.\n");
  else if (exit_code==2)
  {
    printf("(rk_adaptive::set_and_solve_time) ERROR - max count of integration steps reached. (t = %.2e of %.2e, iteration %d of %d)\n", t, tend_, nstep, nmax);
    printf("last value\n");
    for (size_t i = 0; i < ndof_ODE; i++) printf("%.2e ", wkspc_[csn][i]);
    printf("\nDerivative value:\n");
  }
  else if (exit_code==3) printf("(rk_adaptive::set_and_solve_time) ERROR - time step has vanished below tolerance.\n");
}

Sun5::Sun5(ode_system &ode_): rk_fixed(ode_,"Sun5",3), r6(sqrt(6.0)),
dq(new double[ndof_ODE]), kb_wkspc(Tmatrix<double>(3,ndof_ODE)) {}
Sun5::~Sun5() {delete [] dq; free_Tmatrix<double>(kb_wkspc);}

DoPri5::DoPri5(ode_system &ode_): rk_adaptive(ode_,"DoPri5",6,5), ysti(new double[ndof_ODE]) {}
DoPri5::~DoPri5() {delete [] ysti;}

DoP853::DoP853(ode_system &ode_): rk_adaptive(ode_,"DoP853",10,8) {}
DoP853::~DoP853() {}

void Sun5::step(double dt)
{
  int iter=0;
  double  delsq,
          del,
          *c;

  // Clear steps
  for(int i=0;i<ndof_ODE;i++) k1[i]=k2[i]=k3[i]=0.0;

  do
  {
      // Check for too many iterations
      if(++iter>1000)
      {
          fputs("Too many iterations in IRK\n",stderr);
          exit(1);
      }
      // Perform update
      for(int i=0;i<ndof_ODE;i++) dq[i]=q[i]+dt*(a11*k1[i]+a12*k2[i]+a13*k3[i]);
      ff(t+(4-r6)/10*dt,dq,k1b);
      for(int i=0;i<ndof_ODE;i++) dq[i]=q[i]+dt*(a21*k1[i]+a22*k2[i]+a23*k3[i]);
      ff(t+(4+r6)/10*dt,dq,k2b);
      for(int i=0;i<ndof_ODE;i++) dq[i]=q[i]+dt*(a31*k1[i]+a32*k2[i]+a33*k3[i]);
      ff(t+dt,dq,k3b);
      nff+=3;
      // Find size of step from previous iteration
      delsq=0;
      for(int i=0;i<ndof_ODE;i++)
      {
          del=k1[i]-k1b[i]; delsq+=del*del;
          del=k2[i]-k2b[i]; delsq+=del*del;
          del=k3[i]-k3b[i]; delsq+=del*del;
      }
      // Switch k1<->k1b and k2<->k2b array pointers. This will make k1 & k2
      // be used as the new values on the next iteration.
      c=k1b;k1b=k1;k1=c;
      c=k2b;k2b=k2;k2=c;
      c=k3b;k3b=k3;k3=c;
  } while(delsq>1e-25);

  // Complete solution
  for(int i=0;i<ndof_ODE;i++) q[i]+=dt*((16-r6)/36*k1[i]+(16+r6)/36*k2[i]+1./9*k3[i]);
  t+=dt;
}

int DoPri5::solve(double tstart_,double tend_,int snaps_, double **wkspc_)
{
  const double  facc1=1.0/fac1,
                facc2=1.0/fac2,
                expo1=0.2-beta*0.75;
                // expo1=1.0/8.0-beta*0.2;
  double  posneg=(tend_>tstart_)?(1.0):(-1.0),
          facold=1e-4,
          hmax=fabs(tend_-tstart_),
          sf = (snaps_>1)?((tend_-tstart_)/(snaps_-1)):(0.0),
          isf = (snaps_>1)?(1.0/sf):(0.0),
          hlamb = 0.0;
  bool  last=false,
        reject=false;
  t=tstart_; csn=0;
  int nsn;

  init_counters();

  if(snaps_>0) memcpy(wkspc_[0],u_state,ndof_ODE*sizeof(double));
  ff(t,u_state,k1); nff++;
  h=hinit(hmax,posneg);

  while (true)
  {
    if (nstep>nmax)
    {
      if (snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
      return 2;
    }
    if (0.1*fabs(h)<=fabs(t)*uround)
    {
      if(snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
      return 3;
    }
    if ((t+1.01*h-tend_)*posneg>0.0)
    {
      h=tend_-t;
      last=true;
    }

    nstep++;

    // the first 6 stages
    for (int i = 0; i < ndof_ODE; i++) yy1[i] = y[i] + h*a21*k1[i];
    ff(t+c2*h, yy1, k2);

    for (int i = 0; i < ndof_ODE; i++) yy1[i] = y[i] + h*(a31*k1[i] + a32*k2[i]);
    ff(t+c3*h, yy1, k3);

    for (int i = 0; i < ndof_ODE; i++) yy1[i] = y[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
    ff(t+c4*h, yy1, k4);

    for (int i = 0; i < ndof_ODE; i++) yy1[i] = y[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
    ff(t+c5*h, yy1, k5);

    for (int i = 0; i < ndof_ODE; i++) ysti[i] = y[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
    ff(t+h, ysti, k6);

    for (int i = 0; i < ndof_ODE; i++) yy1[i] = y[i] + h*(a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
    ff(t+h, yy1, k2);

    // dense output prep
    if (snaps_>1)
      for (int i = 0; i < ndof_ODE; i++) rc5[i] = h*(d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*k2[i]);

    for (int i = 0; i < ndof_ODE; i++) k4[i] = h*(e1*k1[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*k2[i]);
    nff += 6;

    double err = 0.0;
    for (int i = 0; i < ndof_ODE; i++)
    {
      double sqr = k4[i]/(atoli + rtoli*max_d(fabs(y[i]), fabs(yy1[i])));
      err += sqr*sqr;
    }
    err = sqrt(err/((double)ndof_ODE));

    double  fac11 = pow(err,expo1),
            fac = fac11*pow(facold,-beta);
    fac = max_d(facc2,min_d(facc1,fac/safe));
    double hnew = h/fac;

    if (err <= 1.0)
    {
      facold = max_d(err, 1.0e-4);
      naccpt++;

      // stiffness detection
			if (!(naccpt % nstiff) || (iasti > 0))
      {
				double stnum = 0.0, stden = 0.0;
				for (int i = 0; i < ndof_ODE; i++)
        {
					double sqr = k2[i] - k6[i];
					stnum += sqr*sqr;
					sqr = yy1[i] - ysti[i];
					stden += sqr*sqr;
				}
				if (stden > 0.0) hlamb = h*sqrt(stnum/stden);
				if (hlamb > 3.25)
        {
					nonsti = 0;
					iasti++;
					if (iasti == 15)
          {
            if(snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
						return 1;
					}
				}
				else
        {
					nonsti++;
					if (nonsti == 6) iasti = 0;
				}
			}

      if (snaps_>1)
      {
        nsn=int(((t + h)-tstart_)*isf);
        if(nsn>=snaps_-1) nsn=snaps_-2;
        if(nsn>csn)
          for (int i = 0; i < ndof_ODE; i++)
          {
            double  yd0 = y[i],
                    ydiff = yy1[i] - yd0,
                    bspl = h*k1[i] - ydiff;
            rc1[i] = y[i];
            rc2[i] = ydiff;
            rc3[i] = bspl;
            rc4[i] = -h*k2[i] + ydiff - bspl;
          }
      }
      memcpy(k1, k2, ndof_ODE*sizeof(double));
      memcpy(u_state, yy1, ndof_ODE*sizeof(double));
      told = t;
      t += h;

      if (snaps_>1)
      {
        while (csn<nsn)
        {
          csn++;
          double  ti = tstart_ + csn*sf,
                  theta = (ti-told)/h,
                  theta1 = 1.0 - theta;
          for (size_t i = 0; i < ndof_ODE; i++)
            wkspc_[csn][i] = rc1[i] + theta*(rc2[i] + theta1*(rc3[i] + theta*(rc4[i] + theta1*rc5[i])));
        }
        if (last)
        {
          memcpy(wkspc_[++csn],y,ndof_ODE*sizeof(double));
          return 0;
        }
      }

      if (fabs(hnew)>hmax) hnew = posneg*hmax;
      if (reject) hnew = posneg*min_d(fabs(hnew),fabs(h));
      reject = false;
    }
    else
    {
      hnew = h/min_d(facc1,fac11/safe);
      reject = true;
      if (naccpt>=1) nrejct;
      last = false;
    }
    h = hnew;
  }
}
double DoPri5::hinit(double hmax_, double posneg_)
{
  double  dnf = 0.0,
          dny = 0.0;
  for (int i = 0; i < ndof_ODE; i++)
  {
    double  sk = atoli + rtoli * fabs(y[i]),
            sqr = k1[i]/sk;
    dnf += sqr*sqr;
    sqr = y[i]/sk;
    dny += sqr*sqr;
  }

  double h0=min_d(dnf<=1E-10||dny<=1E-10?1.0E-6:sqrt(dny/dnf)*0.01,hmax_)*posneg_;
	for (int i = 0; i < ndof_ODE; i++) k3[i] = y[i] + h0 * k1[i];

  ff(t+h0,k3,k2); nff++;

  double der2 = 0.0;
  for (int i = 0; i < ndof_ODE; i++)
  {
    double  sk = atoli + rtoli * fabs(y[i]),
	          sqr = (k2[i] - k1[i])/sk;
    der2 += sqr*sqr;
  }
  der2 = sqrt(der2)/h0;

  double  der12 = max_d(fabs(der2),sqrt(dnf)),
          h1 = (der12<=1.0E-15)?max_d(1.0E-6,fabs(h0)*1.0E-3):pow(0.01/der12,0.2);
  return min_d(100.0*fabs(h0),min_d(h1,hmax_))*posneg_;
}

int DoP853::solve(double tstart_,double tend_,int snaps_,double **wkspc_)
{
  const double  facc1=1.0/fac1,
                facc2=1.0/fac2,
                expo1=1.0/8-beta*0.2;
  double  posneg=(tend_>tstart_)?(1.0):(-1.0),
          facold=1e-4,
          hmax=fabs(tend_-tstart_),
          hnew,
          err,
          sf = (snaps_>1)?((tend_-tstart_)/(snaps_-1)):(0.0),
          isf = (snaps_>1)?(1.0/sf):(0.0);
  bool  last=false,
        reject=false;
  t=tstart_; csn=0;

  init_counters();

  if(snaps_>0) memcpy(wkspc_[0],u_state,ndof_ODE*sizeof(double));

  ff(t,u_state,k1); nff++;
  h=hinit(hmax,posneg);

  while(true)
  {
    if (nstep>nmax)
    {
      if (snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
      return 2;
    }
    if (0.1*fabs(h)<=fabs(t)*uround)
    {
      if(snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
      return 3;
    }
    if ((t+1.01*h-tend_)*posneg>0.0)
    {
      h=tend_-t;
      last=true;
    }

    // Perform the eighth-order integration step, and estimate the local
    // error
    nstep++;

    step12(k5);
    err=fabs(h)*error_estimation();

    // Estimate new timestep based on the local error
    double  fac11=pow(err,expo1),
            fac=fac11*pow(facold,-beta);
    fac=max_d(facc2,min_d(facc1,fac/safe));
    hnew=h/fac;

    // Check whether the estimated error is within the tolerance
    if(err<=1.0)
    {

      // If the error is within the tolerance, then accept the step
      facold=max_d(err,1.0E-4);
      naccpt++;
      ff(tph,k5,k4);
      nff++;

      // Carry out stiffness detection at periodic intervals, or if
      // stiffness has previously been detected
      if((!(naccpt%nstiff)||(iasti>0))&&detect_stiffness())
      {
        if(snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
        return 1;
      }

      // Check whether snapshots are required over this integration step,
      // and if so, compute the dense output formulae
      int nsn= (int)((tph-tstart_)*isf);
      if(nsn>=snaps_-1) nsn=snaps_-2;
      if(nsn>csn) dense_output();

      // Copy the computed step into the state vector, and the last
      // Runge--Kutta step into the new first Runge--Kutta step. Update
      // the time.
      memcpy(k1,k4,ndof_ODE*sizeof(double));
      memcpy(u_state,k5,ndof_ODE*sizeof(double));
      told=t;
      t=tph;

      // Compute any snapshots using the dense output
      while(csn<nsn) {csn++; dense(wkspc_[csn],tstart_+csn*sf);}

      // If this is the last timestep, then store a snapshot and return
      if (last)
      {
        if(snaps_>1) memcpy(wkspc_[++csn],u_state,ndof_ODE*sizeof(double));
        return 0;
      }

      // Check for special cases for timestep choice
      if (fabs(hnew)>hmax) hnew=posneg*hmax;
      if (reject) hnew=posneg*min_d(fabs(hnew),fabs(h));
      reject=false;
    }
    else
    {
      // If the local error exceeded the tolerance, then compute a new
      // timestep and try again. Note that as in the original DOP853.F,
      // rejected steps at the start of the computation are not counted.
      hnew = h/min_d(facc1,fac11/safe);
      reject=true;
      if(naccpt >=1 ) nrejct++;
      last=false;
    }
    h=hnew;
  }
}
double DoP853::hinit(double hmax_,double posneg_)
{
  double  dnf=0.0,
          dny=0.0;

  // Compute preliminary step size estimate
  for (size_t i = 0; i < ndof_ODE; i++)
  {
      double  sk = atoli + rtoli*fabs(u_state[i]),
              sqr = k1[i]/sk;
      dnf += sqr*sqr;
      sqr = u_state[i]/sk;
      dny += sqr*sqr;
  }
  double h0=min_d(dnf<=1E-10||dny<=1E-10?1.0E-6:sqrt(dny/dnf)*0.01,hmax_)*posneg_;

  // Perform an explicit Euler step
  for (size_t i = 0; i < ndof_ODE; i++) u_state_new[i] = u_state[i] + h0 * k1[i];

  ff(t+h0,u_state_new,k2); nff++;

  double der2 = 0.0;
  // Estimate the second derivative of the solution
  for (size_t i = 0; i < ndof_ODE; i++)
  {
    double sqr = (k2[i]-k1[i])/(atoli + rtoli*fabs(u_state[i]));
    der2 += sqr*sqr;
  }
  der2=sqrt(der2)/h0;

  // The step size is computed such that h**8*max_d(norm(f0),norm(der2))=0.01
  double  der12 = max_d(fabs(der2),sqrt(dnf)),
          h1=(der12<=1.0E-15)?max_d(1.0E-6,fabs(h0)*1.0E-3):pow(0.01/der12,0.125);
  return min_d(100.0*fabs(h0),min_d(h1,hmax_))*posneg_;
}
void DoP853::step12(double *p1_)
{
  const int n = ndof_ODE;
  int i;
  for(i=0;i<n;i++) ww1[i]=w[i]+h*a21*k1[i];
  ff(t+c2*h,ww1,k2);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a31*k1[i]+a32*k2[i]);
  ff(t+c3*h,ww1,k3);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a41*k1[i]+a43*k3[i]);
  ff(t+c4*h,ww1,k4);
  for(i=0;i <n;i++) ww1[i]=w[i]+h*(a51*k1[i]+a53*k3[i]+a54*k4[i]);
  ff(t+c5*h,ww1,k5);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a61*k1[i]+a64*k4[i]+a65*k5[i]);
  ff(t+c6*h,ww1,k6);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a71*k1[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
  ff(t+c7*h,ww1,k7);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a81*k1[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
  ff(t+c8*h,ww1,k8);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);
  ff(t+c9*h,ww1,k9);
  for(i=0;i<n;i++)
    ww1[i]=w[i]+h*(a101*k1[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+a107*k7[i]+a108*k8[i]+a109*k9[i]);
  ff(t+c10*h,ww1,k10);
  for(i=0;i<n;i++)
    ww1[i]=w[i]+h*(a111*k1[i]+a114*k4[i]+a115*k5[i]+a116*k6[i] +a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
  ff(t+c11*h,ww1,k2);
  tph=t+h;
  for(i=0;i<n;i++)
    ww1[i]=w[i]+h*(a121*k1[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k2[i]);
  ff(tph,ww1,k3);
  nff+=11;
  for(i=0;i<n;i++)
  {
    k4[i]=b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+b11*k2[i]+b12*k3[i];
    p1_[i]=w[i]+h*k4[i];
  }
}
double DoP853::error_estimation()
{
  double  err=0.0,
          err2=0.0;
  // Calculate the contribution to the error from each variable
  for(int i=0;i<ndof_ODE;i++)
  {
    double  sk=1.0/(atoli+rtoli*max_d(fabs(w[i]),fabs(k5[i]))),
            sqr=k4[i]-bhh1*k1[i]-bhh2*k9[i]-bhh3*k3[i];
    sqr*=sk;
    err2+=sqr*sqr;
    sqr=er1*k1[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+er10*k10[i]+er11*k2[i]+er12*k3[i];
    sqr*=sk;
    err+=sqr*sqr;
  }

  // Assemble the combination of the third-order and fifth-order estimators
  double deno=err+0.01*err2;
  return err*sqrt(1.0/(deno<=0.0?((double)ndof_ODE):(deno*((double)ndof_ODE))));
}
bool DoP853::detect_stiffness()
{
  double sqr,stnum=0.0,stden=0.0;

  // Calculate the measure of stiffness, using previously computed RK steps
  for (int i=0;i<ndof_ODE;i++)
  {
    sqr=k4[i]-k3[i]; stnum+=sqr*sqr;
    sqr=k5[i]-ww1[i]; stden+=sqr*sqr;
  }

  // If the stiffness criterion is met, and many recent steps have been
  // stiff, then bail out
  if ((stden>0.0)&&(h*h*stnum>37.21*stden))
  {
    nonsti=0;
    iasti++;
    if(iasti==15) return true;
  }

  // If the system is not stiff after several tests, then reset the stiffness
  // counter
  if (++nonsti==6) iasti=0;
  return false;
}
void DoP853::dense_output()
{
  const int n = ndof_ODE;
  int i;

  // Calculate the contributions to the dense output corrections using
  // the previously computed RK steps
  for(i=0;i<n;i++)
  {
    rc1[i]=w[i];
    double  ydiff=rc2[i]=k5[i]-w[i],
            bspl=rc3[i]=h*k1[i]-ydiff;
    rc4[i]=ydiff-h*k4[i]-bspl;
    rc5[i]=d41*k1[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+d49*k9[i]+d410*k10[i]+d411*k2[i]+d412*k3[i];
    rc6[i]=d51*k1[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+d59*k9[i]+d510*k10[i]+d511*k2[i]+d512*k3[i];
    rc7[i]=d61*k1[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+d69*k9[i]+d610*k10[i]+d611*k2[i]+d612*k3[i];
    rc8[i]=d71*k1[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+d79*k9[i]+d710*k10[i]+d711*k2[i]+d712*k3[i];
  }

  // Carry out the next three function evaluations (steps 14 to 16)
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a141*k1[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]
                 +a1410*k10[i]+a1411*k2[i]+a1412*k3[i]+a1413*k4[i]);
  ff(t+c14*h,ww1,k10);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a151*k1[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]
                 +a1511*k2[i]+a1512*k3[i]+a1513*k4[i]+a1514*k10[i]);
  ff(t+c15*h,ww1,k2);
  for(i=0;i<n;i++) ww1[i]=w[i]+h*(a161*k1[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]
                 +a169*k9[i]+a1613*k4[i]+a1614*k10[i]+a1615*k2[i]);
  ff(t+c16*h,ww1,k3);
  nff+=3;

  // Use the newly computed steps to complete the calculation of the dense
  // output corrections
  for(i=0;i<n;i++)
  {
    rc5[i]=h*(rc5[i]+d413*k4[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
    rc6[i]=h*(rc6[i]+d513*k4[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
    rc7[i]=h*(rc7[i]+d613*k4[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
    rc8[i]=h*(rc8[i]+d713*k4[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
  }
}
