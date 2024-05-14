#ifndef LD_INTEGR_HH
#define LD_INTEGR_HH

#include "LD_ode.hh"

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
  // DoPri5_settings(  int nmax_=100000,int nstiff_=1000, // default by Ashby
  DoPri5_settings(  int nmax_=16777216,int nstiff_=4,
                    // double atoli_=1e-7,double rtoli_=1e-7, // default by Ashby
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

struct DoP853_settings: public rk_adaptive_settings
{
  DoP853_settings(  int nmax_=16777216,int nstiff_=4,
                    double atoli_=1e-14,double rtoli_=0e0,
                    double fac1_=1.0/3.0,double fac2_=6.0,
                    double safe_=0.9,double beta_=0.0,double uround_=2.3e-16):
      rk_adaptive_settings(nmax_,nstiff_,atoli_,rtoli_,fac1_,fac2_,safe_,beta_,uround_) {}
  DoP853_settings(rk_adaptive_settings &s_): rk_adaptive_settings(s_) {}
  ~DoP853_settings() {}
};

class DoP853: public rk_adaptive
{
  public:

    DoP853(ode_system &ode_);
    DoP853(ode_system &ode_,rk_adaptive_settings &s_): DoP853(ode_) {init_settings(s_);}
    ~DoP853();

  protected:

    double  * w = u_state,
            * ww1 = u_state_new;

    double  * k1 = k_wkspc[0],
            * k2 = k_wkspc[1],
            * k3 = k_wkspc[2],
            * k4 = k_wkspc[3],
            * k5 = k_wkspc[4],
            * k6 = k_wkspc[5],
            * k7 = k_wkspc[6],
            * k8 = k_wkspc[7],
            * k9 = k_wkspc[8],
            * k10 = k_wkspc[9];

    double  * rc1 = rc_wkspc[0],
            * rc2 = rc_wkspc[1],
            * rc3 = rc_wkspc[2],
            * rc4 = rc_wkspc[3],
            * rc5 = rc_wkspc[4],
            * rc6 = rc_wkspc[5],
            * rc7 = rc_wkspc[6],
            * rc8 = rc_wkspc[7];

  private:

    double tph;

    int solve(double tstart_,double tend_,int snaps_, double **wkspc_);
    double hinit(double hmax_, double posneg_);
    void step12(double * p1_);
    double error_estimation();
    bool detect_stiffness();
    void dense_output();
    inline void dense(double *wkspc_,double ti_)
    {
      double  s=(ti_-told)/h,
              s1=1.0-s;
      for (int i=0;i<ndof_ODE;i++)
        wkspc_[i]=rc1[i]+s*(rc2[i]+s1*(rc3[i]+s*(rc4[i]+s1*(rc5[i]+s*(rc6[i]+s1*(rc7[i]+s*rc8[i]))))));
    }

    // rk coefficients
    const double c2=0.526001519587677318785587544488E-01, c3=0.789002279381515978178381316732E-01,
                 c4=0.118350341907227396726757197510E+00, c5=0.281649658092772603273242802490E+00,
                 c6=0.333333333333333333333333333333E+00, c7=0.25E+00,
                 c8=0.307692307692307692307692307692E+00, c9=0.651282051282051282051282051282E+00,
                 c10=0.6E+00, c11=0.857142857142857142857142857142E+00,

                 b1=5.42937341165687622380535766363E-2, b6=4.45031289275240888144113950566E0,
                 b7=1.89151789931450038304281599044E0, b8=-5.8012039600105847814672114227E0,
                 b9=3.1116436695781989440891606237E-1, b10=-1.52160949662516078556178806805E-1,
                 b11=2.01365400804030348374776537501E-1, b12=4.47106157277725905176885569043E-2,
                 a21=5.26001519587677318785587544488E-2, a31=1.97250569845378994544595329183E-2,
                 a32=5.91751709536136983633785987549E-2, a41=2.95875854768068491816892993775E-2,
                 a43=8.87627564304205475450678981324E-2, a51=2.41365134159266685502369798665E-1,
                 a53=-8.84549479328286085344864962717E-1, a54=9.24834003261792003115737966543E-1,
                 a61=3.7037037037037037037037037037E-2, a64=1.70828608729473871279604482173E-1,
                 a65=1.25467687566822425016691814123E-1, a71=3.7109375E-2,
                 a74=1.70252211019544039314978060272E-1, a75=6.02165389804559606850219397283E-2,
                 a76=-1.7578125E-2,

                 a81=3.70920001185047927108779319836E-2, a84=1.70383925712239993810214054705E-1,
                 a85=1.07262030446373284651809199168E-1, a86=-1.53194377486244017527936158236E-2,
                 a87=8.27378916381402288758473766002E-3, a91=6.24110958716075717114429577812E-1,
                 a94=-3.36089262944694129406857109825E0, a95=-8.68219346841726006818189891453E-1,
                 a96=2.75920996994467083049415600797E1, a97=2.01540675504778934086186788979E1,
                 a98=-4.34898841810699588477366255144E1, a101=4.77662536438264365890433908527E-1,
                 a104=-2.48811461997166764192642586468E0, a105=-5.90290826836842996371446475743E-1,
                 a106=2.12300514481811942347288949897E1, a107=1.52792336328824235832596922938E1,
                 a108=-3.32882109689848629194453265587E1, a109=-2.03312017085086261358222928593E-2,

                 a111=-9.3714243008598732571704021658E-1, a114=5.18637242884406370830023853209E0,
                 a115=1.09143734899672957818500254654E0, a116=-8.14978701074692612513997267357E0,
                 a117=-1.85200656599969598641566180701E1, a118=2.27394870993505042818970056734E1,
                 a119=2.49360555267965238987089396762E0, a1110=-3.0467644718982195003823669022E0,
                 a121=2.27331014751653820792359768449E0, a124=-1.05344954667372501984066689879E1,
                 a125=-2.00087205822486249909675718444E0, a126=-1.79589318631187989172765950534E1,
                 a127=2.79488845294199600508499808837E1, a128=-2.85899827713502369474065508674E0,
                 a129=-8.87285693353062954433549289258E0, a1210=1.23605671757943030647266201528E1,
                 a1211=6.43392746015763530355970484046E-1;

   // error estimation
    const double bhh1=0.244094488188976377952755905512E+00, bhh2=0.733846688281611857341361741547E+00,
                 bhh3=0.220588235294117647058823529412E-01,
                 er1=0.1312004499419488073250102996E-01, er6=-0.1225156446376204440720569753E+01,
                 er7=-0.4957589496572501915214079952E+00, er8=0.1664377182454986536961530415E+01,
                 er9=-0.3503288487499736816886487290E+00, er10=0.3341791187130174790297318841E+00,
                 er11=0.8192320648511571246570742613E-01, er12=-0.2235530786388629525884427845E-01;

    // dense output
    const double c14 = 0.1E+00, c15 = 0.2E+00, c16 = 0.777777777777777777777777777778E+00,

                 a141 =  5.61675022830479523392909219681E-2, a147 =  2.53500210216624811088794765333E-1,
                 a148 = -2.46239037470802489917441475441E-1, a149 = -1.24191423263816360469010140626E-1,
                 a1410 =  1.5329179827876569731206322685E-1, a1411 =  8.20105229563468988491666602057E-3,
                 a1412 =  7.56789766054569976138603589584E-3, a1413 = -8.298E-3,

                 a151 =  3.18346481635021405060768473261E-2, a156 =  2.83009096723667755288322961402E-2,
                 a157 =  5.35419883074385676223797384372E-2, a158 = -5.49237485713909884646569340306E-2,
                 a1511 = -1.08347328697249322858509316994E-4, a1512 =  3.82571090835658412954920192323E-4,
                 a1513 = -3.40465008687404560802977114492E-4, a1514 =  1.41312443674632500278074618366E-1,

                 a161 = -4.28896301583791923408573538692E-1, a166 = -4.69762141536116384314449447206E0,
                 a167 =  7.68342119606259904184240953878E0, a168 =  4.06898981839711007970213554331E0,
                 a169 =  3.56727187455281109270669543021E-1, a1613 = -1.39902416515901462129418009734E-3,
                 a1614 =  2.9475147891527723389556272149E0, a1615 = -9.15095847217987001081870187138E0,

                 d41 =-0.84289382761090128651353491142E+01, d46 = 0.56671495351937776962531783590E+00,
                 d47 =-0.30689499459498916912797304727E+01, d48 = 0.23846676565120698287728149680E+01,
                 d49 = 0.21170345824450282767155149946E+01, d410=-0.87139158377797299206789907490E+00,
                 d411= 0.22404374302607882758541771650E+01, d412= 0.63157877876946881815570249290E+00,
                 d413=-0.88990336451333310820698117400E-01, d414= 0.18148505520854727256656404962E+02,
                 d415=-0.91946323924783554000451984436E+01, d416=-0.44360363875948939664310572000E+01,

                 d51 = 0.10427508642579134603413151009E+02, d56 = 0.24228349177525818288430175319E+03,
                 d57 = 0.16520045171727028198505394887E+03, d58 =-0.37454675472269020279518312152E+03,
                 d59 =-0.22113666853125306036270938578E+02, d510= 0.77334326684722638389603898808E+01,
                 d511=-0.30674084731089398182061213626E+02, d512=-0.93321305264302278729567221706E+01,
                 d513= 0.15697238121770843886131091075E+02, d514=-0.31139403219565177677282850411E+02,
                 d515=-0.93529243588444783865713862664E+01, d516= 0.35816841486394083752465898540E+02,

                 d61= 0.19985053242002433820987653617E+02, d66=-0.38703730874935176555105901742E+03,
                 d67=-0.18917813819516756882830838328E+03, d68= 0.52780815920542364900561016686E+03,
                 d69=-0.11573902539959630126141871134E+02, d610= 0.68812326946963000169666922661E+01,
                 d611=-0.10006050966910838403183860980E+01, d612= 0.77771377980534432092869265740E+00,
                 d613=-0.27782057523535084065932004339E+01, d614=-0.60196695231264120758267380846E+02,
                 d615= 0.84320405506677161018159903784E+02, d616= 0.11992291136182789328035130030E+02,

                 d71 =-0.25693933462703749003312586129E+02, d76 =-0.15418974869023643374053993627E+03,
                 d77 =-0.23152937917604549567536039109E+03, d78 = 0.35763911791061412378285349910E+03,
                 d79 = 0.93405324183624310003907691704E+02, d710=-0.37458323136451633156875139351E+02,
                 d711= 0.10409964950896230045147246184E+03, d712= 0.29840293426660503123344363579E+02,
                 d713=-0.43533456590011143754432175058E+02, d714= 0.96324553959188282948394950600E+02,
                 d715=-0.39177261675615439165231486172E+02, d716=-0.14972683625798562581422125276E+03;
};

#endif
