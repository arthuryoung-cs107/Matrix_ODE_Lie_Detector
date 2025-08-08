#ifndef ODE_CRV_OBS_HH
#define ODE_CRV_OBS_HH

// #include "LD_parameter_space.hh"

#include "LD_ode.hh"
#include "LD_aux.hh"

/*
  ode_curve_observations is a data structure consisting of observed trajectories,
  interpreted as integral curves of an ordinary differential equation.
  These observations may or may not be noisy, and are viewed both as snapshots of
  a family of curves parameterized by x, as well as a "cloud" of point solutions
  when considering their union. The numerical values are point solutions to an
  order n = eor ordinary differential equation with q = ndep degrees of freedom
  (i.e. dependent variables).

  This data structure essentially comprises the memory footprint of most program.
*/
class ode_curve_observations
{
  inline int comp_ndim(int eor_, int ndep_) {return 1 + ndep_*(eor_+1);}

  public:

    int eor, // equation order (n)
        ndep, // number of dependent variables (q)
        ncrv, // number of curves (m_c)
        nobs, // total number of observed solutions (i.e. capacity of double data as in below)
        * const ode_meta = &eor,
        * const obs_meta = &ncrv;

    int * npts_per_crv = NULL;
    double  * pts_in = NULL,
            * dnp1xu_in = NULL,
            * JFs_in = NULL;

    inline void print_details(const char preamble_[]="", const char postscript_[]="\n")
    {
      printf(
      "%s"
      "(ode_curve_observations::print_details) "
      "n = %d ode system with q = %d degrees of freedom (b = %d)."
      " Observing m_c = %d trivial flows (M = %d total point solutions). \n"
      "   npts_per_crv : %s\n"
      "   pts_in : %s\n"
      "   dnp1xu_in : %s\n"
      "   JFs_in : %s\n"
      "%s",
      preamble_,
        eor, ndep, comp_ndim(eor,ndep),
        ncrv, nobs,
        (npts_per_crv==NULL)?("unallocated"):("allocated"),
        (pts_in==NULL)?("unallocated"):("allocated"),
        (dnp1xu_in==NULL)?("unallocated"):("allocated"),
        (JFs_in==NULL)?("unallocated"):("allocated"),
      postscript_);
    }

    // default constructor
    ode_curve_observations() {}

    /*
      purely memory allocating constructors
    */
    ode_curve_observations(int eor_, int ndep_, int nobs_) :
      eor(eor_), ndep(ndep_), nobs(nobs_),
      pts_in(new double[(1+ndep_*(eor_+1))*nobs_])
      {}
    ode_curve_observations(int eor_, int ndep_, int ncrv_, int nobs_) :
      eor(eor_), ndep(ndep_), ncrv(ncrv_), nobs(nobs_),
      npts_per_crv(new int[ncrv_]), pts_in(new double[(1+ndep_*(eor_+1))*nobs_])
        { for (int i = 0, np = nobs_/ncrv_; i < ncrv_; i++) npts_per_crv[i] = np; }
    ode_curve_observations(int eor_, int ndep_, int ncrv_, int *npts_per_crv_) :
      eor(eor_), ndep(ndep_), ncrv(ncrv_), nobs(LD_linalg::sum_vec<int>(npts_per_crv_,ncrv_)),
      npts_per_crv(new int[ncrv_]), pts_in(new double[(1+ndep_*(eor_+1))*nobs])
        { for (int i = 0; i < ncrv_; i++) npts_per_crv[i] = npts_per_crv_[i]; }
    ode_curve_observations(ode_solspc_meta &meta_, int ncrv_, int nobs_) :
      ode_curve_observations(meta_.eor,meta_.ndep,ncrv_,nobs_) {}
    ode_curve_observations(ode_solspc_meta &meta_, int ncrv_, int *npts_per_crv_) :
      ode_curve_observations(meta_.eor,meta_.ndep,ncrv_,npts_per_crv_) {}

    /*
      memory allocating and setting constructors, via data files.
    */
    ode_curve_observations(const char name_[])
    {
      FILE * file_in = LD_io::fopen_SAFE(name_,"r");
      LD_io::fread_SAFE(ode_meta,sizeof(int),2,file_in);
      LD_io::fread_SAFE(obs_meta,sizeof(int),2,file_in);
      npts_per_crv = new int[ncrv];
      int ndim = 1+(ndep*(eor+1));
      pts_in = new double[ndim*nobs];
      LD_io::fread_SAFE(npts_per_crv,sizeof(int),ncrv,file_in);
      LD_io::fread_SAFE(pts_in,sizeof(double),ndim*nobs,file_in);
      LD_io::fclose_SAFE(file_in);
      printf("(ode_curve_observations::ode_curve_observations) read %s\n",name_);
    }
    ode_curve_observations(const char name_[], const char name1_[]):
      ode_curve_observations(name_) {read_additional_observations(name1_);}
    ode_curve_observations(const char name_[], const char name1_[], const char name2_[]):
      ode_curve_observations(name_,name1_) {read_additional_observations(name2_);}

    // destructor
    ~ode_curve_observations()
    {
      if (npts_per_crv != NULL) delete [] npts_per_crv;
      if (pts_in != NULL) delete [] pts_in;
      if (JFs_in != NULL) delete [] JFs_in;
      if (dnp1xu_in != NULL) delete [] dnp1xu_in;
    }


    inline int nobs_ioffset(int icrv_)
      {int ioffset = 0; for (size_t i = 0; i < icrv_; i++) ioffset += npts_per_crv[i]; return ioffset;}
    inline double *pts_icrv(int icrv_) {return pts_in + (1+ndep*(eor+1))*nobs_ioffset(icrv_);}

    inline bool palloc() {return dnp1xu_in!=NULL;}
    inline bool Jalloc() {return JFs_in!=NULL;}

    void read_basic_observations(const char name_addtl_[], bool force_overwrite_=false);
    void read_additional_observations(const char name_addtl_[]);
    void write_observed_solutions(const char name_[])
    {
      FILE * file_out = LD_io::fopen_SAFE(name_,"wb");
      fwrite(ode_meta,sizeof(int),2,file_out);
      fwrite(obs_meta,sizeof(int),2,file_out);
      fwrite(npts_per_crv,sizeof(int),ncrv,file_out);
      fwrite(pts_in,sizeof(double),(1 + ndep*(eor+1))*nobs,file_out);
      LD_io::fclose_SAFE(file_out);
      printf("(ode_curve_observations::write_observed_solutions) wrote %s\n",name_);
    }
};

#endif
