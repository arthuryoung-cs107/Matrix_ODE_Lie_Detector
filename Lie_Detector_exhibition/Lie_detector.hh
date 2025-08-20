#ifndef LIE_DETECTOR_HH
#define LIE_DETECTOR_HH

#include "ode_curve_observations.hh"
#include "LD_trivial_flows.hh"

class curve_Lie_detector
{
  bool palloc,
       Jalloc;

  public:

    const int nobs_f,
              nobs_h;

    ode_solcurve &crv, // observed trajectory
                  crv_h; // induced trajectory of colocation points
    ode_trivial_curvejet tcrvjets; // deterministic M'th order Hermite jet approximation

    tjet_chart ** const tjcharts;

    curve_Lie_detector(ode_solcurve &crv_,ode_jetspc_meta &jmeta_trivial_,bool palloc_,bool Jalloc_) :
      palloc(palloc_), Jalloc(Jalloc_),
      nobs_f(crv_.nobs), nobs_h(crv_.nobs-1),
      crv(crv_),
        crv_h(crv_.icrv,crv_.meta,nobs_h,palloc,Jalloc),
        tcrvjets(crv_,jmeta_trivial_),
      tjcharts(new tjet_chart*[nobs_h])
      {
        for (int i = 0; i < nobs_h; i++)
          tjcharts[i] = new tjet_chart(*(crv_h.sols[i]),*(crv.sols[i]),*(crv.sols[i+1]));
      }
    ~curve_Lie_detector()
    {
      for (int i = 0; i < nobs_h; i++) delete tjcharts[i];
      delete [] tjcharts;
    }

    inline double recombine_flows(ode_solcurve &crv_alt_,ode_solcurve &crv_)
    {
      double res_out = 0.0;
      res_out += combine_two_sols( *(crv_alt_.sol_0()), *(crv_.sol_0()), (tjchart_0())->sol0_alt );
      for (int isol = 1; isol < nobs_h; isol++)
        res_out += combine_three_sols( *(crv_alt_.sols[isol]),
                              *(crv_.sols[isol]),
                              tjcharts[isol-1]->sol1_alt,
                              tjcharts[isol]->sol0_alt );
      res_out += combine_two_sols( *(crv_alt_.sol_f()), *(crv_.sol_f()), (tjchart_f())->sol1_alt );
      return res_out;
    }

    inline double combine_three_sols(ode_solution &s0_,ode_solution &s1_,ode_solution &s2_,ode_solution &s3_)
    {
      double res_out = 0.0;
      for (int i = 0; i < s0_.ndim; i++)
      {
        double si_old = s0_.pts[i];

                  // aggressive update. Reinforces dominant components.
        si_old -= (s0_.pts[i] = ( s2_.pts[i] + s3_.pts[i] )/2.0);
                  // lazy update, perturbs old solution.
        // si_old -= (s0_.pts[i] = ( s1_.pts[i] + s2_.pts[i] + s3_.pts[i] )/3.0);

        res_out += si_old*si_old;
      }
      return res_out;
    }
    inline double combine_two_sols(ode_solution &s0_,ode_solution &s1_,ode_solution &s2_)
    {
      double res_out = 0.0;
      for (int i = 0; i < s0_.ndim; i++)
      {
        double si_old = s0_.pts[i];

                  // aggressive update. Reinforces dominant components, still one sided.
        si_old -= (s0_.pts[i] = s2_.pts[i]);
                  // lazy update, perturbs old solution.
        // si_old -= (s0_.pts[i] = ( s1_.pts[i] + s2_.pts[i] )/2.0);

        res_out += si_old*si_old;
      }
      return res_out;
    }
    inline tjet_chart * tjchart_0() {return tjcharts[0];}
    inline tjet_chart * tjchart_f() {return tjcharts[nobs_h-1];}
};

class Lie_detector
{
  int palloc,
      Jalloc;

  ode_curve_observations &obs; // raw observational data, read from binary data files.

  public:

    ode_solspc_meta meta0; // specifications of ode system induced

    const int eor,
              ndep,
              ndim,
              ncrv,
              nobs,
              kor;
    int * const npts_per_crv;
    double * const pts,
           * const dnp1xu,
           ** const JFs_mat;

    ode_solspc_subset Sobs; // raw data set via obs, indexed data structure from file data
    curve_Lie_detector ** const cdets;

    ode_solution ** const sols,
                 ** const sols_h;
    ode_solcurve ** const curves,
                 ** const curves_h;

    ode_jetspc_meta jmeta_trivial;
    ode_trivial_curvejet ** const tcrvjets;

    Lie_detector(ode_curve_observations &obs_,bool verbose_=true) :
      palloc(obs_.dnp1xu_in!=NULL), Jalloc(obs_.JFs_in!=NULL),
      obs(obs_), // assume that obs is initialized.
      meta0(obs_.eor,obs_.ndep),
      eor(obs_.eor), ndep(obs_.ndep), ndim(meta0.ndim),
      ncrv(obs_.ncrv), nobs(obs_.nobs),
        kor( ((palloc)?(obs_.eor+1):(obs_.eor)) ),
      npts_per_crv(obs_.npts_per_crv),
      pts(obs_.pts_in), dnp1xu(obs_.dnp1xu_in),
      JFs_mat((Jalloc)?(Tmatrix_ptrs<double>(obs_.JFs_in,obs_.ndep*obs_.nobs,meta0.ndim)):(NULL)),
      Sobs(meta0,obs_.nobs,pts,dnp1xu,JFs_mat),
        sols(Sobs.sols),
        sols_h(new ode_solution*[nobs-ncrv]),
        curves(new ode_solcurve*[ncrv]),
        curves_h(new ode_solcurve*[ncrv]),
      cdets(new curve_Lie_detector*[ncrv]),
      jmeta_trivial(meta0, 2*(kor+1)-1 ),
        tcrvjets(new ode_trivial_curvejet*[ncrv])
      {
        int net_clcp = 0,
            net_Lmat_rows = 0,
            net_Rmat_rows = 0;
        double t0 = LD_threads::tic();
        #pragma omp parallel reduction(+:net_clcp,net_Lmat_rows,net_Rmat_rows)
        {
          LD_lu lu_t(jmeta_trivial.jor+1,jmeta_trivial.jor+1);
          double  ** const LUmat_t = lu_t.LUmat;
          #pragma omp for
          for (int icrv = 0; icrv < ncrv; icrv++)
          {
            const int ipts = sum_Tvec<int>(npts_per_crv,icrv);
            curves[icrv] = new ode_solcurve(icrv,meta0,npts_per_crv[icrv],
                                              Sobs.pts_mat+ipts,
                                              Sobs.sols+ipts,
                                              (palloc)?(Sobs.dnp1xu_mat+ipts):(NULL),
                                              (Jalloc)?(Sobs.JFs_tns+ipts):(NULL)
                                            );
            cdets[icrv] = new curve_Lie_detector(*(curves[icrv]),jmeta_trivial,palloc,Jalloc);
              curves_h[icrv] = &(cdets[icrv]->crv_h);
              tcrvjets[icrv] = &(cdets[icrv]->tcrvjets);

            for (int isol = 0, ipts_h=ipts-icrv; isol < cdets[icrv]->nobs_h; isol++)
            {
              sols_h[ipts_h+isol] = curves_h[icrv]->sols[isol];

              tcrvjets[icrv]->set_trivial_Amat(LUmat_t,isol);
              lu_t.decompose_A();

              for (int i = 0; i < ndep; i++)
                lu_t.solve_system(tcrvjets[icrv]->set_trivial_bvec(isol,i));

              tcrvjets[icrv]->tjets[isol]->set_sol_h( *(curves_h[icrv]->sols[isol]) );
                cdets[icrv]->tjcharts[isol]->reload_sol_h(); // sync sol_h and sol_h_alt
            }

            net_clcp += tcrvjets[icrv]->nxh;
            net_Lmat_rows += tcrvjets[icrv]->Lsvd_f1jet.Muse;
            net_Rmat_rows += tcrvjets[icrv]->R1svd_f1jet.Muse;
          }
        }
        double work_time = LD_threads::toc(t0);
        if (verbose_)
        {
          printf(
          "(Lie_detector::Lie_detector) "
          "computed 2*%d = %d SVDs (%.1f x %d and %.1f x %d, on average), "
          "determined %d*%d = %d jets (jor = %d) via %d LU decompositions (%d x %d) "
          "in %.3f seconds (%d threads)\n",
          ncrv, 2*ncrv,
            ((double)net_Lmat_rows)/((double)ncrv), tcrvjets[0]->fspc_1jet.perm_len,
            ((double)net_Rmat_rows)/((double)ncrv), tcrvjets[0]->fspc_1jet.ndof_full,
          ndep, net_clcp, ndep*net_clcp, jmeta_trivial.jor, net_clcp,
            jmeta_trivial.jor+1, jmeta_trivial.jor+1,
          work_time, LD_threads::numthreads()
          );
          print_details();
        }
      }
    ~Lie_detector()
    {
      if (JFs_mat!=NULL) delete [] JFs_mat;
      for (int i = 0; i < ncrv; i++)
      {
        delete cdets[i];
        delete curves[i];
      }
      delete [] cdets;
      delete [] sols_h;
      delete [] curves;
        delete [] curves_h;
        delete [] tcrvjets;
    }

    void print_details(bool detailed_=false, const char preamble_[]="", const char postscript_[]="\n")
    {
      printf(
      "%s"
      "(Lie_detector::print_details) "
      "n = %d ode system with q = %d degrees of freedom (b = %d), q*(%d+1) = %d observed derivative coordinates (including order 0).\n"
      " ncrv = %d sets of order jor = %d trivial flow models generated (%d total flow charts)"
      "%s",
      preamble_,
        eor, ndep, ndim, kor, ndep*(kor+1),
        ncrv, jmeta_trivial.jor, nobs-ncrv,
      postscript_);

      if (detailed_)
        for (int i = 0; i < ncrv; i++)
          printf(
          "   curve %d : L nulldim = %d, R1 nulldim = %d"
          "\n",
          i,
          tcrvjets[i]->Lsvd_f1jet.nulldim(),
          tcrvjets[i]->R1svd_f1jet.nulldim()
          );
    }

    inline void set_trivial_jets(ode_solcurve ** crvs_)
    {
      #pragma omp parallel
      {
        LD_lu lu_t(jmeta_trivial.jor+1,jmeta_trivial.jor+1);
        double  ** const LUmat_t = lu_t.LUmat;
        #pragma omp for
        for (int icrv = 0; icrv < ncrv; icrv++)
        {
          ode_solcurve &crv_i = *(crvs_[icrv]);
          ode_trivial_curvejet &tcjet_i = *(tcrvjets[icrv]);
          for (int isol = 0; isol < cdets[icrv]->nobs_h; isol++)
          {
            ode_solution &obs0_i = *(crv_i.sols[isol]),
                         &obs1_i = *(crv_i.sols[isol+1]);
            tcjet_i.tjets[isol]->set_trivial_Amat(LUmat_t,obs0_i,obs1_i);
            lu_t.decompose_A();

            for (int i = 0; i < ndep; i++)
              lu_t.solve_system(tcjet_i.tjets[isol]->set_trivial_bvec(i,obs0_i,obs1_i));

            tcrvjets[icrv]->tjets[isol]->set_sol_h( *(curves_h[icrv]->sols[isol]) );
              cdets[icrv]->tjcharts[isol]->reload_sol_h(); // sync sol_h and sol_h_alt
          }
        }
      }
    }

    inline void get_icrvisol_h_given_iobs_h(int &icrv_,int &isol_,int iobs_)
    {
      int iobs_local = iobs_,
          icrv_local = 0;
      do
      {
        iobs_local -= (npts_per_crv[icrv_local]-1);
        if (iobs_local < 0) // the point belongs to icrv
        {
          isol_ = iobs_local + (npts_per_crv[icrv_=icrv_local]-1); // set icrv and isol as well
          break;
        }
        else icrv_local++;
      } while (true);
    }
    inline ode_solution * get_sol_h_i(int &icrv_, int &isol_, int iobs_)
    {
      get_icrvisol_h_given_iobs_h(icrv_,isol_,iobs_);
      return curves_h[icrv_]->sols[isol_];
    }
    inline int comp_nobs_h() {return nobs-ncrv;}


};


#endif
