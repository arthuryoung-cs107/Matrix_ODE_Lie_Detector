#ifndef LIE_DETECTOR_HH
#define LIE_DETECTOR_HH

#include "ode_curve_observations.hh"
#include "LD_trivial_flows.hh"

// struct LD_observation_meta
// {
//   const int &kor_obs; // highest order derivative observation available for analysis
//
//   LD_observation_meta(int &kor_obs_) : kor_obs(kor_obs_) {}
//   ~LD_observation_meta() {}
//
// };


class curve_Lie_detector
{
  const int kor,
            alen;
  bool palloc,
       Jalloc;
  double * const avec_crv;

  ode_trivial_soljet tsj_crv;
  jet_chart tjc_crv;

  const static int combine_flag;
  const static bool truncate_Hermite_exp;

  inline double combine_2coords(double c1_,double c2_)
  {
   // return (c1_+c2_)/2.0; // lazy update, perturbs old solution.
   // return c2_; // aggressive update. Reinforces dominant components, still one sided.
   return (combine_flag)?(c2_):((c1_+c2_)/2.0);
  }
  inline double combine_3coords(double c1_,double c2_,double c3_)
  {
   // return (c1_+c2_+c3_)/3.0; // lazy update, perturbs old solution.
   // return (c2_+c3_)/2.0; // aggressive update. Reinforces dominant components.
   return (combine_flag)?((c2_+c3_)/2.0):((c1_+c2_+c3_)/3.0);
  }

  public:

    // static int kor_upd;

    const int nobs_f,
              nobs_h;

    ode_solcurve &crv, // observed trajectory
                  crv_h; // induced trajectory of colocation points
      ode_solution ** const sols,
                   ** const sols_h;

    ode_trivial_curvejet tcrvjet; // deterministic Hermite jet approximation, organized over whole curve.
      ode_trivial_soljet ** const tsoljets; // array of pointers to individual jets

    jet_chart ** const tjcharts; // charts over which we compute trivial Hermite jets

    /*
      [CONSTRUCTOR] on instantiation, this structure receives an ode_solcurve (defined in LD_ode.hh)
      and proceeds by allocating space for a trivial flow model of the input data, which consists of
      observed point solutions to an ordinary differential equation belonging to a common integral
      curve paramaterized by the independent variable x. This involves the curve_Lie_detector
      instantiating corresponding structures crv_h (ode_solcurve, defined in LD_ode.hh) and tcrvjet
      (ode_trivial_curvejet, defined in LD_trivial_flows.hh); the former consisting of intermediate
      point solutions, sols_h (ode_solution, defined in LD_ode.hh), lying between sequential point
      solutions sols on the given integral curve; the latter comprised by a system of Hermite
      interpolants, tsoljets (ode_trivial_curvejet, defined in LD_trivial_flows.hh), each centered
      at one of the intermediate points, sols_h.

      Together, the points, sols and sols_h, and the interpolants, tsoljets, induce tjcharts
      (jet_chart, defined in LD_ode.hh), a workspace to be interpreted as a neighborhood centered
      at a given sol_h such that the observations sol_0,sol_1 \in sols lie on the closure of the
      neighborhood
    */
    curve_Lie_detector(ode_solcurve &crv_,ode_jetspc_meta &jmeta_trivial_,bool palloc_,bool Jalloc_) :
      kor(ode_jetspc_meta::comp_kor(jmeta_trivial_.jor)), alen( ode_soljet::alen_trivial_jet(jmeta_trivial_.jor,crv_.ndep) ),
      palloc(palloc_), Jalloc(Jalloc_),
        avec_crv( new double[alen] ),
        tsj_crv( jmeta_trivial_,avec_crv,*(crv_.sol_0()),*(crv_.sol_f()) ),
        tjc_crv( tsj_crv.solh,tsj_crv.sol0,tsj_crv.sol1, false ),
      nobs_f(crv_.nobs), nobs_h(crv_.nobs-1),
      crv(crv_),
        sols(crv_.sols),
      crv_h(crv_.icrv,crv_.meta,nobs_h,palloc,Jalloc),
        sols_h(crv_h.sols),
      tcrvjet(crv_,jmeta_trivial_),
        tsoljets(tcrvjet.tjets),
      tjcharts(new jet_chart*[nobs_h])
      {
        for (int i = 0; i < nobs_h; i++)
          tjcharts[i] = new jet_chart(*(crv_h.sols[i]),*(crv.sols[i]),*(crv.sols[i+1]));
      }
    ~curve_Lie_detector()
    {
      for (int i = 0; i < nobs_h; i++) delete tjcharts[i];
      delete [] tjcharts;
      delete [] avec_crv;
    }

    inline double compute_staggered_Hermite_flow_curve(ode_solution **solso_, LD_lu &lu_,LD_spectral_tvfield &tvf_ ,ode_solution **solsi_, int kor_upd_)
    {
      double res_out = 0.0;

      // compute deterministic Hermite jet directly from input point solutions
      (tsoljets[0])->set_and_solve_Hermite_jets( *(sols_h[0]) , lu_, *(solsi_[0]), *(solsi_[1]) );
      // correct 0'th collocation point using input trivial vector field
      tvf_.set_sol_dkxu((tjcharts[0])->solh_alt,kor,*(sols_h[0]));

      for (int iobs = 1; iobs < nobs_h; iobs++) // loop over interior points
      {
        // compute deterministic Hermite jet directly from input point solutions
        (tsoljets[iobs])->set_and_solve_Hermite_jets( *(sols_h[iobs]) , lu_, *(solsi_[iobs]), *(solsi_[iobs+1]) );
        // correct i'th collocation point using input trivial vector field
        tvf_.set_sol_dkxu((tjcharts[iobs])->solh_alt, kor, *(sols_h[iobs])); // correct i'th collocation point using trivial vector field

        ((tjcharts[iobs])->sol1_alt).x = (solsi_[iobs])->x;
        compute_staggered_Hermite_flow(  (tjcharts[iobs])->sol1_alt,
          lu_, tvf_,
          *(tsoljets[iobs]),
          (tjcharts[iobs-1])->solh_alt, // left hand knot, already corrected
          (tjcharts[iobs])->sol0_alt, // induced collocation point
          (tjcharts[iobs])->solh_alt ); // right hand knot, freshly corrected

        // res_out += (solso_[iobs])->copy_comp_sol_residual((tjcharts[iobs])->sol1_alt, kor);
        res_out += combine_two_sols_dxu(*(solso_[iobs]),*(solsi_[iobs]),(tjcharts[iobs])->sol1_alt,kor_upd_);
        // res_out += (combine_flag)?( (solso_[iobs])->copy_comp_sol_residual((tjcharts[iobs])->sol1_alt, kor) )
        //               : ( combine_two_sols_dxu(*(solso_[iobs]),*(solsi_[iobs]),(tjcharts[iobs])->sol1_alt) ) ;
      }
      // We now have updated values for the interior points. Use these to generate updates at edges

      // crvo_.sols[0]->copy_sol( *(crvi_.sols[0]) );
      // generate LEFT edge update, uses smoothened solh_0 and UPDATED sol_1
      ((tjcharts[0])->sol1_alt).x = (solsi_[0])->x;
      compute_staggered_Hermite_flow(  (tjcharts[0])->sol1_alt,
       lu_, tvf_,
       *(tsjet_start()),
       (tjcharts[0])->solh_alt, // left hand knot, already corrected
       (tjcharts[0])->sol0_alt, // induced collocation point
       (tjcharts[1])->sol1_alt ); // right hand knot, already corrected (and updated)

      // res_out += (solso_[0])->copy_comp_sol_residual((tjcharts[0])->sol1_alt, kor);
      res_out += combine_two_sols_dxu(*(solso_[0]),*(solsi_[0]),(tjcharts[0])->sol1_alt,kor_upd_);
      // res_out += (combine_flag)?( (solso_[0])->copy_comp_sol_residual((tjcharts[0])->sol1_alt, kor) )
      //               : ( combine_two_sols_dxu(*(solso_[0]),*(solsi_[0]),(tjcharts[0])->sol1_alt) ) ;

      // crvo_.sols[nobs_h]->copy_sol( *(crvi_.sols[nobs_h]) );
      // generate RIGHT edge update, uses smoothened solh_f and UPDATED sol_f-1
      ((tjchart_final())->sol1_alt).x = (solsi_[nobs_h])->x;
      compute_staggered_Hermite_flow(  (tjchart_final())->sol1_alt,
       lu_, tvf_,
       *(tsjet_final()),
       (tjcharts[nobs_h-1])->sol1_alt, // left hand knot, already corrected (and updated)
       (tjchart_final())->sol0_alt, // induced collocation point
       (tjcharts[nobs_h-1])->solh_alt ); // right hand knot, already corrected

      // res_out += (solso_[nobs_h])->copy_comp_sol_residual((tjchart_final())->sol1_alt, kor);
      res_out += combine_two_sols_dxu(*(solso_[nobs_h]),*(solsi_[nobs_h]),(tjchart_final())->sol1_alt,kor_upd_);
      // res_out += (combine_flag)?( (solso_[nobs_h])->copy_comp_sol_residual((tjchart_final())->sol1_alt, kor) )
      //               : ( combine_two_sols_dxu(*(solso_[nobs_h]),*(solsi_[nobs_h]),(tjchart_final())->sol1_alt) ) ;

      return res_out;
    }
    inline void compute_staggered_Hermite_flow(ode_solution &sol_out_, // output. Assume that x is set
      LD_lu &lu_,LD_spectral_tvfield &tvf_,
      ode_trivial_soljet &tsj_,
      ode_solution &sol0_,
      ode_solution &solh_,
      ode_solution &sol1_)
    {
      // construct Hermite interpolant between given knots
      tsj_.set_and_solve_Hermite_jets(lu_,sol0_,sol1_);

      // record u coordinate values of estimated intermediate solution induced by Hermite jet over knots.
      solh_.x = 0.5*( sol0_.x + sol1_.x ); // set x value of new colocation point
      tsj_.set_sol_h_u(solh_); // set u values of new colocation point using raw Hermite jet
        tvf_.set_sol_dkxu(solh_, kor); // improved estimate of derivatives at new collocation point
      tsj_.set_jet_given_solh(solh_); // use improved derivative estimates to update corresponding jet coeffficients

      // exponentiate resultant Hermite jet
      // if (truncate_Hermite_exp) tsj_.exp_u_trivial( sol_out_.u, sol_out_.x - solh_.x, kor);
      // else tsj_.exp_u_trivial( sol_out_.u, sol_out_.x - solh_.x );

      tsj_.exp_u_trivial( sol_out_.u, sol_out_.x - solh_.x, kor);

      tvf_.set_sol_dkxu(sol_out_, kor); // apply vfield to resultant solution for improved derivatives
    }
    inline double combine_two_sols_dxu(ode_solution &snew_,ode_solution &sold_,ode_solution &sexp_, int kor_upd_)
    {
      snew_.x = sexp_.x;
      double res_out = 0.0;
      for (int i = 1; i < snew_.nvar; i++)
      {
        double si_diff = sold_.pts[i];
        // si_diff -= (
        //   snew_.pts[i] = (combine_flag)?(sexp_.pts[i]):(0.5*( sexp_.pts[i]+sold_.pts[i] ))
        //   );
        si_diff -= ( snew_.pts[i] = sexp_.pts[i] );
        // si_diff -= (snew_.pts[i] = 0.5*( sexp_.pts[i]+sold_.pts[i] ));
        res_out += si_diff*si_diff;
      }
      int ivar = sold_.nvar,
          kor_cap = LD_linalg::min_T_3way<int>(kor_upd_,kor,snew_.eor);
      for (int k = 1 ; k <= kor_cap; k++)
        for (int i = 0; i < snew_.ndep; i++,ivar++)
        {
          double si_diff = sold_.pts[ivar];
          si_diff -= (
            snew_.pts[ivar] = (combine_flag)?(sexp_.pts[ivar]):(0.5*( sexp_.pts[ivar]+sold_.pts[ivar] ))
            );
          // si_diff -= (snew_.pts[ivar] = sexp_.pts[ivar]);
          // si_diff -= (snew_.pts[ivar] = 0.5*( sexp_.pts[ivar]+sold_.pts[ivar] ));
          res_out += si_diff*si_diff;
        }
      for (int k = kor_cap+1; k <= snew_.eor; k++)
        for (int i = 0; i < snew_.ndep; i++,ivar++)
        {
          snew_.pts[ivar] = sold_.pts[ivar];

          // double si_diff = sold_.pts[ivar];
          // si_diff -= (snew_.pts[ivar] = sold_.pts[ivar] );
          // si_diff -= (snew_.pts[ivar] = 0.5*( sexp_.pts[ivar]+sold_.pts[ivar] ));
          // res_out += si_diff*si_diff;
        }
      if (kor>snew_.eor)
      {
        if (kor_upd_>snew_.eor)
        {
          for (int i = 0; i < snew_.ndep; i++)
          {
            double si_diff = sold_.dnp1xu[i];
            si_diff -= (
              snew_.dnp1xu[i] = (combine_flag)?(sexp_.dnp1xu[i]):(0.5*( sexp_.dnp1xu[i]+sold_.dnp1xu[i] ))
              );
            // si_diff -= (snew_.dnp1xu[i] = sexp_.dnp1xu[i]);
            // si_diff -= (snew_.dnp1xu[i] = 0.5*( sexp_.dnp1xu[i]+sold_.dnp1xu[i] ));
            res_out += si_diff*si_diff;
          }
        }
        else for (int i = 0; i < snew_.ndep; i++) snew_.dnp1xu[i] = sold_.dnp1xu[i];
      }
      return res_out;
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
        si_old -= (s0_.pts[i] = combine_3coords(s1_.pts[i],s2_.pts[i],s3_.pts[i]));
        res_out += si_old*si_old;
      }
      if (kor>s0_.eor)
        for (int i = 0; i < s0_.ndep; i++)
        {
          double si_old = s0_.dnp1xu[i];
          si_old -= (s0_.dnp1xu[i] = combine_3coords(s1_.dnp1xu[i],s2_.dnp1xu[i],s3_.dnp1xu[i]));
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
        si_old -= (s0_.pts[i] = combine_2coords(s1_.pts[i],s2_.pts[i]));
        res_out += si_old*si_old;
      }
      if (kor>s0_.eor)
        for (int i = 0; i < s0_.ndep; i++)
        {
          double si_old = s0_.dnp1xu[i];
          si_old -= (s0_.dnp1xu[i] = combine_2coords(s1_.dnp1xu[i],s2_.dnp1xu[i]));
          res_out += si_old*si_old;
        }
      return res_out;
    }

    inline jet_chart * tjchart_start() {return tjcharts[0];}
    inline jet_chart * tjchart_final() {return &(tjc_crv);}
    inline ode_trivial_soljet * tsjet_start() {return tsoljets[0];}
    inline ode_trivial_soljet * tsjet_final() {return &(tsj_crv);}

    inline jet_chart * tjchart_0() {return tjcharts[0];}
    inline jet_chart * tjchart_f() {return tjcharts[nobs_h-1];}
    inline ode_trivial_soljet * tsjet_0() {return tsoljets[0];}
    inline ode_trivial_soljet * tsjet_f() {return tsoljets[nobs_h-1];}

    inline ode_solution * sol_wrk0_i(int i_)
    { return ( i_<nobs_h )?( &(tjcharts[i_]->sol0_alt) ):( (i_==nobs_h)?(&(tjc_crv.sol0_alt)):(NULL) ); }
    inline ode_solution * sol_wrk1_i(int i_)
    { return ( i_<nobs_h )?( &(tjcharts[i_]->sol1_alt) ):( (i_==nobs_h)?(&(tjc_crv.sol1_alt)):(NULL) ); }
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
    curve_Lie_detector ** const cdets; // curve Lie detectors

    ode_solution ** const sols,
                 ** const sols_h,
                 ** const sols_h_alt;
    jet_chart ** const tjcharts;
    ode_solcurve ** const curves,
                 ** const curves_h;

    ode_jetspc_meta jmeta_trivial;
    ode_trivial_curvejet ** const tcrvjets;
    ode_trivial_soljet ** const tsoljets;

    /*
      [CONSTRUCTOR] on instantiation, this structure receives obs a set of ode_curve_observations
      (defined in ode_curve_observations.hh), containing raw data in the form of observed point solutions
      to an ordinary differential equation, each of which belong to curves over the independent variable,
      x, on the ODE's solution space.

      The Lie_detector proceeds by building Sobs, an ode_solspc_subset (defined in LD_ode.hh)
      indexing the observational data into a collection of sols, which are ode_solution (defined in
      LD_ode.hh) data structures storing point solution data. Each curve's set of ode_solutions
      are organized into an ode_solcurve (defined in LD_ode.hh), corresponding to a single flow over the
      ODE's solution space, passing through each of assigned point solutions and generated by the trivial
      vector field, dx s = [ 1 ; dx u^(n) ] = [ 1 ; dx u ; ... ].

      For each curve, the Lie_detector builds a curve_Lie_detector, a data structure receiving the
      corresponding ode_solcurve as input. Each curve_Lie_detector allocates space for a system of trivial
      flow models comprised of deterministic Hermite interpolants.

      The Lie_detector gets the data structures and workspaces produced by each curve_Lie_detector and
      solves the associated linear system determining each Hermite interpolant, centered at the midpoint
      x_hat between sequential point solutions
        s_0 = [ x_0 ; u^(n)_0 ] , s_1 = [ x_1 ; u^(n)_1 ] : x_1 > x_0
      so that the resultant Hermite interpolant induces an estimate of the intermediate point solution
        s_hat = [ x_hat ; u^(n)_hat ],
      where x_hat = 0.5*(x_0 + x_1), and u^(n)_hat = [ u_hat ; dx u_hat ; ... ]. The resultant trivial flow
      model is n--smooth, and consists of smoothly compatible (up to order n) local parameterizations of
      the observed integral curves.
    */

    Lie_detector(ode_curve_observations &obs_,bool verbose_=true) :
      palloc(obs_.dnp1xu_in!=NULL), Jalloc(obs_.JFs_in!=NULL),
      obs(obs_),
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
        sols_h_alt(new ode_solution*[nobs-ncrv]),
          tjcharts(new jet_chart*[nobs-ncrv]),
        curves(new ode_solcurve*[ncrv]),
        curves_h(new ode_solcurve*[ncrv]),
      cdets(new curve_Lie_detector*[ncrv]),
      jmeta_trivial(meta0, 2*(kor+1)-1 ),
        tcrvjets(new ode_trivial_curvejet*[ncrv]),
          tsoljets(new ode_trivial_soljet*[nobs-ncrv])
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
              tcrvjets[icrv] = &(cdets[icrv]->tcrvjet);

            for (int isol = 0, ipts_h=ipts-icrv; isol < cdets[icrv]->nobs_h; isol++,ipts_h++)
            {
              // initializing global point cloud colocation point data
              sols_h[ipts_h] = cdets[icrv]->sols_h[isol];
              tjcharts[ipts_h] = cdets[icrv]->tjcharts[isol];
                sols_h_alt[ipts_h] = &(tjcharts[ipts_h]->solh_alt);
              tsoljets[ipts_h] = cdets[icrv]->tsoljets[isol];

              // solving Hermite jet via LU decomposition, then sets solh accordingly
              tsoljets[ipts_h]->set_and_solve_Hermite_jets( tjcharts[ipts_h]->solh , lu_t, tjcharts[ipts_h]->sol0, tjcharts[ipts_h]->sol1 );
              tjcharts[ipts_h]->reload_sol_h(); // sync sol_h and sol_h_alt
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
        delete [] sols_h_alt;
          delete [] tjcharts;
      delete [] curves;
        delete [] curves_h;
        delete [] tcrvjets;
          delete [] tsoljets;
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
    double combine_trivial_flows(ode_solcurve **crvs_alt_,ode_solcurve **crvs_)
    {
      double res_out = 0.0;
      #pragma omp parallel for reduction(+:res_out)
      for (int icrv = 0; icrv < ncrv; icrv++)
      {
        res_out += cdets[icrv]->recombine_flows(*(crvs_alt_[icrv]),*(crvs_[icrv]));
      }
      return res_out;
    }
    inline void set_trivial_jets(ode_solution ** sols_h_,ode_solcurve ** crvs_)
    {
      const int nobs_h = nobs-ncrv;
      #pragma omp parallel
      {
        LD_lu lu_t(jmeta_trivial.jor+1,jmeta_trivial.jor+1);
        double  ** const LUmat_t = lu_t.LUmat;
        int icrv_t,
            isol_t;
        #pragma omp for
        for (int iobs = 0; iobs < nobs_h; iobs++)
        {
          get_icrvisol_h_given_iobs_h(icrv_t,isol_t,iobs);
          ode_solution &obs0_i = *(crvs_[icrv_t]->sols[isol_t]),
                       &obs1_i = *(crvs_[icrv_t]->sols[isol_t+1]);
          ode_trivial_soljet &tjet_i = *(tcrvjets[icrv_t]->tjets[isol_t]);

          tjet_i.set_and_solve_Hermite_jets(*( sols_h_[iobs]), lu_t, obs0_i,obs1_i );

        }
      }
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

            tcjet_i.tjets[isol]->set_and_solve_Hermite_jets(*(curves_h[icrv]->sols[isol]), lu_t, obs0_i,obs1_i );
            cdets[icrv]->tjcharts[isol]->reload_sol_h();

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

struct Lie_detector_experiment
{

  Lie_detector &det;

  const int nobs,
            ncrv,
            nobs_h;

  // S (solution space) members, points and sets (curves)
  ode_solcurve ** const curves,
               ** const curves_h;
  ode_solution ** const sols,
               ** const sols_h,
               ** const sols_h_alt;

  // trivial jet space jets, curves, and charts
  ode_jetspc_meta &tjmeta;
  jet_chart ** const tjcharts;
  ode_trivial_curvejet ** const tcjets;
  ode_trivial_soljet ** const tsjets;

  // Lie detectors on a single trivial flow
  curve_Lie_detector ** const cdets;

  Lie_detector_experiment(Lie_detector &det_) :
    det(det_),
    nobs(det_.nobs), ncrv(det_.ncrv), nobs_h(det_.nobs - det_.ncrv),
    curves(det.curves),
    curves_h(det.curves_h),
      sols(det.sols),
      sols_h(det.sols_h),
      sols_h_alt(det.sols_h_alt),
    tjmeta(det.jmeta_trivial),
    tjcharts(det.tjcharts),
      tcjets(det.tcrvjets),
      tsjets(det.tsoljets),
    cdets(det.cdets)
    {}
  ~Lie_detector_experiment() {}

  inline void set_sol_trivial_vfield(ode_solution &sol_,LD_spectral_tvfield &tvf_)
  {
    // use 0'th order Hermite jet info to compute local theta via pseudoinverse of Rk matrix
    tvf_.comp_spectral_theta_local(sol_.x,sol_.u);
      tvf_.eval_prn_theta_image(sol_,det.kor); // the result is an improved estimate of the derivatives at interior colocation point
  }
  inline void smooth_trivial_soljet(ode_trivial_soljet &tsjet_, ode_solution &sol_out_, LD_spectral_tvfield  &tvf_)
  {
    set_sol_trivial_vfield(sol_out_,tvf_);
    // now modify the Hermite jet coefficients to reflect R pseudoinverse estimates
    tsjet_.set_jet_given_solh(sol_out_);
  }
  inline void smooth_trivial_soljet(ode_trivial_soljet &tsjet_, ode_solution &sol_out_, LD_spectral_tvfield  &tvf_, ode_solution &sol_in_)
  {
    sol_out_.copy_xu(sol_in_);
    smooth_trivial_soljet(tsjet_,sol_out_,tvf_);
  }

};

class Lie_detector_multinomial_solmodel
{
  ode_solution &sol_alt;
  orthopolynomial_space &fspc;

  public:

    Lie_detector_multinomial_solmodel(ode_solution &sol_alt_, orthopolynomial_space &fspc_) :
      sol_alt(sol_alt_), fspc(fspc_) {}
    ~Lie_detector_multinomial_solmodel() {}

};

class Lie_detector_multinomial_crvmodel
{
  curve_Lie_detector &cdet;
  orthopolynomial_space &fspc;

  double * const svec,
         ** const VTmat;
  public:

    Lie_detector_multinomial_crvmodel(curve_Lie_detector &cdet_, orthopolynomial_space &fspc_, double * svec_,double ** VTmat_) :
      cdet(cdet_), fspc(fspc_),
      svec(svec_), VTmat(VTmat_)
    {

    }
    ~Lie_detector_multinomial_crvmodel() {}

};

class multinomial_experiment : public Lie_detector_experiment, public function_space_element
{
  protected:

    double  ** const svecs_crvs,
            *** const VTmats_crvs;

  public:

    orthopolynomial_space &fspace0;

    Lie_detector_multinomial_crvmodel ** const cmdls;

    multinomial_experiment(Lie_detector &det_, orthopolynomial_space &fspace0_) :
      Lie_detector_experiment(det_), function_space_element(fspace0_),
      svecs_crvs(Tmatrix<double>(ncrv,ndof_full)),
      VTmats_crvs(T3tensor<double>(ncrv,ndof_full,ndof_full)),
      fspace0(fspace0_),
      cmdls(new Lie_detector_multinomial_crvmodel*[ncrv])
    {
      for (int i = 0; i < ncrv; i++)
        cmdls[i] = new Lie_detector_multinomial_crvmodel( *(det_.cdets[i]), fspace0, svecs_crvs[i], VTmats_crvs[i]);
    }

    ~multinomial_experiment()
    {
      free_Tmatrix<double>(svecs_crvs);
      free_T3tensor<double>(VTmats_crvs);
      for (int i = 0; i < ncrv; i++) delete cmdls[i];
      delete [] cmdls;
    }

};

struct global_multinomial_experiment : public multinomial_experiment
{
  protected:

    const int ncod; // number of linear constraints per point solution in nullspace problem
    LD_svd Asvd_global; // SVD of globally encoded matrix

    inline double * Umat_jsol_global(int jsol_) {return Asvd_global.Umat[jsol_*ncod];}

  public:

    // theta (parameter) space variables
    double * const svec_global,
           ** const VTmat_global;

    global_multinomial_experiment(Lie_detector &det_, orthopolynomial_space &fspace0_, int ncod_) :
      multinomial_experiment(det_,fspace0_),
      ncod(ncod_), Asvd_global(ncod_*nobs,ndof_full),
      svec_global(Asvd_global.svec),
      VTmat_global(Tmatrix<double>(ndof_full,ndof_full))
    {

    }

    ~global_multinomial_experiment()
    {
      free_Tmatrix<double>(VTmat_global);
    }

    inline int telescope_decompose_global_matrix(double **VTmg_,LD_svd &Asvdg_,LD_encoder &Aenc_,int nobs_)
    {
      #pragma omp parallel
      {
        double ** const Umatg_t = Asvdg_.Umat,
                        wvec_t[ndof_full];

        #pragma omp for
        for (int icrv = 0; icrv < ncrv; icrv++)
        {
          int iobs_i=0;
          for (int i = 0; i < icrv; i++) iobs_i += det.npts_per_crv[i];

          LD_svd Asvd_icrv( Aenc_.ncod*det.npts_per_crv[icrv], Asvdg_.Nuse ,
            Umatg_t+(Aenc_.ncod*iobs_i), VTmats_crvs[icrv], svecs_crvs[icrv], wvec_t );

          Asvd_icrv.decompose_U();
          Asvd_icrv.transpose_V();
        }

        #pragma omp for
        for (int icrv = 0; icrv < ncrv; icrv++)
        {
          double ** const Umat_i = Umatg_t + (Asvdg_.Nuse*icrv);
          for (int i = 0; i < Asvdg_.Nuse; i++)
            for (int j = 0; j < Asvdg_.Nuse; j++)
              Umat_i[i][j] = svecs_crvs[icrv][i]*VTmats_crvs[icrv][i][j];
        }
      }

      const int Muse_old = Asvdg_.Muse;

      Asvdg_.decompose_U(ncrv*Asvdg_.Nuse,Asvdg_.Nuse);
      Asvdg_.unpack_VTmat(VTmat_global);

      int rank_out = Asvdg_.rank();
      Asvdg_.set_use_dims(Muse_old,Asvdg_.Nuse);

      return rank_out;
    }

};

struct Lie_detector_digital_twin
{
  Lie_detector &det;
  ode_solcurve_chunk &Stwn;

  Lie_detector_digital_twin(Lie_detector &det_, ode_solcurve_chunk &Stwn_) :
    det(det_), Stwn(Stwn)
    {}

  ~Lie_detector_digital_twin()
    {}



};


#endif
