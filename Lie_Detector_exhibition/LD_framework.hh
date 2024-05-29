#ifndef LD_EXHIB_HH
#define LD_EXHIB_HH

#include "LD_function_space.hh"
#include "LD_aux.hh"
#include "LD_io.hh"

struct ode_curve_observations
{
  ode_curve_observations() {}
  ode_curve_observations(int eor_, int ndep_, int nobs_);
  ode_curve_observations(int eor_, int ndep_, int ncrv_, int nobs_);
  ode_curve_observations(ode_solspc_meta &meta_, int ncrv_, int nobs_): ode_curve_observations(meta_.eor,meta_.ndep,ncrv_,nobs_) {}
  ode_curve_observations(const char name_[]);
  ode_curve_observations(const char name_[], const char name1_[]):
    ode_curve_observations(name_) {read_additional_observations(name1_);}
  ode_curve_observations(const char name_[], const char name1_[], const char name2_[]):
    ode_curve_observations(name_,name1_) {read_additional_observations(name2_);}
  ~ode_curve_observations();

  int eor,
      ndep,
      ncrv,
      nobs,
      * const ode_meta = &eor,
      * const obs_meta = &ncrv;

  int * npts_per_crv = NULL;
  double  * pts_in = NULL,
          * dnp1xu_in = NULL,
          * JFs_in = NULL;

  inline int nobs_ioffset(int icrv_)
    {int ioffset = 0; for (size_t i = 0; i < icrv_; i++) ioffset += npts_per_crv[i]; return ioffset;}
  inline double *pts_icrv(int icrv_) {return pts_in + (1+ndep*(eor+1))*nobs_ioffset(icrv_);}

  inline bool palloc() {return dnp1xu_in!=NULL;}
  inline bool Jalloc() {return JFs_in!=NULL;}

  void read_basic_observations(const char name_addtl_[], bool force_overwrite_=false);
  void read_additional_observations(const char name_addtl_[]);

  void write_observed_solutions(const char name_[]);

  private:

    inline int comp_ndim(int eor_, int ndep_) {return 1 + ndep_*(eor_+1);}

};

struct generated_ode_observations: public ode_curve_observations
{
  generated_ode_observations(ode_system &ode_, int nc_, int np_);
  ~generated_ode_observations();

  ode_system &ode;

  const int npts,
            ndim = ode.ndim,
            ndof_ODE = ode.ndep*ode.eor;
  double ** const pts_IC;

  template <class INFGN, class INTGR> void parallel_generate_solution_curves(INFGN &infgn_, INTGR &intgr_, const double *indep_range_)
  {
    double t0 = LD_threads::tic();
    #pragma omp parallel
    {
      INFGN infgn_t(infgn_);
      INTGR intgr_t(infgn_t,intgr_);

      intgr_t.del_t = (indep_range_[1]-indep_range_[0])/((double)(npts-1));

      double  * const u_state_t = intgr_t.get_u_state(),
              ** const integr_wkspc_t = Tmatrix<double>(npts,ndof_ODE);
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrv; icrv++)
      {
        for (size_t i_dof = 0; i_dof < ndof_ODE; i_dof++) u_state_t[i_dof] = pts_IC[icrv][i_dof+1];
        intgr_t.init_curve_integration(icrv);
        intgr_t.set_and_solve_time(indep_range_[0],indep_range_[1],npts,integr_wkspc_t);
        intgr_t.unpack_time_sol(indep_range_[0],npts,integr_wkspc_t,pts_IC[icrv]);
      }
      free_Tmatrix<double>(integr_wkspc_t);
    }
    double work_time = LD_threads::toc(t0);
    printf("(generated_ode_observations::parallel_generate_solution_curves) exponentiated %d integral curves (%d net snaps, %d degrees of freedom) in %.4f seconds (%d threads)\n",
    ncrv, nobs, intgr_.ndof_ODE, work_time, LD_threads::numthreads());
  }
  inline void set_solcurve_ICs(ode_solcurve **crvs_)
  {
    for (size_t i = 0; i < ncrv; i++)
      for (size_t idim = 0; idim <= ndof_ODE; idim++)
        pts_IC[i][idim] = crvs_[i]->pts0[idim];
  }
  void set_random_ICs(LD_rng rng_, const double *IC_range_);
  void generate_solution_curves(ode_integrator &integrator_, const double *indep_range_);
  // void write_solution_curves(const char name_[]);

  void generate_JFs();
  void write_JFs(const char name_[]);

  void generate_dnp1xu();
  void write_dnp1xu(const char name_[]);
};

struct LD_observations_set: public solspc_data_chunk
{
  LD_observations_set(ode_solspc_meta &meta_, ode_curve_observations input_);
  LD_observations_set(ode_solspc_meta &meta_, int ncrv_, int npts_, bool palloc_=false, bool Jalloc_=false);
  ~LD_observations_set();

  const int ncrvs_tot;

  int * const npts_per_crv;
  double  *** const pts_tns;
  ode_solcurve ** const curves;

  double  *** dnp1xu_tns = NULL,
          **** JFs_crv = NULL;

  void load_additional_inputs(ode_curve_observations input_, bool overwrite_basics_=false);

  inline double * get_default_IC_indep_range(int icrv_=0)
  {
    ode_solcurve &crv_i = *(curves[icrv_]);
    indep_range[0] = crv_i.pts0[0];
    indep_range[1] = crv_i.pts_mat[crv_i.nobs-1][0];
    return indep_range;
  }
  inline void configure_centered_domain(orthopolynomial_space &fspace_, double h_min_=-0.99, double h_max_=0.99)
  {
    double **sve_g = Tmatrix<double>(2,ndim);
    get_solspace_val_extrema(sve_g);
    fspace_.set_centered_domain(sve_g,h_min_,h_max_);
    free_Tmatrix<double>(sve_g);
  }
  inline void configure_center_mass_domain(orthopolynomial_space &fspace_, double h_min_=-0.99, double h_max_=0.99)
  {
    double **sve_g = Tmatrix<double>(2,ndim),
            scm_g[ndim];
    get_solspace_val_extrema(sve_g);
    get_solspace_center_mass(scm_g);
    fspace_.set_center_mass_domain(sve_g,scm_g,h_min_,h_max_);
    free_Tmatrix<double>(sve_g);
  }
  inline void configure_0maxmag_0pi05_domain(orthopolynomial_space &fspace_, double h_max_tru_=1.0, double h_max_scl_=0.99)
  {
    double **sme_g = Tmatrix<double>(2,ndim);
    get_solspace_mag_extrema(sme_g);
    fspace_.set_0maxmag_0pi05_domain(sme_g,h_max_tru_,h_max_scl_);
    free_Tmatrix<double>(sme_g);
  }

  void get_solspace_val_extrema(double **sve_g_);
  void get_solspace_mag_extrema(double **sme_g_);
  void get_solspace_center_mass(double *scm_g_);

  inline ode_solution * get_icrv_jsol(int i_, int j_) {return curves[i_]->sols[j_];}
  inline void print_curve_i(int i_) {curves[i_]->print_curve();}
  inline void print_icrv_jsol(int i_, int j_) {curves[i_]->print_jsol(j_);}

  inline int min_npts_curve() {return LD_linalg::min_val<int>(npts_per_crv,ncrvs_tot);}
  inline int max_npts_curve() {return LD_linalg::max_val<int>(npts_per_crv,ncrvs_tot);}

  private:

    double indep_range[2];
};

struct LD_experiment
{
  LD_experiment(LD_observations_set &Sset_): Sset(Sset_) {}
  ~LD_experiment() {}

  LD_observations_set &Sset;

  const int ncrvs_tot = Sset.ncrvs_tot,
            nobs_full = Sset.nobs;
  int * const npts_per_crv = Sset.npts_per_crv;
  double  * const pts_chunk = Sset.pts_chunk,
          ** const pts_mat = Sset.pts_mat,
          *** const pts_tns = Sset.pts_tns;

  ode_solcurve ** const curves = Sset.curves;
  ode_solution ** const sols_full = Sset.sols;

  inline int min_npts_curve() {return Sset.min_npts_curve();}
  inline int max_npts_curve() {return Sset.max_npts_curve();}
};

struct LD_matrix: public function_space_element, public LD_experiment
{
  LD_matrix(function_space &fspc_, LD_observations_set &Sset_, int dim_cnstr_, int net_cols_);
  LD_matrix(function_space &fspc_, LD_observations_set &Sset_, int dim_cnstr_):
    LD_matrix(fspc_,Sset_,dim_cnstr_,fspc_.ndof_full) {}
  ~LD_matrix();

  const int dim_cnstr,
            net_rows = dim_cnstr*nobs_full,
            net_cols,
            net_eles = net_rows*net_cols;
  double  **** const Attns,
          *** const Atns,
          ** const Amat = Atns[0],
          * const Avec = Amat[0];

  inline int nrows_mat_i(int i_) {return dim_cnstr*npts_per_crv[i_];}
  inline double ** Amat_crv_i(int i_) {return Attns[i_][0];}
  inline void print_matrix_i(int i_)
    {LD_linalg::print_A("curve i sub matrix",Attns[i_][0],nrows_mat_i(i_),net_cols);}
  inline void print_matrix_i_submat_j(int i_,int j_)
    {LD_linalg::print_A("curve i, pts j sub matrix",Attns[i_][j_],dim_cnstr,net_cols);}

  void write_matrix(const char name_[]);
  void read_matrix(const char name_[]);

  template <class TBSIS> void populate_matrix(TBSIS **bases_)
  {
    LD_linalg::fill_vec<double>(Avec,net_eles,0.0);
    #pragma omp parallel
    {
      int tid = LD_threads::thread_id();
      TBSIS &basis_t = *(bases_[tid]);
      partial_chunk &chunk_t = basis_t.partials;
      #pragma omp for
      for (size_t iobs = 0; iobs < nobs_full; iobs++)
      {
        ode_solution &sol_i = *(sols_full[iobs]);
        basis_t.fill_partial_chunk(sol_i.pts);
        fill_A_rows(chunk_t,sol_i,Atns[iobs]);
      }
    }
  }

  inline int min_nrow_curve() {return dim_cnstr*min_npts_curve();}
  inline int max_nrow_curve() {return dim_cnstr*max_npts_curve();}

  static void concatenate_matrices(LD_matrix &mat_cat_, LD_matrix &mat1_, LD_matrix &mat2_, bool block_cat_=false)
  {
    if (LD_matrix::net_cols_eq(mat_cat_,mat1_,mat2_))
    {
      if ((!block_cat_)&&(LD_matrix::net_dim_constr_eq(mat_cat_,mat1_,mat2_)))
        LD_matrix::enmesh_curve_matrices(mat_cat_,mat1_,mat2_);
      else
        LD_matrix::block_concat_matrices(mat_cat_,mat1_,mat2_);
    }
    else printf("(LD_matrix::concatenate_matrices) ERROR - incompatible matrices\n");
  }

  protected:

    virtual void fill_A_rows(partial_chunk &chunk_, ode_solution &sol_, double **Amat_i_) {}

    static void enmesh_curve_matrices(LD_matrix &mat0_, LD_matrix &mat1_, LD_matrix &mat2_)
    {
      const size_t  ncrvs_tot_ = mat0_.ncrvs_tot,
                    net_cols_ = mat0_.net_cols,
                    sze_ele_ = sizeof(double),
                    submat1_len = mat1_.dim_cnstr*net_cols_,
                    submat2_len = mat2_.dim_cnstr*net_cols_,
                    submat1_sze = submat1_len*sze_ele_,
                    submat2_sze = submat2_len*sze_ele_;
      int * const nppc0 = mat0_.npts_per_crv,
          * const nppc1 = mat1_.npts_per_crv,
          * const nppc2 = mat2_.npts_per_crv;
      double  **** const A0ttns = mat0_.Attns,
              **** const A1ttns = mat1_.Attns,
              **** const A2ttns = mat2_.Attns;
      bool trunc_flag = false;
      #pragma omp parallel for reduction( || : trunc_flag)
      for (size_t icrv = 0; icrv < ncrvs_tot_; icrv++)
      {
        int nobs_icrv = nppc0[icrv];
        if (!( ((nppc0[icrv])==(nppc1[icrv])) && ((nppc1[icrv])==(nppc2[icrv])) ))
        {
          nobs_icrv = LD_linalg::min_T_3way<int>(nppc0[icrv],nppc1[icrv],nppc2[icrv]);
          if (nobs_icrv == nppc0[icrv]) trunc_flag = (trunc_flag)||(true);
        }
        for (size_t iobs = 0; iobs < nobs_icrv; iobs++)
        {
          memcpy(A0ttns[icrv][iobs][0],A1ttns[icrv][iobs][0],submat1_sze);
          memcpy(A0ttns[icrv][iobs][0]+submat1_len,A2ttns[icrv][iobs][0],submat2_sze);
        }
      }
      if (trunc_flag)
        printf("(LD_matrix::enmesh_curve_matrices) WARNING - snapshots discarded from concatenated matrix.\n");
    }

    static void block_concat_matrices(LD_matrix &mat0_, LD_matrix &mat1_, LD_matrix &mat2_)
    {
      const size_t  ncrvs_tot_ = mat0_.ncrvs_tot,
                    net_cols_ = mat0_.net_cols,
                    sze_ele_ = sizeof(double),
                    submat0_len = mat0_.dim_cnstr*net_cols_,
                    submat1_len = mat1_.dim_cnstr*net_cols_,
                    submat2_len = mat2_.dim_cnstr*net_cols_;
      int * const nppc0 = mat0_.npts_per_crv,
          * const nppc1 = mat1_.npts_per_crv,
          * const nppc2 = mat2_.npts_per_crv;
      double  **** const A0ttns = mat0_.Attns,
              **** const A1ttns = mat1_.Attns,
              **** const A2ttns = mat2_.Attns;
      bool trunc_flag = false;
      #pragma omp parallel for reduction( || : trunc_flag)
      for (size_t icrv = 0; icrv < ncrvs_tot_; icrv++)
      {
        size_t  len0_i = submat0_len*nppc0[icrv],
                len1_i = submat1_len*nppc1[icrv],
                len2_i = submat2_len*nppc2[icrv],
                sze0_i = sze_ele_*len0_i,
                sze1_i = sze_ele_*len1_i,
                sze2_i = sze_ele_*len2_i;
        if (sze0_i >= sze1_i)
        {
          memcpy(A0ttns[icrv][0][0],A1ttns[icrv][0][0],sze1_i);
          if (sze0_i >= sze1_i+sze2_i) memcpy(A0ttns[icrv][0][0]+len1_i,A2ttns[icrv][0][0],sze2_i);
          else
          {
            memcpy(A0ttns[icrv][0][0]+len1_i,A2ttns[icrv][0][0],sze0_i-sze1_i);
            trunc_flag = (trunc_flag)||(true);
          }
        }
        else
        {
          memcpy(A0ttns[icrv][0][0],A2ttns[icrv][0][0],sze0_i);
          trunc_flag = (trunc_flag)||(true);
        }
      }
      if (trunc_flag)
        printf("(LD_matrix::block_concat_matrices) WARNING - inadequately large concatenation matrix. Rows truncated.\n");
    }

    static bool compatible_curve_concat(LD_matrix &m0_, LD_matrix &m1_, LD_matrix &m2_)
      {return (LD_matrix::net_cols_eq(m0_,m1_,m2_))&&(LD_matrix::net_dim_constr_eq(m0_,m1_,m2_));}
    static bool net_dim_constr_eq(LD_matrix &m0_, LD_matrix &m1_, LD_matrix &m2_)
      {return m0_.dim_cnstr == (m1_.dim_cnstr+m2_.dim_cnstr);}
    static bool net_cols_eq(LD_matrix &m0_, LD_matrix &m1_, LD_matrix &m2_)
      {return ((m0_.net_cols==m1_.net_cols)&&(m1_.net_cols==m2_.net_cols));}
};

#endif
