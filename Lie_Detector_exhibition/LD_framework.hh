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

  void write_observed_solutions(const char name_[]);
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

  inline void set_solcurve_ICs(ode_solcurve **crvs_)
  {
    for (size_t i = 0; i < ncrv; i++)
      for (size_t idim = 0; idim <= ndof_ODE; idim++)
        pts_IC[i][idim] = crvs_[i]->pts0[idim];
  }
  void set_random_ICs(LD_rng rng_, const double *IC_range_);
  void generate_solution_curves(ode_integrator &integrator_, const double *indep_range_);
  void write_solution_curves(const char name_[]);

  void generate_JFs();
  void write_JFs(const char name_[]);

  void generate_dnp1xu();
  void write_dnp1xu(const char name_[]);
};

struct input_ode_observations: public ode_curve_observations
{
  input_ode_observations(const char name_[]);
  input_ode_observations(const char name_[], const char name1_[]):
    input_ode_observations(name_) {read_additional_observations(name1_);}
  input_ode_observations(const char name_[], const char name1_[], const char name2_[]):
    input_ode_observations(name_,name1_) {read_additional_observations(name2_);}
  ~input_ode_observations();

  char * const name;
  void print_details();
  void read_additional_observations(const char name_addtl_[]);
};

struct LD_observations_set: public solspc_data_chunk
{
  LD_observations_set(ode_solspc_meta &meta_, input_ode_observations input_);
  ~LD_observations_set();

  const int ncrvs_tot;

  int * const npts_per_crv;
  double  *** const pts_tns;
  ode_solcurve ** const curves;

  double  *** dnp1xu_tns = NULL,
          **** JFs_crv = NULL;

  void load_additional_inputs(input_ode_observations input_, bool overwrite_basics_=false);

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
  inline void configure_0maxmag_0pi05_domain(orthopolynomial_space &fspace_, double h_max_tru_=1.0, double h_max_scl_=0.99)
  {
    double **sme_g = Tmatrix<double>(2,ndim);
    get_solspace_mag_extrema(sme_g);
    fspace_.set_0maxmag_0pi05_domain(sme_g,h_max_tru_,h_max_scl_);
    free_Tmatrix<double>(sme_g);
  }

  void get_solspace_val_extrema(double **sve_g_);
  void get_solspace_mag_extrema(double **sme_g_);

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
  LD_matrix(function_space &fspc_, LD_observations_set &Sset_, int dim_cnstr_);
  ~LD_matrix();

  const int dim_cnstr,
            net_rows = dim_cnstr*nobs_full,
            net_eles = net_rows*ndof_full;
  double  **** const Attns,
          *** const Atns,
          ** const Amat = Atns[0],
          * const Avec = Amat[0];

  inline int nrows_mat_i(int i_) {return dim_cnstr*npts_per_crv[i_];}
  inline double ** Amat_crv_i(int i_) {return Attns[i_][0];}
  inline void print_matrix_i(int i_)
    {LD_linalg::print_A("curve i sub matrix",Attns[i_][0],nrows_mat_i(i_),ndof_full);}
  inline void print_matrix_i_submat_j(int i_,int j_)
    {LD_linalg::print_A("curve i, pts j sub matrix",Attns[i_][j_],dim_cnstr,ndof_full);}

  void write_matrix(const char name_[]);
  void read_matrix(const char name_[]);

  inline int min_nrow_curve() {return dim_cnstr*min_npts_curve();}
  inline int max_nrow_curve() {return dim_cnstr*max_npts_curve();}
};

class LD_R_matrix: public LD_matrix
{
  public:
    LD_R_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.eor*fspc_.ndep) {}
    ~LD_R_matrix() {}

    template <class TBSIS> void populate_R_matrix(TBSIS **bases_)
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
          fill_Rn_rows(chunk_t,sol_i,Atns[iobs]);
        }
      }
    }

  protected:

    inline void fill_Rn_rows(partial_chunk &chunk_, ode_solution &sol_, double **Rmat_i_)
    {
      int i_dof = 0;
      double  * const u = sol_.u,
              * const dxu = sol_.dxu,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
        {
          fill_x_Rn_columns(i_dof,Rmat_i_,dxu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L]);
          i_dof++;
        }

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
          {
            fill_u_Rn_columns(i_dof,idep,Rmat_i_,Lambda_xu_mat[idep+1][i_L],Jac_utheta_vdxu_mat[i_L]);
            i_dof++;
          }
    }
    inline void fill_x_Rn_columns(int i_dof, double **Rmat_, double *dxu_, double &lambda_, double **Jac_xtheta_i_vdxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++) Rmat_[idep][i_dof] = dxu_[idep]*lambda_; // 1st order contribution
      for (size_t k = 1, i_R = ndep; k < eor; k++)
        for (size_t idep = 0; idep < ndep; idep++, i_R++)
          Rmat_[i_R][i_dof] = dxu_[i_R]*lambda_ - Jac_xtheta_i_vdxu_[k-1][idep]; // subtract contribution to v1n
    }
    inline void fill_u_Rn_columns(int i_dof, int idep, double **Rmat_, double &lambda_, double *del_utheta_i)
    {
      Rmat_[idep][i_dof] = -(lambda_);
      for (size_t k = 1, i_R = idep + ndep; k < eor; k++, i_R += ndep)
        Rmat_[i_R][i_dof] = -(del_utheta_i[k-1]);
    }
};

class LD_P_matrix: public LD_matrix
{
  public:
    LD_P_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep) {}
    ~LD_P_matrix() {}

    template <class TBSIS> void populate_P_matrix(TBSIS **bases_)
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
          fill_P_rows(chunk_t,sol_i,Atns[iobs]);
        }
      }
    }

  protected:

    inline void fill_P_rows(partial_chunk &chunk_, ode_solution &sol_, double **Pmat_i_)
    {
      int i_dof = 0;
      double  * const dnp1xu = sol_.dnp1xu,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
        {
          fill_x_P_columns(i_dof,Pmat_i_,dnp1xu,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L][eor-1]);
          i_dof++;
        }

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
          {
            fill_u_P_columns(Pmat_i_[idep][i_dof],Jac_utheta_vdxu_mat[i_L][eor-1]);
            i_dof++;
          }
    }
    inline void fill_x_P_columns(int i_dof, double **Pmat_, double *dnp1xu_, double &lambda_, double *Jac_xtheta_i_vdnxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++)
        Pmat_[idep][i_dof] = dnp1xu_[idep]*lambda_ - Jac_xtheta_i_vdnxu_[idep];
    }
    inline void fill_u_P_columns(double &Pmat_ij_,double par_utheta_i) {Pmat_ij_ = -(par_utheta_i);}
};

class LD_G_matrix: public LD_matrix
{
  public:
    LD_G_matrix(function_space &fspc_,LD_observations_set &Sset_): LD_matrix(fspc_,Sset_,fspc_.ndep) {}
    ~LD_G_matrix() {}

    template <class TBSIS> void populate_G_matrix(TBSIS **bases_)
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
          fill_G_rows(chunk_t,sol_i,Atns[iobs]);
        }
      }
    }

  protected:

    inline void fill_G_rows(partial_chunk &chunk_, ode_solution &sol_, double **Gmat_i_)
    {
      int i_dof = 0;
      double  ** const JFs = sol_.JFs,
              * const lambda_x_vec = chunk_.Jac_mat[0],
              ** const Lambda_xu_mat = chunk_.Jac_mat,
              *** const Jac_xtheta_vdxu_tns = chunk_.C_x,
              ** const Jac_utheta_vdxu_mat = chunk_.C_u;

      for (size_t i_L = 0; i_L < perm_len; i_L++)
        if (dof_tun_flags[i_L])
        {
          fill_x_G_columns(i_dof,Gmat_i_,JFs,lambda_x_vec[i_L],Jac_xtheta_vdxu_tns[i_L]);
          i_dof++;
        }

      for (size_t idep = 0; idep < ndep; idep++)
        for (size_t i_L = 0; i_L < perm_len; i_L++)
          if (dof_tun_flags[i_L])
          {
            fill_u_G_columns(i_dof,idep,Gmat_i_,JFs,Lambda_xu_mat[idep+1][i_L],Jac_utheta_vdxu_mat[i_L]);
            i_dof++;
          }
    }
    inline void fill_x_G_columns(int i_dof,double **Gmat_,double **JFs_,double &lambda_,double **Jac_xtheta_i_vdxu_)
    {
      for (size_t idep = 0; idep < ndep; idep++) // tangent space constraint
      {
        Gmat_[idep][i_dof] = JFs_[idep][0]*lambda_;
        for (size_t k = 1, idim = nvar; k <= eor; k++)
          for (size_t jdep = 0; jdep < ndep; jdep++, idim++)
            Gmat_[idep][i_dof] += JFs_[idep][idim]*Jac_xtheta_i_vdxu_[k-1][jdep];
      }
    }
    inline void fill_u_G_columns(int i_dof,int idep,double **Gmat_,double **JFs_,double &lambda_,double *del_utheta_i)
    {
      for (size_t jdep = 0; jdep < ndep; jdep++) // tangent space contraint
      {
        Gmat_[jdep][i_dof] = JFs_[jdep][1+idep]*lambda_;
        for (size_t k = 1, idim = nvar+idep; k <= eor; k++, idim+=ndep)
          Gmat_[jdep][i_dof] += JFs_[jdep][idim]*del_utheta_i[k-1];
      }
    }
};



#endif
