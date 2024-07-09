#ifndef MAT_LD_HH
#define MAT_LD_HH

#include "LD_framework.hh"

struct LD_SVD_space // assumes M>N
{
  LD_SVD_space(int M_, int N_, bool fast_flag_=false):
  U_gsl(gsl_matrix_alloc(M_,N_)), V_gsl(gsl_matrix_alloc(N_,N_)),
  X_gsl((fast_flag_)?(gsl_matrix_alloc(N_,N_)):(NULL)),
  s_gsl(gsl_vector_alloc(N_)),  w_gsl(gsl_vector_alloc(N_)) {}
  ~LD_SVD_space()
  {
    gsl_matrix_free(U_gsl); gsl_matrix_free(V_gsl);
    if (X_gsl != NULL) gsl_matrix_free(X_gsl);
    gsl_vector_free(s_gsl); gsl_vector_free(w_gsl);
  }

  gsl_matrix  * const U_gsl,
              * const V_gsl,
              * X_gsl = NULL;
  gsl_vector  * const s_gsl,
              * const w_gsl;

  inline void load_and_decompose(double **mat_)
  {
    LD_gsl::load_gsl_mat(U_gsl,mat_);
    decompose_loaded_U();
  }
  inline void decompose_loaded_U()
  {
    int status_decomp = (X_gsl == NULL)?
    (gsl_linalg_SV_decomp(U_gsl,V_gsl,s_gsl,w_gsl))
    :(gsl_linalg_SV_decomp_mod(U_gsl,X_gsl,V_gsl,s_gsl,w_gsl));
  }
  inline void decompose_loaded_U_fast()
  {
    int status_decomp = (X_gsl == NULL)?
    (gsl_linalg_SV_decomp_mod(U_gsl,X_gsl=gsl_matrix_alloc(V_gsl->size1,V_gsl->size2),V_gsl,s_gsl,w_gsl))
    :(gsl_linalg_SV_decomp_mod(U_gsl,X_gsl,V_gsl,s_gsl,w_gsl));
  }
  inline void unpack_s_VT(double *s_, double **VT_)
  {
    LD_gsl::unpack_gsl_vec(s_,s_gsl);
    LD_gsl::unpack_gsl_matT(VT_,V_gsl);
  }
  inline int unpack_rank_s_VT(double *s_, double **VT_)
  {
    unpack_s_VT(s_,VT_);
    double tol_eps = ((double) U_gsl->size1)*(LD_linalg::eps(s_[0]));
    int rank_out = 0;
    for (size_t i = 0; i < U_gsl->size2; i++)
      if (s_[i]>tol_eps) rank_out++;
      else break;
    return rank_out;
  }
  inline void unpack_s(double *s_, int offset_=0) {LD_gsl::unpack_gsl_vec(s_,s_gsl,offset_);}
  inline void unpack_VT(double ** VT_, int offset_=0)
  {
    for (size_t jcol = offset_, jjcol = 0; jcol < V_gsl->size2; jcol++, jjcol++)
      LD_gsl::unpack_gsl_mat_col_j(VT_[jjcol],V_gsl,jcol);
  }

  inline void set_Uij(int i_,int j_, double val_) {gsl_matrix_set(U_gsl,i_,j_,val_);}
  inline void load_Uirow(int i_,double *row_)
    {for (size_t j = 0; j < U_gsl->size2; j++) gsl_matrix_set(U_gsl,i_,j,row_[j]);}
  inline double get_Vij(int i_,int j_) {return gsl_matrix_get(V_gsl,i_,j_);}

  inline int nrows() {return U_gsl->size1;}
  inline int ncols() {return U_gsl->size2;}
};

struct LD_svd_file
{
  LD_svd_file(const char name_[]);
  ~LD_svd_file()
  {
    if (rank_vec_in != NULL) delete [] rank_vec_in;
    if (Smat_in != NULL) free_Tmatrix<double>(Smat_in);
    if (VTtns_in != NULL) free_T3tensor<double>(VTtns_in);
  }

  int hlen_in,
      ncrvs_in,
      ncols_in,
      ncol_use_in,
      * const header_in = &ncrvs_in;

  int * rank_vec_in = NULL;
  double  ** Smat_in = NULL,
          *** VTtns_in = NULL;
};

struct LD_matrix_svd_result
{
  LD_matrix_svd_result(int ncrvs_,int ncols_,int ncol_use_=0);
  LD_matrix_svd_result(LD_matrix &mat_): LD_matrix_svd_result(mat_.ncrvs_tot,mat_.net_cols) {}
  LD_matrix_svd_result(LD_svd_file svdfile_);

  ~LD_matrix_svd_result() {delete [] rank_vec; free_Tmatrix<double>(Smat);}

  LD_vector_bundle V_bndle;

  const int ncrvs,
            ncols;
  int ncol_use,
      * const rank_vec;
  double  ** const Smat,
          *** const VTtns = V_bndle.Vtns;

  void write_svd_results(const char name_[]);
  void read_svd_results(const char name_[]);

  inline void print_details(const char name_[] =  "Amat_SVD")
  {
    const char preamble[] = "(LD_matrix_svd_result::print_details)";
    char name_buf[strlen(preamble) + strlen(name_) + 20];
    sprintf(name_buf,"%s %s rank_vec", preamble, name_);
    LD_linalg::print_xT(name_buf,rank_vec,ncrvs);
    printf("ncols = %d, ncol_use = %d, min_nulldim = %d, max_nulldim = %d \n",
              ncols, ncol_use, min_nulldim(), max_nulldim());
  }

  inline int nulldim_i(int i_) {return ncol_use-rank_vec[i_];}
  inline int min_rank() {return LD_linalg::min_val<int>(rank_vec,ncrvs);}
  inline int max_rank() {return LD_linalg::max_val<int>(rank_vec,ncrvs);}
  inline int min_nulldim() {return ncol_use - max_rank();}
  inline int max_nulldim() {return ncol_use - min_rank();}

  inline int kappa_def() {int min_nulldim_val; return (min_nulldim_val=min_nulldim())?(min_nulldim_val):(1);}
  virtual double ** Kmat_crvi(int icrv_,int kappa_=0)
    {return VTtns[icrv_]+( ncol_use-((kappa_)?(kappa_):(kappa_def())) );}
};

struct LD_alternate_svd_result: public LD_matrix_svd_result
{
  LD_alternate_svd_result(int ncrvs_, int ncols_): LD_matrix_svd_result(ncrvs_,ncols_),
  V_bndle_alt(LD_vector_bundle(ncrvs_,ncols_)) {}
  LD_alternate_svd_result(LD_matrix &mat_): LD_alternate_svd_result(mat_.ncrvs_tot,mat_.ndof_full) {}
  ~LD_alternate_svd_result() {}

  LD_vector_bundle V_bndle_alt;

  double *** const VTtns_alt = V_bndle_alt.Vtns;

  void compute_restricted_curve_svds(LD_matrix &Amat_,LD_matrix_svd_result &Lsvd_,int nrows_,int rho_L_=0,bool verbose_=true);
  LD_matrix * make_AYmat_compute_restricted_curve_svds(LD_matrix &Amat_,LD_matrix_svd_result &Lsvd_,int nrows_,int rho_L_=0,bool verbose_=true);
  inline void set_VTtns_alt_Y_VAY(LD_matrix_svd_result &Lsvd_)
  {
    // assumes that ncol_use is already set, and that VTtns, VTtns_L have been set
    rho_L = ncol_use/( nvar = ( ncols/( ncol_L = Lsvd_.ncols ) ) );
    set_VTtns_alt_Y_VAY(Lsvd_.VTtns);
  }
  inline void set_VTtns_alt_Y_VAY(double ***VTtns_L_)
  {
    // assumes all parameters have been initialized, and computes VTtns_alt
    #pragma omp parallel for
    for (size_t icrv = 0; icrv < ncrvs; icrv++)
    {
      fill_Y_VAY_mat_icrv(VTtns_alt[icrv],VTtns_L_[icrv],VTtns[icrv]);
    }
  }

  void compute_regularized_curve_svds(LD_SVD_space &Bglb_svd,LD_matrix_svd_result &Asvd_,int nrows_B_,int kappa_A_=0,bool verbose_=true);
  inline void fill_KA_T_KB(double **KAi_T_KB_,double **KAi_,double **KBglb_,int ncol_B_)
  {
    // compute inner products between KA and KB columns
    for (size_t j = 0; j < ncol_use; j++)
      for (size_t l = 0; l < ncol_use; l++)
      {
        double  inner_KAicolj_KBcoll = 0.0;
        for (size_t k = 0; k < ncol_B_; k++) inner_KAicolj_KBcoll += KAi_[j][k]*KBglb_[l][k];
        KAi_T_KB_[l][j] = inner_KAicolj_KBcoll; // storing cols of KAi_T_KB adjacent in memory
      }
  }
  inline int comp_eff_rank_fill_sAB_KBA(double *sAi_B_,double **KB_Ai_,double **KAi_,double **KAi_T_KB_,double *skAi_,int ncol_B_)
  {
    for (size_t j = 0; j < ncol_use; j++)
      for (size_t k = 0; k < ncol_B_; k++)
      {
        double KB_Ai_jk = 0.0;
        for (size_t l = 0; l < ncol_use; l++) KB_Ai_jk += KAi_T_KB_[j][l]*KAi_[l][k];
        KB_Ai_[j][k] = KB_Ai_jk;
      }

    double mag_acc = 0.0;
    for (size_t j = 0; j < ncol_use; j++)
    {
      double mag_acc_j = 0.0;
      for (size_t k = 0; k < ncol_B_; k++) mag_acc_j += (KB_Ai_[j][k])*(KB_Ai_[j][k]);
      mag_acc += sqrt(mag_acc_j); // each column of KB_Ai_ has magnitude between 0 and 1

      // "singular values" of Ai under B regularization is weighted average of inner product contributions
      double w_acc = sAi_B_[j] = 0.0;
      for (size_t l = 0; l < ncol_use; l++)
      {
        double w_l = fabs(KAi_T_KB_[j][l]);
        w_acc += w_l;
        sAi_B_[j] += skAi_[l]*w_l;
      }
      sAi_B_[j] /= w_acc;
    }
    // effective rank is small if KB_Ai_ is close to orthogonal, i.e. if every column of KB_Ai_ is nearly unitary
    return ncol_use - ((int) mag_acc);
  }
  inline void fill_KBA(double **KB_Ai_,double **KAi_,double **KAi_T_KB_)
  {
    for (size_t j = 0; j < ncol_use; j++)
      for (size_t k = 0; k < ncols; k++)
      {
        double KB_Ai_jk = 0.0;
        for (size_t l = 0; l < ncol_use; l++) KB_Ai_jk += KAi_T_KB_[j][l]*KAi_[l][k];
        KB_Ai_[j][k] = KB_Ai_jk;
      }
  }
  inline void set_VTtns_alt_KBA(LD_matrix_svd_result &Asvd_)
  {
    // assumes that ncol_use is already set, and that VTtns, VTtns_L have been set
    #pragma omp parallel for
    for (size_t icrv = 0; icrv < ncrvs; icrv++)
    {
      fill_KBA(VTtns_alt[icrv],Asvd_.Kmat_crvi(icrv,ncol_use),VTtns[icrv]);
    }
  }

  virtual double ** Kmat_crvi(int icrv_,int kappa_=0)
    {return VTtns_alt[icrv_]+(ncol_use-(kappa_)?(kappa_):(kappa_def()));}

  protected:

    int nvar,
        ncol_L,
        rho_L;

    inline void fill_AY_mat_icrv(gsl_matrix *AY_t_,double **Amat_rows_,double **Ycols_L_,int nrows_)
    {
      for (size_t irow = 0; irow < nrows_; irow++)
        for (size_t ivar = 0, jcolA_start = 0, jcolAY = 0; ivar < nvar; ivar++, jcolA_start+=ncol_L)
          for (size_t kcolY_L = 0; kcolY_L < rho_L; kcolY_L++, jcolAY++)
          {
            double AY_irowjcol = 0.0;
            for (size_t pL = 0, jcolA = jcolA_start; pL < ncol_L; pL++, jcolA++)
              AY_irowjcol += Amat_rows_[irow][jcolA]*Ycols_L_[kcolY_L][pL];
            gsl_matrix_set(AY_t_,irow,jcolAY,AY_irowjcol);
          }
    }
    inline void fill_AY_mat_icrv(gsl_matrix *AY_t_,double **Amat_rows_,double **Ycols_L_,int nrows_,double **AY_mat_)
    {
      for (size_t irow = 0; irow < nrows_; irow++)
        for (size_t ivar = 0, jcolA_start = 0, jcolAY = 0; ivar < nvar; ivar++, jcolA_start+=ncol_L)
          for (size_t kcolY_L = 0; kcolY_L < rho_L; kcolY_L++, jcolAY++)
          {
            double AY_irowjcol = 0.0;
            for (size_t pL = 0, jcolA = jcolA_start; pL < ncol_L; pL++, jcolA++)
              AY_irowjcol += Amat_rows_[irow][jcolA]*Ycols_L_[kcolY_L][pL];
            gsl_matrix_set(AY_t_,irow,jcolAY,AY_mat_[irow][jcolAY]=AY_irowjcol);
          }
    }
    inline void fill_Y_VAY_mat_icrv(double **Y_VAY_T_mat_,double **Ycols_L_,double **VTrows_AY_)
    {
      for (size_t irow = 0; irow < ncol_use; irow++)
        for (size_t ivar = 0, jcolYV_T = 0, jcolVT_AY_start = 0; ivar < nvar; ivar++, jcolVT_AY_start+=rho_L)
          for (size_t pL = 0; pL < ncol_L; pL++, jcolYV_T++)
          {
            double Y_VAY_irowjcol = 0.0;
            for (size_t kcolY_L = 0, jcolVT_AY = jcolVT_AY_start; kcolY_L < rho_L; kcolY_L++, jcolVT_AY++)
              Y_VAY_irowjcol += VTrows_AY_[irow][jcolVT_AY]*Ycols_L_[kcolY_L][pL];
            Y_VAY_T_mat_[irow][jcolYV_T] = Y_VAY_irowjcol;
          }
    }

};

struct infinitesimal_generator: public ode_system
{
  infinitesimal_generator(function_space &fspace_, int eor_, int ndep_): ode_system(eor_,ndep_), fspace(fspace_) {}
  ~infinitesimal_generator() {}

  function_space &fspace;
  const int perm_len = fspace.perm_len,
            ndof_full = fspace.ndof_full,
            ndof_ODE = eor*ndep;
};

struct rspace_infinitesimal_generator: public infinitesimal_generator
{
  rspace_infinitesimal_generator(function_space &fspace_,double ***V_crv_tns_,int eor_,int kappa_,int ncol_use_):
  infinitesimal_generator(fspace_,eor_,fspace_.ndep),
  kappa(kappa_), ncol_use((ncol_use_)?(ncol_use_):(ndof_full)),
  V_crv_tns(V_crv_tns_), Kmat(V_crv_tns[0] + (ncol_use-kappa)),
  lamvec(new double[perm_len]), theta_wkspc(new double[ndof_full]),
  vxu_wkspc(vxu_workspace(nvar,fspace_.comp_ord_len())) {}

  ~rspace_infinitesimal_generator() {delete [] lamvec; delete [] theta_wkspc;}

  int kappa,
      ncol_use;
  double  *** const V_crv_tns,
          ** Kmat;

  void init_dudx_eval(int icrv_) {Kmat = V_crv_tns[icrv_] + (ncol_use-kappa);}

  void JacF_eval(double x_, double *u_, double **dls_out_) {} // do later
  void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) {} // do later

  inline void init_svd_default(LD_matrix_svd_result &svd_, int kappa_=0)
    {ncol_use = svd_.ncol_use; kappa = (kappa_)?(kappa_):(svd_.kappa_def()); init_dudx_eval(0);}

  protected:

    double  * const lamvec,
            * const theta_wkspc;

    vxu_workspace vxu_wkspc;
};

class r1space_infinitesimal_generator: public rspace_infinitesimal_generator
{

  double  * const xu,
          * const vx_vec = theta_wkspc;

  public:

    r1space_infinitesimal_generator(function_space &fspace_,double ***V_crv_tns_, int kappa_=1, int ncol_use_=0):
    rspace_infinitesimal_generator(fspace_,V_crv_tns_,1,kappa_,ncol_use_), xu(new double[nvar]) {}
    r1space_infinitesimal_generator(r1space_infinitesimal_generator &r1infgen_):
    r1space_infinitesimal_generator(r1infgen_.fspace,r1infgen_.V_crv_tns,r1infgen_.kappa,r1infgen_.ncol_use) {}

    ~r1space_infinitesimal_generator() {delete [] xu;}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      xu[0] = x_;
      for (size_t i = 0; i < ndep; i++) xu[i+1] = u_[i] + (dudx_[i] = 0.0);
      fspace.lamvec_eval(xu,lamvec,vxu_wkspc);
      double Vx2 = 0.0;
      for (size_t i_k = 0; i_k < kappa; i_k++)
      {
        double &vx_ik = vx_vec[i_k] = 0.0;
        for (size_t i = 0; i < perm_len; i++) vx_ik += Kmat[i_k][i]*lamvec[i];
        for (size_t idep = 0,i_theta = perm_len; idep < ndep; idep++)
          for (size_t iperm = 0; iperm < perm_len; iperm++,i_theta++)
            dudx_[idep] += vx_ik*Kmat[i_k][i_theta]*lamvec[iperm];
        Vx2 += vx_ik*vx_ik;
      }
      for (size_t i = 0; i < ndep; i++) dudx_[i] /= Vx2;
    }

    static void init_extended_observations(ode_curve_observations &obs_out_,ode_curve_observations &obs_in_)
    {
      // clear output observation buffers
      if (obs_out_.npts_per_crv != NULL) delete [] obs_out_.npts_per_crv;
      if (obs_out_.pts_in != NULL) delete [] obs_out_.pts_in;
      int ncrv = obs_out_.ncrv = obs_in_.ncrv,
          ndep = (obs_out_.ndep!=obs_in_.ndep)?(obs_out_.ndep=obs_in_.ndep):(obs_out_.ndep),
          nobs = (obs_out_.nobs!=obs_in_.nobs)?(obs_out_.nobs=obs_in_.nobs):(obs_out_.nobs),
          ndim_in = 1 + ndep*(obs_in_.eor + 1),
          ndim_out =  1 + ndep*(obs_out_.eor + 1);

      double  * const pts_chunk_in = obs_in_.pts_in,
              * const pts_chunk_out = obs_out_.pts_in = new double[ndim_out*nobs];

      obs_out_.npts_per_crv = new int[ncrv];
      LD_linalg::copy_vec<int>(obs_out_.npts_per_crv,obs_in_.npts_per_crv,ncrv);

      #pragma omp parallel for
      for (size_t iobs = 0; iobs < nobs; iobs++)
      {
        double  *pts_i_in = pts_chunk_in + ndim_in*iobs,
                *pts_i_out = pts_chunk_out + ndim_out*iobs;
        for (size_t idim = 0; idim < ndim_in; idim++) pts_i_out[idim] = pts_i_in[idim];
        for (size_t idim = ndim_in; idim < ndim_out; idim++) pts_i_out[idim] = 0.0;
      }
    }

    inline void extend_curve_observations(function_space_basis &basis_, double *theta_vec_, double *v_, double *pts_, int npts_)
    {
      const int ndim_full = basis_.ndim,
                eor_full = basis_.eor;
      for (size_t iobs = 0; iobs < npts_; iobs++)
      {
        double * const pts_i = pts_ + ndim_full*iobs;
        fspace.lamvec_eval(pts_i,lamvec,vxu_wkspc);
        LD_linalg::fill_vec<double>(theta_vec_,ndof_full,0.0);
        double Vx2 = 0.0;
        for (size_t i_k = 0; i_k < kappa; i_k++)
        {
          double &vx_ik = vx_vec[i_k] = 0.0;
          for (size_t i = 0; i < perm_len; i++) vx_ik += Kmat[i_k][i]*lamvec[i];
          for (size_t i = 0; i < ndof_full; i++) theta_vec_[i] += Kmat[i_k][i]*vx_ik;
          Vx2 += vx_ik*vx_ik;
        }
        for (size_t i = 0; i < ndof_full; i++) theta_vec_[i] /= Vx2;
        for (size_t k = 2, ioffset = ndep*k + 1; k <= eor_full; k++, ioffset+=ndep)
        {
          basis_.v_eval(pts_i,v_,theta_vec_,k-1);
          for (size_t i = 0, idkxu = i+ioffset; i < ndep; i++, idkxu++) pts_i[idkxu] = v_[idkxu-ndep]; // never using n'th order terms
        }
      }
    }
};

class rnspace_infinitesimal_generator: public rspace_infinitesimal_generator
{

  double  * const s_local,
          * const v_local,
          * const theta_local = theta_wkspc;

  function_space_basis &fbasis;

  public:

    rnspace_infinitesimal_generator(function_space_basis &fbasis_,double ***V_crv_tns_,int kappa_=1, int ncol_use_=0):
    rspace_infinitesimal_generator(fbasis_.fspc,V_crv_tns_,fbasis_.eor,kappa_,ncol_use_),
    s_local(new double[ndim]), v_local(new double[ndim]),
    fbasis(fbasis_) {}

    rnspace_infinitesimal_generator(rnspace_infinitesimal_generator &rninfgen_,function_space_basis &fbasis_):
    rnspace_infinitesimal_generator(fbasis_,rninfgen_.V_crv_tns,rninfgen_.kappa,rninfgen_.ncol_use) {}

    ~rnspace_infinitesimal_generator() {delete [] s_local; delete [] v_local;}

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      s_local[0] = x_;
      for (size_t i = 0; i < ndep; i++) s_local[i+1] = u_[i];
      fspace.lamvec_eval(s_local,lamvec,vxu_wkspc);
      LD_linalg::fill_vec<double>(theta_local,ndof_full,0.0);
      double Vx2 = 0.0;
      for (size_t i_k = 0; i_k < kappa; i_k++) // compute local parameter values
      {
        double  vx_ik = 0.0;
        for (size_t i = 0; i < perm_len; i++) vx_ik += Kmat[i_k][i]*lamvec[i];
        for (size_t i = 0; i < ndof_full; i++) theta_local[i] += vx_ik*Kmat[i_k][i];
        Vx2 += vx_ik*vx_ik;
      }
      for (size_t i = 0; i < ndof_full; i++) theta_local[i] /= Vx2; // normalize local parameters
      for (size_t i = ndep; i < ndof_ODE; i++) s_local[i+1] = dudx_[i-ndep] = u_[i]; // load 1 : n-1 derivative terms
      for (size_t i = ndof_ODE+1; i < ndim; i++) s_local[i] = 0.0; // pad n'th order terms with zeroes
      fbasis.v_eval(s_local,v_local,theta_local,eor-1);
      for (size_t i = ndof_ODE-ndep; i < ndof_ODE; i++) dudx_[i] = v_local[i+1];
    }

};

struct matrix_Lie_detector
{
  matrix_Lie_detector() {}
  ~matrix_Lie_detector() {}

  static int compare_relaxed_subspaces(LD_matrix_svd_result &svd_, LD_vspace_record &rec_, const char letter_[], double tol_=0.0)
  {
    const char fcn_name[] = "matrix_Lie_detector::compare_relaxed_subspaces";
    char name_buf[20];
    const int ncrvs = svd_.ncrvs;
    bool nsat_sf[ncrvs];
    int nnsat_s = 0,
        nnsat_diff[ncrvs];
    double diff_acc = 0.0;
    printf("(%s)\n  %s rank (1 x %d)\n", fcn_name,letter_,ncrvs);
    for (size_t i = 0; i < ncrvs; i++)
    {
      printf("%d ", svd_.rank_vec[i]);
      diff_acc += (double)(nnsat_diff[i] = svd_.rank_vec[i]-rec_.nV_spcvec[i]);
      nnsat_s+=(int)(nsat_sf[i]=(nnsat_diff[i] >= 0));
    }
    sprintf(name_buf,"\n  nsat_%s",letter_); LD_linalg::print_xT(name_buf,rec_.nV_spcvec,ncrvs);
    sprintf(name_buf,"  nsat_diff_%s",letter_); LD_linalg::print_xT(name_buf,nnsat_diff,ncrvs);
    if (tol_) printf("  nsat_smaller_%s = %d (of %d, tol = %.1e). ",letter_,nnsat_s,ncrvs,tol_);
    else printf(" nsat_smaller_%s = %d (of %d). ",letter_,nnsat_s,ncrvs);
    printf("Average diff: %.2f \n", diff_acc/((double) ncrvs));
    return nnsat_s;
  }

  static void compute_AYmat_curve_svds(LD_matrix_svd_result &svd_,LD_matrix &Amat_,LD_Theta_bundle &Tbndle_,bool verbose_=true)
  {
    const int ncrvs = Amat_.ncrvs_tot,
              nrows_max = Amat_.max_nrow_curve(),
              ncol_full = Amat_.net_cols;
    int * const rank_vec = svd_.rank_vec;
    double  ** const Smat_out = svd_.Smat,
            *** const VTtns_out = svd_.VTtns,
            **** const Attns = Amat_.Attns;
    LD_Theta_space ** const Tspaces = Tbndle_.Tspaces;

    double  mrows_acc = 0.0,
            ncols_acc = 0.0,
            t0 = LD_threads::tic();
    #pragma omp parallel reduction(+:mrows_acc,ncols_acc)
    {
      LD_svd svd_t(nrows_max,ncol_full);
      double ** const Umat_t = svd_t.Umat;
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrvs; icrv++)
      {
        LD_Theta_space &Tspc_i = *(Tspaces[icrv]);
        const int mrows_Ai = Amat_.nrows_mat_i(icrv),
                  ncols_Ai = Tspc_i.ndof_use;
        Tspc_i.post_multiply_bases(Umat_t,Amat_.Amat_crv_i(icrv),mrows_Ai);
        svd_t.decompose_U(mrows_Ai,ncols_Ai);
        rank_vec[icrv] = svd_t.unpack_rank_svec_VTmat(Smat_out[icrv],VTtns_out[icrv]);

        mrows_acc += (double) mrows_Ai; ncols_acc += (double) ncols_Ai;
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
      printf("(matrix_Lie_detector::compute_relaxed_curve_svds) computed %d svds (%.1f x %.1f, on avg.) in %.4f seconds (%d threads)\n",
        ncrvs,mrows_acc/((double)ncrvs),ncols_acc/((double)ncrvs),work_time, LD_threads::numthreads());
  }

  static void compute_global_svd(LD_SVD_space &svd_, LD_matrix &mat_, int ncrvs_cmp_=1, bool verbose_=true)
  {
    const char fcn_name[] = "(matrix_Lie_detector::compute_global_svd)";
    bool dim_mismatch = false;
    const int nrows_net = svd_.nrows(),
              ncrvs_use = ( dim_mismatch = (nrows_net % ncrvs_cmp_) )?(1):(ncrvs_cmp_),
              nrows_per_crv = nrows_net / ncrvs_use;
    if (dim_mismatch)
      printf("%s WARNING - svd space rows (%d) indivisable by input curves (%d). Defaulting to first %d rows of input matrix\n",
        fcn_name, nrows_net, ncrvs_cmp_, nrows_net);

    double  **** const Attns = mat_.Attns,
            t0 = LD_threads::tic();
    #pragma omp parallel for
    for (size_t icrv = 0; icrv < ncrvs_use; icrv++)
    {
      for (size_t irow = 0, irow_glb = icrv*nrows_per_crv; irow < nrows_per_crv; irow++, irow_glb++)
        svd_.load_Uirow(irow_glb,Attns[icrv][0][irow]);
    }

    svd_.decompose_loaded_U_fast();
    double work_time = LD_threads::toc(t0);

    if (verbose_)
      printf("%s loaded and computed %d x %d svd in %.4f seconds (%d threads)\n",
        fcn_name,nrows_net,svd_.ncols(),work_time,LD_threads::numthreads());
  }

  static void compute_curve_svds(LD_matrix &mat_,LD_matrix_svd_result &mat_svd_,int nrows_,int ncol_use_=0,bool verbose_=true)
  {
    const int ncrvs = mat_.ncrvs_tot,
              ncols = mat_svd_.ncol_use = (ncol_use_)?(ncol_use_):(LD_linalg::min_T<int>(mat_.net_cols,mat_svd_.ncols));
    int * const rank_vec_out = mat_svd_.rank_vec;
    double  ** const Smat_out = mat_svd_.Smat,
            *** const VTtns_out = mat_svd_.VTtns,
            t0 = LD_threads::tic();
    #pragma omp parallel
    {
      LD_SVD_space svd_t(nrows_,ncols);
      #pragma omp for
      for (size_t i = 0; i < ncrvs; i++)
      {
        svd_t.load_and_decompose(mat_.Amat_crv_i(i));
        rank_vec_out[i] = svd_t.unpack_rank_s_VT(Smat_out[i],VTtns_out[i]);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
      printf("(matrix_Lie_detector::compute_curve_svds) computed %d svds (%d x %d) in %.4f seconds (%d threads)\n",
        ncrvs,nrows_,ncols,work_time, LD_threads::numthreads());
  }

  template <class INFGN, class BSIS> static void extend_ode_observations(ode_curve_observations &obs_out_, INFGN &infgen_, BSIS **bases_, bool verbose_=true)
  {
    const int ncrv = obs_out_.ncrv,
              ndof_full = infgen_.ndof_full;
    int * const npts_per_crv = obs_out_.npts_per_crv;
    double t0 = LD_threads::tic();
    #pragma omp parallel
    {
      INFGN infgen_t(infgen_);
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      double  theta_vec_t[ndof_full],
              v_t[basis_t.ndim];
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrv; icrv++)
      {
        infgen_t.init_dudx_eval(icrv);
        infgen_t.extend_curve_observations(basis_t,theta_vec_t,v_t,obs_out_.pts_icrv(icrv),npts_per_crv[icrv]);
      }
    }
    double work_time = LD_threads::toc(t0);
    if (verbose_)
      printf("(matrix_Lie_detector::extend_ode_observations) extended %d integral curves (%d tot snaps, dim1 = %d -> dim2 = %d) in %.4f seconds (%d threads)\n",
        ncrv,obs_out_.nobs,infgen_.ndim,1+obs_out_.ndep*(obs_out_.eor+1),work_time, LD_threads::numthreads());
  }

};

#endif
