#ifndef MAT_LD_HH
#define MAT_LD_HH

#include "LD_framework.hh"

#include <gsl/gsl_linalg.h>

struct LD_SVD_space // assumes M>N
{
  LD_SVD_space(int M_, int N_):
  U_gsl(gsl_matrix_alloc(M_,N_)), V_gsl(gsl_matrix_alloc(N_,N_)), s_gsl(gsl_vector_alloc(N_)), w_gsl(gsl_vector_alloc(N_)) {}
  ~LD_SVD_space() {gsl_matrix_free(U_gsl); gsl_matrix_free(V_gsl); gsl_vector_free(s_gsl); gsl_vector_free(w_gsl);}

  gsl_matrix  * const U_gsl,
              * const V_gsl;
  gsl_vector  * const s_gsl,
              * const w_gsl;
  inline void load_and_decompose(double **mat_)
  {
    LD_gsl::load_gsl_mat(U_gsl,mat_);
    decompose_loaded_U();
  }
  inline void decompose_loaded_U()
    {int status_decomp = gsl_linalg_SV_decomp(U_gsl,V_gsl,s_gsl,w_gsl);}
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
};

struct LD_matrix_svd_result
{
  LD_matrix_svd_result(int ncrvs_,int ncols_): ncrvs(ncrvs_), ncols(ncols_),
  rank_vec(new int[ncrvs_]), Smat(Tmatrix<double>(ncrvs_,ncols_)), VTtns(T3tensor<double>(ncrvs_,ncols_,ncols_)) {}
  LD_matrix_svd_result(LD_matrix &mat_): LD_matrix_svd_result(mat_.ncrvs_tot,mat_.net_cols) {}
  ~LD_matrix_svd_result() {delete [] rank_vec; free_Tmatrix<double>(Smat); free_T3tensor<double>(VTtns);}

  const int ncrvs,
            ncols;
  int ncol_use = ncols,
      * const rank_vec;
  double  ** const Smat,
          *** const VTtns;

  void write_svd_results(const char name_[]);
  void read_svd_results(const char name_[]);

  inline void print_details()
  {
    LD_linalg::print_xT("(LD_matrix_svd_result::print_details) rank_vec",rank_vec,ncrvs);
    printf("ncols = %d, min_nulldim = %d, max_nulldim = %d \n", ncols, min_nulldim(), max_nulldim());
  }

  inline int nulldim_i(int i_) {return ncol_use-rank_vec[i_];}
  inline int min_rank() {return LD_linalg::min_val<int>(rank_vec,ncrvs);}
  inline int max_rank() {return LD_linalg::max_val<int>(rank_vec,ncrvs);}
  inline int min_nulldim() {return ncol_use - max_rank();}
  inline int max_nulldim() {return ncol_use - min_rank();}

  inline int kappa_def() {int min_nulldim_val; return (min_nulldim_val=min_nulldim())?(min_nulldim_val):(1);}
  inline double ** Kmat_crvi(int icrv_,int kappa_=0) {return VTtns[icrv_]+(ncol_use-(kappa_)?(kappa_):(kappa_def()));}
};

struct infinitesimal_generator: public ode_system
{
  infinitesimal_generator(function_space &fspace_): ode_system(1,fspace_.ndep), fspace(fspace_) {}
  ~infinitesimal_generator() {}

  function_space &fspace;
  const int perm_len = fspace.perm_len,
            ndof_full = fspace.ndof_full;
};

struct rspace_infinitesimal_generator: public infinitesimal_generator
{
  rspace_infinitesimal_generator(function_space &fspace_,double ***V_crv_tns_,int kappa_=1):
  infinitesimal_generator(fspace_), V_crv_tns(V_crv_tns_),
  xu(new double[nvar]), lamvec(new double[perm_len]), vx_vec(new double[ndof_full]),
  vxu_wkspc(vxu_workspace(nvar,fspace_.comp_ord_len())), kappa(kappa_) {}

  rspace_infinitesimal_generator(rspace_infinitesimal_generator &rinfgen_):
  rspace_infinitesimal_generator(rinfgen_.fspace,rinfgen_.V_crv_tns,rinfgen_.kappa) {ncol_use = rinfgen_.ncol_use;}

  ~rspace_infinitesimal_generator() {delete [] xu; delete [] lamvec; delete [] vx_vec;}

  int kappa = 1;
  double  ** Kmat;

  void init_dudx_eval(int icrv_) {Kmat = V_crv_tns[icrv_] + (ncol_use-kappa);}

  void JacF_eval(double x_, double *u_, double **dls_out_) {} // do later
  void dnp1xu_eval(double x_, double *u_, double *dnp1xu_) {} // do later

  static void init_extended_observations(ode_curve_observations &obs_out_,ode_curve_observations &obs_in_)
  {
    int ncrv = obs_out_.ncrv = obs_in_.ncrv,
        ndep = (obs_out_.ndep!=obs_in_.ndep)?(obs_out_.ndep=obs_in_.ndep):(obs_out_.ndep),
        nobs = (obs_out_.nobs!=obs_in_.nobs)?(obs_out_.nobs=obs_in_.nobs):(obs_out_.nobs),
        ndim_in = 1 + ndep*(obs_in_.eor + 1),
        ndim_out =  1 + ndep*(obs_out_.eor + 1);

    double  * const pts_chunk_in = obs_in_.pts_in,
            * const pts_chunk_out = obs_out_.pts_in = new double[ndim_out*nobs];

    if (obs_out_.npts_per_crv == NULL) obs_out_.npts_per_crv = new int[ncrv];
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
        basis_.v_eval(pts_i,v_,theta_vec_);
        for (size_t i = 0, idkxu = i+ioffset; i < ndep; i++, idkxu++) pts_i[idkxu] = v_[idkxu-ndep];
      }
    }
  }

  inline void init_svd_default(LD_matrix_svd_result &svd_) {kappa = svd_.kappa_def(); init_dudx_eval(0);}

  protected:

    int ncol_use = ndof_full;
    double  *** const V_crv_tns,
            * const xu,
            * const lamvec,
            * const vx_vec;

    vxu_workspace vxu_wkspc;
};

struct restricted_rspace_infgen: public rspace_infinitesimal_generator
{
  restricted_rspace_infgen(function_space &fspace_,int ncrvs_,LD_matrix_svd_result &AY_L_svd_):
    rspace_infinitesimal_generator(fspace_,T3tensor<double>(ncrvs_,fspace_.ndof_full,fspace_.ndof_full)),
    AY_L_svd(AY_L_svd_), V_crv_tns_owner(true) {}
  restricted_rspace_infgen(restricted_rspace_infgen &rrinfgen_):
    rspace_infinitesimal_generator(rrinfgen_),
    AY_L_svd(rrinfgen_.AY_L_svd), V_crv_tns_owner(false) {rho_L = rrinfgen_.rho_L; ncols_PL = rrinfgen_.ncols_PL;}
  ~restricted_rspace_infgen() {if (V_crv_tns_owner) free_T3tensor<double>(V_crv_tns);}

  LD_matrix_svd_result &AY_L_svd;

  const bool V_crv_tns_owner;
  const int ncrvs = AY_L_svd.ncrvs,
            ncols = AY_L_svd.ncols;
  int * const rank_vec = AY_L_svd.rank_vec;
  double  ** const Smat = AY_L_svd.Smat,
          *** const VTtns = AY_L_svd.VTtns;

  void compute_restricted_curve_svds(LD_matrix &Amat_, LD_matrix_svd_result &Lsvd_, int nrows_)
  {
    rho_L = Lsvd_.min_rank(),
    ncols_PL = AY_L_svd.ncol_use = nvar*rho_L;

    double  **** const Attns = Amat_.Attns,
            *** const VTtns_L = Lsvd_.VTtns;
    #pragma omp parallel
    {
      LD_SVD_space svd_t(nrows_,ncols_PL);
      gsl_matrix * const U_gsl_t = svd_t.U_gsl;
      #pragma omp for
      for (size_t icrv = 0; icrv < ncrvs; icrv++)
      {
        fill_AY_mat_icrv(U_gsl_t,Attns[icrv][0],VTtns_L[icrv],nrows_);
        svd_t.decompose_loaded_U();
        rank_vec[icrv] = svd_t.unpack_rank_s_VT(Smat[icrv],VTtns[icrv]);
        fill_YV_mat_icrv(V_crv_tns[icrv],VTtns_L[icrv],VTtns[icrv]);
      }
    }
  }

  inline void init_svd_default() {kappa = AY_L_svd.kappa_def(); init_dudx_eval(0);}

  protected:

    int rho_L = perm_len - 1,
        &ncols_PL = ncol_use = nvar*rho_L;

    inline void fill_AY_mat_icrv(gsl_matrix *AY_t_,double **Amat_rows_,double **Ycols_L_,int nrows_)
    {
      for (size_t irow = 0; irow < nrows_; irow++)
        for (size_t ivar = 0, jcolA_start = 0, jcolAY = 0; ivar < nvar; ivar++, jcolA_start+=perm_len)
          for (size_t kcolY_L = 0; kcolY_L < rho_L; kcolY_L++, jcolAY++)
          {
            double AY_irowjcol = 0.0;
            for (size_t pL = 0, jcolA = jcolA_start; pL < perm_len; pL++, jcolA++)
              AY_irowjcol += Amat_rows_[irow][jcolA]*Ycols_L_[kcolY_L][pL];
            gsl_matrix_set(AY_t_,irow,jcolAY,AY_irowjcol);
          }
    }
    inline void fill_YV_mat_icrv(double **YV_T_rows_,double **Ycols_L_,double **VTrows_AY_)
    {
      for (size_t irow = 0; irow < ncols_PL; irow++)
        for (size_t ivar = 0, jcolYV_T = 0, jcolVT_AY_start = 0; ivar < nvar; ivar++, jcolVT_AY_start+=rho_L)
          for (size_t pL = 0; pL < perm_len; pL++, jcolYV_T++)
          {
            YV_T_rows_[irow][jcolYV_T] = 0.0;
            for (size_t kcolY_L = 0, jcolVT_AY = jcolVT_AY_start; kcolY_L < rho_L; kcolY_L++, jcolVT_AY++)
              YV_T_rows_[irow][jcolYV_T] += VTrows_AY_[irow][jcolVT_AY]*Ycols_L_[kcolY_L][pL];
          }
    }
};

struct matrix_Lie_detector
{
  matrix_Lie_detector() {}
  ~matrix_Lie_detector() {}

  static void compute_curve_svds(LD_matrix &mat_,LD_matrix_svd_result &mat_svd_,int nrows_)
  {
    const int ncols = mat_svd_.ncols,
              ncrvs = mat_.ncrvs_tot;
    int * const rank_vec_out = mat_svd_.rank_vec;
    double  ** const Smat_out = mat_svd_.Smat,
            *** const VTtns_out = mat_svd_.VTtns;
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
  }

  template <class INFGN, class BSIS> static void extend_ode_observations(ode_curve_observations &obs_out_, INFGN &infgen_, BSIS **bases_)
  {
    const int ncrv = obs_out_.ncrv,
              ndof_full = infgen_.ndof_full;
    int * const npts_per_crv = obs_out_.npts_per_crv;
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
  }

};

#endif
