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
    int status_decomp = gsl_linalg_SV_decomp(U_gsl,V_gsl,s_gsl,w_gsl);
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
};

struct LD_matrix_svd_result
{
  LD_matrix_svd_result(int ncrvs_,int ncols_): ncrvs(ncrvs_), ncols(ncols_),
  rank_vec(new int[ncrvs_]), Smat(Tmatrix<double>(ncrvs_,ncols_)), VTtns(T3tensor<double>(ncrvs_,ncols_,ncols_)) {}
  ~LD_matrix_svd_result() {delete [] rank_vec; free_Tmatrix<double>(Smat); free_T3tensor<double>(VTtns);}

  const int ncrvs,
            ncols;
  int * const rank_vec;
  double  ** const Smat,
          *** const VTtns;

  void write_svd_results(const char name_[]);
  void read_svd_results(const char name_[]);

  inline void print_details(const char name_[]="Amat_svd")
  {
    LD_linalg::print_xT("(LD_matrix_svd_result::print_details) rank_vec",rank_vec,ncrvs);
    printf("min_nulldim = %d, max_nulldim = %d \n", min_nulldim(), max_nulldim());
  }

  inline int nulldim_i(int i_) {return ncols-rank_vec[i_];}
  inline int min_rank() {return LD_linalg::min_val<int>(rank_vec,ncrvs);}
  inline int max_rank() {return LD_linalg::max_val<int>(rank_vec,ncrvs);}
  inline int min_nulldim() {return ncols - max_rank();}
  inline int max_nulldim() {return ncols - min_rank();}

  inline int kappa_def() {int min_nulldim_val; return (min_nulldim_val=min_nulldim())?(min_nulldim_val):(1);}
  inline double ** Kmat_crvi(int icrv_,int kappa_=0) {return VTtns[icrv_]+(ncols-(kappa_)?(kappa_):(kappa_def()));}
};

struct infinitesimal_generator: public ode_system
{
  infinitesimal_generator(function_space &fspace_): ode_system(1,fspace_.ndep), fspace(fspace_) {}
  ~infinitesimal_generator() {}

  function_space &fspace;
  const int perm_len = fspace.perm_len,
            ndof = fspace.ndof_full;
};

class rspace_infinitesimal_generator: public infinitesimal_generator
{
  double  *** const Ktns,
          * const xu,
          * const lamvec,
          * const vx_vec;
  vxu_workspace vxu_wkspc;

  public:

    rspace_infinitesimal_generator(function_space &fspace_,double ***Ktns_,int kappa_=1):
    infinitesimal_generator(fspace_), Ktns(Ktns_),
    xu(new double[nvar]), lamvec(new double[perm_len]), vx_vec(new double[ndof]),
    vxu_wkspc(vxu_workspace(nvar,fspace_.comp_ord_len())), kappa(kappa_) {}

    rspace_infinitesimal_generator(rspace_infinitesimal_generator &rinfgen_):
    rspace_infinitesimal_generator(rinfgen_.fspace,rinfgen_.Ktns,rinfgen_.kappa) {}

    ~rspace_infinitesimal_generator() {delete [] xu; delete [] lamvec; delete [] vx_vec;}

    int kappa = 1;
    double  ** Kmat;

    void init_dudx_eval(int icrv_) {Kmat = Ktns[icrv_] + (ndof-kappa);}

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
        LD_linalg::fill_vec<double>(theta_vec_,ndof,0.0);
        double Vx2 = 0.0;
        for (size_t i_k = 0; i_k < kappa; i_k++)
        {
          double &vx_ik = vx_vec[i_k] = 0.0;
          for (size_t i = 0; i < perm_len; i++) vx_ik += Kmat[i_k][i]*lamvec[i];
          for (size_t i = 0; i < ndof; i++) theta_vec_[i] += Kmat[i_k][i]*vx_ik;
          Vx2 += vx_ik*vx_ik;
        }
        for (size_t i = 0; i < ndof; i++) theta_vec_[i] /= Vx2;
        for (size_t k = 2, ioffset = ndep*k + 1; k <= eor_full; k++, ioffset+=ndep)
        {
          basis_.v_eval(pts_i,v_,theta_vec_);
          for (size_t i = 0, idkxu = i+ioffset; i < ndep; i++, idkxu++) pts_i[idkxu] = v_[idkxu-ndep];
        }
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
              ndof = infgen_.ndof;
    int * const npts_per_crv = obs_out_.npts_per_crv;
    #pragma omp parallel
    {
      INFGN infgen_t(infgen_);
      BSIS &basis_t = *(bases_[LD_threads::thread_id()]);
      double  theta_vec_t[ndof],
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
