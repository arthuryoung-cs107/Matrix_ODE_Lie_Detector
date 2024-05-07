#ifndef MAT_LD_HH
#define MAT_LD_HH

#include "LD_framework.hh"

// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_vector.h>
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
  LD_matrix_svd_result(int ncrvs_,int ndof_): ncrvs(ncrvs_), ndof(ndof_),
  rank_vec(new int[ncrvs_]), Smat(Tmatrix<double>(ncrvs_,ndof_)), VTtns(T3tensor<double>(ncrvs_,ndof_,ndof_)) {}
  ~LD_matrix_svd_result() {delete [] rank_vec; free_Tmatrix<double>(Smat); free_T3tensor<double>(VTtns);}

  const int ncrvs,
            ndof;
  int * const rank_vec;
  double  ** const Smat,
          *** const VTtns;

  void write_svd_results(const char name_[]);
  void read_svd_results(const char name_[]);

  inline void print_details(const char name_[]="Amat_svd")
  {
    LD_linalg::print_xT("(LD_matrix_svd_result::print_details) rank_vec",rank_vec,ncrvs);
    printf("min_nulldim = %d, max_nulldim = %d \n",
      min_nulldim(),
      max_nulldim());
  }

  inline int nulldim_i(int i_) {return ndof-rank_vec[i_];}
  inline int min_rank() {return LD_linalg::min_val<int>(rank_vec,ncrvs);}
  inline int max_rank() {return LD_linalg::max_val<int>(rank_vec,ncrvs);}
  inline int min_nulldim() {return ndof - max_rank();}
  inline int max_nulldim() {return ndof - min_rank();}

  inline int kappa_def() {int min_nulldim_val; return (min_nulldim_val=min_nulldim())?(min_nulldim_val):(1);}
  inline double ** Kmat_crvi(int icrv_,int kappa_=0) {return VTtns[icrv_]+(ndof-(kappa_)?(kappa_):(kappa_def()));}
};

struct infinitesimal_generator: public ode_system
{
  infinitesimal_generator(function_space &fspace_): ode_system(1,fspace_.ndep), fspace(fspace_) {}
  ~infinitesimal_generator() {}

  function_space &fspace;

  const int nvar = fspace.nvar,
            perm_len = fspace.perm_len;
};

class rspace_infinitesimal_generator: public infinitesimal_generator
{
  public:
    rspace_infinitesimal_generator(function_space &fspace_,int kappa_,double **Kmat_):
    infinitesimal_generator(fspace_),
    kappa(kappa_), Kmat(Kmat_),
    xu(new double[nvar]), lamvec(new double[perm_len]),
    vx_vec(new double[kappa_]),
    vxu_wkspc(vxu_workspace(nvar,fspace_.comp_ord_len()))
    {}
    ~rspace_infinitesimal_generator()
    {
      delete [] xu; delete [] lamvec;
      delete [] vx_vec;
    }

    const int kappa;
    double  ** const Kmat,
            * const xu,
            * const lamvec,
            * const vx_vec;

    vxu_workspace vxu_wkspc;

    void dudx_eval(double x_, double *u_, double *dudx_)
    {
      xu[0] = x_;
      for (size_t i = 0; i < ndep; i++) xu[i+1] = u_[i];
      fspace.lamvec_eval(xu,lamvec,vxu_wkspc);
      double Vx2 = 0.0;
      for (size_t i_k = 0; i_k < kappa; i_k++)
      {
        vx_vec[i_k] = 0.0;
        for (size_t i = 0; i < perm_len; i++) vx_vec[i_k] += Kmat[i_k][i]*lamvec[i];
        Vx2 += vx_vec[i_k]*vx_vec[i_k];
      }
    }
    void JacF_eval(double x_, double *u_, double **dls_out_) {}

};

struct matrix_Lie_detector
{
  matrix_Lie_detector() {}
  ~matrix_Lie_detector() {}

  static void compute_curve_svds(LD_matrix &mat_,LD_matrix_svd_result &mat_svd_,int nrows_)
  {
    const int ndof = mat_.ndof_full,
              ncrvs = mat_.ncrvs_tot;
    int * const rank_vec_out = mat_svd_.rank_vec;
    double  ** const Smat_out = mat_svd_.Smat,
            *** const VTtns_out = mat_svd_.VTtns;
    #pragma omp parallel
    {
      LD_SVD_space svd_t(nrows_,ndof);
      #pragma omp for
      for (size_t i = 0; i < ncrvs; i++)
      {
        svd_t.load_and_decompose(mat_.Amat_crv_i(i));
        rank_vec_out[i] = svd_t.unpack_rank_s_VT(Smat_out[i],VTtns_out[i]);
      }
    }
  }

};

#endif
