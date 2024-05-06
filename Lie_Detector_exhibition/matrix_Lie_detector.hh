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
};

struct LD_matrix_svd_result
{
  LD_matrix_svd_result(int ncrvs_,int ndof_): ncrvs(ncrvs_), ndof(ndof_),
  Smat(Tmatrix<double>(ncrvs_,ndof_)), VTtns(T3tensor<double>(ncrvs_,ndof_,ndof_)) {}
  ~LD_matrix_svd_result() {free_Tmatrix<double>(Smat); free_T3tensor<double>(VTtns);}
  const int ncrvs,
            ndof;
  double  ** const Smat,
          *** const VTtns;
};

struct matrix_Lie_detector
{
  matrix_Lie_detector(LD_matrix &mat_, function_space &fspace_);
  ~matrix_Lie_detector();

  LD_matrix &mat;
  function_space &fspace;
  const int ndof_full = fspace.ndof_full;

  LD_matrix_svd_result * compute_curve_svds(int nrows_);
};

#endif
