#include "matrix_Lie_detector.hh"
// #include <gsl/gsl_linalg.h>

matrix_Lie_detector::matrix_Lie_detector(LD_matrix &mat_, function_space &fspace_): mat(mat_), fspace(fspace_)
{}
matrix_Lie_detector::~matrix_Lie_detector()
{}
LD_matrix_svd_result * matrix_Lie_detector::compute_curve_svds(int nrows_)
{
  const int ncrvs = mat.ncrvs_tot,
            ndof = mat.ndof_full;
  LD_matrix_svd_result * svd_result = new LD_matrix_svd_result(ncrvs,ndof);
  double  ** const Smat_out = svd_result->Smat,
          *** const VTtns_out = svd_result->VTtns;
  #pragma omp parallel
  {
    LD_SVD_space svd_t(nrows_,ndof);
    #pragma omp for
    for (size_t i = 0; i < ncrvs; i++)
    {
      svd_t.load_and_decompose(mat.Amat_crv_i(i));
      svd_t.unpack_s_VT(Smat_out[i],VTtns_out[i]);
    }
  }
  return svd_result;
}
