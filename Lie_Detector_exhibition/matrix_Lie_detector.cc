#include "matrix_Lie_detector.hh"

void LD_matrix_svd_result::write_svd_results(const char name_[])
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  int hlen = 2,
      header[] = {hlen,ncrvs,ndof};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(rank_vec, sizeof(int), ncrvs, file);
  fwrite(Smat[0], sizeof(double), ncrvs*ndof, file);
  fwrite(VTtns[0][0], sizeof(double), ncrvs*ndof*ndof, file);
  fclose(file);
  printf("(LD_matrix_svd_result::write_svd_results) wrote %s\n",name_);
}
