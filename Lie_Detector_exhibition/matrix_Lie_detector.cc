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

void LD_matrix_svd_result::read_svd_results(const char name_[])
{
  const int hlen_check = 2,
            Slen = ncrvs*ndof,
            VTlen = ncrvs*ndof*ndof;
  int hlen_in,
      header_in[hlen_check],
      &ncrvs_in = header_in[0],
      &ndof_in = header_in[1];
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if (hlen_in==hlen_check)
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_in,file_in);
    if ((ncrvs_in==ncrvs)&&(ndof_in==ndof))
    {
      LD_io::fread_SAFE(rank_vec,sizeof(int),ncrvs,file_in);
      LD_io::fread_SAFE(Smat[0],sizeof(double),Slen,file_in);
      LD_io::fread_SAFE(VTtns[0][0],sizeof(double),VTlen,file_in);
      LD_io::fclose_SAFE(file_in);
      printf("(LD_matrix_svd_result::read_svd_results) read %s\n",name_);
    }
    else
    {
      printf("(LD_matrix_svd_result::read_svd_results) ERROR: ncrvs_in, ndof_in incorrect (%d vs. %d, %d vs. %d)\n",ncrvs_in,ncrvs,ndof_in,ndof);
      LD_io::fclose_SAFE(file_in);
      return;
    }
  }
  else
  {
    printf("(LD_matrix_svd_result::read_svd_results) ERROR: hlen_in != hlen_check (%d vs. %d)\n", hlen_in, hlen_check);
    LD_io::fclose_SAFE(file_in);
    return;
  }
}
