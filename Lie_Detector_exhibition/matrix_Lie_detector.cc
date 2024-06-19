#include "matrix_Lie_detector.hh"

LD_svd_file::LD_svd_file(const char name_[])
{
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if ((hlen_in == 2)||(hlen_in == 3))
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_in,file_in);
    if (hlen_in == 2) ncol_use_in = ncols_in;
    rank_vec_in = new int[ncrvs_in];
    Smat_in = Tmatrix<double>(ncrvs_in,ncol_use_in);
    VTtns_in = T3tensor<double>(ncrvs_in,ncol_use_in,ncol_use_in);
    LD_io::fread_SAFE(rank_vec_in,sizeof(int),ncrvs_in,file_in);
    LD_io::fread_SAFE(Smat_in[0],sizeof(double),ncrvs_in*ncols_in,file_in);
    LD_io::fread_SAFE(VTtns_in[0][0],sizeof(double),ncrvs_in*ncols_in*ncols_in,file_in);
    printf("(LD_svd_file::LD_svd_file) read %s \n", name_);
  }
  else printf("(LD_svd_file::LD_svd_file) ERROR: hlen_in unsupported (%d vs. %d or %d)\n", hlen_in,2,3);
  LD_io::fclose_SAFE(file_in);
}

LD_matrix_svd_result::LD_matrix_svd_result(int ncrvs_,int ncols_,int ncol_use_):
  ncrvs(ncrvs_), ncols(ncols_), ncol_use((ncol_use_)?(ncol_use_):(ncols_)),
  rank_vec(new int[ncrvs_]), Smat(Tmatrix<double>(ncrvs_,ncols_)), VTtns(T3tensor<double>(ncrvs_,ncols_,ncols_)) {}
LD_matrix_svd_result::LD_matrix_svd_result(LD_svd_file svdfile_):
  LD_matrix_svd_result(svdfile_.ncrvs_in,svdfile_.ncols_in,svdfile_.ncol_use_in)
{
  memcpy(rank_vec,svdfile_.rank_vec_in,ncrvs*sizeof(int));
  const int col_chunk_len = ncol_use*sizeof(double);
  if (ncol_use==ncols)
  {
    const int mat_chunk_len = ncrvs*col_chunk_len;
    memcpy(Smat[0],svdfile_.Smat_in[0],mat_chunk_len);
    memcpy(VTtns[0][0],svdfile_.VTtns_in[0][0],ncols*mat_chunk_len);
  }
  else
    for (size_t icrv = 0; icrv < ncrvs; icrv++)
    {
      memcpy(Smat[icrv],svdfile_.Smat_in[icrv],col_chunk_len);
      for (size_t icol = 0; icol < ncol_use; icol++) memcpy(VTtns[icrv][icol],svdfile_.VTtns_in[icrv][icol],col_chunk_len);
    }
}

void LD_matrix_svd_result::write_svd_results(const char name_[])
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  if (ncol_use == ncols)
  {
    int hlen = 2,
        header[] = {hlen,ncrvs,ncols};
    fwrite(header, sizeof(int), hlen+1, file);
    fwrite(rank_vec, sizeof(int), ncrvs, file);
    fwrite(Smat[0], sizeof(double), ncrvs*ncols, file);
    fwrite(VTtns[0][0], sizeof(double), ncrvs*ncols*ncols, file);
  }
  else
  {
    int hlen = 3,
        header[] = {hlen,ncrvs,ncols,ncol_use};
    fwrite(header, sizeof(int), hlen+1, file);
    fwrite(rank_vec, sizeof(int), ncrvs, file);
    for (size_t icrv = 0; icrv < ncrvs; icrv++) fwrite(Smat[icrv], sizeof(double), ncol_use,file);
    for (size_t icrv = 0; icrv < ncrvs; icrv++)
      for (size_t icol = 0; icol < ncol_use; icol++)
        fwrite(VTtns[icrv][icol], sizeof(double), ncol_use, file);
  }
  fclose(file);
  printf("(LD_matrix_svd_result::write_svd_results) wrote %s\n",name_);
}

void LD_matrix_svd_result::read_svd_results(const char name_[])
{
  const int hlen_max = 3,
            Slen = ncrvs*ncols,
            VTlen = ncrvs*ncols*ncols;
  int hlen_in,
      header_in[hlen_max],
      &ncrvs_in = header_in[0],
      &ncols_in = header_in[1],
      &ncol_use_in = header_in[2];
  FILE * file_in = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen_in,sizeof(int),1,file_in);
  if (hlen_in<=hlen_max)
  {
    LD_io::fread_SAFE(header_in,sizeof(int),hlen_in,file_in);
    if (hlen_in==2)
    {
      if ((ncrvs_in==ncrvs)&&(ncols_in==ncols))
      {
        ncol_use = ncols;
        LD_io::fread_SAFE(rank_vec,sizeof(int),ncrvs,file_in);
        LD_io::fread_SAFE(Smat[0],sizeof(double),Slen,file_in);
        LD_io::fread_SAFE(VTtns[0][0],sizeof(double),VTlen,file_in);
        LD_io::fclose_SAFE(file_in);
        printf("(LD_matrix_svd_result::read_svd_results) read %s\n",name_);
      }
      else
      {
        printf("(LD_matrix_svd_result::read_svd_results) ERROR: ncrvs_in, ncols_in incorrect (%d vs. %d, %d vs. %d)\n",ncrvs_in,ncrvs,ncols_in,ncols);
        LD_io::fclose_SAFE(file_in);
        return;
      }
    }
    else if (hlen_in==3)
    {
      if ((ncrvs_in==ncrvs)&&(ncol_use_in<=ncols))
      {
        ncol_use = ncol_use_in;
        LD_io::fread_SAFE(rank_vec,sizeof(int),ncrvs,file_in);
        for (size_t icrv = 0; icrv < ncrvs; icrv++) LD_io::fread_SAFE(Smat[icrv],sizeof(double),ncol_use_in,file_in);
        for (size_t icrv = 0; icrv < ncrvs; icrv++)
          for (size_t icol = 0; icol < ncol_use_in; icol++)
            LD_io::fread_SAFE(VTtns[icrv][icol], sizeof(double), ncol_use_in, file_in);
        LD_io::fclose_SAFE(file_in);
        printf("(LD_matrix_svd_result::read_svd_results) read %s\n",name_);
      }
      else
      {
        printf("(LD_matrix_svd_result::read_svd_results) ERROR: ncrvs_in, ncol_use_in incorrect (%d vs. %d, %d vs. %d)\n",ncrvs_in,ncrvs,ncol_use_in,ncols);
        LD_io::fclose_SAFE(file_in);
        return;
      }
    }
    else
    {
      printf("(LD_matrix_svd_result::read_svd_results) ERROR: hlen_in != 2 or 3 (%d)\n", hlen_in);
      LD_io::fclose_SAFE(file_in);
      return;
    }
  }
  else
  {
    printf("(LD_matrix_svd_result::read_svd_results) ERROR: hlen_in > hlen_max (%d vs. %d)\n", hlen_in, hlen_max);
    LD_io::fclose_SAFE(file_in);
    return;
  }
}
