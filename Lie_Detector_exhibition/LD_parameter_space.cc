#include "LD_parameter_space.hh"

// #include <cstring>
// #include <cstdio>
#include "LD_io.hh"

/*
  constructors
*/

LD_vector_space::LD_vector_space(int vlen_): data_owner(true),
  Vmat_data(Tmatrix<double>(vlen_,vlen_)),
  vlen_full(vlen_), nvec_use(vlen_), vlen_use(vlen_),
  inds_V(new int[vlen_]), Vmat(new double*[vlen_])
  // ,vdata(LD_vspace_data(vlen_full,nvec_use,vlen_use,inds_V,Vmat_data))
  {for (size_t i = 0; i < vlen_; i++) Vmat[i] = Vmat_data[inds_V[i] = i];}
LD_vector_space::LD_vector_space(int vlen_,int *inds_V_,double **Vmat_data_,double **Vmat_): data_owner(false),
  Vmat_data(Vmat_data_),
  vlen_full(vlen_), nvec_use(vlen_), vlen_use(vlen_),
  inds_V(inds_V_), Vmat(Vmat_)
  // ,vdata(LD_vspace_data(vlen_full,nvec_use,vlen_use,inds_V,Vmat_data))
  {}
LD_vector_space::LD_vector_space(LD_vector_space &Vspc_,bool deepcopy_): data_owner(deepcopy_),
  Vmat_data((deepcopy_)?(Tmatrix<double>(Vspc_.vlen_full,Vspc_.vlen_full)):(Vspc_.get_Vmat_data())),
  vlen_full(Vspc_.vlen_full), nvec_use(Vspc_.nvec_use), vlen_use(Vspc_.vlen_use),
  inds_V((deepcopy_)?(new int[vlen_full]):(Vspc_.inds_V)),
  Vmat((deepcopy_)?(new double*[vlen_full]):(Vspc_.Vmat))
  // ,vdata(LD_vspace_data(vlen_full,nvec_use,vlen_use,inds_V,Vmat_data))
  {
    if (deepcopy_)
    {
      double ** const Vdata_in = Vspc_.get_Vmat_data();
      for (size_t i = 0; i < vlen_full; i++)
        for (size_t j = 0; j < vlen_full; j++) Vmat_data[i][j] = Vdata_in[i][j];
      for (size_t i = 0; i < vlen_full; i++) Vmat[i] = Vmat_data[inds_V[i] = Vspc_.inds_V[i]];
    }
  }
LD_vector_space::~LD_vector_space()
  {if (data_owner) {delete [] inds_V; free_Tmatrix<double>(Vmat_data); delete [] Vmat;}}

LD_vector_bundle::LD_vector_bundle(int nspc_, int vlen_): data_owner(true),
  Vtns_data(T3tensor<double>(nspc_,vlen_,vlen_)), Vrows(new double*[nspc_*vlen_]),
  rec_ptr(new LD_vspace_record(nspc_,vlen_,vlen_,Vtns_data)),
  nspc(nspc_), vlen_full(vlen_),
  Vtns(new double**[nspc_]), Vspaces(new LD_vector_space*[nspc_])
  {
    for (size_t i = 0; i < nspc_; i++)
    {
      Vspaces[i] = new LD_vector_space(nV_spcvec[i]=vlen_,iV_spcmat[i],Vtns_data[i],Vtns[i] = Vrows + (i*vlen_));
      for (size_t j = 0; j < vlen_; j++) Vtns[i][j] = Vtns_data[i][iV_spcmat[i][j] = (int)(j)];
    }
  }
LD_vector_bundle::LD_vector_bundle(LD_vector_bundle &bndle_): data_owner(false),
  Vtns_data(bndle_.Vtns_data), Vrows(bndle_.Vrows),
  rec_ptr(bndle_.rec_ptr),
  nspc(bndle_.nspc), vlen_full(bndle_.vlen_full),
  Vtns(bndle_.Vtns), Vspaces(bndle_.Vspaces) {}
LD_vector_bundle::~LD_vector_bundle()
  {
    if (data_owner)
    {
      for (size_t i = 0; i < nspc; i++) delete Vspaces[i];
      delete [] Vspaces; delete [] Vtns;
      delete rec_ptr;
      delete [] Vrows; free_T3tensor<double>(Vtns_data);
    }
  }


LD_Theta_space::LD_Theta_space(ode_solspc_meta &meta_, int ndof_):
  ode_solspc_element(meta_), vspc_owner(true), data_owner(true),
  perm_len(ndof_/nvar), spc_ptr(new LD_vector_space(ndof_)), bse_ptrs(new LD_vector_space*[nvar]),
  ndof_use(ndof_), pbse_nvar_use(new int[nvar])
  {for (size_t i = 0; i < nvar; i++) {pbse_nvar_use[i] = perm_len; bse_ptrs[i] = NULL;}}
LD_Theta_space::LD_Theta_space(ode_solspc_meta &meta_, int ndof_, LD_vector_space *spc_ptr_, LD_vector_space **bse_ptrs_, int *pbse_):
  ode_solspc_element(meta_), vspc_owner(false), data_owner(false),
  perm_len(ndof_/nvar), spc_ptr(spc_ptr_), bse_ptrs(bse_ptrs_),
  ndof_use(ndof_), pbse_nvar_use(pbse_) {}
LD_Theta_space::~LD_Theta_space()
  {
    if (data_owner) {delete [] pbse_nvar_use; delete [] bse_ptrs;}
    if (vspc_owner) delete spc_ptr;
  }

LD_Theta_bundle::LD_Theta_bundle(ode_solspc_meta &meta_,int nspc_,int ndof_): ode_solspc_element(meta_),
  Vbndl_owner(true), data_owner(true),
  Vbndle_ptr(new LD_vector_bundle(nspc_,ndof_)), Yspaces_nvar(new LD_vector_space*[nvar*nspc_]), perm_len(ndof_/nvar),
  pbse_nvar_spcmat(Tmatrix<int>(nspc_,nvar)),
  Tspaces(new LD_Theta_space*[nspc_])
  {
    for (size_t i = 0; i < nspc_; i++)
    {
      Tspaces[i] = new LD_Theta_space(meta_,ndof_,Vspaces[i],Yspaces_nvar+(nvar*i),pbse_nvar_spcmat[i]);
      for (size_t j = 0; j < nvar; j++) {*(Yspaces_nvar+(nvar*i) + j) = NULL; pbse_nvar_spcmat[i][j] = perm_len;}
    }
  }
LD_Theta_bundle::LD_Theta_bundle(LD_Theta_bundle &Tbndle_): ode_solspc_element(Tbndle_.meta),
  Vbndl_owner(false), data_owner(false),
  Vbndle_ptr(Tbndle_.Vbndle_ptr), Yspaces_nvar(Tbndle_.Yspaces_nvar), perm_len(perm_len),
  pbse_nvar_spcmat(Tbndle_.pbse_nvar_spcmat),
  Tspaces(Tbndle_.Tspaces) {}
LD_Theta_bundle::~LD_Theta_bundle()
  {
    if (data_owner)
    {
      for (size_t i = 0; i < nspc; i++) delete Tspaces[i];
      delete [] Tspaces;
      free_Tmatrix<int>(pbse_nvar_spcmat);
      delete [] Yspaces_nvar;
    }
    if (Vbndl_owner) delete Vbndle_ptr;
  }

LD_vspace_evaluator::LD_vspace_evaluator(ode_solspc_subset &Sset_,int ncon_,int nvec_,int nsol_,double tol_): Sset(Sset_),
  ncon(ncon_), nvec_max(nvec_), nsol_max(nsol_),
  sat_flags_mat(Tmatrix<bool>(nsol_,nvec_)),
  tol(tol_), nvec_evl(nvec_), nsol_evl(nsol_) {}

/*
  definitions
*/
void LD_vspace_record::init_combined_records(LD_vspace_record &rec1_,LD_vspace_record &rec2_)
{
  bool fV_cmpmat[2][nvec];
  const int min_nspc = LD_linalg::min_T<int>(rec1_.nspc,rec2_.nspc);
  for (size_t ispc = 0; ispc < min_nspc; ispc++)
  {
    for (size_t iV = 0; iV < nvec; iV++) fV_cmpmat[0][iV] = false;
    for (size_t iV = 0; iV < rec1_.nV_spcvec[ispc]; iV++) fV_cmpmat[0][rec1_.iV_spcmat[ispc][iV]] = true;
    for (size_t iV = 0; iV < nvec; iV++) fV_cmpmat[1][iV] = false;
    for (size_t iV = 0; iV < rec2_.nV_spcvec[ispc]; iV++) fV_cmpmat[1][rec2_.iV_spcmat[ispc][iV]] = true;
    int &nV_i = nV_spcvec[ispc] = 0;
    for (size_t iV = 0; iV < nvec; iV++) if ((fV_cmpmat[0][iV])&&(fV_cmpmat[1][iV])) iV_spcmat[ispc][nV_i++] = iV;
  }
}
void LD_Theta_space::post_multiply_evry_basis(double **AYmat_, double **Amat_, int mrows_)
{
  for (size_t irow = 0; irow < mrows_; irow++)
    for (size_t ivar = 0, jcolAY = 0, jcolA = 0; ivar < nvar; jcolAY+=pbse_nvar_use[ivar], jcolA+=perm_len, ivar++)
    {
      if (bse_ptrs[ivar]==NULL)
        for (size_t jbse = 0; jbse < perm_len; jbse++) AYmat_[irow][jcolAY+jbse] = Amat_[irow][jcolA+jbse];
      else bse_ptrs[ivar]->post_multiply_rowvec(AYmat_[irow]+jcolAY,Amat_[irow]+jcolA);
    }
}
void LD_Theta_space::post_multiply_ivar_basis(double **AYmat_, double **Amat_, int mrows_, int ivar_)
{
  for (size_t irow = 0; irow < mrows_; irow++)
  {
    int jcolA = 0;
    for (size_t ivar = 0; ivar < ivar_; ivar++, jcolA+=perm_len)
      for (size_t jbse = 0; jbse < perm_len; jbse++) AYmat_[irow][jcolA+jbse] = Amat_[irow][jcolA+jbse];
    bse_ptrs[ivar_]->post_multiply_rowvec(AYmat_[irow]+jcolA,Amat_[irow]+jcolA);
    int jcolAY = jcolA+pbse_nvar_use[ivar_]; jcolA += perm_len;
    for (size_t ivar = ivar_+1; ivar < nvar; jcolAY+=perm_len, jcolA+=perm_len, ivar++)
      for (size_t jbse = 0; jbse < perm_len; jbse++) AYmat_[irow][jcolAY+jbse] = Amat_[irow][jcolA+jbse];
  }
}
int LD_Theta_space::init_Vspce_premult(double **Wmat_)
{
  const int nV = verify_Yspaces();
  for (size_t irow = 0; irow < nV; irow++)
  {
    double * const vi = Tmat[irow] = spc.Vrowi_data(spc.inds_V[irow] = irow);
    for (size_t ivar = 0, jcolV = 0, jcolW = 0; ivar < nvar; jcolV+=perm_len, jcolW+=pbse_nvar_use[ivar], ivar++)
    {
      if (bse_ptrs[ivar]==NULL) for (size_t jbse = 0; jbse < perm_len; jbse++) vi[jcolV+jbse] = Wmat_[irow][jcolW+jbse];
      else bse_ptrs[ivar]->pre_multiply_colvec(vi+jcolV,Wmat_[irow]+jcolW);
    }
  }
  return spc.nvec_use = nV;
}

void LD_svd_bundle::compute_Acode_curve_svds(LD_encoding_bundle &Abndle_, bool verbose_)
{
  const int mrows_max = Abndle_.max_nrows();

  LD_encoded_matrix ** const Acodes = Abndle_.Acodes;

  int mrows_acc = 0,
      ncols_acc = 0;
  double t0 = LD_threads::tic();

  printf("(LD_svd_bundle::compute_Acode_curve_svds) mrows_max = %d, vlen_full = %d, nspc = %d\n", mrows_max, vlen_full, nspc);
  #pragma omp parallel reduction(+:mrows_acc,ncols_acc)
  {
    LD_svd svd_t(mrows_max,vlen_full);
    int mrows_Ai, ncols_Ai;
    #pragma omp for
    for (size_t ispc = 0; ispc < nspc; ispc++)
    {
      mrows_acc += (mrows_Ai = Acodes[ispc]->verify_nrow());
      ncols_acc += (ncols_Ai = nV_spcvec[ispc] = Vspaces[ispc]->vlen_use = Acodes[ispc]->ncol);
      svd_t.load_and_decompose_U(Acodes[ispc]->Amat,mrows_Ai,ncols_Ai);
      rank_vec[ispc] = svd_t.unpack_rank_svec_VTmat(Smat[ispc],VTtns[ispc]);
    }
  }
  double work_time = LD_threads::toc(t0);
  if (verbose_)
    printf("(LD_svd_bundle::compute_Acode_curve_svds) computed %d svds (%.1f x %.1f, on avg.) in %.4f seconds (%d threads)\n",
      nspc,((double)mrows_acc)/((double)nspc),((double)ncols_acc)/((double)nspc),
      work_time, LD_threads::numthreads());

}

void LD_svd_bundle::compute_AYmat_curve_svds(LD_encoding_bundle &Abndle_,LD_Theta_space ** const Tspaces_,bool verbose_)
{
  const int mrows_max = Abndle_.max_nrows();
  int mrows_acc = 0,
      ncols_acc = 0;
  double  *** const Amats = Abndle_.Amats,
          t0 = LD_threads::tic();

  #pragma omp parallel reduction(+:mrows_acc,ncols_acc)
  {
    LD_svd svd_t(mrows_max,vlen_full);
    int mrows_Ai, ncols_Ai;
    double ** const Umat_t = svd_t.Umat;
    #pragma omp for
    for (size_t ispc = 0; ispc < nspc; ispc++)
    {
      mrows_acc += (mrows_Ai = Abndle_.nrows_mat_i(ispc));
      ncols_acc += (ncols_Ai = nV_spcvec[ispc] = Vspaces[ispc]->vlen_use = Tspaces_[ispc]->ndof_use);
      Tspaces_[ispc]->post_multiply_bases(Umat_t,Amats[ispc],mrows_Ai);
      svd_t.decompose_U(mrows_Ai,ncols_Ai);
      rank_vec[ispc] = svd_t.unpack_rank_svec_VTmat(Smat[ispc],VTtns[ispc]);
    }
  }
  double work_time = LD_threads::toc(t0);
  if (verbose_)
    printf("(LD_svd_bundle::compute_AYmat_curve_svds) computed %d svds (%.1f x %.1f, on avg.) in %.4f seconds (%d threads)\n",
      nspc,((double)mrows_acc)/((double)nspc),((double)ncols_acc)/((double)nspc),
      work_time, LD_threads::numthreads());
}

/*
  diagnostics
*/

void LD_vspace_record::print_subspace_details(const char name_[], bool longwinded_, bool print_comp_)
{
  printf("\n(LD_vspace_record::print_subspace_details) %s ", name_);
  if (longwinded_)
  {
    printf("isat_vec\n", name_);
    for (size_t ispc = 0; ispc < nspc; ispc++)
    {
      printf("spc %d (nsat = %d): ", ispc, nV_spcvec[ispc]);
      for (size_t isat = 0; isat < nV_spcvec[ispc]; isat++) printf("%d ", iV_spcmat[ispc][isat]);
      if (print_comp_)
      {
        printf("|| ");
        for (size_t insat = nV_spcvec[ispc]; insat < nvec; insat++) printf("%d ", iV_spcmat[ispc][insat]);
      }
      printf("\n");
    }
  }
  LD_linalg::print_xT("nsat_vec",nV_spcvec,nspc);
  double nV_stats[3]; LD_linalg::comp_Tvec_maM_stats<int>(nV_stats,nV_spcvec,nspc);
  printf("  nvec_use = %d, vlen_use = %d, min nsat = %d, avg. nsat = %.2f, max nsat = %d \n",
            nvec, vlen, (int)(nV_stats[0]), nV_stats[1], (int)(nV_stats[2]));
  printf("\n");
}

void LD_vspace_record::print_subspace_details(LD_vspace_record &rec1_,const char name_[], bool longwinded_, bool print_comp_)
{
  printf("\n");
  if (longwinded_)
  {
    printf("(LD_vspace_record::print_subspace_details) V1: %s isat_vec\n", name_);
    for (size_t ispc = 0; ispc < nspc; ispc++)
    {
      printf("spc %d (nsat = %d): ", ispc, rec1_.nV_spcvec[ispc]);
      for (size_t isat = 0; isat < rec1_.nV_spcvec[ispc]; isat++) printf("%d ", rec1_.iV_spcmat[ispc][isat]);
      if (print_comp_)
      {
        printf("|| ");
        for (size_t insat = rec1_.nV_spcvec[ispc]; insat < nV_spcvec[ispc]; insat++) printf("%d ", rec1_.iV_spcmat[ispc][insat]);
      }
      printf("(of %d)\n", nV_spcvec[ispc]);
    }
  }
  int nvec0_acc = 0,
      nvec1_acc = 0;
  for (size_t i = 0; i < nspc; i++)
  {
    nvec0_acc += nV_spcvec[i];
    nvec1_acc += rec1_.nV_spcvec[i];
  }
  printf("(LD_vspace_record::print_subspace_details) V1: %s (V0 nspc = %d). Avg. nV0 = %.2f, nV1 = %.2f, (embedded in %d dimensions)\n nV0: --v\n  ",
  name_,nspc,((double)nvec0_acc)/((double)nspc),((double)nvec1_acc)/((double)nspc),vlen);
  for (size_t i = 0; i < nspc; i++) printf("%d ", nV_spcvec[i]); printf("\n  ");
  for (size_t i = 0; i < nspc; i++) printf("%d ", rec1_.nV_spcvec[i]); printf("<-- %s nV\n",name_);
  printf("\n");
}

void LD_vspace_record::compare_subspaces(LD_vspace_record &rec1_,const char name1_[],LD_vspace_record &rec2_,const char name2_[])
{
  printf("\n");
  int nvec0_acc = 0,
      nvec1_acc = 0,
      nvec2_acc = 0,
      diff_acc = 0,
      nV_diff[nspc];
  for (size_t i = 0; i < nspc; i++)
  {
    diff_acc += (nV_diff[i] = rec1_.nV_spcvec[i]-rec2_.nV_spcvec[i]);
    nvec0_acc += nV_spcvec[i];
    nvec1_acc += rec1_.nV_spcvec[i];
    nvec2_acc += rec2_.nV_spcvec[i];
  }
  printf("(LD_vspace_record::compare_subspaces) comparing V1: %s vs. V2: %s (V0 nspc = %d). Avg. nV0 = %.2f, nV1 = %.2f, nV2 = %.2f, nV1-nV2 = %.2f (embedded in %d dimensions)\n nV0: --v\n  ",
  name1_,name2_,nspc, ((double)nvec0_acc)/((double)nspc),
  ((double)nvec1_acc)/((double)nspc),((double)nvec2_acc)/((double)nspc),
  ((double)diff_acc)/((double)nspc), vlen);
  for (size_t i = 0; i < nspc; i++) printf("%d ", nV_spcvec[i]); printf("\n  ");
  for (size_t i = 0; i < nspc; i++) printf("%d ", rec1_.nV_spcvec[i]); printf("<-- %s nV\n  ",name1_);
  for (size_t i = 0; i < nspc; i++) printf("%d ", rec2_.nV_spcvec[i]); printf("<-- %s nV\n  ",name2_);
  for (size_t i = 0; i < nspc; i++) printf("%d ", nV_diff[i]); printf("\n %s nV - %s nV --^\n",name1_,name2_);
  printf("\n");
}
void LD_svd_bundle::print_details(const char name_[])
{
  printf("\n(LD_matrix_svd_result::print_details) %s ", name_);
  LD_linalg::print_xT("rank_vec",rank_vec,nspc);
  for (size_t i = 0; i < nspc; i++) printf("%d ", nV_spcvec[i]); printf("\n  ^-- (out of n columns)\n");
  int ncol_acc = 0, max_ncol = 0, min_ncol = vlen_full,
      rank_acc = 0, max_rank = 0, min_rank = vlen_full,
      nuld_acc = 0, max_nuld = 0, min_nuld = vlen_full,
      ncol_i, rank_i, nuld_i;
  for (size_t i = 0; i < nspc; i++)
  {
    ncol_acc += (ncol_i = nV_spcvec[i]);
    rank_acc += (rank_i = rank_vec[i]);
    nuld_acc += (nuld_i = ncol_i-rank_i);
    if (min_ncol>ncol_i) min_ncol=ncol_i;
    if (max_ncol<ncol_i) max_ncol=ncol_i;
    if (min_rank>rank_i) min_rank=rank_i;
    if (max_rank<rank_i) max_rank=rank_i;
    if (min_nuld>nuld_i) min_nuld=nuld_i;
    if (max_nuld<nuld_i) max_nuld=nuld_i;
    printf("%d ", nuld_i);
  }
  printf("\n  ^-- (nullspace dimension)\n  ");
  printf("min. ncol = %d, avg. ncol = %.2f, max. ncol = %d; ", min_ncol,((double)(ncol_acc))/((double)(nspc)),max_ncol);
  printf("min. rank = %d, avg. rank = %.2f, max. rank = %d; ", min_rank,((double)(rank_acc))/((double)(nspc)),max_rank);
  printf("min. nuld = %d, avg. nuld = %.2f, max. nuld = %d\n", min_nuld,((double)(nuld_acc))/((double)(nspc)),max_nuld);
}

/*
  io
*/

void LD_vspace_record::write_vspace_record(const char name_[],bool write_Vtns_data_)
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  int hlen = 5,
      iV_spcmat_len = comp_iV_spcmat_len(),
      header[] = {hlen,nspc,nvec,vlen,nspc+iV_spcmat_len,(write_Vtns_data_)?(iV_spcmat_len*vlen):(0)};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(nV_spcvec, sizeof(int), nspc, file);
  if (iV_spcmat_len != (vlen*nspc))
    for (size_t i = 0; i < nspc; i++) fwrite(iV_spcmat[i], sizeof(int), nV_spcvec[i], file);
  if (write_Vtns_data_)
    for (size_t i = 0; i < nspc; i++)
      for (size_t j = 0; j < nV_spcvec[i]; j++)
        fwrite(Vtns_data[i][iV_spcmat[i][j]], sizeof(double), vlen, file);
  fclose(file);
  printf("(LD_vspace_record::write_vspace_record) wrote %s\n",name_);
}

void LD_vector_bundle::write_LD_vector_bundle(const char name_[],bool write_Vtns_)
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  int hlen = 4,
      iV_spcmat_len = rec.comp_iV_spcmat_len(),
      header[] = {hlen,nspc,vlen_full,nspc+iV_spcmat_len,(write_Vtns_)?(comp_Vtns_len()):(0)};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(nV_spcvec, sizeof(int), nspc, file);
  if (iV_spcmat_len != (vlen_full*nspc))
    for (size_t i = 0; i < nspc; i++) fwrite(iV_spcmat[i], sizeof(int), nV_spcvec[i], file);
  if (write_Vtns_)
    for (size_t i = 0; i < nspc; i++)
      for (size_t j = 0; j < nV_spcvec[i]; j++)
        fwrite(Vtns[i][j], sizeof(double), Vspaces[i]->vlen_use, file);
  fclose(file);
  printf("(LD_vector_bundle::write_LD_vector_bundle) wrote %s\n",name_);
}

void LD_svd_bundle::write_LD_svd_bundle(const char name_[])
{
  FILE * file = LD_io::fopen_SAFE(name_,"wb");
  int hlen = 4,
      spcmat_len = rec.comp_iV_spcmat_len(),
      header[] = {hlen,nspc,vlen_full,(2*nspc)+spcmat_len,comp_Vtns_len() + spcmat_len};
  fwrite(header, sizeof(int), hlen+1, file);
  fwrite(nV_spcvec, sizeof(int), nspc, file);
  if (spcmat_len != (vlen_full*nspc))
    for (size_t i = 0; i < nspc; i++) fwrite(iV_spcmat[i], sizeof(int), nV_spcvec[i], file);
  fwrite(rank_vec, sizeof(int), nspc, file);
  for (size_t i = 0; i < nspc; i++)
    for (size_t j = 0; j < nV_spcvec[i]; j++) fwrite(Vtns[i][j], sizeof(double), Vspaces[i]->vlen_use, file);
  for (size_t i = 0; i < nspc; i++) fwrite(Smat[i], sizeof(int), nV_spcvec[i], file);
  fclose(file);
  printf("(LD_svd_bundle::write_LD_svd_bundle) wrote %s\n",name_);
}
