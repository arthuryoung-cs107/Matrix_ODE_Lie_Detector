#include "LD_parameter_space.hh"

#include <cstring>

/*
  constructors
*/

LD_vector_space::LD_vector_space(int vlen_): data_owner(true),
  Vmat_data(Tmatrix<double>(vlen_,vlen_)),
  vlen_full(vlen_), nvec_use(vlen_), vlen_use(vlen_),
  inds_V(new int[vlen_]), Vmat(new double*[vlen_])
  {for (size_t i = 0; i < vlen_; i++) Vmat[i] = Vmat_data[inds_V[i] = i];}
LD_vector_space::LD_vector_space(int vlen_,int *inds_V_,double **Vmat_data_,double **Vmat_): data_owner(false),
  Vmat_data(Vmat_data_),
  vlen_full(vlen_), nvec_use(vlen_), vlen_use(vlen_),
  inds_V(inds_V_), Vmat(Vmat_) {}
LD_vector_space::~LD_vector_space()
  {if (data_owner) {delete [] inds_V; free_Tmatrix<double>(Vmat_data); delete [] Vmat;}}

LD_vector_bundle::LD_vector_bundle(int nspc_, int vlen_): data_owner(true),
  Vtns_data(T3tensor<double>(nspc_,vlen_,vlen_)), Vrows(new double*[nspc_*vlen_]),
  rec_ptr(new LD_vspace_record(nspc_,vlen_,vlen_,Vtns_data)),
  nspc(nspc_), vlen_full(vlen_), nvec_use(vlen_), vlen_use(vlen_),
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
  nspc(bndle_.nspc), vlen_full(bndle_.vlen_full), nvec_use(bndle_.nvec_use), vlen_use(bndle_.vlen_use),
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
  ndof_spcvec(new int[nspc_]), pbse_nvar_spcmat(Tmatrix<int>(nspc_,nvar)),
  Tspaces(new LD_Theta_space*[nspc_])
{
  for (size_t i = 0; i < nspc_; i++)
  {
    Tspaces[i] = new LD_Theta_space(meta_,ndof_spcvec[i]=ndof_,Vspaces[i],Yspaces_nvar+(nvar*i),pbse_nvar_spcmat[i]);
    for (size_t j = 0; j < nvar; j++) {*(Yspaces_nvar+(nvar*i) + j) = NULL; pbse_nvar_spcmat[i][j] = perm_len;}
  }
}
LD_Theta_bundle::LD_Theta_bundle(LD_Theta_bundle &Tbndle_): ode_solspc_element(Tbndle_.meta),
  Vbndl_owner(false), data_owner(false),
  Vbndle_ptr(Tbndle_.Vbndle_ptr), Yspaces_nvar(Tbndle_.Yspaces_nvar), perm_len(perm_len),
  ndof_spcvec(Tbndle_.ndof_spcvec), pbse_nvar_spcmat(Tbndle_.pbse_nvar_spcmat),
  Tspaces(Tbndle_.Tspaces) {}

LD_Theta_bundle::~LD_Theta_bundle()
{
  if (data_owner)
  {
    for (size_t i = 0; i < nspc; i++) delete Tspaces[i];
    delete [] Tspaces;
    delete [] ndof_spcvec;
    free_Tmatrix<int>(pbse_nvar_spcmat);
    delete [] Yspaces_nvar;
  }
  if (Vbndl_owner) delete Vbndle_ptr;
}

/*
  definitions
*/

void LD_vspace_record::print_selected_details(const char name_[], bool longwinded_)
{
  if (longwinded_)
  {
    printf("(LD_vspace_record::print_selected_details) %s isat_vec\n", name_);
    for (size_t ispc = 0; ispc < nspc; ispc++)
    {
      printf("spc %d (nsat = %d): ", ispc, nV_spcvec[ispc]);
      for (size_t isat = 0; isat < nV_spcvec[ispc]; isat++) printf("%d ", iV_spcmat[ispc][isat]);
      printf("\n");
    }
  }
  char name_buf[strlen(name_) + 20];
  sprintf(name_buf,"  %s nsat_vec", name_);
  LD_linalg::print_xT(name_buf,nV_spcvec,nspc);
  printf("  nvec_use = %d, vlen_use = %d, min nsat = %d, max nsat = %d \n",
            nvec, vlen,
            LD_linalg::min_val<int>(nV_spcvec,nspc), LD_linalg::max_val<int>(nV_spcvec,nspc));
}

void LD_Theta_space::post_multiply_evry_basis(double **AYmat_, double **Amat_, int mrows_)
{
  for (size_t irow = 0; irow < mrows_; irow++)
    for (size_t ivar = 0, jcolAY = 0, jcolA = 0; ivar < nvar; ivar++, jcolAY+=pbse_nvar_use[ivar], jcolA+=perm_len)
    {
      if (bse_ptrs[ivar]==NULL) for (size_t jbse = 0; jbse < perm_len; jbse++) AYmat_[irow][jcolAY+jbse] = Amat_[irow][jcolA+jbse];
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
    for (size_t ivar = ivar_+1; ivar < nvar; ivar++, jcolAY+=pbse_nvar_use[ivar], jcolA+=perm_len)
      for (size_t jbse = 0; jbse < perm_len; jbse++) AYmat_[irow][jcolAY+jbse] = Amat_[irow][jcolA+jbse];
  }
}
int LD_Theta_space::init_Vspce_premult(double **Wmat_)
{
  const int nV = verify_Yspaces();
  for (size_t irow = 0; irow < nV; irow++)
  {
    double * const vi = spc.Vmat[irow] = spc.Vrowi_data(spc.inds_V[irow] = irow);
    for (size_t ivar = 0, jcolV = 0, jcolW = 0; ivar < nvar; ivar++, jcolV+=perm_len, jcolW+=pbse_nvar_use[ivar])
    {
      if (bse_ptrs[ivar]==NULL) for (size_t jbse = 0; jbse < perm_len; jbse++) vi[jcolV+jbse] = Wmat_[irow][jcolW+jbse];
      else bse_ptrs[ivar]->pre_multiply_colvec(vi+jcolV,Wmat_[irow]+jcolW);
    }
  }
  return spc.nvec_use = ndof_use;
}
