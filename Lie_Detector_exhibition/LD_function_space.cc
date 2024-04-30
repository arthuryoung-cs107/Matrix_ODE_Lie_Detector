#include "LD_function_space.hh"
#include "LD_aux.hh"
#include "LD_io.hh"

vxu_workspace::vxu_workspace(int nvar_, int ord_len_):
xu(new double[nvar_]), xu_vals(Tmatrix<double>(nvar_,ord_len_)) {}
vxu_workspace::~vxu_workspace() {delete xu; free_Tmatrix<double>(xu_vals);}

function_space::function_space(ode_solspc_meta &meta_, int bor_, int perm_len_):
ode_solspc_element(meta_), bor(bor_), perm_len(perm_len_),
dof_flags_mat(Tmatrix<bool>(2,ndof_full)),
order_mat(Tmatrix<int>(perm_len,nvar)), dof_indices_mat(Tmatrix<int>(2,ndof_full)),
theta_full(new double[ndof_full]),
dof_hot_flags_mat(new bool*[nvar]), dof_tun_flags_mat(new bool*[nvar]),
dof_hot_indices_mat(new int*[nvar]), dof_tun_indices_mat(new int*[nvar]),
theta_mat(new double*[nvar])
{
  for (size_t i = 0, itheta= 0; i < nvar; i++, itheta += perm_len)
  {
    dof_hot_flags_mat[i] = dof_hot_flags + itheta;
    dof_tun_flags_mat[i] = dof_tun_flags + itheta;
    dof_hot_indices_mat[i] = dof_hot_indices + itheta;
    dof_tun_indices_mat[i] = dof_tun_indices + itheta;
    theta_mat[i] = theta_full + itheta;
  }
  LD_linalg::fill_x(dof_hot_flags,ndof_full,true);
  LD_linalg::fill_x(dof_tun_flags,ndof_full,true);
  LD_linalg::fill_x_012(dof_hot_indices,ndof_full);
  LD_linalg::fill_x_012(dof_tun_indices,ndof_full);
}
function_space::~function_space()
{
  free_Tmatrix<bool>(dof_flags_mat);
  free_Tmatrix<int>(order_mat); free_Tmatrix<int>(dof_indices_mat);
  delete [] theta_full;
  delete [] dof_hot_flags_mat; delete [] dof_tun_flags_mat;
  delete [] dof_hot_indices_mat; delete [] dof_tun_indices_mat;
  delete [] theta_mat;
}
void function_space::vxu_eval(double *xu_,double *vxu_,vxu_workspace &wkspc_)
{
  for (int i = 0; i < ndim; i ++) vxu_[i] = 0.0; // clear output v eval space
  init_vxu(xu_,wkspc_); // precompute powers of variables and derivatives
  for (size_t ilam = 0; ilam < perm_len; ilam++)
  {
    double  lambda_i = comp_lambda_i(ilam,wkspc_.xu_vals);
    for (int ivar = 0; ivar < nvar; ivar++)
      if (dof_hot_flags_mat[ivar][ilam])
        vxu_[ivar] += theta_mat[ivar][ilam]*lambda_i;
  }
}

orthopolynomial_space::orthopolynomial_space(ode_solspc_meta &meta_, int bor_):
power_space(meta_,bor_), icoeff_mat(Tsym<int>(bor+1)),
poly_coeffs(Tsym_lower<double>(bor+1)), map_params(Tmatrix<double>(4,nvar))
{
  // precompute exponent multipliers for polynomial derivatives
  for (int j = 0; j <= bor; j++)
  {
    icoeff_mat[0][j] = 1;
    for (int i = 1, jpow=j, jshift=j-1; i <= j; i++, jpow--, jshift--)
      icoeff_mat[i][jshift] = icoeff_mat[i-1][jshift+1]*jpow;
  }
  // initialize linear mappings to identity
  for (size_t i = 0; i < nvar; i++)
  {
    fmap_m[i] = bmap_m[i] = 1.0;
    fmap_b[i] = bmap_b[i] = 0.0;
  }
}
orthopolynomial_space::~orthopolynomial_space()
{
  free_Tmatrix<int>(icoeff_mat);
  free_Tmatrix<double>(poly_coeffs); free_Tmatrix<double>(map_params);
}

partial_chunk::partial_chunk(int perm_len_, int eor_, int ndep_):
chunk(new double[perm_len_*(eor_+1)*(ndep_+1)]),
Jac_mat(new double*[ndep_+1]),
C_x(new double**[perm_len_]), C_u(new double*[perm_len_]),
C_x_ptrs(new double*[perm_len_*eor_]), par_chunk(chunk+(perm_len_*(ndep_+1)))
{
  for (int ivar = 0, iperm = 0; ivar <= ndep_; ivar++, iperm+=perm_len_) Jac_mat[ivar] = chunk + iperm;
  for (int ic = 0, i_C_x = 0, i_skip_x = 0, i_skip_Cx = 0, i_skip_u = perm_len_*eor_*ndep_; ic < perm_len_; ic++, i_skip_Cx+=eor_, i_skip_u+=eor_)
  {
    C_x[ic] = C_x_ptrs + i_skip_Cx;
    C_u[ic] = par_chunk + i_skip_u;
    for (int k_der = 0; k_der < eor_; k_der++, i_C_x++, i_skip_x += ndep_)
      C_x_ptrs[i_C_x] = par_chunk + i_skip_x;
  }
}
partial_chunk::~partial_chunk()
{
  delete [] chunk; delete [] Jac_mat;
  delete [] C_x; delete [] C_x_ptrs; delete [] C_u;
}

function_space_basis::function_space_basis(function_space &fspc_):
function_space_element(fspc_),
partials(partial_chunk(perm_len,eor,ndep)),
ord_len(fspc_.comp_ord_len()), ord_i_len(fspc_.comp_ord_i_len()),
dof_hot_flags_i(new bool[nvar]),
ord_xu(new int[nvar]), theta_xu(new double[nvar]),
s_j(new double[ndim]), v_j(new double[ndim]),
xu_vals(Tmatrix<double>(nvar,ord_len)),
dxu_vals(Tmatrix<double>(ndep,eor)),
dxu_pow_chunk((dxu_pows_allocate)?(new double[ndep*(((eorm2)*(eorm2+1))/2)]):(NULL)),
dnxu_pow_ptrs((dxu_pows_allocate)?(new double*[ndep*(eorm2)]):(NULL)),
dnxu_val((dxu_pows_allocate)?(new double**[ndep]):(NULL)),
couples(new coupling_term*[eor+1])
{
  for (size_t i = 0; i <= eor; i++) couples[i] = new coupling_term(meta,i);
  if (dxu_pows_allocate) // if we need extra space for high order differential equations
  {
    int delchunk = ((eorm2)*(eorm2+1))/2;
    for (int idep = 0; idep < ndep; idep++)
    {
      dnxu_val[idep] = dnxu_pow_ptrs+(idep*(eorm2));
      dnxu_val[idep][0] = dxu_pow_chunk + (idep*delchunk);
      for (int ixp = 1, iskip = eorm2; ixp < eorm2; ixp++, iskip--) dnxu_val[idep][ixp] = dnxu_val[idep][ixp-1] + iskip;
    }
  }
}
function_space_basis::~function_space_basis()
{
  delete [] dof_hot_flags_i;
  delete [] ord_xu; delete [] theta_xu;
  delete [] s_j; delete [] v_j;
  free_Tmatrix<double>(xu_vals);
  free_Tmatrix<double>(dxu_vals);
  if (dxu_pows_allocate) {delete [] dxu_pow_chunk; delete [] dnxu_pow_ptrs; delete [] dnxu_val;}
  for (int i = 0; i <= eor; i++) delete couples[i];
  delete [] couples;
}

void function_space_basis::v_eval(double *s_,double *v_)
{
  for (int i = 0; i < ndim; i ++) v_[i] = 0.0; // clear output v eval space
  init_workspace(s_); // precompute powers of variables and derivatives
  coupling_term &c0 = *(couples[0]);
  for (int ixord = 0, i_L_start = 0; ixord < ord_len; i_L_start+=ord_i_len, ixord++)
  {
    double Li_x = stage_indep_var(ixord);
    for (int iu_perm = 0, i_L = i_L_start; iu_perm < ord_i_len; iu_perm++, i_L++)
    {
      double  Li_u = stage_dep_var(i_L),
              Li = Li_x*Li_u;
      for (int ivar = 0; ivar < nvar; ivar++)
        if (dof_hot_flags_i[ivar] = dof_hot_flags_mat[ivar][i_L])
          v_[ivar] += ((theta_xu[ivar] = theta_mat[ivar][i_L])*Li);

      if (dof_hot_flags_i[0])
      {
        stage_coupling(-1.0); compute_x_coupling(c0);
        for (int idnxu = nvar; idnxu < ndim; idnxu++) v_[idnxu] += theta_xu[0]*v_j[idnxu];
      }
      stage_coupling(1.0); compute_u_coupling(1,c0);
      for (int ider = 1, idnxu = nvar, idnxu_skip = idnxu; ider <= eor; ider++, idnxu_skip+=ndep)
        for (int iu = 1; iu < nvar; iu++, idnxu++)
          if (dof_hot_flags_i[iu]) // accumulate in each dimension, scaled by corresponding coeff
            v_[idnxu] += theta_xu[iu]*v_j[idnxu_skip];
    }
  }
}
void function_space_basis::fill_partial_chunk(double *s_)
{
  init_workspace(s_);
  coupling_term &c0 = *(couples[0]);
  for (int ixord = 0, i_L_start = 0; ixord < ord_len; i_L_start+=ord_i_len, ixord++)
  {
    double Li_x = stage_indep_var(ixord);
    for (int iu_perm = 0, i_L = i_L_start; iu_perm < ord_i_len; iu_perm++, i_L++)
    {
      double  Li_u = stage_dep_var(i_L),
              Li = Li_x*Li_u;
      // store Ltxi data for later processing
      bool  indep_i_tun = false,
            any_dep_i_tun = false;
      if (indep_i_tun = dof_tun_flags_mat[0][i_L]) Jac_mat[0][i_L] = Li;
      for (int ivar = 1; ivar <= ndep; ivar++)
      {
        any_dep_i_tun = (any_dep_i_tun)||(dof_tun_flags_mat[ivar][i_L]);
        if (dof_tun_flags_mat[ivar][i_L]) Jac_mat[ivar][i_L] = Li;
      }
      if (indep_i_tun)
      {
        stage_coupling(-1.0); compute_x_coupling(c0); // perform dtde couplings
        // accumulate dt/de contribution De, store partial derivative
        for (int idnxu = nvar, iparx = 0; idnxu < ndim; idnxu++, iparx++)
          C_x[i_L][0][iparx] = v_j[idnxu];
      }
      // perform dxde couplings
      if (any_dep_i_tun)
      {
        stage_coupling(1.0); compute_u_coupling(1,c0);
        for (int ider = 1, idnxu = nvar, idnxu_skip = idnxu; ider <= eor; ider++, idnxu_skip+=ndep)
          C_u[i_L][ider-1] = v_j[idnxu_skip];
      }
    }
  }
}
void function_space_basis::compute_x_coupling(coupling_term &c0_)
{
  double  ncoeff,
          dLsi,
          dLy;
  check_coupling_derivatives(c0_);
  if (c0_.xderiv_flag)
  {
    double acc =  comp_dx_dnLi(c0_,ncoeff,dLsi);
    for (int k = 1, iv = nvar; k <= eor; k++) // apply resultant dt L to all orders
    {
      // apply resultant dt L to each dim (idep) of each order (k)
      for (int idep = 0; idep < ndep; idep++, iv++) // premultiply by corresponding derivative term (k'th derivative of idep'th dep. var.)
        v_j[iv] += acc*s_j[iv]; // accumulate in corresponding dp/de dimension
      if (k<eor) // if term contributes to k+1 prolongation
      {
        stage_dx_couple(c0_,*(couples[k]),ncoeff,dLsi); // copy c0_ information regarding 1st derivative wrt t to ck_ coupling
        compute_u_coupling(k+1,*(couples[k]),k); // pass ck to next order of prolongation
      }
    }
  }
  for (int iu = 1; iu <= ndep; iu++)
    if (c0_.deriv_flags[iu])
    {
      double acc = comp_dui_dnLi(c0_,iu,ncoeff,dLsi,dLy)*s_j[iu+ndep];
      for (int k = 1, iv = nvar; k <= eor; k++) // apply resultant dt L to all orders
      {
        // apply resultant dt L to each dim (idep) of each order (k)
        for (int idep = 0; idep < ndep; idep++, iv++) // premultiply by corresponding derivative term (k'th derivative of idep'th dep. var.)
          v_j[iv] += acc*s_j[iv]; // accumulate in corresponding dp/de dimension
        if (k<eor) // if term contributes to k+1 prolongation
        {
          stage_du_couple(iu,c0_,*(couples[k]),ncoeff,dLsi,dLy); // copy ckm1_ information regarding freshly computed derivative wrt i'th dep. var. to ck coupling
          compute_u_coupling(k+1,*(couples[k]),k); // pass ck to next order of prolongation
        }
      }
    }
}
void function_space_basis::compute_u_coupling(int k_, coupling_term &ckm1_)
{
  bool continue_coupling = k_<eor;
  int k_m1 = k_-1;
  double  ckm1_dnxu_power_product = comp_dnxu_power_product(ckm1_),
          &v_j_dkui = v_j[comp_v_k_start_index(k_)],
          ncoeff,
          dLsi,
          dLy;
  coupling_term &ck = *(couples[k_]);
  check_coupling_derivatives(ckm1_);
  if (ckm1_.xderiv_flag)
  {
    v_j_dkui += comp_dx_dnLi(ckm1_,ncoeff,dLsi)*ckm1_dnxu_power_product;
    if (continue_coupling)
    {
      stage_dx_couple(ckm1_,ck,ncoeff,dLsi); // copy ckm1_ information regarding freshly computed derivative wrt t to ck coupling
      compute_u_coupling(k_+1,ck); // pass ck to next order of prolongation
    }
  }
  for (int iu = 1; iu <= ndep; iu++)
    if (ckm1_.deriv_flags[iu])
    {
      v_j_dkui += comp_dui_dnLi(ckm1_,iu,ncoeff,dLsi,dLy)*ckm1_dnxu_power_product*s_j[iu+ndep];
      if (continue_coupling)
      {
        stage_du_couple(iu,ckm1_,ck,ncoeff,dLsi,dLy);
        compute_u_coupling(k_+1,ck);
      }
    }
  // for xp powers, we strictly have polynomial characteristics
  double Lxu = ckm1_.Lxh*ckm1_.Luh;
  for (int kk = 1, ind_dnxuj = nvar; kk <= k_m1; kk++) // consider dxp input L contribution for each order kk up to kk = k-1
    for (int jdep = 0; jdep < ndep; jdep++, ind_dnxuj++) // consider dkk xpj (partial wrt kk'th derivative of j'th dep. var.) for input L
      if (ckm1_.Pvec[ind_dnxuj]>0) // if input L has derivative in dkk xpj
      {
        ncoeff = ((double) ckm1_.Pvec[ind_dnxuj])*ckm1_.coeff, // update scalar coefficient according to power rule
        v_j_dkui += comp_ddnxui_power_product(ind_dnxuj,ckm1_,ncoeff*Lxu)*s_j[ind_dnxuj+ndep]; // accumulate products forming L evaluation, chaining j'th dep. var. derivative
        if (continue_coupling)
        {
          stage_ddnxu_couple(ind_dnxuj,ckm1_,ck,ncoeff);
          compute_u_coupling(k_+1,ck);
        }
      }
}
void function_space_basis::compute_u_coupling(int k_, coupling_term &ckm1_, int indep_source_)
{
  bool continue_coupling = k_<eor;
  int k_m1 = k_-1,
      iv_start = comp_v_k_start_index(k_), // start index of k'th prolongation dimension
      ipre_dnxu_start = 1+(indep_source_*ndep); // start index of premultiplied derivative dep. var. terms
  double  ckm1_dnxu_power_product = comp_dnxu_power_product(ckm1_),
          ncoeff,
          dLsi,
          dLy;
  coupling_term &ck = *(couples[k_]);
  check_coupling_derivatives(ckm1_);
  if (ckm1_.xderiv_flag)
  {
    double acc = comp_dx_dnLi(ckm1_,ncoeff,dLsi)*ckm1_dnxu_power_product;
    for (int idep = 0, iv = iv_start, ipre_dnxu = ipre_dnxu_start; idep < ndep; idep++, iv++, ipre_dnxu++)
      v_j[iv] += acc*s_j[ipre_dnxu]; // accumulate in corresponding dp/de dimension
    if (continue_coupling)
    {
      stage_dx_couple(ckm1_,ck,ncoeff,dLsi); // copy ckm1_ information regarding freshly computed derivative wrt t to ck coupling
      compute_u_coupling(k_+1,ck,indep_source_); // pass ck to next order of prolongation
    }
  }
  for (int iu = 1; iu <= ndep; iu++)
    if (ckm1_.deriv_flags[iu])
    {
      double acc = comp_dui_dnLi(ckm1_,iu,ncoeff,dLsi,dLy)*ckm1_dnxu_power_product*s_j[iu+ndep];
      for (int idep = 0, iv = iv_start, ipre_dnxu = ipre_dnxu_start; idep < ndep; idep++, iv++, ipre_dnxu++)
        v_j[iv] += acc*s_j[ipre_dnxu]; // accumulate in corresponding dp/de dimension
      if (continue_coupling)
      {
        stage_du_couple(iu,ckm1_,ck,ncoeff,dLsi,dLy); // copy ckm1_ information regarding freshly computed derivative wrt i'th dep. var. to ck coupling
        compute_u_coupling(k_+1,ck,indep_source_); // pass ck to next order of prolongation
      }
    }
  // for xp powers, we strictly have polynomial characteristics
  double Lxu = ckm1_.Lxh*ckm1_.Luh;
  for (int kk = 1, ind_dnxuj = nvar; kk <= k_m1; kk++) // consider dxp input L contribution for each order kk up to kk = k-1
    for (int jdep = 0; jdep < ndep; jdep++, ind_dnxuj++) // consider dkk xpj (partial wrt kk'th derivative of j'th dep. var.) for input L
      if (ckm1_.Pvec[ind_dnxuj]>0) // if input L has derivative in dkk xpj
      {
        ncoeff = ((double) ckm1_.Pvec[ind_dnxuj])*ckm1_.coeff; // update scalar coefficient according to power rule
        double acc = comp_ddnxui_power_product(ind_dnxuj,ckm1_,ncoeff*Lxu)*s_j[ind_dnxuj+ndep]; // accumulate products forming L evaluation, chaining j'th dep. var. derivative
        for (int idep = 0, iv = iv_start, ipre_dnxu = ipre_dnxu_start; idep < ndep; idep++, iv++, ipre_dnxu++)
          v_j[iv] += acc*s_j[ipre_dnxu]; // accumulate in corresponding dp/de dimension
        if (continue_coupling)
        {
          stage_ddnxu_couple(ind_dnxuj,ckm1_,ck,ncoeff); // copy ckm1_ information regarding freshly computed derivative wrt kk'th derivative of j'th dep. var. to ck coupling
          compute_u_coupling(k_+1,ck,indep_source_); // pass ck to next order of prolongation
        }
      }
  /*
    Finally, add contribution of product rule term for premultiplied vector
    Initiate accumulation of products forming L evaluation.
    Note that this is equivalent to input L evaluation, but premultiplied by
    derivatives of premultiplied vector
  */
  double acc = ckm1_.coeff*Lxu*ckm1_dnxu_power_product;
  // compute accumulation to each dep. var. dimension
  for (int jjdep = 0, iv = iv_start, ipre_dnxu = ipre_dnxu_start+ndep; jjdep < ndep; jjdep++, iv++, ipre_dnxu++)
    v_j[iv] += acc*(s_j[ipre_dnxu]); // accumulate in corresponding dp/de dimension
  if (continue_coupling)
  {
    spawn_ddnxu_couple(k_,ckm1_,ck); // copy ckm1_ information regarding freshly computed derivative wrt kk'th derivative of j'th dep. var. to ck coupling
    compute_u_coupling(k_+1,ck,k_); // pass ck to next order of prolongation
  }
}

void power_basis::init_workspace(double *s_)
{
  for (int i = 0; i < ndim; i++) s_j[i] = s_[i]; // copy the input position vector
  for (int ivar = 0; ivar < nvar; ivar++) // precompute powers of indep. and dep. vars.
  {
    xu_vals[ivar][0] = 1.0;
    for (int ipow = 1; ipow <= bor; ipow++) xu_vals[ivar][ipow] = xu_vals[ivar][ipow-1]*s_j[ivar];
  }
  init_workspace_dnxu_vals();
}
void orthopolynomial_basis::debugging_description()
{

  printf("(orthopolynomial_basis::debugging_description) order %d, %d degrees of freedom\n", bor, ndof_full);

  printf("forward mappings: x_hat = %.2e x + %.2e, ", map_params[0][0], map_params[1][0]);
  for (int iu= 1; iu <= ndep; iu++) printf("u%d_hat = %.2e u%d + %.2e, ", iu, map_params[0][iu], iu, map_params[1][iu]);
  printf("\ninverse mappings: x = %.2e x_hat + %.2e, ", map_params[2][0], map_params[3][0]);
  for (int iu = 1; iu <= ndep; iu++) printf("u%d = %.2e u%d_hat + %.2e, ", iu, map_params[2][iu], iu, map_params[3][iu]);

  printf("\nx  orders: ");
  for (int i = 0; i < perm_len; i++) printf("%d ", order_mat[i][0]);
  for (int iu = 1; iu <= ndep; iu++)
  {
    printf("\nu%d orders: ", iu);
    for (int i = 0; i < perm_len; i++) printf("%d ", order_mat[i][iu]);
  }
  printf("\npolynomials:\n");
  for (int i = 0; i <= bor; i++)
  {
    printf("L%d: ",i);
    for (int j = 0; j <= i; j++) printf(" + (%.2e)u^%d", poly_coeffs[i][j], j);
    printf("\n");
  }
}

void orthopolynomial_space::debugging_description()
{
    // printf("orthogonal domain: (%.2e, %.2e)", h_min, h_max);
  // if (solspc_val_extrema!=NULL)
  // {
  //   printf("\n%.2e < t < %.2e, ", solspc_val_extrema[0][0], solspc_val_extrema[1][0]);
  //   for (int ix = 1; ix <= ndep; ix++) printf("%.2e < x%d < %.2e, ", solspc_val_extrema[0][ix], ix, solspc_val_extrema[1][ix]);
  // }
  printf("(orthopolynomial_space::debugging_description) order %d, %d degrees of freedom\n", bor, ndof_full);

  printf("forward mappings: x_hat = %.2e x + %.2e, ", fmap_m[0], fmap_b[0]);
  for (int iu= 1; iu <= ndep; iu++) printf("u%d_hat = %.2e u%d + %.2e, ", iu, fmap_m[iu], iu, fmap_b[iu]);
  printf("\ninverse mappings: x = %.2e x_hat + %.2e, ", bmap_m[0], bmap_b[0]);
  for (int iu = 1; iu <= ndep; iu++) printf("u%d = %.2e u%d_hat + %.2e, ", iu, bmap_m[iu], iu, bmap_b[iu]);

  printf("\nx  orders: ");
  for (int i = 0; i < perm_len; i++) printf("%d ", order_mat[i][0]);
  for (int iu = 1; iu <= ndep; iu++)
  {
    printf("\nu%d orders: ", iu);
    for (int i = 0; i < perm_len; i++) printf("%d ", order_mat[i][iu]);
  }
  printf("\npolynomials:\n");
  for (int i = 0; i <= bor; i++)
  {
    printf("L%d: ",i);
    for (int j = 0; j <= i; j++) printf(" + (%.2e)u^%d", poly_coeffs[i][j], j);
    printf("\n");
  }
}
void orthopolynomial_space::configure_self(const char name_[])
{
  const int hlen_check = 3,
            omlen_check = (perm_len*nvar),
            mlen_check = 4*nvar,
            slen_check = (((bor+1)*(bor+2))/2);
  int hlen,
      header[hlen_check],
      &ord_mat_len = header[0],
      &map_len = header[1],
      &len_sym = header[2];
  FILE * file = LD_io::fopen_SAFE(name_,"r");
  LD_io::fread_SAFE(&hlen,sizeof(int),1,file);
  if (hlen == hlen_check)
    LD_io::fread_SAFE(header,sizeof(int),hlen,file);
    else
    {
      printf("(orthopolynomial_space::configure_self) skipped header (len = %d)\n", hlen);
      LD_io::fseek_SAFE(file,hlen*sizeof(int), SEEK_CUR);
    }
  if (ord_mat_len == omlen_check) LD_io::fread_SAFE(order_mat[0],sizeof(int),ord_mat_len,file);
    else
    {
      printf("(orthopolynomial_space::configure_self) skipped order mat (len = %d)\n", ord_mat_len);
      LD_io::fseek_SAFE(file,ord_mat_len*sizeof(int), SEEK_CUR);
    }
  if (map_len == mlen_check) LD_io::fread_SAFE(map_params[0],sizeof(double),map_len,file);
    else
    {
        printf("(orthopolynomial_space::configure_self) skipped mapping mat (len = %d)\n", map_len);
        LD_io::fseek_SAFE(file,map_len*sizeof(double), SEEK_CUR);
    }
  if (len_sym == slen_check) LD_io::fread_SAFE(poly_coeffs[0],sizeof(double),len_sym,file);
    else
    {
      printf("(orthopolynomial_space::configure_self) skipped mapping mat (len = %d)\n", len_sym);
      LD_io::fseek_SAFE(file,len_sym*sizeof(double), SEEK_CUR);
    }
  LD_io::fclose_SAFE(file);
  printf("(orthopolynomial_space::configure_self) read %s\n", name_);
}
void orthopolynomial_space::write_configuration_file(const char name_[])
{
  int hlen = 3,
      ord_mat_len = perm_len*nvar,
      map_len = 4*nvar,
      len_sym = (((bor+1)*(bor+2))/2),
      header[] = {hlen, ord_mat_len, map_len, len_sym};
  FILE * file = LD_io::fopen_SAFE(name_,"wb");

  fwrite(header,sizeof(int),hlen+1,file);
  fwrite(order_mat[0],sizeof(int),ord_mat_len,file);
  fwrite(map_params[0],sizeof(double),map_len,file);
  fwrite(poly_coeffs[0],sizeof(double),len_sym,file);
  LD_io::fclose_SAFE(file);
  printf("(orthopolynomial_space::write_configuration_file) wrote %s\n", name_);
}
