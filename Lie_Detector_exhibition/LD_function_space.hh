#ifndef LD_FCN_SPC_HH
#define LD_FCN_SPC_HH

// #define M_PI 3.14159265358979323846

#include "LD_ode.hh"

// #include <cstdio>

#ifdef _OPENMP
  #include "omp.h"
#endif

struct permutation_computer
{
  permutation_computer() {}
  ~permutation_computer() {}

  static int comp_full_power_perm_len(int bor_, int nvar_)
    {int out = 1; for (size_t i = 0; i < nvar_; i++) out*=(bor_+1); return out;}
  static void set_order_mat_full_power(int ** omat_, int bor_, int nvar_)
  {
    int k_perm = 0;
    for (int i = 0; i <= bor_; i++)
    {
      int delk = permutation_computer::set_order_mat_full_power_recursive(omat_,bor_,nvar_,1,k_perm);
      for (int k = k_perm; k < k_perm+delk; k++) omat_[k][0] = i;
      k_perm+=delk;
    }
  }
  static int set_order_mat_full_power_recursive(int **omat_, int o_, int d_, int ilevel_, int k_perm_)
  {
    if (ilevel_>=d_) return 0;
    else
    {
      int delk_level = 0;
      if (ilevel_ == d_-1) for (int i = 0; i <= o_; i++, delk_level++) omat_[k_perm_+delk_level][ilevel_] = i; // last
      else for (int i = 0; i <= o_; i++)
        {
          int delk = permutation_computer::set_order_mat_full_power_recursive(omat_,o_,d_,ilevel_+1,delk_level+k_perm_);
          for (int k = 0; k < delk; k++, delk_level++) omat_[k_perm_+delk_level][ilevel_] = i;
        }
      return delk_level;
    }
  }
};

struct vxu_workspace
{
  vxu_workspace(int nvar_, int ord_len_);
  ~vxu_workspace();

  double  * const xu,
          ** const xu_vals;
};

struct function_space: public ode_solspc_element
{
  function_space(ode_solspc_meta &meta_, int bor_, int perm_len_);
  ~function_space();

  const int bor,
            perm_len,
            ndof_full = perm_len*(ndep+1);

  bool  ** const dof_flags_mat,
        * const dof_hot_flags = dof_flags_mat[0],
        * const dof_tun_flags = dof_flags_mat[1],
        ** const dof_hot_flags_mat,
        ** const dof_tun_flags_mat;

  int ** const order_mat,
      ** const dof_indices_mat,
      * const dof_hot_indices = dof_indices_mat[0],
      * const dof_tun_indices = dof_indices_mat[1],
      ** const dof_hot_indices_mat,
      ** const dof_tun_indices_mat,
      ndof_hot,
      ndof_tun;

  double  * const theta_full,
          ** const theta_mat;

  virtual int comp_ord_len() = 0;
  virtual int comp_ord_i_len() = 0;

  inline void lamvec_eval(double *xu_, double *lamvec_, vxu_workspace &wkspc_)
  {
    init_vxu(xu_,wkspc_); // precompute powers of variables
    for (size_t ilam = 0; ilam < perm_len; ilam++)
      lamvec_[ilam] = comp_lambda_i(ilam,wkspc_.xu_vals);
  }
  void vxu_eval(double *xu_,double *v_,vxu_workspace &wkspc_);
  void vx_spc_eval(double *xu_, double *vx_, vxu_workspace &wkspc_, double ** Kmat_, int kappa_);
  virtual void init_vxu(double *xu_in_, vxu_workspace &wkspc_) = 0;
  inline double comp_lambda_i(int ilam_, double **xu_vals_)
  {
    int * const Ovec_i = order_mat[ilam_];
    double  out = 1.0;
    for (size_t ivar = 0; ivar < nvar; ivar++) out *= xu_vals_[ivar][Ovec_i[ivar]];
    return out;
  }

  virtual double eval_Li(double *v_, int O_) = 0;
  virtual double eval_dnLi(int n_, double *v_, int O_) = 0;
  virtual double get_dcoeff(int idim_, int O_) = 0;

  #ifdef _OPENMP
  inline int thread_num() {return omp_get_thread_num();}
  inline int get_nt() {return omp_get_max_threads();}
  #else
  inline int thread_num() {return 0;}
  inline int get_nt() {return 1;}
  #endif

};

struct power_space: public function_space
{
  power_space(ode_solspc_meta &meta_, int bor_):
  function_space(meta_,bor_,permutation_computer::comp_full_power_perm_len(bor_,meta_.nvar))
    {permutation_computer::set_order_mat_full_power(order_mat,bor,nvar);}
  ~power_space() {}

  int comp_ord_len() {return bor+1;}
  int comp_ord_i_len() {return perm_len/comp_ord_len();}

  inline void init_xu_vals(double *xu_,double **xu_vals_)
  {
    for (int ivar = 0; ivar < nvar; ivar++) // precompute powers of indep. and dep. vars.
    {
      xu_vals_[ivar][0] = 1.0;
      for (int ipow = 1; ipow <= bor; ipow++) xu_vals_[ivar][ipow] = xu_vals_[ivar][ipow-1]*xu_[ivar];
    }
  }
};

struct orthopolynomial_space: public power_space
{
  orthopolynomial_space(ode_solspc_meta &meta_, int bor_);
  ~orthopolynomial_space();

  int ** const icoeff_mat;
  double  ** const poly_coeffs,
          ** const map_params,
          * const fmap_m = map_params[0],
          * const fmap_b = map_params[1],
          * const bmap_m = map_params[2],
          * const bmap_b = map_params[3];

  void debugging_description();
  void configure_self(const char name_[]);
  void write_configuration_file(const char name_[]);

  inline void set_centered_domain(double **sve_, double h_min_, double h_max_)
  {
    for (int ivar = 0; ivar <= ndep; ivar++)
    {
      double  v_min = sve_[0][ivar],
              v_max = sve_[1][ivar],
              fvmap_m = fmap_m[ivar] = (h_max_-h_min_)/(v_max-v_min),
              fvmap_b = fmap_b[ivar] = 0.5*(h_min_+h_max_+((h_min_-h_max_)*(v_max+v_min)/(v_max-v_min))),
              bvmap_m = bmap_m[ivar] = 1.0/fvmap_m,
              bvmap_b = bmap_b[ivar] = (-1.0*fvmap_b)/(fvmap_m);
    }
  }
  inline void set_0maxmag_0pi05_domain(double **sme_, double h_max_tru_, double h_max_scl_)
  {
    // map the time domain to be be 0 - max magnitude, always,
    fmap_m[0]=h_max_tru_/sme_[1][0];
    bmap_m[0]=1.0/fmap_m[0];
    fmap_b[0]=bmap_b[0]=0.0;
    double  pi05 = 3.14159265358979323846/2.0,
            hmax_use = h_max_tru_*h_max_scl_;
    for (int ivar = 1; ivar <= ndep; ivar++)
    {
      fmap_m[ivar] = hmax_use/((sme_[1][ivar]>pi05)?(sme_[1][ivar]):(pi05));
      bmap_m[ivar] = 1.0/fmap_m[ivar];
      fmap_b[ivar] = bmap_b[ivar] = 0.0;
    }
  }

  void init_vxu(double *xu_in_, vxu_workspace &wkspc_)
  {
    for (size_t ivar = 0; ivar < nvar; ivar++) wkspc_.xu[ivar] = (fmap_m[ivar]*xu_in_[ivar]) + fmap_b[ivar];
    init_xu_vals(wkspc_.xu,wkspc_.xu_vals);
    for (int ivar = ndep; ivar >= 0; ivar--)
      for (int iord = bor; iord >= 0; iord--)
        wkspc_.xu_vals[ivar][iord] = eval_Li(wkspc_.xu_vals[ivar],iord);
  }

  double eval_Li(double *v_, int O_)
    {double L_acc = 0.0; for (int i = 0; i <= O_ ; i++) L_acc += poly_coeffs[O_][i]*v_[i]; return L_acc;}
  double eval_dnLi(int n_, double *v_, int O_)
    {double dnL_acc=0.0; for (int i = n_, ishift=0; i <= O_; i++, ishift++) dnL_acc += ((double)icoeff_mat[n_][ishift])*poly_coeffs[O_][i]*v_[i-n_]; return dnL_acc;}
  double get_dcoeff(int idim_, int O_) {return fmap_m[idim_];}

  inline void set_Legendre_coeffs()
  {
    poly_coeffs[0][0] = 1.0;
    if (bor>0)
    {
      poly_coeffs[1][0] = 0.0; poly_coeffs[1][1] = 1.0;
      for (int k = 2; k <= bor; k++)
      {
        double k_d=(double)k, k2_m1=(double)(2*k-1), k_m1=(double)(k-1);
        for (int j = 0; j <= k; j++) poly_coeffs[k][j] = (k2_m1*(get_coeff(k-1, j-1)) - k_m1*(get_coeff(k-2,j)))/k_d;
      }
    }
  }
  inline void set_Chebyshev1_coeffs()
  {
    poly_coeffs[0][0] = 1.0;
    if (bor>0)
    {
      poly_coeffs[1][0] = 0.0; poly_coeffs[1][1] = 1.0;
      for (int k = 2; k <= bor; k++)
      {
        poly_coeffs[k][0] = -1.0*poly_coeffs[k-2][0];
        for (int j = 1; j <= k; j++) poly_coeffs[k][j] = 2.0*get_coeff(k-1,j-1) - get_coeff(k-2,j);
      }
    }
  }
  inline void set_Chebyshev2_coeffs()
  {
    poly_coeffs[0][0] = 1.0;
    if (bor>0)
    {
      poly_coeffs[1][0] = 0.0; poly_coeffs[1][1] = 2.0;
      for (int k = 2; k <= bor; k++)
      {
        poly_coeffs[k][0] = -1.0*poly_coeffs[k-2][0];
        for (int j = 1; j <= k; j++) poly_coeffs[k][j] = 2.0*get_coeff(k-1,j-1) - get_coeff(k-2,j);
      }
    }
  }
  inline void set_Hermite_coeffs()
  {
    poly_coeffs[0][0] = 1.0;
    if (bor>0)
    {
      poly_coeffs[1][0] = 0.0; poly_coeffs[1][1] = 2.0;
      for (int k = 2; k <= bor; k++)
      {
        poly_coeffs[k][0] = -1.0*poly_coeffs[k-1][1];
        for (int j = 1; j <= k; j++) poly_coeffs[k][j] = 2.0*get_coeff(k-1,j-1) - ((double)(j+1))*get_coeff(k-1,j+1);
      }
    }
  }

  inline double eval_idim_forward_map(int idim_, double in_) {return fmap_m[idim_]*in_ + fmap_b[idim_];}
  inline double get_coeff(int k_, int j_) {return ((j_>=0)&&(j_<=k_))?(poly_coeffs[k_][j_]):(0.0);}
};

struct function_space_element: public ode_solspc_element
{
  function_space_element(function_space &fspc_): ode_solspc_element(fspc_.meta), fspc(fspc_) {}
  ~function_space_element() {}

  function_space &fspc;

  const int &bor = fspc.bor,
            &perm_len = fspc.perm_len,
            &ndof_full = fspc.ndof_full;

  bool  * const dof_hot_flags = fspc.dof_hot_flags,
        * const dof_tun_flags = fspc.dof_tun_flags,
        ** const dof_hot_flags_mat = fspc.dof_hot_flags_mat,
        ** const dof_tun_flags_mat = fspc.dof_tun_flags_mat;

  int ** const order_mat = fspc.order_mat,
      * const dof_hot_indices = fspc.dof_hot_indices,
      * const dof_tun_indices = fspc.dof_tun_indices,
      ** const dof_hot_indices_mat = fspc.dof_hot_indices_mat,
      ** const dof_tun_indices_mat = fspc.dof_tun_indices_mat,
      &ndof_hot = fspc.ndof_hot,
      &ndof_tun = fspc.ndof_tun;

  double  * const theta_full = fspc.theta_full,
          ** const theta_mat = fspc.theta_mat;

  protected:

    inline int thread_num() {return fspc.thread_num();}
    inline int get_nt() {return fspc.get_nt();}
};

struct coupling_term: public ode_solspc_element
{
  coupling_term(ode_solspc_meta &meta_, int cpl_): ode_solspc_element(meta_),
    cpl(cpl_), Pvec(new int[len_Pvec]), deriv_flags(new bool[nvar]), Lxuh_vec(new double[len_Lxuh_vec]) {}
  ~coupling_term() {delete [] deriv_flags; delete [] Pvec; delete [] Lxuh_vec;}

  const int cpl;
  const size_t  len_Pvec = 1 + (ndep*(cpl+1)),
                len_Lxuh_vec = nvar + 1;
  bool * const  deriv_flags,
                &xderiv_flag = deriv_flags[0];
  int * const Pvec,
      &Px = Pvec[0];
  double  coeff,
          * const Lxuh_vec,
          &Lxh = Lxuh_vec[0],
          &Luh = Lxuh_vec[nvar];
};

struct partial_chunk
{
  partial_chunk(int perm_len_, int eor_, int ndep_);
  ~partial_chunk();

  double  * const chunk,
          ** const Jac_mat,
          *** const C_x,
          ** const C_u,
          ** const C_x_ptrs,
          * const par_chunk;

  inline double ** Jac_xtheta_i_vdxu(int i_) {return C_x[i_];}
  inline double * grad_utheta_i_vdxu(int i_) {return C_u[i_];}
};

class function_space_basis: public function_space_element
{
  public:
    function_space_basis(function_space &fspc_);
    ~function_space_basis();

    partial_chunk partials;
    double  ** const Jac_mat = partials.Jac_mat,
            *** const C_x = partials.C_x,
            ** const C_u = partials.C_u;

  void v_eval(double *s_, double *v_);
  void fill_partial_chunk(double *s_);
  // void fill_partial_chunk(double *s_, double *v_);

  protected:
    const bool dxu_pows_allocate = eor > 2;

    const int eorm2 = eor-2,
              ord_len,
              ord_i_len;

    bool  * const dof_hot_flags_i;

    int * const ord_xu;

    double  * const theta_xu,
            * const s_j,
            * const v_j,
            ** const xu_vals,
            ** const dxu_vals;

    // conditionally defined data chunks for storing powers of high order derivatives
    double  * const dxu_pow_chunk,
            ** const dnxu_pow_ptrs,
            *** const dnxu_val;

    coupling_term ** const couples;

    void compute_x_coupling(coupling_term &c0_);
    void compute_u_coupling(int k_, coupling_term &ckm1_);
    void compute_u_coupling(int k_, coupling_term &ckm1_, int indep_source_);

    virtual void init_workspace(double *s_) = 0;
    virtual double stage_indep_var(int tord_) = 0;
    virtual double stage_dep_var(int i_L_) = 0;

    virtual void check_coupling_derivatives(coupling_term &ckm1_) = 0;
    virtual double comp_dx_dnLi(coupling_term &ckm1_, double &ncoeff_, double &dLx_) = 0;
    virtual double comp_dui_dnLi(coupling_term &ckm1_, int iu_, double &ncoeff_, double &dLui_, double &dLu_) = 0;

    inline void init_workspace_dnxu_vals()
    {
      for (int idxu = 0, idim=nvar; idxu < ndep; idxu++, idim++) // precompute powers of first derivative terms
      {
        dxu_vals[idxu][0] = s_j[idim]*s_j[idim];
        for (int i = 1; i < eor; i++) dxu_vals[idxu][i] = dxu_vals[idxu][i-1]*s_j[idim];
      }
      if (dxu_pows_allocate) // the diff eq. is higher than second order, we need additional precomputed components
      {
        for (int idep = 0; idep < ndep; idep++)
        {
          double ** dnxu_dep = dnxu_val[idep];
          for (int ider = 2, iskip = eorm2; ider < eor; ider++, iskip--)
          {
            double  *dnxu_der = dnxu_dep[ider-2],
                    dnxuk = s_j[1 + idep + (ider*ndep)];
            dnxu_der[0] = dnxuk*dnxuk;
            for (int ipow = 1; ipow < iskip; ipow++) dnxu_der[ipow] = dnxuk*dnxu_der[ipow-1];
          }
        }
      }
    }
    inline void stage_coupling(double coeff_)
    {
      for (int idim = nvar; idim < ndim; idim++) v_j[idim] = 0.0;
      couples[0]->coeff = coeff_;
    }
    inline double comp_dnxu_power_product(coupling_term &c_)
    {
      double acc = 1.0;
      for (int ider = 1, idim = nvar; ider <= c_.cpl; ider++)
        for (int idep = 0; idep < ndep; idep++, idim++)
          acc*=(get_dnxu_pow(idep,ider,c_.Pvec[idim]));
      return acc;
    }
    inline double comp_ddnxui_power_product(int jdnxu_,coupling_term &c_, double acc_start_)
    {
      double acc = acc_start_;
      for (int ider = 1, idnxu = nvar; ider <= c_.cpl; ider++)
        for (int idep = 0; idep < ndep; idep++, idnxu++)
          acc*=(get_dnxu_pow(idep,ider,(idnxu!=jdnxu_)?(c_.Pvec[idnxu]):(c_.Pvec[jdnxu_]-1)));
      return acc;
    }
    inline void stage_dx_couple(coupling_term &ckm1_, coupling_term &ck_, double ncoeff_, double dLx_)
    {
      update_vec_set_rest(0,ckm1_.len_Pvec,ckm1_.Px-1,ckm1_.Pvec,ck_.Pvec,ck_.len_Pvec); // update power vector, set remaining dims to zero (def. arg.)
      ck_.coeff = ncoeff_; // update coefficient
      update_vec(0,ckm1_.len_Lxuh_vec,dLx_,ckm1_.Lxuh_vec,ck_.Lxuh_vec); // update Ltx history vector
    }
    inline void stage_du_couple(int iu_, coupling_term &ckm1_, coupling_term &ck_, double ncoeff_, double dLui_, double dLu_)
    {
      update_vec_set_rest(iu_,ckm1_.len_Pvec,ckm1_.Pvec[iu_]-1,ckm1_.Pvec,ck_.Pvec,ck_.len_Pvec); // update power vector, set remaining dims to zero (def. arg.)
      ck_.Pvec[iu_+ndep]++; // increment power of derivative term
      ck_.coeff = ncoeff_; // update coefficient
      update_vec_set_rest(iu_,ckm1_.len_Lxuh_vec-1,dLui_,ckm1_.Lxuh_vec,ck_.Lxuh_vec,ckm1_.len_Lxuh_vec,dLu_); // update Ltx history vector, set last term to new product
    }
    inline void stage_ddnxu_couple(int jdnxu_, coupling_term &ckm1_, coupling_term &ck_, double ncoeff_)
    {
      update_vec_set_rest(jdnxu_,ckm1_.len_Pvec,ckm1_.Pvec[jdnxu_]-1,ckm1_.Pvec,ck_.Pvec,ck_.len_Pvec); // update power vector, set remaining dims to zero (def. arg.)
      ck_.Pvec[jdnxu_+ndep]++; // increment power of derivative term
      ck_.coeff = ncoeff_; // update coefficient
      copy_vec(ck_.len_Lxuh_vec, ckm1_.Lxuh_vec, ck_.Lxuh_vec); // copy Ltx history vector
    }
    inline void spawn_ddnxu_couple(int k_, coupling_term &ckm1_, coupling_term &ck_)
    {
      copy_vec_set_rest(ckm1_.len_Pvec,ckm1_.Pvec,ck_.Pvec,ck_.len_Pvec); // copy power vector, set remaining dims to zero (def. arg.)
      ck_.coeff = ckm1_.coeff; // copy coefficient
      copy_vec(ck_.len_Lxuh_vec, ckm1_.Lxuh_vec, ck_.Lxuh_vec); // copy Ltx history vector
    }

        //  dim         deriv     power
    inline double get_dnxu_pow(int idep_, int ider_, int Pdnxu_)
    {
      if (Pdnxu_<=1) return (Pdnxu_<=0)?(1.0):(s_j[1+idep_+(ider_*ndep)]);
      else return (ider_==1)?(dxu_vals[idep_][Pdnxu_-2]):(dnxu_val[idep_][ider_-2][Pdnxu_-2]);
    }

    inline void update_vec(int i_, int len_, int val_, int *vec_in_, int *vec_out_)
    {
      for (size_t j = 0; j < i_; j++) vec_out_[j] = vec_in_[j];
      vec_out_[i_] = val_;
      for (size_t j = i_+1; j < len_; j++) vec_out_[j] = vec_in_[j];
    }
    inline void update_vec(int i_, int len_, double val_, double *vec_in_, double *vec_out_)
    {
      for (size_t j = 0; j < i_; j++) vec_out_[j] = vec_in_[j];
      vec_out_[i_] = val_;
      for (size_t j = i_+1; j < len_; j++) vec_out_[j] = vec_in_[j];
    }
    inline void update_vec_set_rest(int i_, int len_short_, int val_, int *vec_in_, int *vec_out_, int len_long_, int val_rest_=0)
    {
      update_vec(i_,len_short_,val_,vec_in_,vec_out_);
      for (size_t j = len_short_; j < len_long_; j++) vec_out_[j] = val_rest_;
    }
    inline void update_vec_set_rest(int i_, int len_short_, double val_, double *vec_in_, double *vec_out_, int len_long_, double val_rest_=0.0)
    {
      update_vec(i_,len_short_,val_,vec_in_,vec_out_);
      for (size_t j = len_short_; j < len_long_; j++) vec_out_[j] = val_rest_;
    }
    inline void copy_vec(int len_, bool *in_, bool *out_)
      {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}
    inline void copy_vec(int len_, int *in_, int *out_)
      {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}
    inline void copy_vec(int len_, double *in_, double *out_)
      {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}
    inline void copy_vec_set_rest(int len_short_, int *vec_in_, int *vec_out_, int len_long_, int val_rest_=0)
    {
      copy_vec(len_short_,vec_in_,vec_out_);
      for (size_t j = len_short_; j < len_long_; j++) vec_out_[j] = val_rest_;
    }
    inline void copy_vec_set_rest(int len_short_, double *vec_in_, double *vec_out_, double len_long_, double val_rest_=0.0)
    {
      copy_vec(len_short_,vec_in_,vec_out_);
      for (size_t j = len_short_; j < len_long_; j++) vec_out_[j] = val_rest_;
    }

    inline int comp_v_k_start_index(int k_) {return 1 + (k_*ndep);}
    inline double update_product(double prod_old_, double old_, double new_) {return (prod_old_/old_)*new_;}

    inline bool is_zero(double in_) {return ((in_+0.0)==0.0);}
};

class power_basis: public function_space_basis
{
  public:
    power_basis(power_space &fspc_): function_space_basis(fspc_) {}
    ~power_basis() {}

  protected:

    virtual void init_workspace(double *s_);

    double stage_indep_var(int xord_) {return couples[0]->Lxh = comp_Li(xu_vals[0],couples[0]->Px = ord_xu[0] = xord_);}
    double stage_dep_var(int i_L_)
    {
      double Liu_acc = 1.0;
      for (int iu = 1; iu <= ndep; iu++)
        Liu_acc *= (couples[0]->Lxuh_vec[iu] = comp_Li(xu_vals[iu],couples[0]->Pvec[iu] = ord_xu[iu] = order_mat[i_L_][iu]));
      return couples[0]->Luh = Liu_acc;
    }
    void check_coupling_derivatives(coupling_term &ckm1_)
      {for (size_t i = 0; i <= ndep; i++) ckm1_.deriv_flags[i] = (ckm1_.Pvec[i] > 0);}
    double comp_dx_dnLi(coupling_term &ckm1_, double &ncoeff_, double &dLx_)
    {
      ncoeff_ = comp_dcoeff(0,ckm1_.Px)*ckm1_.coeff; // update scalar coefficient according to basis ruleset
      dLx_ = comp_dnLi(ord_xu[0]-ckm1_.Px+1,xu_vals[0],ord_xu[0]); // compute NEXT (pow_tx[0]-Pt+1) derivative of L wrt to t (0'th txval row, 0'th pow_tx row)
      return ncoeff_*dLx_*ckm1_.Luh;
    }
    double comp_dui_dnLi(coupling_term &ckm1_, int iu_, double &ncoeff_, double &dLui_, double &dLu_)
    {
      ncoeff_ = comp_dcoeff(iu_,ckm1_.Pvec[iu_])*ckm1_.coeff; // update scalar coefficient according to basis ruleset
      dLui_ = comp_dnLi(ord_xu[iu_]-ckm1_.Pvec[iu_]+1,xu_vals[iu_],ord_xu[iu_]); // compute NEXT (pow_tx[ix]-Pt+1) derivative of L wrt to i'th dep. var. (i'th txval row, i'th pow_tx row)
      dLu_ = (is_zero(ckm1_.Lxuh_vec[iu_]))?(0.0):(update_product(ckm1_.Luh,ckm1_.Lxuh_vec[iu_],dLui_)); // update product of dep. var. terms, ensuring we don't divide by zero
      return ncoeff_*dLu_*ckm1_.Lxh;
    }

    inline double comp_Li(double *v_, int i_) {return fspc.eval_Li(v_,i_);}
    inline double comp_dnLi(int n_, double *v_, int i_) {return fspc.eval_dnLi(n_,v_,i_);}
    inline double comp_dcoeff(int idim_, int i_) {return fspc.get_dcoeff(idim_,i_);}
};

class orthopolynomial_basis: public power_basis
{
  public:
    orthopolynomial_basis(orthopolynomial_space &fspc_): power_basis(fspc_),
    map_params(fspc_.map_params), poly_coeffs(fspc_.poly_coeffs),
    s_j_raw(new double[ndim]), s_j_hat(new double[ndim]) {}
    ~orthopolynomial_basis() {delete [] s_j_raw; delete [] s_j_hat;}

    double  ** const map_params,
            ** const poly_coeffs;

    void debugging_description();

    protected:

      double  * const s_j_raw,
              * const s_j_hat;

      void init_workspace(double *s_)
      {
        for (size_t idim = 0; idim < ndim; idim++) s_j_raw[idim] = s_j_hat[idim] = s_[idim];
        for (size_t ivar = 0; ivar < nvar; ivar++) s_j_hat[ivar] = (map_params[0][ivar]*s_j_hat[ivar]) + map_params[1][ivar];
        power_basis::init_workspace(s_j_hat);
      }
};


#endif
