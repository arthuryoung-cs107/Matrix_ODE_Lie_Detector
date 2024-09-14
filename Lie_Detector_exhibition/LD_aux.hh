#ifndef LD_AUX_HH
#define LD_AUX_HH

#include "LD_util.hh"

#include <cstdio>
#include <cmath>

#include <cfloat>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

template <class T, class Y> T ** make_evaluation_bases(Y &fspace_, int nbases_=-1)
{
  const int nbases = (nbases_<0)?(LD_threads::numthreads()):(nbases_);
  T ** bases_out = new T*[nbases];
  for (size_t i = 0; i < nbases; i++) bases_out[i] = new T(fspace_);
  return bases_out;
}
template <class T> void free_evaluation_bases(T **bases_, int nbases_=-1)
{
  const int nbases = (nbases_<0)?(LD_threads::numthreads()):(nbases_);
  for (size_t i = 0; i < nbases; i++) delete bases_[i];
  delete [] bases_;
}

struct LD_linalg
{
  LD_linalg() {}
  ~LD_linalg() {}

  template <typename T> static void abs_vec(T *vec_, int len_)
    {for (size_t i = 0; i < len_; i++) if (vec_[i] < 0) vec_[i] = -vec_[i];}

  static void normalize_vec_l2(double *vec_, int len_)
  {
    double wrk = 0.0;
    for (size_t i = 0; i < len_; i++) wrk += vec_[i]*vec_[i];
    wrk = sqrt(wrk);
    for (size_t i = 0; i < len_; i++) vec_[i] /= wrk;
  }

  static double norm_l2(double *x_, int n_)
  {
    double r2_acc = 0.0;
    for (size_t i = 0; i < n_; i++) r2_acc += x_[i]*x_[i];
    return sqrt(r2_acc);
  }

  template <typename T> static T min_T(T a_, T b_) {return (a_<b_)?(a_):b_;}
  template <typename T> static T max_T(T a_, T b_) {return (a_>b_)?(a_):b_;}

  template <typename T> static T min_T_3way(T a_, T b_, T c_) {return min_T<T>(min_T<T>(a_,b_),c_);}
  template <typename T> static T max_T_3way(T a_, T b_, T c_) {return max_T<T>(max_T<T>(a_,b_),c_);}

  template <typename T> static void comp_Tvec_stats(double stats_[], T *vec_, int len_)
    {T wrk[len_]; LD_linalg::comp_Tvec_stats<T>(stats_,vec_,wrk,len_);}
  template <typename T> static void comp_Tvec_stats(double stats_[], T *vec_, T *wrk_, int len_)
  {
    T acc = 0;
    for (size_t i = 0; i < len_; i++) acc += (wrk_[i] = vec_[i]);
    LD_linalg::sort_vec_inc<T>(wrk_,len_);
    stats_[0] = (double) wrk_[0]; // min
    stats_[1] = (len_%2)?((double) wrk_[len_/2]):(0.5*((double) (wrk_[(len_/2) - 1] + wrk_[len_/2]) )); // med
    stats_[2] = ((double) acc)/((double) len_);
    stats_[3] = wrk_[len_-1];
  }

  template <typename T> static void comp_Tvec_maM_stats(double stats_[], T *vec_, int len_)
  {
    T vec_acc = vec_[0], max_vec = vec_[0], min_vec = vec_[0], vec_i;
    for (size_t i = 1; i < len_; i++)
    {
      vec_acc += (vec_i = vec_[i]);
      if (min_vec > vec_i) min_vec = vec_i;
      if (max_vec < vec_i) max_vec = vec_i;
    }
    stats_[0] = (double)(min_vec);
    stats_[1] = ((double)(vec_acc))/((double)(len_));
    stats_[2] = (double)(max_vec);
  }
  template <typename T> static T min_val(T *vec_, int len_, int &ind_, int ind_offset_=0)
  {
    T val_out = vec_[ind_=0];
    for (size_t i = 1; i < len_; i++)
      if (val_out > vec_[i])
        val_out = vec_[ind_ = i];
    ind_ += ind_offset_;
    return val_out;
  }
  template <typename T> static T max_val(T *vec_, int len_, int &ind_, int ind_offset_=0)
  {
    T val_out = vec_[ind_=0];
    for (size_t i = 1; i < len_; i++)
      if (val_out < vec_[i])
        val_out = vec_[ind_ = i];
    ind_ += ind_offset_;
    return val_out;
  }
  template <typename T> static T min_val(T *vec_, int len_)
  {
    T val_out = vec_[0];
    for (size_t i = 1; i < len_; i++) if (val_out > vec_[i]) val_out = vec_[i];
    return val_out;
  }
  template <typename T> static T max_val(T *vec_, int len_)
  {
    T val_out = vec_[0];
    for (size_t i = 1; i < len_; i++) if (val_out < vec_[i]) val_out = vec_[i];
    return val_out;
  }

  template <typename T> static T sum_vec(T *vec_,int len_)
  {
    T acc = 0;
    for (size_t i = 0; i < len_; i++) acc += vec_[i];
    return acc;
  }

  template <typename T> static void scale_vec(T *vec_,int len_,T val_)
    {for (size_t i = 0; i < len_; i++) vec_[i] *= val_;}

  template <typename T> static void copy_vec(T *out_, T *in_, int len_)
    {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}

  template <typename T> static void fill_vec(T *out_, int len_, T val_)
    {for (size_t i = 0; i < len_; i++) out_[i] = val_;}
  template <typename T> static void fill_vec_012(T *vec_, int len_, T offset_=0)
    {for (int i = 0; i < len_; i++) vec_[i] = ((T)(i)) + offset_;}

  template <typename T> static void sort_vec_inc(T *vec_, int len_)
    {T m_; int im_; LD_linalg::sort_vec_inc<T>(vec_,len_,m_,im_);}
  template <typename T> static void sort_vec_inc(T *vec_, int len_, T &m_, int &im_)
  {
    if (len_>1)
    {
      m_ = vec_[im_=0];
      for (size_t i = 1; i < len_; i++) if (m_>vec_[i]) m_ = vec_[im_ = i];
      vec_[im_] = vec_[0]; vec_[0] = m_;
      if (len_>2) LD_linalg::sort_vec_inc<T>(vec_+1,len_-1,m_,im_);
    }
  }
  template <typename T> static void sort_vec_dec(T *vec_, int len_)
    {T m_; int im_; LD_linalg::sort_vec_dec<T>(vec_,len_,m_,im_);}
  template <typename T> static void sort_vec_dec(T *vec_, int len_, T &m_, int &im_)
  {
    if (len_>1)
    {
      m_ = vec_[im_=0];
      for (size_t i = 1; i < len_; i++) if (m_<vec_[i]) m_ = vec_[im_ = i];
      vec_[im_] = vec_[0]; vec_[0] = m_;
      if (len_>2) LD_linalg::sort_vec_dec<T>(vec_+1,len_-1,m_,im_);
    }
  }

  template <typename T> static void sort_vec_inc(T *vec_, int len_, int &ind_m_)
  {
    if (len_>1)
    {
      T vm = LD_linalg::min_val<T>(vec_,len_,ind_m_);
      vec_[ind_m_] = vec_[0]; vec_[0] = vm;
    }
    if (len_>2) LD_linalg::sort_vec_inc<T>(vec_+1,len_-1,ind_m_);
  }
  template <typename T> static void sort_vec_dec(T *vec_, int len_, int &ind_m_)
  {
    if (len_>1)
    {
      T vm = LD_linalg::max_val<T>(vec_,len_,ind_m_);
      vec_[ind_m_] = vec_[0]; vec_[0] = vm;
    }
    if (len_>2) LD_linalg::sort_vec_dec<T>(vec_+1,len_-1,ind_m_);
  }

  static void A__B_CT(double **Amat_,double **Bmat_,double **Cmat_,int m_,int l_,int n_,bool init0_=false)
  {
    if (init0_) LD_linalg::fill_vec<double>(Amat_[0],m_*n_,0.0);
    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        for (size_t ii = 0; ii < l_; ii++)
          Amat_[i][j] += Bmat_[i][ii]*Cmat_[j][ii];
  }
  static void A_x_b(double **A_, double *x_, double *b_, int m_, int n_)
  {
    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        b_[i] += A_[i][j]*x_[j];
  }

  static void AT_x_b(double ** A_, double *x_, double *b_, int m_, int n_)
  {
    for (size_t i = 0; i < m_; i++)
      for (size_t j = 0; j < n_; j++)
        b_[j] += A_[i][j]*x_[i];
  }

  static void set_Asym_ATA(double **Asym_, double **A_, int m_, int n_)
  {
    for (size_t i = 0; i < n_; i++)
      for (size_t j = 0; j < (n_-i); j++)
      {
        Asym_[i][j] = 0.0;
        for (size_t k = 0; k < m_; k++) Asym_[i][j] += A_[k][i]*A_[k][j+i]; // not very cache efficient
      }
  }

  // computes vT A v, where A is symmetric and upper stored. w_ is work vec, equal to Ax on output
  static double xT_Asym_x(double **A_, double *x_, double *w_, int n_)
  {
    double inner_acc = 0.0;
    for (size_t i = 0; i < n_; i++)
    {
      w_[i] += A_[i][0]*x_[i]; // diagonal element
      for (size_t j = 1; j < (n_ - i); j++)
        {w_[i] += A_[i][j]*x_[i+j]; w_[i+j] += A_[i][j]*x_[i];}
      inner_acc += (w_[i])*(x_[i]);
    }
    return inner_acc;
  }

  // Computes matrix product, assumes b_ is initialized, added onto by c_*Ax
  static void Asym_x_b(double ** A_, double *x_, double *b_, int n_, double c_=1.0)
  {
    for (size_t i = 0; i < n_; i++)
    {
      b_[i] += c_*A_[i][0]*x_[i]; // diagonal element
      for (size_t j = 1; j < (n_ - i); j++)
        {b_[i] += c_*A_[i][j]*x_[i+j]; b_[i+j] += c_*A_[i][j]*x_[i];}
    }
  }

  static double xTy(double *x_, double *y_, int n_)
  {
    double acc = 0.0;
    for (size_t i = 0; i < n_; i++) acc += x_[i]*y_[i];
    return acc;
  }

  static void fill_x_012(double *x_, int n_, double offset_=0.0)
    {for (size_t i = 0; i < n_; i++) x_[i] = ((double)(i)) + offset_;}
  static void fill_x_012(int *x_, int n_, int offset_=0)
    {for (size_t i = 0; i < n_; i++) x_[i] = i + offset_;}

  static void copy_x(double *x_in_, double *x_out_, int n_)
    {for (size_t i = 0; i < n_; i++) x_out_[i] = x_in_[i];}
  static void copy_x(int *x_in_, int *x_out_, int n_)
    {for (size_t i = 0; i < n_; i++) x_out_[i] = x_in_[i];}
  static void copy_x(bool *x_in_, bool *x_out_, int n_)
    {for (size_t i = 0; i < n_; i++) x_out_[i] = x_in_[i];}

  static double min_val_mag(double *x_, int n_, int &i_min_, int i_start_=0)
  {
    if (i_start_<n_)
    {
      double  out = fabs(x_[i_min_=i_start_]),
              x_i_mag;
      for (size_t i = i_start_+1; i < n_; i++)
        if ((x_i_mag=fabs(x_[i])) < out)
          {out = x_i_mag; i_min_ = i;}
      return out;
    }
    else return fabs(x_[i_min_=(n_-1)]); // return last element
  }

  static int max_mag_elements(double *x_, double *v_, int *inds_, int n_, int k_)
  {
    for (int i = 0; i < k_; i++)
    {
      v_[i] = x_[i]; // initializing the magnitudes to be the first few
      inds_[i] = i; // setting first few indices
    }
    int i_it;
    double  min_large = min_val_mag(v_, k_, i_it); // determine the smallest magnitude on this list
    for (size_t i = k_; i < n_; i++) // go through the rest of them
      if (fabs(x_[i]) > min_large) // if we find somethng bigger than the current smallest element on the top list
      {
        v_[i_it] = x_[i]; // kick off the smallest 'large' element
        inds_[i_it] = i;
        min_large = min_val_mag(v_,k_,i_it); // determine the new smallest 'large' element
      }
    return i_it;
  }

  static double eps(double x_=1.0) {return nextafter(x_,DBL_MAX)-x_;}

  static void print_x(const char vec_name_[], double *x_, int n_)
  {
    printf("%s (%d x 1)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%.2e\n", x_[i]);
  }
  static void print_x(const char vec_name_[], int *x_, int n_)
  {
    printf("%s (%d x 1)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%d\n", x_[i]);
  }
  static void print_x(const char vec_name_[], bool *x_, int n_)
  {
    printf("%s (%d x 1)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%d\n", (int) x_[i]);
  }
  static void print_xT(const char vec_name_[], double *x_, int n_)
  {
    printf("%s (1 x %d)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%.2e ", x_[i]);
    printf("\n");
  }
  static void print_xT(const char vec_name_[], int *x_, int n_)
  {
    printf("%s (1 x %d)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%d ", x_[i]);
    printf("\n");
  }
  static void print_xT(const char vec_name_[], bool *x_, int n_)
  {
    printf("%s (1 x %d)\n", vec_name_, n_);
    for (size_t i = 0; i < n_; i++) printf("%d ", (int) x_[i]);
    printf("\n");
  }
  static void print_A(const char mat_name_[], double **A_, int m_, int n_)
  {
    printf("%s (%d x %d)\n", mat_name_, m_, n_);
    for (size_t i = 0; i < m_; i++)
      {for (size_t j = 0; j < n_; j++) printf("%.2e ", A_[i][j]); printf("\n");}
  }
  static void print_A(const char mat_name_[], int **A_, int m_, int n_)
  {
    printf("%s (%d x %d)\n", mat_name_, m_, n_);
    for (size_t i = 0; i < m_; i++)
      {for (size_t j = 0; j < n_; j++) printf("%d ", A_[i][j]); printf("\n");}
  }
  static void print_Asym(const char mat_name_[], double **Asym_, int n_)
  {
    printf("%s (%d x %d)\n", mat_name_, n_, n_);
    for (size_t i = 0, j_cap = n_; i < n_; i++, j_cap--)
      {for (size_t j = 0; j < j_cap; j++) printf("%.2e ", Asym_[i][j]); printf("\n");}
  }
};

struct LD_gsl
{
  LD_gsl() {}
  ~LD_gsl() {}

  static void load_gsl_vec(gsl_vector *gsl_vec_, double *vec_, int offset_=0)
    {for (int i = offset_, ii=0; i < gsl_vec_->size; i++, ii++) gsl_vector_set(gsl_vec_,i,vec_[ii]);}
  static void unpack_gsl_vec(double *vec_, gsl_vector *gsl_vec_, int offset_=0)
    {for (int i = offset_, ii=0; i < gsl_vec_->size; i++, ii++) vec_[ii] = gsl_vector_get(gsl_vec_,i);}

  static void load_gsl_mat_row_i(gsl_matrix *gsl_mat_, double *vec_, int i_, int offset_=0)
    {for (int j = offset_, jj = 0; j < gsl_mat_->size2; j++, jj++) gsl_matrix_set(gsl_mat_,i_,j,vec_[jj]);}
  static void load_gsl_mat_col_j(gsl_matrix *gsl_mat_, double *vec_, int j_, int offset_=0)
    {for (int i = offset_, ii = 0; i < gsl_mat_->size1; i++, ii++) gsl_matrix_set(gsl_mat_,i,j_,vec_[ii]);}

  static void unpack_gsl_mat_row_i(double *vec_, gsl_matrix *gsl_mat_, int i_, int offset_=0)
    {for (int j = offset_, jj = 0; j < gsl_mat_->size2; j++, jj++) vec_[jj] = gsl_matrix_get(gsl_mat_,i_,j);}
  static void unpack_gsl_mat_col_j(double *vec_, gsl_matrix *gsl_mat_, int j_, int offset_=0)
    {for (int i = offset_, ii = 0; i < gsl_mat_->size1; i++, ii++) vec_[ii] = gsl_matrix_get(gsl_mat_,i,j_);}

  static void load_gsl_mat(gsl_matrix *gsl_mat_, double **mat_)
  {
    for (int i = 0; i < gsl_mat_->size1; i++)
      for (int j = 0; j < gsl_mat_->size2; j++)
        gsl_matrix_set(gsl_mat_,i,j,mat_[i][j]);
  }
  static void load_gsl_matT(gsl_matrix *gsl_mat_, double **mat_)
  {
    for (int i = 0; i < gsl_mat_->size1; i++)
      for (int j = 0; j < gsl_mat_->size2; j++)
        gsl_matrix_set(gsl_mat_,j,i,mat_[i][j]);
  }
  static void unpack_gsl_mat(double **mat_, gsl_matrix *gsl_mat_)
  {
    for (int i = 0; i < gsl_mat_->size1; i++)
      for (int j = 0; j < gsl_mat_->size2; j++)
        mat_[i][j] = gsl_matrix_get(gsl_mat_,i,j);
  }
  static void unpack_gsl_matT(double **mat_, gsl_matrix *gsl_mat_)
  {
    for (int i = 0; i < gsl_mat_->size1; i++)
      for (int j = 0; j < gsl_mat_->size2; j++)
        mat_[i][j] = gsl_matrix_get(gsl_mat_,j,i);
  }

  static void print_A_gsl(const char mat_name_[], gsl_matrix *mat_)
  {
    printf("%s (%d x %d)\n", mat_name_, mat_->size1, mat_->size2);
    for (size_t i = 0; i < mat_->size1; i++)
      {for (size_t j = 0; j < mat_->size2; j++) printf("%.2e ", gsl_matrix_get(mat_,i,j)); printf("\n");}
  }

};

struct LD_rng
{
  LD_rng(int seed_=1): seed(seed_), gsl_gen(gsl_rng_alloc(gsl_rng_taus2))
    {gsl_rng_set(gsl_gen,seed); current_seed=seed;}
  ~LD_rng() {gsl_rng_free(gsl_gen);}

  const int seed;
  gsl_rng * const gsl_gen;

  int current_seed,
      old_seed=0;

  inline void set_seed(int seed_)
    {old_seed=current_seed; gsl_rng_set(gsl_gen,seed_); current_seed=seed_;}

  inline double rand_uni(double low_=0.0,double high_=1.0)
    {return (gsl_rng_uniform(gsl_gen)-0.5)*(high_-low_) + 0.5*(high_+low_);}
  inline double rand_gau(double mu_=0.0, double sigma_=1.0)
    {return mu_+gsl_ran_gaussian(gsl_gen,sigma_);}

  inline void set_vec_rand_uni(double *v_, int n_, double low_=0.0, double high_=1.0)
    {for (size_t i = 0; i < n_; i++) v_[i] = rand_uni(low_,high_);}
  inline void set_vec_rand_gau(double *v_, int n_, double mu_=0.0, double sigma_=1.0)
    {for (size_t i = 0; i < n_; i++) v_[i] = rand_gau(mu_,sigma_);}
};

struct LD_svd
{
  LD_svd(size_t Mmax_, size_t Nmax_): Mmax(LD_linalg::max_T<int>(Mmax_,2)), Nmax(LD_linalg::min_T<int>(Mmax,Nmax_)),
    Muse(Mmax), Nuse(Nmax),
    Umat(Tmatrix<double>(Mmax,Nmax)), Vmat(Tmatrix<double>(Nmax,Nmax)), svec(new double[Nmax]), wvec(new double[Nmax]) {}
  ~LD_svd()
  {
    delete [] wvec; delete [] svec;
    free_Tmatrix<double>(Vmat); free_Tmatrix<double>(Umat);
  }

  const int Mmax,
            Nmax;
  int Muse,
      Nuse,
      Iuse,
      Juse;
  double  ** const Umat,
          ** const Vmat,
          * const svec,
          * const wvec;
  double &s0 = svec[0];

  inline double norm_oprt() {return svec[0];}
  inline double norm_nuke()
  {
    double acc = 0.0;
    for (size_t i = 0; i < Nuse; i++) acc += svec[i];
    return acc;
  }

  inline void print_result(const char name_[])
  {
    printf("(LD_svd::print_result) %s (%d x %d) SVD (rank = %d):\n", name_, Muse, Nuse, rank());
    LD_linalg::print_xT("  s",svec,Nuse);
  }

  inline void load_and_decompose_U(double **Amat_, int Muse_=0, int Nuse_=0)
  {
    set_use_dims(Muse_,Nuse_);
    for (size_t i = 0; i < Muse; i++) for (size_t j = 0; j < Nuse; j++) Umat[i][j] = Amat_[i][j];
    gsl_matrix_view U_gsl = gsl_matrix_view_array_with_tda(Umat[0],Muse,Nuse,Nmax),
                    V_gsl = gsl_matrix_view_array_with_tda(Vmat[0],Nuse,Nuse,Nmax);
    gsl_vector_view s_gsl = gsl_vector_view_array(svec,Nuse),
                    w_gsl = gsl_vector_view_array(wvec,Nuse);
    int status_decomp = gsl_linalg_SV_decomp(&(U_gsl.matrix),&(V_gsl.matrix),&(s_gsl.vector),&(w_gsl.vector));
  }

  inline void decompose_U(int Muse_=0, int Nuse_=0)
  {
    set_use_dims(Muse_,Nuse_);
    gsl_matrix_view U_gsl = gsl_matrix_view_array_with_tda(Umat[0],Muse,Nuse,Nmax),
                    V_gsl = gsl_matrix_view_array_with_tda(Vmat[0],Nuse,Nuse,Nmax);
    gsl_vector_view s_gsl = gsl_vector_view_array(svec,Nuse),
                    w_gsl = gsl_vector_view_array(wvec,Nuse);
    int status_decomp = gsl_linalg_SV_decomp(&(U_gsl.matrix),&(V_gsl.matrix),&(s_gsl.vector),&(w_gsl.vector));
  }

  inline void unpack_svec_VTmat(double *svec_, double **VTmat_)
  {
    for (size_t i = 0; i < Nuse; i++) svec_[i] = svec[i];
    for (size_t i = 0; i < Nuse; i++)
      for (size_t j = 0; j < Nuse; j++)
        VTmat_[i][j] = Vmat[j][i];
  }
  inline int unpack_rank_svec_VTmat(double *svec_, double **VTmat_)
  {
    unpack_svec_VTmat(svec_,VTmat_);
    return comp_rank(((double) Muse)*(LD_linalg::eps(svec_[0])));
  }
  inline int rank(int iscl_=0) {return comp_rank( ((double) ((iscl_)?(iscl_):(Muse)))*(LD_linalg::eps(svec[0])) );}
  inline int comp_rank(double tol_eps_)
  {
    int rank_out = 0;
    for (size_t i = 0; i < Nuse; i++)
      if (svec[i]>tol_eps_) rank_out++;
      else break;
    return rank_out;
  }

  static int rank(double *s_, int len_, int iscl_)
  {
    const double tol_eps = ((double) iscl_) * LD_linalg::eps(s_[0]);
    int rank_out = 0;
    for (size_t i = 0; i < len_; i++)
      if (s_[i]>tol_eps) rank_out++;
      else break;
    return rank_out;
  }

  inline void init_use_dims() { Muse = Mmax; Nuse = Nmax; }
  inline void set_use_dims(int Muse_, int Nuse_) {if (Muse_) Muse = Muse_; if (Nuse_) Nuse = Nuse_; verify_use_dims();}

  private:

    inline void verify_use_dims()
    {
      // verify adequate space
      if (!( (0<Muse)&&(Muse<Mmax) )) Muse = Mmax; if (!( (0<Nuse)&&(Nuse<Nmax) )) Nuse = Nmax;

      // verify sub matrix
      if (Muse>Mmax) Muse = Mmax; if (Nuse>Nmax) Nuse = Nmax;

      // verify thin matrix
      if (Nuse>Muse) Nuse = Muse;

      if (Nuse==Muse)
        printf("(LD_svd::verify_use_dims) WARNING - intialized as eigendecomposition (square matrix, %d x %d)\n",Muse,Nuse);
    }

};

struct k_medoids_package
{
  k_medoids_package(double **dsym_,int npts_): dsym(dsym_), npts(npts_),
    i_meds(new int[2*npts_]), i_nmeds(i_meds+npts_) {}
  k_medoids_package(int kclst_,double **dsym_,int npts_): k_medoids_package(dsym_,npts_)
    {comp_k_medoids(kclst_);}
  ~k_medoids_package() {delete [] i_meds;}

  // input data
  const int npts;
  double  ** const dsym;

  int * const i_meds,
      * const i_nmeds;
  double  dclst,
          silh_coeff;

  inline void comp_sort_cluster_distances(double dpairs_[], int ipairs_[], int ijpairs_[],int imems_[],int nmem_)
  {
    const int nmem_m1 = nmem_-1;
    int inn0,iinn00,iinn01;
    double  dij_min = DBL_MAX,
            dij_try;
    for (size_t i = 0, iprs = 0, i_ij = 0; i < nmem_m1; i++)
      for (size_t j = i+1; j < nmem_; j++, iprs++, i_ij+=2)
        if (dij_min > (dpairs_[ipairs_[iprs] = iprs] = (dij_try = dij(ijpairs_[i_ij] = imems_[i],ijpairs_[i_ij+1] = imems_[j]))))
          {dij_min = dij_try; inn0=iprs; iinn00 = imems_[i]; iinn01 = imems_[j];}
    dpairs_[inn0] = dpairs_[0]; ipairs_[inn0] = ipairs_[0]; ijpairs_[2*inn0] = ijpairs_[0]; ijpairs_[(2*inn0)+1] = ijpairs_[1];
    dpairs_[0] = dij_min; ipairs_[0] = inn0; ijpairs_[0] = iinn00; ijpairs_[1] = iinn01;
    if (nmem_>2) sort_cluster_distances(dpairs_+1,ipairs_+1,ijpairs_+2,(((nmem_)*(nmem_-1))/2)-1);
  }
  inline double get_nearest_neighbors(int inn_[],int imems_[],int nmem_)
  {
    const int nmem_m1 = nmem_-1;
    double  dij_min = DBL_MAX,
            dij_try;
    for (size_t i = 0; i < nmem_m1; i++)
      for (size_t j = i+1; j < nmem_; j++)
        if (dij_min > (dij_try = dij(imems_[i],imems_[j])))
          {dij_min = dij_try; inn_[0] = imems_[i]; inn_[1] = imems_[j];}
    return dij_min;
  }

  void assign_clusters(int nmem_v_[],int imem_v_[],int *imem_m_[],int kclst_,int inmed_[],bool fnmed_[],bool verbose_=false)
  {
    populate_nmed(inmed_,fnmed_,kclst_);
    const int n_nmed = npts-kclst_;
    bool fmem_clst[kclst_][npts];

    for (size_t k = 0; k < kclst_; k++) // assign cluster centres to selves
    {
      for (size_t kk = 0; kk < k; kk++) fmem_clst[kk][i_meds[k]] = 0;
      fmem_clst[k][i_meds[k]] = (bool)(nmem_v_[k] = 1);
      for (size_t kk = k+1; kk < kclst_; kk++) fmem_clst[kk][i_meds[k]] = 0;
    }
    for (size_t inm = 0; inm < n_nmed; inm++) // assign non medoid members
    {
      const int kclst_inm = assign_ipt_kcluster(inmed_[inm],kclst_);
      for (size_t k = 0; k < kclst_inm; k++) fmem_clst[k][inmed_[inm]] = 0;
      nmem_v_[kclst_inm] += (int)(fmem_clst[kclst_inm][inmed_[inm]] = 1);
      for (size_t k = kclst_inm+1; k < kclst_; k++) fmem_clst[k][inmed_[inm]] = 0;
    }
    for (size_t k = 0, iv = 0; k < kclst_; k++)
    {
      imem_m_[k] = imem_v_ + iv;
      for (size_t i = 0; i < npts; i++) if (fmem_clst[k][i]) imem_v_[iv++] = i;
    }
    if (verbose_)
    {
      printf("(k_medoids_package::assign_clusters) k = %d cluster assignments:\n", kclst_);
      for (size_t k = 0; k < kclst_; k++)
      {
        printf("  (k=%d) imed=%d, nmem=%d - ", k,i_meds[k],nmem_v_[k]);
        for (size_t i = 0; i < nmem_v_[k]; i++) printf("%d ", imem_m_[k][i]);
        printf("\n");
      }
    }
  }
  inline void assign_clusters(int nmem_v_[],int imem_v_[],int *imem_m_[],int kclst_,bool verbose_=false)
    {bool  f_nmeds[npts]; assign_clusters(nmem_v_,imem_v_,imem_m_,kclst_,i_nmeds,f_nmeds,verbose_);}

  int comp_kSC_medoids(int klow_=2,bool verbose_=true)
  {
    bool f_nmeds[npts];
    int k_SC;
    double SC_max = comp_k_medoids(k_SC=find_kSC(klow_,f_nmeds),f_nmeds);
    if (verbose_)
      printf("(k_medoids_package::comp_kSC_medoids) optimal k=%d clusters identified for %d points (SC= %e)\n",
        k_SC,npts,SC_max);
    return k_SC;
  }
  inline int find_kSC(int klow_=2) {bool f_nmeds[npts]; return find_kSC(klow_,f_nmeds);}
  int find_kSC(int klow_,bool fnmed_[])
  {
    int k_SC = klow_;
    double  SC_max = comp_k_medoids(k_SC,fnmed_),
            SC_try;
    while ((k_SC < npts)&&((SC_try=comp_k_medoids(++k_SC,fnmed_))>SC_max)) SC_max=SC_try;
    silh_coeff = SC_max;
    return k_SC-1;
  }

  void comp_SC_krange_medoids(double SCvec_[],int klow_,int khigh_)
  {
    bool f_nmeds[npts];
    for (size_t k = klow_, ik=0; k <= khigh_; k++, ik++)
      SCvec_[ik] = comp_k_medoids(k,f_nmeds);
  }
  int comp_kSC_krange_medoids(int klow_,int khigh_,bool verbose_=true)
  {
    bool f_nmeds[npts];
    int k_SC;
    double SC_max = comp_k_medoids(k_SC=find_kSC_krange(klow_,khigh_,f_nmeds),f_nmeds);
    if (verbose_)
      printf("(k_medoids_package::comp_kSC_medoids) optimal k=%d clusters identified of %d points (SC= %e, considered %d <= k <= %d)\n",
        k_SC,npts,SC_max,klow_,khigh_);
    return k_SC;
  }

  inline int find_kSC_krange(int klow_,int khigh_)
    {bool f_nmeds[npts]; return find_kSC_krange(klow_,khigh_,f_nmeds);}
  int find_kSC_krange(int klow_,int khigh_,bool f_nmeds_[])
  {
    int k_SC;
    double  SC_max = 0.0, SC_try;
    for (size_t k = klow_; k <= khigh_; k++)
      if (SC_max < (SC_try = comp_k_medoids(k,f_nmeds_)))
        {SC_max = SC_try; k_SC = k;}
    silh_coeff = SC_max;
    return k_SC;
  }

  double comp_sc_k(int kclst_,int inmed_[],bool fnmed_[])
  {
    const int n_nmed = npts-kclst_;
    bool fmem_clst[kclst_][npts];
    populate_nmed(inmed_,fnmed_,kclst_);
    assign_cluster_members(fmem_clst[0],kclst_,inmed_);

    double s_acc = 0.0;
    for (size_t k = 0; k < kclst_; k++)
    {
      for (size_t ipt = 0; ipt < i_meds[k]; ipt++)
        if (fmem_clst[k][ipt]) s_acc += (1.0 - (dij(ipt,i_meds[k])/get_min_dclst_i_kskp(ipt,kclst_,k)));
      for (size_t ipt = i_meds[k]+1; ipt < npts; ipt++)
        if (fmem_clst[k][ipt]) s_acc += (1.0 - (dij(ipt,i_meds[k])/get_min_dclst_i_kskp(ipt,kclst_,k)));
    }
    return s_acc/((double) n_nmed);
  }

  inline double comp_k_medoids(int kclst_) {bool  f_nmeds[npts]; return comp_k_medoids(kclst_,f_nmeds);}
  double comp_k_medoids(int kclst_, bool f_nmeds_[])
  {
    double  dclst_greedy = comp_greedy_k_medoids(kclst_,i_nmeds,f_nmeds_); // perform greedy k medoids to start
    dclst = comp_k_medoids(kclst_,i_nmeds,f_nmeds_,dclst_greedy);
    return silh_coeff = comp_sc_k(kclst_,i_nmeds,f_nmeds_);
  }
  double comp_k_medoids(int kclst_,int inmed_[],bool fnmed_[],double dclst_greedy_=0.0)
  {
    const int n_nmed = npts-kclst_;
    double dclst_out = (dclst_greedy_)?(dclst_greedy_):(comp_greedy_k_medoids(kclst_,inmed_,fnmed_));
    while (true)
    {
      populate_nmed(inmed_,fnmed_,kclst_);
      int mswap_best = -1,
          oswap_best;
      double  dswap_best = dclst_out,
              dswap_try;
      for (size_t k = 0; k < kclst_; k++)
        for (size_t ik = 0; ik < n_nmed; ik++)
          if (dswap_best > (dswap_try = comp_dclst_itry_kskp(inmed_,kclst_,inmed_[ik],k)))
          {
            mswap_best = k;
            oswap_best = inmed_[ik];
            dswap_best = dswap_try;
          }

      if (mswap_best == -1) break;
      else
      {
        i_meds[mswap_best] = oswap_best;
        dclst_out = dswap_best;
      }
    }
    return dclst_out;
  }

  inline void set_one_medoid() {i_meds[0] = imin_net_d();}
  double comp_greedy_k_medoids(int kclst_,int inmed_[],bool fnmed_[])
  {
    i_meds[0] = imin_net_d();
    double dclst_out;
    for (size_t km1 = 1, n_nmed = npts-1; km1 < kclst_; km1++, n_nmed--)
    {
      populate_nmed(inmed_,fnmed_,km1);
      int &ik = i_meds[km1] = inmed_[0];
      double  dik = dclst_out = comp_itry_dclst(inmed_,km1,i_meds[km1]);
      for (size_t iik = 1; iik < n_nmed; iik++)
        if (dclst_out > (dik = comp_itry_dclst(inmed_,km1,inmed_[iik])))
          {ik = inmed_[iik]; dclst_out = dik;}
    }
    return dclst_out;
  }

  private:

    inline void sort_cluster_distances(double dpairs_[],int ipairs_[],int ijpairs_[], int npairs_)
    {
      if (npairs_ >= 2)
      {
        double d_m = dpairs_[0];
        int i_m = 0,
            inn_m = ipairs_[0],
            iinn_m0 = ijpairs_[0],
            iinn_m1 = ijpairs_[1];
        for (size_t i = 1, i_ij = 2; i < npairs_; i++, i_ij+=2)
          if (d_m > dpairs_[i]) {d_m = dpairs_[i]; i_m = i; inn_m = ipairs_[i]; iinn_m0 = ijpairs_[i_ij]; iinn_m1 = ijpairs_[i_ij+1];}
        dpairs_[i_m] = dpairs_[0]; ipairs_[i_m] = ipairs_[0]; ijpairs_[2*i_m] = ijpairs_[0]; ijpairs_[(2*i_m)+1] = ijpairs_[1];
        dpairs_[0] = d_m; ipairs_[0] = inn_m; ijpairs_[0] = iinn_m0; ijpairs_[1] = iinn_m1;
      }
      if (npairs_ > 2) sort_cluster_distances(dpairs_+1,ipairs_+1,ijpairs_+2,npairs_-1);
    }
    inline double comp_dclst_itry_kskp(int inmed_[],int k_,int ik_try_,int k_skp_)
    {
      const int n_nmed = npts-k_;
      double dclst_acc = get_min_dclst_i_kskp(i_meds[k_skp_],k_,ik_try_,k_skp_); // initialize w replaced k distance
      for (size_t i = 0; i < n_nmed; i++)
        if (inmed_[i] != ik_try_) dclst_acc += get_min_dclst_i_kskp(inmed_[i],k_,ik_try_,k_skp_);
      return dclst_acc;
    }
    inline double get_min_dclst_i_kskp(int i_,int k_,int ik_try_,int k_skp_)
    {
      double  dclst_min = dij(i_,ik_try_),
              dclst_kk;
      for (size_t kk = 0; kk < k_skp_; kk++)
        if (dclst_min > (dclst_kk = dij(i_,i_meds[kk]))) dclst_min = dclst_kk;
      for (size_t kk = k_skp_+1; kk < k_; kk++)
        if (dclst_min > (dclst_kk = dij(i_,i_meds[kk]))) dclst_min = dclst_kk;
      return dclst_min;
    }

    inline double comp_itry_dclst(int inmed_[],int k_,int ik_try_)
    {
      const int n_nmed = npts-k_;
      double dclst_acc = 0.0;
      for (size_t i = 0; i < n_nmed; i++)
        if (inmed_[i] != ik_try_) dclst_acc += get_min_dclst_i(inmed_[i],k_,ik_try_);
      return dclst_acc;
    }
    inline double get_min_dclst_i(int i_,int k_,int ik_try_)
    {
      double  dclst_min = dij(i_,ik_try_),
              dclst_kk;
      for (size_t kk = 0; kk < k_; kk++)
        if (dclst_min > (dclst_kk = dij(i_,i_meds[kk]))) dclst_min = dclst_kk;
      return dclst_min;
    }

    inline double dij(int i_,int j_) {return (i_<j_)?(dsym[i_][j_-i_-1]):((i_>j_)?(dsym[j_][i_-j_-1]):(0.0));}
    inline void populate_nmed(int inmed_[],bool fnmed_[],int k_)
    {
      for (size_t i = 0; i < npts; i++) fnmed_[i] = true;
      for (size_t i = 0; i < k_; i++) fnmed_[i_meds[i]] = false;
      for (size_t i = 0, ii = 0; i < npts; i++) if (fnmed_[i]) inmed_[ii++] = i;
    }
    inline double comp_ipts_net_d(int i_)
    {
      double acc = 0.0;
      for (size_t j = 0; j < npts; j++) acc += dij(i_,j);
      return acc;
    }
    inline int imin_net_d()
    {
      int imin_out = 0;
      double  dmin = comp_ipts_net_d(0),
              dmin_comp;
      for (size_t i = 1; i < npts; i++)
        if (dmin > ( dmin_comp = comp_ipts_net_d(i) ))
          {imin_out = i; dmin = dmin_comp;}
      return imin_out;
    }
    inline void assign_cluster_members(bool fmem_clst_[],int k_,int inmed_[])
    {
      const int n_nmed = npts-k_;
      for (size_t k = 0; k < k_; k++) // assign cluster centres to selves
      {
        int ifm = i_meds[k];
        for (size_t kk = 0; kk < k; kk++, ifm+=npts) fmem_clst_[ifm] = 0;
        fmem_clst_[ifm] = 1; ifm+=npts;
        for (size_t kk = k+1; kk < k_; kk++, ifm+=npts) fmem_clst_[ifm] = 0;
      }
      for (size_t inm = 0; inm < n_nmed; inm++) // assign non medoid members
      {
        const int kclst_inm = assign_ipt_kcluster(inmed_[inm],k_);
        int ifnm = inmed_[inm];
        for (size_t k = 0; k < kclst_inm; k++, ifnm+=npts) fmem_clst_[ifnm] = 0;
        fmem_clst_[ifnm] = 1; ifnm+=npts;
        for (size_t k = kclst_inm+1; k < k_; k++, ifnm+=npts) fmem_clst_[ifnm] = 0;
      }
    }
    inline int assign_ipt_kcluster(int i_, int k_)
    {
      int kclst_out = 0;
      double  dclst_min = dij(i_,i_meds[kclst_out]),
              dclst_try;
      for (size_t k = 1; k < k_; k++)
        if (dclst_min > (dclst_try=dij(i_,i_meds[k])))
          {dclst_min = dclst_try; kclst_out = k;}

      return kclst_out;
    }
    inline double get_min_dclst_i_kskp(int i_,int k_,int k_skp_)
    {
      double  dclst_min = DBL_MAX,
              dclst_kk;
      for (size_t kk = 0; kk < k_skp_; kk++)
        if (dclst_min > (dclst_kk = dij(i_,i_meds[kk]))) dclst_min = dclst_kk;
      for (size_t kk = k_skp_+1; kk < k_; kk++)
        if (dclst_min > (dclst_kk = dij(i_,i_meds[kk]))) dclst_min = dclst_kk;
      return dclst_min;
    }
};

#endif
