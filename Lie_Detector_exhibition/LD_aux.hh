#ifndef LD_AUX_HH
#define LD_AUX_HH

#include "LD_util.hh"

#include <cstddef>
#include <cstdio>
#include <cmath>

#include <cfloat>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

  static double eps(double x_=1.0) {return nextafter(x_,DBL_MAX)-x_;}

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

  template <typename T> static void copy_vec(T *out_, T *in_, int len_)
    {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}

  template <typename T> static void fill_vec(T *out_, int len_, T val_=0)
    {for (size_t i = 0; i < len_; i++) out_[i] = val_;}

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

  static double l2_x(double *x_, int n_)
  {
    double r2_acc = 0.0;
    for (size_t i = 0; i < n_; i++) r2_acc += x_[i]*x_[i];
    return sqrt(r2_acc);
  }

  static void fill_x_012(double *x_, int n_, double offset_=0.0)
    {for (size_t i = 0; i < n_; i++) x_[i] = ((double)(i)) + offset_;}
  static void fill_x_012(int *x_, int n_, int offset_=0)
    {for (size_t i = 0; i < n_; i++) x_[i] = i + offset_;}
  static void fill_x(double *x_, int n_, double val_=0.0)
    {for (size_t i = 0; i < n_; i++) x_[i] = val_;}
  static void fill_x(int *x_, int n_, int val_=0)
    {for (size_t i = 0; i < n_; i++) x_[i] = val_;}
  static void fill_x(bool *x_, int n_, bool val_=false)
    {for (size_t i = 0; i < n_; i++) x_[i] = val_;}

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

  static void load_gsl_vec(gsl_vector *gsl_vec_, double *vec_)
    {for (size_t i = 0; i < gsl_vec_->size; i++) gsl_vector_set(gsl_vec_,i,vec_[i]);}
  static void unpack_gsl_vec(double *vec_, gsl_vector *gsl_vec_)
    {for (size_t i = 0; i < gsl_vec_->size; i++) vec_[i] = gsl_vector_get(gsl_vec_,i);}

  static void load_gsl_matT(gsl_matrix * gsl_mat_, double ** mat_)
  {
    for (size_t i = 0; i < gsl_mat_->size1; i++)
      for (size_t j = 0; j < gsl_mat_->size2; j++)
        gsl_matrix_set(gsl_mat_,j,i,mat_[i][j]);
  }
  static void load_gsl_mat(gsl_matrix * gsl_mat_, double ** mat_)
  {
    for (size_t i = 0; i < gsl_mat_->size1; i++)
      for (size_t j = 0; j < gsl_mat_->size2; j++)
        gsl_matrix_set(gsl_mat_,i,j,mat_[i][j]);
  }
  static void unpack_gsl_mat(double ** mat_, gsl_matrix * gsl_mat_)
  {
    for (size_t i = 0; i < gsl_mat_->size1; i++)
      for (size_t j = 0; j < gsl_mat_->size2; j++)
        mat_[i][j] = gsl_matrix_get(gsl_mat_,i,j);
  }
  static void unpack_gsl_matT(double ** mat_, gsl_matrix * gsl_mat_)
  {
    for (size_t i = 0; i < gsl_mat_->size1; i++)
      for (size_t j = 0; j < gsl_mat_->size2; j++)
        mat_[i][j] = gsl_matrix_get(gsl_mat_,j,i);
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

#endif
