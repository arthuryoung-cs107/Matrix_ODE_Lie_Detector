#ifndef LD_UTIL_HH
#define LD_UTIL_HH

#include <cstddef>
#include <cstdio>

#ifdef _OPENMP
  #include "omp.h"
#endif

struct LD_program_meta
{
  LD_program_meta() {}
  ~LD_program_meta() {}

  static const char * dir_name;
  static const char * dat_suff;

  static const int eor_meta0;
  static const int ndep_meta0;
  static const int bor_fspace0;
};

struct LD_threads
{
  LD_threads() {}
  ~LD_threads() {}

#ifdef _OPENMP
  static int thread_id() {return omp_get_thread_num();}
  static int numthreads() {return omp_get_max_threads();}
  static double tic() {return omp_get_wtime();}
  static double toc(double t0_) {return omp_get_wtime()-t0_;}
#else
  static int thread_id() {return 0;}
  static int numthreads() {return 1;}
  static double tic() {return 0.0;}
  static double toc(double t0_) {return 0.0;}
#endif
};

template <typename T> T sum_Tvec(T *in_, size_t len_)
{
  T Tacc = 0;
  for (size_t i = 0; i < len_; i++) Tacc += in_[i];
  return Tacc;
}
template <typename T> T * Tvec_copy(T *in_, size_t len_)
{
  T * Tout = new T[len_];
  for (size_t i = 0; i < len_; i++) Tout[i] = in_[i];
  return Tout;
}
template <typename T> void copy_Tvec(T *out_, T *in_, size_t len_)
  {for (size_t i = 0; i < len_; i++) out_[i] = in_[i];}

template <typename T> T ** Tmatrix_ptrs(T *in_, size_t nptrs_, size_t nskip_)
{
  T ** out_ = new T*[nptrs_];
  for (size_t i = 0; i < nptrs_; i++) out_[i] = in_+(i*nskip_);
  return out_;
}
template <typename T> T ** Tmatrix(int M_, int N_)
{
  T * chunk = new T[M_*N_],
    ** rows = new T*[M_];
  for (int i = 0,j=0; i < M_; i++,j+=N_)
    rows[i] = chunk+j;
  return rows;
}
template <typename T> void free_Tmatrix(T ** Tmat_)
{
  delete [] Tmat_[0];
  delete [] Tmat_;
}
template <typename T> T *** T3tensor(int L_, int M_, int N_)
{
  T * chunk = new T[L_*M_*N_],
    ** rows = new T*[L_*M_],
    *** mats = new T**[L_];
  for (size_t l = 0; l < L_; l++)
  {
    mats[l] = rows + l*M_;
    for (size_t m = 0; m < M_; m++)
      mats[l][m] = chunk + (l*M_*N_) + (m*N_);
  }
  return mats;
}
template <typename T> void free_T3tensor(T *** Ttens_)
{
  delete [] Ttens_[0][0];
  delete [] Ttens_[0];
  delete [] Ttens_;
}
template <typename T> T **** T4tensor(int K_, int L_, int M_, int N_)
{
  T * chunk = new T[K_*L_*M_*N_],
    ** rows = new T*[K_*L_*M_],
    *** mats = new T**[K_*L_],
    **** tsrs = new T***[K_];
  for (size_t k = 0; k < K_; k++)
  {
    tsrs[k] = mats + (k*L_);
    for (size_t l = 0; l < L_; l++)
    {
      tsrs[k][l] = rows + (k*L_*M_) + (l*M_);
      for (size_t m = 0; m < M_; m++)
        tsrs[k][l][m] = chunk + (k*L_*M_*N_) + (l*M_*N_) + (m*N_);
    }
  }
  return tsrs;
}
template <typename T> void free_T4tensor(T **** Ttens_)
{
  delete [] Ttens_[0][0][0];
  delete [] Ttens_[0][0];
  delete [] Ttens_[0];
  delete [] Ttens_;
}
template <typename T> T ** Tsym(int M_)
{
  T * chunk = new T[((M_+1)*(M_))/2],
    ** rows = new T*[M_];
  rows[0] = chunk;
  for (int i = 1, j=M_; i < M_; i++,j--) rows[i] = rows[i-1]+j;
  return rows;
}
template <typename T> T ** Tsym_lower(int M_)
{
  T * chunk = new T[((M_+1)*(M_))/2],
    ** rows = new T*[M_];
  rows[0] = chunk;
  for (int i = 1; i < M_; i++) rows[i] = rows[i-1]+i;
  return rows;
}

class LD_util
{
  public:

    LD_util() {}
    ~LD_util() {}

    template <typename T> static void flip_Tdata(T *dat_, int len_)
    {
      for (int i = 0, len_o2 = ((int)(len_/2)), iend = len_-1; i < len_o2; i++, iend--)
      {
        double di = dat_[i];
        dat_[i] = dat_[iend];
        dat_[iend] = di;
      }
    }
    template <typename T> static void flip_Tmatrix_rows(T **dat_, int len1_, int len2_)
    {
      const int len2_o2 = ((int)(len2_/2));
      for (int i = 0; i < len1_; i++)
        for (int j = 0, jend = len2_-1; j < len2_o2; j++, jend--)
        {
          // LD_util::flip_Tdata<T>(dat_[i],len2_);
          double dij = dat_[i][j];
          dat_[i][j] = dat_[i][jend];
          dat_[i][jend] = dij;
        }
    }
    template <typename T> static void transpose_Tmatrix(T **dat_, int len_)
    {
      for (int i = 0, lenm1 = len_-1; i < lenm1; i++)
        for (int j = i+1; j < len_; j++)
        {
          T dij = dat_[i][j];
          dat_[i][j] = dat_[j][i];
          dat_[j][i] = dij;
        }
    }
    template <typename T> static T min_Tdata(T *dat_, int len_, int &ind_)
    {
      T val_out = dat_[ind_=0];
      for (size_t i = 1; i < len_; i++)
        if (val_out > dat_[i])
          val_out = dat_[ind_ = i];
      return val_out;
    }
    template <typename T> static T max_Tdata(T *dat_, int len_, int &ind_)
    {
      T val_out = dat_[ind_=0];
      for (size_t i = 1; i < len_; i++)
        if (val_out < dat_[i])
          val_out = dat_[ind_ = i];
      return val_out;
    }
    template <typename T> static T min_Tdata(T *dat_,int len_)
    {
      T val_out = dat_[0];
      for (size_t i = 1; i < len_; i++) if (val_out > dat_[i]) val_out = dat_[i];
      return val_out;
    }
    template <typename T> static T max_Tdata(T *dat_,int len_)
    {
      T val_out = dat_[0];
      for (size_t i = 1; i < len_; i++) if (val_out < dat_[i]) val_out = dat_[i];
      return val_out;
    }
    template <typename T> static T min_T(T a_, T b_) {return (a_<b_)?(a_):b_;}
    template <typename T> static T max_T(T a_, T b_) {return (a_>b_)?(a_):b_;}
    template <typename T> static T min_T(T a_, T b_, T c_) {return min_T<T>(min_T<T>(a_,b_),c_);}
    template <typename T> static T max_T(T a_, T b_, T c_) {return max_T<T>(max_T<T>(a_,b_),c_);}
    template <typename T> static void copy_Tdata(T *dato_,int len_,T *dati_)
      { for (int i = 0; i < len_; i++) dato_[i] = dati_[i]; }
    template <typename T> static T sum_Tdata(T *dat_,int len_)
    {
      T acc = 0;
      for (int i = 0; i < len_; i++) acc += dat_[i];
      return acc;
    }

};

#endif
