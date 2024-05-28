#ifndef LD_UTIL_HH
#define LD_UTIL_HH

#include <cstddef>

#ifdef _OPENMP
  #include "omp.h"
// #else
//   #include <ctime.h>
#endif

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

#endif
