#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <iostream>
#include <limits>
#include <type_traits>

#include <pthread.h>
#include <sys/time.h>
#include <time.h>

  
double timestamp() {
  struct timeval tv;
  gettimeofday(&tv, nullptr);
  return tv.tv_sec + 1e-6*static_cast<double>(tv.tv_usec);
}


//template <typename T>

extern "C" {
  void igemm(int *Cg, int *Ag, int *Bg, int n, int m, int _k);
}


template <typename T> 
void naive_gemm(T *Cg, T *Ag, T *Bg, int n, int m, int _k) {
  //the whole matrix
  for(int i = 0; i < n; i++) {
    for(int k = 0; k < _k; k++) {
      for(int j = 0; j < m; j++) {
	  Cg[i*n+j] += Ag[i*_k+k]*Bg[k*m+j];
      }
    }
  }

  
}


template <typename T>
T max(T a, T b) {
  return (a<b) ? b : a;
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, T>::type* = nullptr>
T rand() {
  static const size_t M = static_cast<size_t>(1) << std::numeric_limits<T>::digits;
  T r = static_cast<T>((rand() & (M-1)));
  return (r / static_cast<T>(M/16));
}

#define print_var(X) std::cout << #X << " = " << X << "\n";

void *__dso_handle = nullptr;

int main() {
  uint32_t l = 1800;
  assert((l % 12) == 0);
  uint32_t n = l ,m= l ,k= l;

  print_var(n);
  print_var(m);
  print_var(k);
    
  int *C0, *C1, *A, *B;

  
  C0 = new int[n*m];
  C1 = new int[n*m];
  A = new int[n*k];
  B = new int[k*m];

  memset(C0,0,sizeof(int)*n*m);
  memset(C1,0,sizeof(int)*n*m);

  uint64_t flops = 2UL * n * m * k;
  double gflops = static_cast<double>(flops)*1e-9;
  std::cout << "working set = " << (sizeof(float)*((n*m) + (n*k) + (k*m))) / 1024 << " kbytes\n";
  std::cout << "total gflops = " << gflops  << "\n";
  
  for(int i = 0; i < (n*k); i++) 
    A[i] = rand();
  for(int i = 0; i < (m*k); i++) 
    B[i] = rand();

  double t0 = timestamp();
  igemm(C0,A,B,n,m,k);
  double t1 = timestamp();

  std::cout << "optimized gemm " << gflops/(t1-t0) << " gflops/s\n";

  // while(1) {
  t0 = timestamp();
  naive_gemm<int>(C1,A,B,n,m,k);
  t1 = timestamp();
  std::cout << "naive gemm " << gflops/(t1-t0) << " gflops/s\n";
  //}

  
  double max_e = 0.0f;
  for(int i = 0; i < (n*m); i++) {
    double e = C0[i]-C1[i];
    e *= e;
    max_e = max(max_e, e);
    if(e > 1.0) {
      std::cout << "naive = " << C0[i] << "\n";
      std::cout << "opt   = " << C1[i] << "\n";
      break;
    }
  }
  std::cout << "Error = " << max_e << std::endl;

  delete [] (C0);
  delete [] (C1);
  delete [] (A);
  delete [] (B);
  return 0;
}
