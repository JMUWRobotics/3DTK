#include "sc_fixed/sc_fixed_math.h"


f_float sc_fixed_heron_sqrt(f_float s) {

  f_float x = (s+1)/2;
  f_float x_n;
  
  for(int i = HERON_ITERATIONS; i > 0; i--) {
    x_n = (x + s / x) / 2;
    x = x_n;
  }

  return x_n;
  
}

static inline bool sc_choldc(f_float A[3][3], f_float diag[3])                                    
{
  
  unsigned int N = 3;
  const f_float epsilon = f_float(1e-7);                                            
  
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = i; j < N; j++) {
      f_float sum = A[i][j];
      for (int k=i-1; k >= 0; k--)
        sum -= A[i][k] * A[j][k];
      if (i == j) {
        if (sum < epsilon)
          return false;
        diag[i] = sc_fixed_heron_sqrt(sum);
      } else {
        A[j][i] = sum / diag[i];
      }
    }
  }
  return true;
}

static inline void sc_cholsl(f_float A[3][3],
                          f_float diag[3],
                          f_float B[3],
                          f_float x[3])
{
  int N = 3;
  for (int i=0; i < N; i++) {
    f_float sum = B[i];
    for (int k=i-1; k >= 0; k--)
      sum -= A[i][k] * x[k];
    x[i] = sum / diag[i];
  }
  for (int i=N-1; i >= 0; i--) {
    f_float sum = x[i];
    for (int k=i+1; k < N; k++)
      sum -= A[k][i] * x[k];
    x[i] = sum / diag[i];
  }
}


