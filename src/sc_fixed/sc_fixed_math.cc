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
