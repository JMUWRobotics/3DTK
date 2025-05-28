#include "sc_fixed_math/sc_fixed_math.h"


f_float sc_fixed_heron_sqrt(f_float s) {

  f_float x = (s+1)/2;
  f_float x_n;
  
  for(int i = HERON_ITERATIONS; i > 0; i--) {
    x_n = (x + s / x) / 2;
    x = x_n;
  }

  return x_n;
  
}

// For testing
int sc_main(int argc, char* argv[]){

  f_float s = 163;
  f_float res = sc_fixed_heron_sqrt(s);

  cout << res << endl;
  
  return 0;

}
// END testing
