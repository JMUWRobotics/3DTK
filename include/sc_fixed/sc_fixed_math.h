#ifndef SC_FIXED_MATH_H
#define SC_FIXED_MATH_H

#define SC_INCLUDE_FX

#include <systemc.h>
#include <array>
#include <cstring>
#include <iostream>

#ifndef FIXED_WORD_LENGTH
#define FIXED_WORD_LENGTH 48
#endif

#ifndef FIXED_INT_WORD_LENGTH
#define FIXED_INT_WORD_LENGTH 36
#endif

// Additional datatype for small, precise values
#ifndef FIXED_SMALL_WORD_LENGTH
#define FIXED_SMALL_WORD_LENGTH 48
#endif

#ifndef FIXED_SMALL_INT_WORD_LENGTH
#define FIXED_SMALL_INT_WORD_LENGTH 10
#endif

#ifndef FIXED_HERON_ITERATIONS
#define FIXED_HERON_ITERATIONS 10
#endif

using namespace sc_dt;
using fixed_val = sc_fixed<FIXED_WORD_LENGTH, FIXED_INT_WORD_LENGTH>;
using fixed_small_val = sc_fixed<FIXED_SMALL_WORD_LENGTH, FIXED_SMALL_INT_WORD_LENGTH>;

// function to calculate the square root via heron method
template <typename T>

T sc_fixed_heron_sqrt(T s) {

  if (s == 0) {
    return 0;
  }
  if (s < 0) {
    std::cerr << "warning: square root from negative numbers not allowed --> will return 0" << std::endl;
    return 0;
  }

  T x = s / (2 + 1);
  T x_n;

  for(int i = FIXED_HERON_ITERATIONS; i > 0; i--) {
    x_n = (x + s / x) / 2;
    x = x_n;
  }

  return x_n;

}

bool sc_choldc(fixed_val A[3][3], fixed_val diag[3]);
void sc_cholsl(fixed_val A[3][3], fixed_val diag[3], fixed_val B[3], fixed_small_val x[3]);

void transform(std::vector<std::array<fixed_val, 3>>& scan, fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::array<fixed_small_val, 9>& transMat_R, std::array<fixed_val, 3>& transMat_t , std::array<fixed_small_val, 9>& dalignxf_R, std::array<fixed_val, 3>& dalignxf_t, std::ofstream& frame, int islum);
void transformPoints(const fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::vector<std::array<fixed_val, 3>>& scan);
void transformMatrix(const fixed_small_val alignxf_R[9], const fixed_val alignxf_t[3], std::array<fixed_small_val, 9>& transMat_R, std::array<fixed_val, 3>& transMat_t, std::array<fixed_small_val, 9>& dalignxf_R, std::array<fixed_val, 3>& dalignxf_t);
void transform3(const fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::array<fixed_val, 3>& point);
void MMult(const fixed_small_val M1[9], const std::array<fixed_small_val, 9>& M2, std::array<fixed_small_val, 9>& Mout);
fixed_small_val sc_abs(fixed_small_val);

#endif // SC_FIXED_MATH_H
