#ifndef SC_FIXED_MATH_H
#define SC_FIXED_MATH_H

#define SC_INCLUDE_FX

#include <systemc.h>
#include <array>
#include <cstring>

#ifndef FIXED_WORD_LENGTH
#define FIXED_WORD_LENGTH 48
#endif

#ifndef FIXED_INT_WORD_LENGTH
#define FIXED_INT_WORD_LENGTH 36
#endif

#ifndef FIXED_HERON_ITERATIONS
#define FIXED_HERON_ITERATIONS 3
#endif

using namespace sc_dt;
using f_float = sc_fixed<FIXED_WORD_LENGTH, FIXED_INT_WORD_LENGTH>;

f_float sc_fixed_heron_sqrt(f_float s);
bool sc_choldc(f_float A[3][3], f_float diag[3]);
void sc_cholsl(f_float A[3][3], f_float diag[3], f_float B[3], f_float x[3]);

void transform(std::vector<std::array<f_float, 3>>& scan, f_float alignxf[16], std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf, int islum);
void transformReduced(const f_float alignxf[16], std::vector<std::array<f_float, 3>>& scan);
void transformMatrix(const f_float alignxf[16], std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf);
void transform3(const f_float alignxf[16], std::array<f_float, 3>& point);
void MMult(const f_float M1[16], const std::array<f_float, 16>& M2, std::array<f_float, 16>& Mout);
f_float sc_abs(f_float);

#endif // SC_FIXED_MATH_H
