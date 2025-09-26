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

#ifndef FIXED_HERON_ITERATIONS
#define FIXED_HERON_ITERATIONS 3
#endif

using namespace sc_dt;
using fixed_val = sc_fixed<FIXED_WORD_LENGTH, FIXED_INT_WORD_LENGTH>;

fixed_val sc_fixed_heron_sqrt(fixed_val s);
bool sc_choldc(fixed_val A[3][3], fixed_val diag[3]);
void sc_cholsl(fixed_val A[3][3], fixed_val diag[3], fixed_val B[3], fixed_val x[3]);

void transform(std::vector<std::array<fixed_val, 3>>& scan, fixed_val alignxf[16], std::array<fixed_val, 16>& transMat, std::array<fixed_val, 16>& dalignxf, std::ofstream& frame, int islum);
void transformPoints(const fixed_val alignxf[16], std::vector<std::array<fixed_val, 3>>& scan);
void transformMatrix(const fixed_val alignxf[16], std::array<fixed_val, 16>& transMat, std::array<fixed_val, 16>& dalignxf);
void transform3(const fixed_val alignxf[16], std::array<fixed_val, 3>& point);
void MMult(const fixed_val M1[16], const std::array<fixed_val, 16>& M2, std::array<fixed_val, 16>& Mout);
fixed_val sc_abs(fixed_val);

#endif // SC_FIXED_MATH_H
