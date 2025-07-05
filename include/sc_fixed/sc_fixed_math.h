#ifndef SC_FIXED_MATH_H
#define SC_FIXED_MATH_H

#define SC_INCLUDE_FX

#include <systemc.h>
#include <array>
#include <cstring>

const int WORD_LENGTH = 24;
const int INT_WORD_LENGTH = 8;

const int HERON_ITERATIONS = 5;

using namespace sc_dt;
using f_float = sc_fixed<WORD_LENGTH, INT_WORD_LENGTH>;

f_float sc_fixed_heron_sqrt(f_float s);
bool sc_choldc(f_float A[3][3], f_float diag[3]);
void sc_cholsl(f_float A[3][3], f_float diag[3], f_float B[3], f_float x[3]);

void transform(std::vector<std::array<f_float, 3>>& scan, const std::array<f_float, 16> alignxf, std::array<f_float, 16> transMat, std::array<f_float, 16> dalignxf, int islum);
void transformReduced(const std::array<f_float, 16> alignxf, std::vector<std::array<f_float, 3>>& scan);
void transformMatrix(const std::array<f_float, 16> alignxf, std::array<f_float, 16> transMat, std::array<f_float, 16> dalignxf);
void transform3(const std::array<f_float, 16>& alignxf, std::array<f_float, 3>& point);
void MMult(const std::array<f_float, 16>& M1, const std::array<f_float, 16>& M2, std::array<f_float, 16>& Mout);

#endif // SC_FIXED_MATH_H
