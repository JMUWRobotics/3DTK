#ifndef SC_FIXED_MATH_H
#define SC_FIXED_MATH_H

#define SC_INCLUDE_FX

#include <systemc.h>

constexpr int WORD_LENGTH = 24;
constexpr int INT_WORD_LENGTH = 16;

constexpr int HERON_ITERATIONS = 5;

using namespace sc_dt;
using f_float = sc_fixed<WORD_LENGTH, INT_WORD_LENGTH>;

f_float sc_fixed_heron_sqrt(f_float s);

#endif // SC_FIXED_MATH_H