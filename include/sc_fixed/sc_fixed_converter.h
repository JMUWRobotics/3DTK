#ifndef SC_FIXED_CONVERTER
#define SC_FIXED_CONVERTER

#include "slam6d/data_types.h"
#include <vector>
#include <array>
#include "sc_fixed/sc_fixed_math.h"


std::vector<std::array<f_float, 3>> array2fixedArray(const DataXYZ &input);

DataXYZ fixedArray2array(TripleArray<f_float>);

#endif //SC_FIXED_CONVERTER