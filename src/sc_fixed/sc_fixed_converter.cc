#include "sc_fixed/sc_fixed_converter.h"

// converts TripleArray<double> (aka DataXYZ) into vactor array holding f_floats (fixed points)

std::vector<std::array<f_float, 3>> array2fixedArray(const DataXYZ &input) {
  std::vector<std::array<f_float, 3>> result;
  result.reserve(input.size());

  for (unsigned int i = 0; i < input.size(); ++i) {
    std::array<f_float, 3> point{{
      f_float(input[i][0]),
      f_float(input[i][1]),
      f_float(input[i][2])}
    };
    result.push_back(point);
  }
  cout << "*******" << result.size() << endl;
  return result;
}

// converts transMat and dalignxf (as double[16]) to fixed transMat and dalignxf (as std::array<f_float, 16>)

std::array<f_float, 16> array2fixedArray16(const double input[16]) {
  std::array<f_float, 16> result;
  for (unsigned int i = 0; i < 16; ++i) {
    result[i] = input[i];
  }
  return result;
}

// prints the coordinates from the points from a point array

void printPoints(const std::vector<std::array<f_float, 3>>& points) {
  for (const auto& point : points) {
    std::cout << "x: " << point[0]
              << ", y: " << point[1]
              << ", z: " << point[2] << std::endl;
  }
}
