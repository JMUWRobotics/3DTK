/*
 * sc_fixed_converter implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of conversion functions between double/DataXYZ and SystemC fixed-point
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer. Institute of Computer Science, University of Wuerzburg, Germany.
 */

#include "sc_fixed/sc_fixed_converter.h"

// converts TripleArray<double> (aka DataXYZ) into vactor array holding fixed_vals (fixed_val points)
std::vector<std::array<fixed_val, 3>> array2fixedArray(const DataXYZ &input) {
  std::vector<std::array<fixed_val, 3>> result;
  result.reserve(input.size());

  for (unsigned int i = 0; i < input.size(); ++i) {
    std::array<fixed_val, 3> point{{
      fixed_val(input[i][0]),
      fixed_val(input[i][1]),
      fixed_val(input[i][2])}
    };
    result.push_back(point);
  }
  cout << "*******" << result.size() << " points"  << endl;
  return result;
}

// converts transMat and dalignxf (as double[16]) to fixed_val transMat and dalignxf (as std::array<fixed_val, 16>)
std::array<fixed_val, 16> array2fixedArray16(const double input[16]) {
  std::array<fixed_val, 16> result;
  for (unsigned int i = 0; i < 16; ++i) {
    result[i] = input[i];
  }
  return result;
}

// prints the coordinates from the points from a point array
void printPoints(const std::vector<std::array<fixed_val, 3>>& points) {
  for (const auto& point : points) {
    std::cout << "x: " << point[0]
              << ", y: " << point[1]
              << ", z: " << point[2] << std::endl;
  }
}

// writes a transformation matrix to a .frame-file
void writeFrame(std::ofstream& frame, std::array<fixed_val, 16>& matrix, int viewFactor) {
  for(unsigned int i = 0; i < matrix.size(); i++) {
      frame << static_cast<double>(matrix[i]);
      frame << " ";
    }
    frame << viewFactor << "\n";
}
