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

// splits transformation matrix in rotation matrix and translation vector and transforms into fixed point values

void splitMat(const double input[16], std::array<fixed_small_val, 9>& R, std::array<fixed_val, 3>& t){
  R[0] = fixed_small_val(input[0]);
  R[1] = fixed_small_val(input[1]);
  R[2] = fixed_small_val(input[2]);
  R[3] = fixed_small_val(input[4]);
  R[4] = fixed_small_val(input[5]);
  R[5] = fixed_small_val(input[6]);
  R[6] = fixed_small_val(input[8]);
  R[7] = fixed_small_val(input[9]);
  R[8] = fixed_small_val(input[10]);
  t[0] = fixed_val(input[12]);
  t[1] = fixed_val(input[13]);
  t[2] = fixed_val(input[14]);
}


// converts TripleArray<double> (aka DataXYZ) into vector array holding fixed_vals (fixed_val points)
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

// converts transMat and dalignxf (as double[16]) to fixed_small_val transMat and dalignxf (as std::array<fixed_small_val, 16>)
std::array<fixed_small_val, 16> array2fixedArray16(const double input[16]) {
  std::array<fixed_small_val, 16> result;
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


void writeFrame(std::ofstream& frame, std::array<fixed_small_val, 9>& R, std::array<fixed_val, 3>& t, int viewFactor) {
  
  for(unsigned int i = 0; i < 3; i++) {
    for(unsigned j = 0;j < 3; ++j){
      frame << static_cast<double>(R[i * 3 + j])<< " ";
    }
      frame << "0 ";
    }
    for(unsigned int i = 0; i < 3; ++i){
      frame << static_cast<double>(t[i]) << " ";
    }
    frame << "1 ";
    frame << viewFactor << "\n";
  }