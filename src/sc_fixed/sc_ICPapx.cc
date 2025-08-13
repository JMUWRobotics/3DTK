/*
 * sc_ICPapx implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */


/**
 *  @file
 *  @brief Implementation of the ICP error function minimization via small angle approximation
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#include "sc_fixed/sc_ICPapx.h"

#include "slam6d/globals.icc"
#include "sc_fixed/sc_fixed_math.h"
#include <iomanip>
#include <cstring>

/**
 * computes the rotation matrix consisting
 * of a rotation and translation that
 * minimizes the root-mean-square error
 * of the point pairs, using the <b>approximation</b>
 * sin(x) = x.
 *
 * @params matchedSource, matchedTarget: point pairs (pairs of corresponding points)
 * @params centerSource, centerTarget: centroids of source and target
 * @param alignxf: resulting transformation matrix
 * @return Error estimation of the matching (rms)
 */
f_float sc_ICPapx::Align(const std::vector<std::array<f_float, 3>>& matchedSource,
                      const std::vector<std::array<f_float, 3>>& matchedTarget,
                      f_float *alignxf,
                      const std::array<f_float, 3> centerSource,
                      const std::array<f_float, 3> centerTarget){
  int n = matchedSource.size();
  
  f_float A[3][3];
  f_float B[3]; 

  for(int i=0; i<3;i++) {
    B[i] = f_float(0.0);
    for(int j=0; j<3;j++){
      A[i][j] = f_float(0.0); 
    }
  }	

  f_float sum = 0.0; 
  f_float p1[3], p2[3]; 

  for(size_t i = 0; i < matchedSource.size(); ++i){
    p1[0] = matchedSource[i][0];
    p1[1] = matchedSource[i][1];
    p1[2] = matchedSource[i][2];
    p2[0] = matchedTarget[i][0];
    p2[1] = matchedTarget[i][1];
    p2[2] = matchedTarget[i][2]; 

    f_float p12[3] = {p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]};
    f_float p2c[3] = {p2[0] - centerTarget[0], p2[1] - centerTarget[1], p2[2] - centerTarget[2]};  

    sum += p12[0]*p12[0] + p12[1]*p12[1] + p12[2]*p12[2];

    B[0] += (p12[2]*p2c[1] - p12[1]*p2c[2]);
    B[1] += (p12[0]*p2c[2] - p12[2]*p2c[0]);
    B[2] += (p12[1]*p2c[0] - p12[0]*p2c[1]);
    A[0][0] += ((p2c[1])*(p2c[1]) + (p2c[2])*(p2c[2]));
    A[0][1] -= p2c[0] * p2c[1];
    A[0][2] -= p2c[0] * p2c[2];
    A[1][1] += ((p2c[0])*(p2c[0]) + (p2c[2])*(p2c[2]));
    A[1][2] -= p2c[1] * p2c[2];
    A[2][2] += ((p2c[0])*(p2c[0]) + (p2c[1])*(p2c[1]));
  }
  
  f_float error = sc_fixed_heron_sqrt(sum / n);
  std::cout << "Error: " << error << std::endl;

  f_float diag[3];
 
  if (!sc_choldc(A, diag)) {
    printf("Couldn't find transform.\n");
    return -1.0;
  }
  f_float x[3];
  sc_cholsl(A, diag, B, x);

  // Interpret results
  f_float sx = x[0];
  f_float cx = sc_fixed_heron_sqrt(1.0 - sx*sx);
  f_float sy = x[1];
  f_float cy = sc_fixed_heron_sqrt(1.0 - sy*sy);
  f_float sz = x[2];
  f_float cz = sc_fixed_heron_sqrt(1.0 - sz*sz);
  
  std::cout << "x[0] sx: " << diag[0] << std::endl;
  std::cout << "x[1] sy: " << diag[1] << std::endl;
  std::cout << "x[2] sz: " << diag[2] << std::endl;

  alignxf[0]  = cy*cz;
  alignxf[1]  = sx*sy*cz + cx*sz;
  alignxf[2]  = -cx*sy*cz + sx*sz;
  alignxf[3]  = 0;
  alignxf[4]  = -cy*sz;
  alignxf[5]  = -sx*sy*sz + cx*cz;
  alignxf[6]  = cx*sy*sz + sx*cz;
  alignxf[7]  = 0;
  alignxf[8]  = sy;
  alignxf[9]  = -sx*cy;
  alignxf[10] = cx*cy;
  alignxf[11] = 0;
  alignxf[12] = centerSource[0] - alignxf[0]*centerTarget[0] - alignxf[4]*centerTarget[1] - alignxf[8]*centerTarget[2];
  alignxf[13] = centerSource[1] - alignxf[1]*centerTarget[0] - alignxf[5]*centerTarget[1] - alignxf[9]*centerTarget[2];
  alignxf[14] = centerSource[2] - alignxf[2]*centerTarget[0] - alignxf[6]*centerTarget[1] - alignxf[10]*centerTarget[2];
  alignxf[15] = 1;
  
  return error; // gib den Fehler zurück
}
