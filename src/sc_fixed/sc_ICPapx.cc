/*
 * sc_ICPapx implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of the ICP error function minimization via small angle approximation
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer. Institute of Computer Science, University of Wuerzburg, Germany.
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
 * of the point pairs, using the approximation sin(x) = x.
 *
 * @params matchedSource, matchedTarget: point pairs (pairs of corresponding points)
 * @params centerSource, centerTarget: centroids of source and target
 * @param alignxf: resulting transformation matrix
 * @return Error estimation of the matching (rms)
 */
fixed_small_val sc_ICPapx::Align(const std::vector<std::array<fixed_val, 3>>& matchedSource,
                      const std::vector<std::array<fixed_val, 3>>& matchedTarget,
                      fixed_small_val *alignxf_R,
                      fixed_val *alignxf_t,
                      const std::array<fixed_val, 3> centerSource,
                      const std::array<fixed_val, 3> centerTarget){
  int n = matchedSource.size();

  fixed_val A[3][3];
  fixed_val B[3];

  for(int i = 0; i < 3; i++) {
    B[i] = fixed_val(0.0);
    for(int j = 0; j < 3; j++){
      A[i][j] = fixed_val(0.0);
    }
  }

  fixed_val sum = 0.0;
  fixed_val p1[3], p2[3];

  for(size_t i = 0; i < matchedSource.size(); ++i){
    p1[0] = matchedSource[i][0];
    p1[1] = matchedSource[i][1];
    p1[2] = matchedSource[i][2];
    p2[0] = matchedTarget[i][0];
    p2[1] = matchedTarget[i][1];
    p2[2] = matchedTarget[i][2];

    fixed_val p12[3] = {p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]};
    fixed_val p2c[3] = {p2[0] - centerTarget[0], p2[1] - centerTarget[1], p2[2] - centerTarget[2]};

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

  fixed_small_val error = sc_fixed_heron_sqrt<fixed_small_val>(sum / n);

  fixed_val diag[3];

  if (!sc_choldc(A, diag)) {
    std::cout << "Couldn't find transform" << std::endl;
    return -1.0;
  }

  fixed_small_val x[3];
  sc_cholsl(A, diag, B, x);

  // Interpret results
  fixed_small_val sx = x[0];
  fixed_small_val cx = sc_fixed_heron_sqrt<fixed_small_val>(1.0 - sx*sx);
  fixed_small_val sy = x[1];
  fixed_small_val cy = sc_fixed_heron_sqrt<fixed_small_val>(1.0 - sy*sy);
  fixed_small_val sz = x[2];
  fixed_small_val cz = sc_fixed_heron_sqrt<fixed_small_val>(1.0 - sz*sz);
/*
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
  alignxf[15] = 1;*/

  alignxf_R[0]  = cy*cz;
  alignxf_R[1]  = sx*sy*cz + cx*sz;
  alignxf_R[2]  = -cx*sy*cz + sx*sz;
  alignxf_R[3]  = cy*sz;
  alignxf_R[4]  = -sx*sy*sz + cx*cz;
  alignxf_R[5]  = cx*sy*sz + sx*cz;
  alignxf_R[6]  = sy;
  alignxf_R[7]  = -sx*cy;
  alignxf_R[8] = cx*cy;
  alignxf_t[0] = centerSource[0] - alignxf_R[0]*centerTarget[0] - alignxf_R[3]*centerTarget[1] - alignxf_R[6]*centerTarget[2];
  alignxf_t[1] = centerSource[1] - alignxf_R[1]*centerTarget[0] - alignxf_R[4]*centerTarget[1] - alignxf_R[7]*centerTarget[2];
  alignxf_t[2] = centerSource[2] - alignxf_R[2]*centerTarget[0] - alignxf_R[5]*centerTarget[1] - alignxf_R[8]*centerTarget[2];

  return error; // aktueller Fehler
}
