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
 *  @brief Implementation of the ICP error function minimization via
           small angle approximation
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#include "sc_fixed/sc_ICPapx.h"

#include "slam6d/globals.icc"		//TODO: durch mathematische Fixpoint-Dateien ersetzen
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
 * @param Pairs Vector of point pairs (pairs of corresponding points)
 * @param alignxf The resulting transformation matrix
 * @return Error estimation of the matching (rms)
 */
void sc_ICPapx::Align(const std::vector<std::array<f_float, 3>>& source,
                      const std::vector<std::array<f_float, 3>>& matchedTarget,
                      f_float *alignxf,
                      const std::array<f_float, 3> centerSource,
                      const std::array<f_float, 3> centerTarget){

  f_float A[3][3];
  f_float B[3]; 

  for(int i=0; i<3;i++) {
    B[i] = f_float(0.0);
    for(int j=0; j<3;j++){
      A[i][j] = f_float(0.0); 
    }
  }	

  f_float sum = 0; 
  f_float p1[3], p2[3]; 

  for(size_t i = 0; i < source.size(); ++i){
    const auto& p1 = source[i];
    const auto& p2 = matchedTarget[i];  

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

  f_float diag[3];
 
  if (!sc_choldc(A, diag)) {
    printf("Couldn't find transform.\n");
    return; //-1.0;
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
  alignxf[12] = centerTarget[0] - alignxf[0]*centerSource[0] - alignxf[4]*centerSource[1] - alignxf[8]*centerSource[2];
  alignxf[13] = centerTarget[1] - alignxf[1]*centerSource[0] - alignxf[5]*centerSource[1] - alignxf[9]*centerSource[2];
  alignxf[14] = centerTarget[2] - alignxf[2]*centerSource[0] - alignxf[6]*centerSource[1] - alignxf[10]*centerSource[2];
  alignxf[15] = 1;

}

double sc_ICPapx::Align(const std::vector<sc_PtPair>& Pairs,
                        f_float *alignxf,
                        const f_float centroid_m[3],
                        const f_float centroid_d[3])
{
  int n = Pairs.size();

  // ?!? <= 3
  if (n <= 3) {
    M4identity(alignxf);
    return 0;
  }

  int i;

  f_float A[3][3];
  f_float B[3];
  memset(&A[0][0], 0, 9 * sizeof(f_float));
  memset(&B[0], 0, 3 * sizeof(f_float));

  f_float sum = 0;
  f_float p1[3], p2[3];

  for (i = 0; i < n; i++) {
    p1[0] = Pairs[i].p1.x;
    p1[1] = Pairs[i].p1.y;
    p1[2] = Pairs[i].p1.z;
    p2[0] = Pairs[i].p2.x;
    p2[1] = Pairs[i].p2.y;
    p2[2] = Pairs[i].p2.z;

    f_float p12[3] = { p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] };
    f_float p2c[3] = { p2[0] - centroid_d[0], p2[1] - centroid_d[1],
                      p2[2] - centroid_d[2] };

    sum += Len2(p12);

    B[0] += (p12[2]*p2c[1] - p12[1]*p2c[2]);
    B[1] += (p12[0]*p2c[2] - p12[2]*p2c[0]);
    B[2] += (p12[1]*p2c[0] - p12[0]*p2c[1]);
    A[0][0] += (sqr(p2c[1]) + sqr(p2c[2]));
    A[0][1] -= p2c[0] * p2c[1];
    A[0][2] -= p2c[0] * p2c[2];
    A[1][1] += (sqr(p2c[0]) + sqr(p2c[2]));
    A[1][2] -= p2c[1] * p2c[2];
    A[2][2] += (sqr(p2c[0]) + sqr(p2c[1]));
  }

  f_float error = sqrt(sum / n);
  if (!quiet) {
    std::cout.setf(std::ios::basefield);
    std::cout << "APX RMS point-to-point error = "
         << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
         << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
         << std::setw(10) << std::setprecision(7)
         << error
         << "  using " << std::setw(6) << (int)Pairs.size()
         << " points" << std::endl;
  }

  // Solve eqns
  f_float diag[3];
  if (!sc_choldc(A, diag)) {
    printf("Couldn't find transform.\n");
    return -1.0;
  }
  f_float x[3];
  sc_cholsl(A, diag, B, x);

  // Interpret results
  f_float sx = x[0];
  f_float cx = sqrt(1.0 - sx*sx);
  f_float sy = x[1];
  f_float cy = sqrt(1.0 - sy*sy);
  f_float sz = x[2];
  f_float cz = sqrt(1.0 - sz*sz);

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
  alignxf[12] = centroid_m[0] - alignxf[0]*centroid_d[0] -
    alignxf[4]*centroid_d[1] - alignxf[8]*centroid_d[2];
  alignxf[13] = centroid_m[1] - alignxf[1]*centroid_d[0] -
    alignxf[5]*centroid_d[1] - alignxf[9]*centroid_d[2];
  alignxf[14] = centroid_m[2] - alignxf[2]*centroid_d[0] -
    alignxf[6]*centroid_d[1] - alignxf[10]*centroid_d[2];
  alignxf[15] = 1;

  return error;
}


void sc_ICPapx::computeRt(const double *x, const double *dx, double *alignxf)
{
  double sx = x[0];
  double cx = sqrt(1.0 - sx*sx);
  double sy = x[1];
  double cy = sqrt(1.0 - sy*sy);
  double sz = x[2];
  double cz = sqrt(1.0 - sz*sz);

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
  alignxf[12] = dx[0];
  alignxf[13] = dx[1];
  alignxf[14] = dx[2];
  alignxf[15] = 1;
}
