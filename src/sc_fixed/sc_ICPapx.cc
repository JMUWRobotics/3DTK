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
double sc_ICPapx::Align(const std::vector<sc_PtPair>& Pairs,
                        double *alignxf,
                        const double centroid_m[3],
                        const double centroid_d[3])
{
  int n = Pairs.size();

  // ?!? <= 3
  if (n <= 3) {
    M4identity(alignxf);
    return 0;
  }

  int i;

  double A[3][3];
  double B[3];
  memset(&A[0][0], 0, 9 * sizeof(double));
  memset(&B[0], 0, 3 * sizeof(double));

  double sum = 0;
  double p1[3], p2[3];

  for (i = 0; i < n; i++) {
    p1[0] = Pairs[i].p1.x;
    p1[1] = Pairs[i].p1.y;
    p1[2] = Pairs[i].p1.z;
    p2[0] = Pairs[i].p2.x;
    p2[1] = Pairs[i].p2.y;
    p2[2] = Pairs[i].p2.z;

    double p12[3] = { p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2] };
    double p2c[3] = { p2[0] - centroid_d[0], p2[1] - centroid_d[1],
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

  double error = sqrt(sum / n);
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
  double diag[3];
  if (!choldc(A, diag)) {
    printf("Couldn't find transform.\n");
    return -1.0;
  }
  double x[3];
  cholsl(A, diag, B, x);

  // Interpret results
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
  alignxf[12] = centroid_m[0] - alignxf[0]*centroid_d[0] -
    alignxf[4]*centroid_d[1] - alignxf[8]*centroid_d[2];
  alignxf[13] = centroid_m[1] - alignxf[1]*centroid_d[0] -
    alignxf[5]*centroid_d[1] - alignxf[9]*centroid_d[2];
  alignxf[14] = centroid_m[2] - alignxf[2]*centroid_d[0] -
    alignxf[6]*centroid_d[1] - alignxf[10]*centroid_d[2];
  alignxf[15] = 1;

  return error;
}


double sc_ICPapx::Align_Parallel(const int openmp_num_threads,
                                 const unsigned int n[OPENMP_NUM_THREADS],
                                 const double sum[OPENMP_NUM_THREADS],
                                 const double centroid_m[OPENMP_NUM_THREADS][3],
                                 const double centroid_d[OPENMP_NUM_THREADS][3],
                                 const std::vector<sc_PtPair> pairs[OPENMP_NUM_THREADS],
                                 double *alignxf)

{

#ifdef _OPENMP

  double At[OPENMP_NUM_THREADS][3][3];
  double Bt[OPENMP_NUM_THREADS][3];

  for (int j=0;j < OPENMP_NUM_THREADS; j++)
    for (int k = 0;k < 3; k++) {
      for (int l = 0; l < 3; l++)
        At[j][k][l] = 0.0;
      Bt[j][k] = 0.0;
    }
  double A[3][3];
  double B[3];
  memset(&A[0][0], 0, 9 * sizeof(double));
  memset(&B[0], 0, 3 * sizeof(double));

  double s = 0.0;
  double error;
  unsigned int pairs_size = 0;
  double cm[3] = {0.0, 0.0, 0.0};  // centroid m
  double cd[3] = {0.0, 0.0, 0.0};  // centroid d

  double cms[3], cds[3];
  for (int i = 0; i < 3; i++) cms[i] = cds[i] = 0.0;

  for (int i = 0; i < openmp_num_threads; i++) {
    s += sum[i];
    pairs_size += n[i];

    // compute centroids for all the pairs
    cm[0] += n[i] * centroid_m[i][0];
    cm[1] += n[i] * centroid_m[i][1];
    cm[2] += n[i] * centroid_m[i][2];
    cd[0] += n[i] * centroid_d[i][0];
    cd[1] += n[i] * centroid_d[i][1];
    cd[2] += n[i] * centroid_d[i][2];

    cms[0] += centroid_m[i][0];
    cms[1] += centroid_m[i][1];
    cms[2] += centroid_m[i][2];
    cds[0] += centroid_d[i][0];
    cds[1] += centroid_d[i][1];
    cds[2] += centroid_d[i][2];
  }

  cm[0] /= pairs_size;
  cm[1] /= pairs_size;
  cm[2] /= pairs_size;
  cd[0] /= pairs_size;
  cd[1] /= pairs_size;
  cd[2] /= pairs_size;

  error = sqrt(s / (double)pairs_size);

#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    for (unsigned int i = 0 ; i < (unsigned int)pairs[thread_num].size() ; i++)
      {
        At[thread_num][0][0] +=
          (pairs[thread_num][i].p2.y - cd[1]) *
          (pairs[thread_num][i].p2.y - cd[1]) +
          (pairs[thread_num][i].p2.z - cd[2]) *
          (pairs[thread_num][i].p2.z - cd[2]);
        At[thread_num][0][1] -=
          (pairs[thread_num][i].p2.x - cd[0]) *
          (pairs[thread_num][i].p2.y - cd[1]);
        At[thread_num][0][2] -=
          (pairs[thread_num][i].p2.x - cd[0]) *
          (pairs[thread_num][i].p2.z - cd[2]);
        At[thread_num][1][1] +=
          (pairs[thread_num][i].p2.x - cd[0]) *
          (pairs[thread_num][i].p2.x - cd[0]) +
          (pairs[thread_num][i].p2.z - cd[2]) *
          (pairs[thread_num][i].p2.z - cd[2]);
        At[thread_num][1][2] -=
          (pairs[thread_num][i].p2.y - cd[1]) *
          (pairs[thread_num][i].p2.z - cd[2]);
        At[thread_num][2][2] +=
          (pairs[thread_num][i].p2.x - cd[0]) *
          (pairs[thread_num][i].p2.x - cd[0]) +
          (pairs[thread_num][i].p2.y - cd[1]) *
          (pairs[thread_num][i].p2.y - cd[1]);

        Bt[thread_num][0] +=
          (pairs[thread_num][i].p1.z - pairs[thread_num][i].p2.z) *
          (pairs[thread_num][i].p2.y - cd[1]) -
          (pairs[thread_num][i].p1.y - pairs[thread_num][i].p2.y) *
          (pairs[thread_num][i].p2.z - cd[2]);
        Bt[thread_num][1] +=
          (pairs[thread_num][i].p1.x - pairs[thread_num][i].p2.x) *
          (pairs[thread_num][i].p2.z - cd[2]) -
          (pairs[thread_num][i].p1.z - pairs[thread_num][i].p2.z) *
          (pairs[thread_num][i].p2.x - cd[0]);
        Bt[thread_num][2] +=
          (pairs[thread_num][i].p1.y - pairs[thread_num][i].p2.y) *
          (pairs[thread_num][i].p2.x - cd[0]) -
          (pairs[thread_num][i].p1.x - pairs[thread_num][i].p2.x) *
          (pairs[thread_num][i].p2.y - cd[1]);
      }
  }

  for (int j = 0;j < OPENMP_NUM_THREADS; j++)
    for (int k = 0; k < 3; k++) {
      for (int l = 0; l < 3; l++)
        A[k][l] += At[j][k][l] ;
      B[k] += Bt[j][k];
    }

  // continue with linear solution

  if (!quiet) {
    std::cout.setf(std::ios::basefield);
    std::cout << "PAPX RMS point-to-point error = "
         << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
         << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
         << std::setw(10) << std::setprecision(7)
         << error
         << "  using " << std::setw(6) << pairs_size << " points" << std::endl;
  }

  // Solve eqns
  double diag[3];
  if (!choldc(A, diag)) {
    printf("Couldn't find transform.\n");
    return -1.0;
  }
  double x[3];
  cholsl(A, diag, B, x);

  // Interpret results
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
  alignxf[12] = cm[0] - alignxf[0]*cd[0] - alignxf[4]*cd[1] - alignxf[8]*cd[2];
  alignxf[13] = cm[1] - alignxf[1]*cd[0] - alignxf[5]*cd[1] - alignxf[9]*cd[2];
  alignxf[14] = cm[2] - alignxf[2]*cd[0] - alignxf[6]*cd[1] - alignxf[10]*cd[2];
  alignxf[15] = 1;

  return error;
#else
  std::cout << "Point_Point_Align_Parallel:"<< std::endl
       << "Please compile with OpenMP support to use this function" << std::endl;
  exit(-1);
#endif
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
