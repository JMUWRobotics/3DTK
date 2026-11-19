/*
 * sc_fixed_math implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of mathematical functions with SystemC fixpoint numbers (fixed_val's)
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer. Institute of Computer Science, University of Wuerzburg, Germany.
 */

#include "sc_fixed/sc_fixed_math.h"
#include "sc_fixed/sc_fixed_converter.h"




// function to compute cholesky decomposition
bool sc_choldc(fixed_val A[3][3], fixed_val diag[3]) {

  unsigned int N = 3;
  const fixed_val epsilon = (fixed_val) float(1e-3);

  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = i; j < N; j++) {
      fixed_val sum = A[i][j];
      for (int k = i-1; k >= 0; k--) {
        sum -= A[i][k] * A[j][k];
      }

      if (i == j) {
        if (sum < epsilon) {
          return false;
        }
        diag[i] = sc_fixed_heron_sqrt<fixed_val>(sum);
      }
      else {
        A[j][i] = sum / diag[i];
      }
    }
  }
  return true;
}

// function to solve cholesky
void sc_cholsl(fixed_val A[3][3], fixed_val diag[3], fixed_val B[3], fixed_small_val x[3]) {
  int N = 3;
  for (int i = 0; i < N; i++) {
    fixed_val sum = B[i];
    for (int k = i-1; k >= 0; k--) {
      sum -= A[i][k] * x[k];
    }
    x[i] = sum / diag[i];
  }
  for (int i = N-1; i >= 0; i--) {
    fixed_val sum = x[i];
    for (int k = i+1; k < N; k++) {
      sum -= A[k][i] * x[k];
    }
    x[i] = sum / diag[i];
  }
}

// transform method without any algorithm type, because ICP algorithm is only with APX and Brute Force
void transform(std::vector<std::array<fixed_val, 3>>& scan, fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::array<fixed_small_val, 9>& transMat_R, std::array<fixed_val, 3>& transMat_t , std::array<fixed_small_val, 9>& dalignxf_R, std::array<fixed_val, 3>& dalignxf_t, std::ofstream& frame, int islum){

  // transform points
  transformPoints(alignxf_R, alignxf_t, scan);

  // update matrices
  transformMatrix(alignxf_R, alignxf_t, transMat_R, transMat_t, dalignxf_R, dalignxf_t);

  // speichere Transformation in .frames-Datei (falls islum == 0 statt -1)
  if (islum == 0) {
    writeFrame(frame, transMat_R, transMat_t, 1);
  }
}

// Internal function of transform which alters the points
void transformPoints(const fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::vector<std::array<fixed_val, 3>>& scan){
  for (size_t i = 0; i < scan.size(); ++i) {
    transform3(alignxf_R, alignxf_t, scan[i]);
  }
}

// Internal function of transformPoints which alters a point
void transform3(const fixed_small_val alignxf_R[9], fixed_val alignxf_t[3], std::array<fixed_val, 3>& point){
  fixed_val x_neu, y_neu, z_neu;
  x_neu = point[0] * alignxf_R[0] + point[1] * alignxf_R[3] + point[2] * alignxf_R[6];
  y_neu = point[0] * alignxf_R[1] + point[1] * alignxf_R[4] + point[2] * alignxf_R[7];
  z_neu = point[0] * alignxf_R[2] + point[1] * alignxf_R[5] + point[2] * alignxf_R[8];
  point[0] = x_neu + alignxf_t[0];
  point[1] = y_neu + alignxf_t[1];
  point[2] = z_neu + alignxf_t[2];
}


// Internal function of transform which handles the matrices
void transformMatrix(const fixed_small_val alignxf_R[9], const fixed_val alignxf_t[3], std::array<fixed_small_val, 9>& transMat_R, std::array<fixed_val, 3>& transMat_t, std::array<fixed_small_val, 9>& dalignxf_R, std::array<fixed_val, 3>& dalignxf_t){
  std::array<fixed_small_val, 9> tempxf_R;
  std::array<fixed_val, 3> tempxf_t;


  // apply alignxf_R to transMat_R and update pose vectors, copy to transMat_R
  MMult(alignxf_R, transMat_R, tempxf_R);
  std::copy(tempxf_R.begin(), tempxf_R.end(), transMat_R.begin());

  // apply alignxf_R to dalignxf_R, copy to dalignxf_R
  MMult(alignxf_R, dalignxf_R, tempxf_R);
  std::copy(tempxf_R.begin(), tempxf_R.end(), dalignxf_R.begin());

  // apply alignxf_t to transMat_t: t_new = R_align * t_old + t_align
  for(unsigned int i = 0; i < 3; ++i){
    tempxf_t[i] = alignxf_R[i] * transMat_t[0] + alignxf_R[i+ 3] * transMat_t[1] + alignxf_R[i + 6] * transMat_t[2] + alignxf_t[i];
  }
  std::copy(tempxf_t.begin(), tempxf_t.end(), transMat_t.begin());


  // apply alignxf_t to dalignxf_t
  for(unsigned int i = 0; i < 3; ++i){
    tempxf_t[i] = alignxf_R[i] * dalignxf_t[0] + alignxf_R[i+ 3] * dalignxf_t[1] + alignxf_R[i + 6] * dalignxf_t[2] + alignxf_t[i];
  }
    std::copy(tempxf_t.begin(), tempxf_t.end(), dalignxf_t.begin());


} 


// Matrixmultiplikation von M1 und M2, Ergebnis wird in Mout gespeichert
void MMult(const fixed_small_val M1[9], const std::array<fixed_small_val, 9>& M2, std::array<fixed_small_val, 9>& Mout){
  Mout[0] = M1[0]*M2[0]+M1[3]*M2[1]+M1[6]*M2[2];
  Mout[1] = M1[1]*M2[0]+M1[4]*M2[1]+M1[7]*M2[2];
  Mout[2] = M1[2]*M2[0]+M1[5]*M2[1]+M1[8]*M2[2];

  Mout[3] = M1[0]*M2[3]+M1[3]*M2[4]+M1[6]*M2[5];
  Mout[4] = M1[1]*M2[3]+M1[4]*M2[4]+M1[7]*M2[5];
  Mout[5] = M1[2]*M2[3]+M1[5]*M2[4]+M1[8]*M2[5];

  Mout[6] = M1[0]*M2[6]+M1[3]*M2[7]+M1[6]*M2[8];
  Mout[7] = M1[1]*M2[6]+M1[4]*M2[7]+M1[7]*M2[8];
  Mout[8] = M1[2]*M2[6]+M1[5]*M2[7]+M1[8]*M2[8];

}

// Absolutbetrag einer fixed-point Zahl
fixed_small_val sc_abs(fixed_small_val x) {
  return (x < fixed_small_val(0)) ? fixed_small_val(-x) : x;
}
