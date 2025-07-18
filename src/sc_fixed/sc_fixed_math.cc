#include "sc_fixed/sc_fixed_math.h"

f_float sc_fixed_heron_sqrt(f_float s) {

  //std::cout << "sqrt(" << s << ")" << std::endl;

  if (s == 0) {
    return 0;
  }
  if (s < 0) {
    std::cerr << "warning: square root from negative numbers not allowed --> will return 0" << std::endl;
    return 0;
  }

  f_float x = s / 2 + 1;
  f_float x_n;
  
  for(int i = FIXED_HERON_ITERATIONS; i > 0; i--) {
    x_n = (x + s / x) / 2;
    //std::cout << "Iteration " << i << ": x = " << x_n.to_double() << std::endl;
    x = x_n;
  }

  return x_n;
  
}

bool sc_choldc(f_float A[3][3], f_float diag[3]) {

  unsigned int N = 3;
  const f_float epsilon = (f_float) float(1e-3);  
  
  //std::cout << "Matrix A" << std::endl;
  //for (int i = 0; i < 3; ++i) {
  //  for (int j = 0; j < 3; ++j) {
  //    std::cout << A[i][j] << " ";
  //  }
  //  std::cout << std::endl;
  //}
  //std::cout << "Diagonalmatrix" << std::endl; //ist hier wohl immer die Nullmatrix
  //for (int i = 0; i < 3; ++i) {
  //  std::cout << diag[i] << " " << std::endl;
  //}
  
  for (unsigned int i = 0; i < N; i++) {
    for (unsigned int j = i; j < N; j++) {
      f_float sum = A[i][j];
      //for (unsigned int k = 0; k < i; k++) {
      for (int k=i-1; k >= 0; k--) {
        sum -= A[i][k] * A[j][k];
      }
      if (i == j) {
        //std::cout << "Summe: " << sum << " / Iteration " << i << std::endl;
        //std::cout << "Epsilon: " << epsilon << std::endl;
        if (sum < epsilon) {
          //std::cout << "Fehler: Summe < Epsilon" << std::endl;
          return false;
        }
        //std::cout << "Wurzel wird berechnet" << std::endl;
        diag[i] = sc_fixed_heron_sqrt(sum);
      }
      else {
        A[j][i] = sum / diag[i];
      }
    }
  }
  return true;
}

void sc_cholsl(f_float A[3][3], f_float diag[3], f_float B[3], f_float x[3]) {
  int N = 3;
  for (int i=0; i < N; i++) {
    f_float sum = B[i];
    for (int k=i-1; k >= 0; k--) {
      sum -= A[i][k] * x[k];
    }
    x[i] = sum / diag[i];
  }
  for (int i=N-1; i >= 0; i--) {
    f_float sum = x[i];
    for (int k=i+1; k < N; k++) {
      sum -= A[k][i] * x[k];
    }
    x[i] = sum / diag[i];
  }
}

/**
 * Transforms the scan by a given transformation and writes a new frame.
 * The idea is to write for every transformation in all files,
 * such that the show program is able to determine,
 * which scans have to be drawn in which color. Hidden scans
 * (or later processed scans) are written with INVALID.
 *
 * @param alignxf Transformation matrix
 * @param islum Is the transformtion part of LUM, i.e., all scans
 *              are transformed?
 *              In this case only LUM transformation is stored, otherwise all
 *              scans are processed
 *        -1  no transformation is stored
 *         0  ICP transformation
 */
void transform(std::vector<std::array<f_float, 3>>& scan, f_float alignxf[16], std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf, int islum){
//für jeden Scan haben wir die f_float transMat[16] = (4x4) transformation matrix
//und dalignxf[16] transformation represents the delta transformation virtually applied to the tree and is used to compute are actual corresponding points.

//algoType wird weggelassen, da wir nur ICP haben
  // transform points
  transformPoints(alignxf, scan);

  // update matrices
  transformMatrix(alignxf, transMat, dalignxf);
  
  std::cout << "Transformationsmatrix:" << std::endl;
  for(int i = 0; i < 16; i++) {
    std::cout << transMat[i] << " ";
  }
  std::cout << std::endl;
  
  std::cout << "dalignxf-Matrix:" << std::endl;
  for(int i = 0; i < 16; i++) {
    std::cout << dalignxf[i] << " ";
  }
  std::cout << std::endl;
  // store transformation in frames - falls islum=0 statt -1
  // islum wird nicht benötigt, da Schreiben nur am Ende?
}

//! Internal function of transform which alters the points
void transformPoints(const f_float alignxf[16], std::vector<std::array<f_float, 3>>& scan){
  
  for (size_t i = 0; i < scan.size(); ++i) {
    transform3(alignxf, scan[i]);
  }

}

void transform3(const f_float alignxf[16], std::array<f_float, 3>& point){
  f_float x_neu, y_neu, z_neu;
  x_neu = point[0] * alignxf[0] + point[1] * alignxf[4] + point[2] * alignxf[8];
  y_neu = point[0] * alignxf[1] + point[1] * alignxf[5] + point[2] * alignxf[9];
  z_neu = point[0] * alignxf[2] + point[1] * alignxf[6] + point[2] * alignxf[10];
  point[0] = x_neu + alignxf[12];
  point[1] = y_neu + alignxf[13];
  point[2] = z_neu + alignxf[14];
}

//! Internal function of transform which handles the matrices
void transformMatrix(const f_float alignxf[16], std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf){
  std::array<f_float, 16> tempxf;

  // apply alignxf to transMat and update pose vectors
  MMult(alignxf, transMat, tempxf);
  std::copy(tempxf.begin(), tempxf.end(), transMat.begin());
  //std::memcpy(transMat.data(), tempxf.data(), 16 * sizeof(f_float));

  // apply alignxf to dalignxf
  MMult(alignxf, dalignxf, tempxf);
  // kopiere - nach dalignfx - von tempxf - anzahlBytes: sizeof(transMat), definiert in <cstring>
  //std::memcpy(dalignxf.data(), tempxf.data(), 16 * sizeof(f_float));
  std::copy(tempxf.begin(), tempxf.end(), dalignxf.begin());
}

void MMult(const f_float M1[16], const std::array<f_float, 16>& M2, std::array<f_float, 16>& Mout){
  Mout[ 0] = M1[ 0]*M2[ 0]+M1[ 4]*M2[ 1]+M1[ 8]*M2[ 2]+M1[12]*M2[ 3];
  Mout[ 1] = M1[ 1]*M2[ 0]+M1[ 5]*M2[ 1]+M1[ 9]*M2[ 2]+M1[13]*M2[ 3];
  Mout[ 2] = M1[ 2]*M2[ 0]+M1[ 6]*M2[ 1]+M1[10]*M2[ 2]+M1[14]*M2[ 3];
  Mout[ 3] = M1[ 3]*M2[ 0]+M1[ 7]*M2[ 1]+M1[11]*M2[ 2]+M1[15]*M2[ 3];
  Mout[ 4] = M1[ 0]*M2[ 4]+M1[ 4]*M2[ 5]+M1[ 8]*M2[ 6]+M1[12]*M2[ 7];
  Mout[ 5] = M1[ 1]*M2[ 4]+M1[ 5]*M2[ 5]+M1[ 9]*M2[ 6]+M1[13]*M2[ 7];
  Mout[ 6] = M1[ 2]*M2[ 4]+M1[ 6]*M2[ 5]+M1[10]*M2[ 6]+M1[14]*M2[ 7];
  Mout[ 7] = M1[ 3]*M2[ 4]+M1[ 7]*M2[ 5]+M1[11]*M2[ 6]+M1[15]*M2[ 7];
  Mout[ 8] = M1[ 0]*M2[ 8]+M1[ 4]*M2[ 9]+M1[ 8]*M2[10]+M1[12]*M2[11];
  Mout[ 9] = M1[ 1]*M2[ 8]+M1[ 5]*M2[ 9]+M1[ 9]*M2[10]+M1[13]*M2[11];
  Mout[10] = M1[ 2]*M2[ 8]+M1[ 6]*M2[ 9]+M1[10]*M2[10]+M1[14]*M2[11];
  Mout[11] = M1[ 3]*M2[ 8]+M1[ 7]*M2[ 9]+M1[11]*M2[10]+M1[15]*M2[11];
  Mout[12] = M1[ 0]*M2[12]+M1[ 4]*M2[13]+M1[ 8]*M2[14]+M1[12]*M2[15];
  Mout[13] = M1[ 1]*M2[12]+M1[ 5]*M2[13]+M1[ 9]*M2[14]+M1[13]*M2[15];
  Mout[14] = M1[ 2]*M2[12]+M1[ 6]*M2[13]+M1[10]*M2[14]+M1[14]*M2[15];
  Mout[15] = M1[ 3]*M2[12]+M1[ 7]*M2[13]+M1[11]*M2[14]+M1[15]*M2[15];
}


f_float sc_abs(f_float x) {
  return (x < f_float(0)) ? f_float(-x) : x;
}
