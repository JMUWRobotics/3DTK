/*
 * icp implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of 3D scan matching with ICP
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of
 * Wuerzburg, Germany
 */

#include "sc_fixed/sc_ICP.h"

// #include "slam6d/metaScan.h"			//brauchen wir nicht
#include <iomanip>

#include "sc_fixed/sc_fixed_math.h"
#include "slam6d/globals.icc"  //brauchen wir für M4identity
using std::cerr;

#include <string.h>

/**
 * Constructor
 *
 * @param my_sc_ICPminimizer Pointer to the ICP-minimizer
 * @param max_dist_match Maximum distance to which point pairs are collected
 * @param max_num_iterations Maximum number of iterations
 * @param quiet Whether to print to the standard output
 * @param meta Match against a meta scan?
 * @param rnd Randomized point selection
 * @param eP Extrapolate odometry?
 * @param anim Animate which frames?
 * @param epsilonICP Termination criterion
 * @param nns_method Selects NNS method to be used
 */
sc_ICP::sc_ICP(sc_ICPminimizer* my_sc_ICPminimizer, double max_dist_match,
               int max_num_iterations, bool quiet, bool meta, int rnd, bool eP,
               int anim, double epsilonICP, int nns_method, bool cuda_enabled,
               bool cad_matching, int max_num_metascans) {
  this->my_sc_ICPminimizer = my_sc_ICPminimizer;
  this->anim = anim;
  this->cuda_enabled = cuda_enabled;
  this->nns_method = nns_method;

  if (!quiet) {
    cout << "Maximal distance match      : " << max_dist_match << endl
         << "Maximal number of iterations: " << max_num_iterations << endl
         << endl;
  }

  // checks
  if (max_dist_match < 0.0) {
    cerr << "ERROR [SC_ICP]: first parameter (max_dist_match) "
         << "has to be >= 0," << endl;
    exit(1);
  }
  if (max_num_iterations < 0) {
    cerr << "ERROR [SC_ICP]: second parameter (max_num_iterations)"
         << "has to be >= 0." << endl;
    exit(1);
  }

  this->max_dist_match2 = sqr(max_dist_match);
  this->max_num_iterations = max_num_iterations;
  this->quiet = quiet;
  this->meta = meta;
  this->rnd = rnd;
  this->eP = eP;
  this->epsilonICP = epsilonICP;

  // Set initial seed (for "real" random numbers)
  //  srand( (unsigned)time( NULL ) );
  this->cad_matching = cad_matching;

  this->max_num_metascans = max_num_metascans;

  // set the number of point pairs to zero
  nr_pointPair = 0;
}

// match Methode mit konvertiertem Datentyp als Übergabeparameter
// fehlende includes Array und vector, Methode
// match überladen sollte funktionieren
int sc_ICP::match(std::vector<std::array<f_float, 3>>& source,  
                  std::vector<std::array<f_float, 3>>& target, std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf) {
  
  std::vector<std::array<f_float, 3>> matchedTarget;
  // TO-DO validierung der übergebenen Listen
  // nns brute forece
  for (const auto& src : source) {
    //std::cout << "source iteration" << std::endl;
    f_float minDist;
    bool first = true;
    std::array<f_float, 3> closest;

    for (const auto& tgt : target) {
      //std::cout << ".";
      f_float dx = src[0] - tgt[0];
      f_float dy = src[1] - tgt[1];
      f_float dz = src[2] - tgt[2];
      f_float dist = dx * dx + dy * dy + dz * dz;

      if (first) {
	//std::cout << "if first" << std::endl;
        minDist = dist;
        closest = tgt;
        first = false;
      } else if (dist < minDist) {
	//std::cout << "else if first" <<std::endl;
        minDist = dist;
        closest = tgt;
      }
    }
    //std::cout << std::endl;
    //std::cout << "before matchedTarget.push" << std::endl;
    matchedTarget.push_back(closest);
    //std::cout << "after matchedTarget.push" << std::endl;
  }

  // TO-DO Test matchedTarget

  // Schwerpunkte bestimmen
  std::cout << "Schwerpunkte" << std::endl;
  
  std::array<f_float, 3> centerSource = {0, 0, 0};
  std::array<f_float, 3> centerTarget = {0, 0, 0};

  // size_T Rückgabewert von source.size(); Include sollte in array dabei sein, sonst include von <cstddef>.

  if (source.size() != matchedTarget.size()) {
    std::cerr << "warning: source size does not match target size" << std::endl;
    if (source.size() > matchedTarget.size()) {
      std::cerr << "source size is greater than target size" << std::endl;
    }
  }

  size_t count = std::min(source.size(), matchedTarget.size());
  
  for (size_t i = 0; i < count; ++i) {
    centerSource[0] += source[i][0];
    centerSource[1] += source[i][1];
    centerSource[2] += source[i][2];

    centerTarget[0] += matchedTarget[i][0];
    centerTarget[1] += matchedTarget[i][1];
    centerTarget[2] += matchedTarget[i][2];
  }
  std::cout << "Schwerpunkte Ende" << std::endl;

  // typeCast wohl nötig, genauigkeit?? Trotzdem mal wegen möglichem Typkonflikt fragen
  std::cout << "type cast" << std::endl;
  f_float srcSize = static_cast<f_float>(source.size());
  f_float trgSize = static_cast<f_float>(matchedTarget.size());
  centerSource[0] /= srcSize;
  centerSource[1] /= srcSize;
  centerSource[2] /= srcSize;

  centerTarget[0] /= trgSize;
  centerTarget[1] /= trgSize;
  centerTarget[2] /= trgSize;

  // Rotation und Translation berechnen
  f_float alignxf[16];
  /*int ret = */ my_sc_ICPminimizer->Align(source, matchedTarget, alignxf, centerSource, centerTarget);
  for(int i = 0; i < 16; i++) {
    std::cout << alignxf[i] << " ";
  }
  std::cout << std::endl;
  
  // TODO Transformation berechnen
  //std::array<f_float, 16> transMat;
  //transMat[0] = transMat[5] = transMat[10] = transMat[15] = 1;
  //std::array<f_float, 16> dalignxf;
  //dalignxf[0] = dalignxf[5] = dalignxf[10] = dalignxf[15] = 1;
  transform(target, alignxf, transMat, dalignxf, 0);
  //  4 x 4 Matix auf Konsole ausgeben
  return 0;  // Rückgabewert int für Iterationen, vielleicht langfristig auf 4x4
             // Matrix ändern
}

/**
 * Matches a 3D Scan against a 3D Scan
 * @param PreviousScan The scan or metascan forming the model
 * @param CurrentScan The current scan thas is to be matched
 * @return The number of iterations done in this matching run
 */
int sc_ICP::match2(std::vector<std::array<f_float, 3>>& PreviousScan, std::vector<std::array<f_float, 3>>& CurrentScan) {
  f_float id[16];
  M4identity(id);
  std::array<f_float, 16> transMat;
  std::array<f_float, 16> dalignxf;
  
  transform(CurrentScan, id, transMat, dalignxf, 0);
  // statt CurrentScan->transform(id, Scan::ICP, 0);  // write end pose
  
  // If ICP shall not be applied, then just write
  // the identity matrix and return
  if (max_num_iterations == 0) {
    return 0;
  }

  // icp main loop
  f_float ret = 0.0, prev_ret = 0.0, prev_prev_ret = 0.0;
  int iter = 0;
  f_float alignxf[16];
  long time = GetCurrentTimeInMilliSec();

  for (iter = 0; iter < max_num_iterations; iter++) {
    prev_prev_ret = prev_ret;
    prev_ret = ret;

    if (iter == 1) time = GetCurrentTimeInMilliSec();

    f_float centroid_m[3] = {0.0, 0.0, 0.0};
    f_float centroid_d[3] = {0.0, 0.0, 0.0};
    vector<sc_PtPair> pairs;

    // TODO: Punktepaare generieren
    // Scan::getPtPairs(&pairs, PreviousScan, CurrentScan, 0, rnd, max_dist_match2, ret, centroid_m, centroid_d, pairing_mode);

    // set the number of point paira
    nr_pointPair = pairs.size();

    // do we have enough point pairs?
    if (pairs.size() > 3) {
      ret = my_sc_ICPminimizer->Align(pairs, alignxf, centroid_m, centroid_d);
    } else {
      break;
    }

    if ((iter == 0 && anim != -2) || ((anim > 0) && (iter % anim == 0))) {
      // transform the current scan
      transform(CurrentScan, id, transMat, dalignxf, 0);
    } else {
      // transform the current scan
      transform(CurrentScan, id, transMat, dalignxf, -1);
    }

    if (((fabs(ret - prev_ret) < epsilonICP) &&
         (fabs(ret - prev_prev_ret) < epsilonICP)) ||
        (iter == max_num_iterations - 1)) {
      f_float id[16];
      M4identity(id);
      if (anim == -2) {
        // write end pose
        transform(CurrentScan, id, transMat, dalignxf, -1);
      } else {
        // write end pose
        transform(CurrentScan, id, transMat, dalignxf, 0);
      }
      break;
    }
  }

  long endtime = GetCurrentTimeInMilliSec() - time;
  cout << "TIME  " << endtime << "   ITER " << iter << endl;
  return iter;
}

/**
 * Computes the point to point error between two scans
 *
 *
 */
double sc_ICP::Point_Point_Error(Scan* PreviousScan, Scan* CurrentScan,
                                 double max_dist_match, unsigned int* np,
                                 double scale_max) {
  double scale = log(scale_max) / (max_dist_match * max_dist_match);
  double error = 0;
  unsigned int nr_ppairs = 0;

  double centroid_m[3] = {0.0, 0.0, 0.0};
  double centroid_d[3] = {0.0, 0.0, 0.0};
  vector<sc_PtPair> pairs;

  // TODO: ersetzen
  // Scan::getPtPairs(&pairs, PreviousScan, CurrentScan, 0, rnd, sqr(max_dist_match), error, centroid_m, centroid_d, CLOSEST_POINT);

  // getPtPairs computes error as sum of squared distances
  error = 0;

  for (unsigned int i = 0; i < pairs.size(); i++) {
    double dist = sqr(pairs[i].p1.x - pairs[i].p2.x) +
                  sqr(pairs[i].p1.y - pairs[i].p2.y) +
                  sqr(pairs[i].p1.z - pairs[i].p2.z);
    error -= 0.39894228 * exp(dist * scale);
  }
  nr_ppairs = pairs.size();

  if (np) *np = nr_ppairs;
  return error / nr_ppairs;
}

/**
 * This function matches the scans only with ICP
 *
 * @param allScans Contains all necessary scans.
 */
void sc_ICP::doICP(std::vector<std::vector<std::array<f_float, 3>>> allScans) {
  f_float id[16];
  M4identity(id);

  for (unsigned int i = 0; i < allScans.size(); i++) {
    cout << i << "*" << endl;

    std::vector<std::array<f_float, 3>> CurrentScan = allScans[i];
    std::vector<std::array<f_float, 3>> PreviousScan;

    if (i > 0) {
      PreviousScan = allScans[i - 1];
    //  if (eP) {  // extrapolate odometry - für uns überflüssig?
    //    CurrentScan->mergeCoordinatesWithRoboterPosition(PreviousScan);
    //  }
    }

    if (i > 0) {
    //TODO falls doICP nötig: Parameter match anpassen
    //  match(PreviousScan, CurrentScan);
    }
  }
}
