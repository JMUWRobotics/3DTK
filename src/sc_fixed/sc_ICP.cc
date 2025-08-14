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

#include <iomanip>
#include <iostream>

#include "sc_fixed/sc_fixed_math.h"
#include "slam6d/globals.icc"  //nötig für M4identity
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
 * @param epsilonICPexp Exponent for termination criterion
 * @param nns_method Selects NNS method to be used
 */
sc_ICP::sc_ICP(sc_ICPminimizer* my_sc_ICPminimizer, double max_dist_match,
               int max_num_iterations, bool quiet, bool meta, int rnd, bool eP,
               int anim, int epsilonICPexp, int nns_method, bool cuda_enabled,
               bool cad_matching, int max_num_metascans) {
  this->my_sc_ICPminimizer = my_sc_ICPminimizer;
  this->anim = anim;
  this->cuda_enabled = cuda_enabled;
  this->nns_method = nns_method;

  if (!quiet) {
    cout << "Maximal distance match      : " << max_dist_match << endl
         << "Maximal number of iterations: " << max_num_iterations << endl;
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

  this->max_dist_match2 = (max_dist_match*max_dist_match);
  this->max_num_iterations = max_num_iterations;
  this->quiet = quiet;
  this->meta = meta;
  this->rnd = rnd;
  this->eP = eP;
  this->epsilonICP = f_float(std::pow(10.0, -epsilonICPexp));
  
  if(!quiet){
    cout << "Epsilon                     : " << epsilonICP << endl << endl;
  }
  
  // Set initial seed (for "real" random numbers)
  //  srand( (unsigned)time( NULL ) );
  this->cad_matching = cad_matching;

  this->max_num_metascans = max_num_metascans;

  // set the number of point pairs to zero
  nr_pointPair = 0;
}

// match-Methode mit konvertiertem Datentyp als Übergabeparameter
int sc_ICP::match(std::vector<std::array<f_float, 3>>& source, std::vector<std::array<f_float, 3>>& target, std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf, std::ofstream& frame) {
  f_float id[16];
  M4identity(id);
  
  transform(target, id, transMat, dalignxf, frame, 0);
  
  // wenn ICP nicht angewendet werden soll, nur die Identitätsmatrix schreiben und return 0
  if (max_num_iterations == 0) {
    return 0;
  }
  
  f_float alignxf[16];
  f_float ret = 0.0, prev_ret = 0.0, prev_prev_ret = 0.0;
  int iter = 0;
  
  //ICP main loop
  for (iter = 0; iter < max_num_iterations; iter++) {
    std::cout << std::endl;
    std::cout << "*** ITERATION " << iter << " ***" << std::endl;
    
    // nns Brute Force: finde zu jedem Punkt aus target (data) den nächstgelegenen Punkt aus source (model)
    std::vector<std::array<f_float, 3>> matchedTarget;
    std::vector<std::array<f_float, 3>> matchedSource;
    for (size_t j = 0; j < target.size(); j++){
      std::array<f_float, 3> tgt = target[j];
      
      f_float minDist;
      bool first = true;
      std::array<f_float, 3> closest;
      
      for (size_t i = 0; i < source.size(); i++){
        std::array<f_float, 3> src = source[i];
      
        f_float dx = tgt[0] - src[0];
        f_float dy = tgt[1] - src[1];
        f_float dz = tgt[2] - src[2];
        f_float dist = (dx * dx) + (dy * dy) + (dz * dz);

        if (first) {
	  //std::cout << "if first" << std::endl;
          minDist = dist;
          closest = src;
          first = false;
        } else if (dist < minDist) {
	  //std::cout << "else if first" <<std::endl;
          minDist = dist;
          closest = src;
        }
      }
      //prüfe, ob das Punktpaar überhaupt in die Wertung eingehen soll (aus Distanzgründen)
      if (minDist <= max_dist_match2) {
        //std::cout << "before matchedTarget.push" << std::endl;
        matchedTarget.push_back(tgt);
        matchedSource.push_back(closest);
        //std::cout << "after matchedTarget.push" << std::endl;
      }
    }
    
    //TODO remove Debug-Prints sizes
    if(false){
    std::cout << "Anzahl in source: " << source.size() << std::endl;
    std::cout << "Anzahl in target: " << target.size() << std::endl;
    std::cout << "Anzahl in matchedTarget: " << matchedTarget.size() << std::endl;
    std::cout << "Anzahl in matchedSource: " << matchedSource.size() << std::endl;
    }
    
    //TODO remove Debug-Prints Punktpaare
    if(false){
    for(size_t i = 0; i < 10; i++){
      for(size_t j = 0; j < 3; j++){
        std::cout << static_cast<double>(matchedSource[i][j]) << " vs " << static_cast<double>(matchedTarget[i][j]) << std::endl;
      }
    }
    }
    
    prev_prev_ret = prev_ret;
    prev_ret = ret;
  
    // Schwerpunkte bestimmen
  
    std::array<f_float, 3> centerSource = {0.0, 0.0, 0.0};
    std::array<f_float, 3> centerTarget = {0.0, 0.0, 0.0};

    size_t count = std::min(matchedSource.size(), matchedTarget.size());
  
    for (size_t i = 0; i < count; ++i) {
      //entspricht centroid_m (model = source = not moving)
      centerSource[0] += matchedSource[i][0];
      centerSource[1] += matchedSource[i][1];
      centerSource[2] += matchedSource[i][2];

      //entspricht centroid_d (data = target)
      centerTarget[0] += matchedTarget[i][0];
      centerTarget[1] += matchedTarget[i][1];
      centerTarget[2] += matchedTarget[i][2];
    }

    f_float srcSize = static_cast<f_float>(matchedSource.size());
    f_float trgSize = static_cast<f_float>(matchedTarget.size());
    centerSource[0] /= srcSize;
    centerSource[1] /= srcSize;
    centerSource[2] /= srcSize;

    centerTarget[0] /= trgSize;
    centerTarget[1] /= trgSize;
    centerTarget[2] /= trgSize;

    // Rotation und Translation berechnen
    ret = my_sc_ICPminimizer->Align(matchedSource, matchedTarget, alignxf, centerSource, centerTarget);
    std::cout << "alignxf after calculation" << std::endl;
    for(int i = 0; i < 16; i++) {
      std::cout << alignxf[i] << " ";
    }
    std::cout << std::endl;

    transform(target, alignxf, transMat, dalignxf, frame, 0);
    
    std::cout << "ret: " << ret << ", prev_ret: " << prev_ret << ", prev_prev_ret: " << prev_prev_ret << " (epsilon: " << epsilonICP << ")" << std::endl;

    // Abbruchbedingung
    if ( ((sc_abs(ret - prev_ret) < epsilonICP) && (sc_abs(ret - prev_prev_ret) < epsilonICP)) || (iter == max_num_iterations -1) ) {
      // write end pose -> Transformation mit Identitätsmatrix, vgl. icp6D.cc Z.289
      transform(target, id, transMat, dalignxf, frame, 0);
      break;
    }
  }
  return iter;  // Anzahl der durchgeführten Iterationen
}
