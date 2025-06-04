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
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#include "sc_fixed/sc_ICP.h"

#include "slam6d/metaScan.h"			//brauchen wir nicht
#include "slam6d/globals.icc"			//TODO: durch mathematische sc-fixed ersetzen
#include "sc_fixed/sc_fixed_math.h"

#include <iomanip>
using std::cerr;

#include <string.h>

#ifdef _MSC_VER
#if !defined _OPENMP && defined OPENMP
#define _OPENMP
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
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
sc_ICP::sc_ICP(sc_ICPminimizer *my_sc_ICPminimizer, double max_dist_match,
	     int max_num_iterations, bool quiet, bool meta, int rnd, bool eP,
	     int anim, double epsilonICP, int nns_method, bool cuda_enabled,
	     bool cad_matching, int max_num_metascans)
{
  this->my_sc_ICPminimizer = my_sc_ICPminimizer;
  this->anim              = anim;
  this->cuda_enabled      = cuda_enabled;
  this->nns_method        = nns_method;

  if (!quiet) {
    cout << "Maximal distance match      : " << max_dist_match << endl
	 << "Maximal number of iterations: " << max_num_iterations << endl << endl;
  }

  // checks
  if (max_dist_match < 0.0) {
    cerr << "ERROR [SC_ICP]: first parameter (max_dist_match) "
	 << "has to be >= 0,"
	 << endl;
    exit(1);
  }
  if (max_num_iterations < 0) {
    cerr << "ERROR [SC_ICP]: second parameter (max_num_iterations)"
	 << "has to be >= 0."
	 << endl;
    exit(1);
  }

  this->max_dist_match2    = sqr(max_dist_match);
  this->max_num_iterations = max_num_iterations;
  this->quiet              = quiet;
  this->meta               = meta;
  this->rnd                = rnd;
  this->eP                 = eP;
  this->epsilonICP         = epsilonICP;

  // Set initial seed (for "real" random numbers)
  //  srand( (unsigned)time( NULL ) );
  this->cad_matching = cad_matching;

  this->max_num_metascans = max_num_metascans;

  //set the number of point pairs to zero
  nr_pointPair = 0;
}

/**
 * Matches a 3D Scan against a 3D Scan
 * @param PreviousScan The scan or metascan forming the model
 * @param CurrentScan The current scan thas is to be matched
 * @return The number of iterations done in this matching run
 */
int sc_ICP::match(Scan* PreviousScan, Scan* CurrentScan,
                 PairingMode pairing_mode)
{
  double id[16];
  M4identity(id);
  CurrentScan->transform(id, Scan::ICP, 0);  // write end pose
  // If ICP shall not be applied, then just write
  // the identity matrix and return
  if (max_num_iterations == 0) {
    return 0;
  }

  // icp main loop
  double ret = 0.0, prev_ret = 0.0, prev_prev_ret = 0.0;
  int iter = 0;
  double alignxf[16];
  long time = GetCurrentTimeInMilliSec();

  for (iter = 0; iter < max_num_iterations; iter++) {

    prev_prev_ret = prev_ret;
    prev_ret = ret;

    if (iter == 1) time = GetCurrentTimeInMilliSec();

#ifdef _OPENMP
    // Implementation according to the paper
    // "The Parallel Iterative Closest Point Algorithm"
    // by Langis / Greenspan / Godin, IEEE 3DIM 2001
    //
    // The same information are given in (ecrm2007.pdf)
    // Andreas Nüchter. Parallelization of Scan Matching
    // for Robotic 3D Mapping. In Proceedings of the 3rd
    // European Conference on Mobile Robots (ECMR '07),
    // Freiburg, Germany, September 2007
    omp_set_num_threads(OPENMP_NUM_THREADS);

    int max = (int)CurrentScan->size<DataXYZ>("xyz reduced");
    int step = ceil(max / (double)OPENMP_NUM_THREADS);

    vector<sc_PtPair> pairs[OPENMP_NUM_THREADS];
    double sum[OPENMP_NUM_THREADS];
    double centroid_m[OPENMP_NUM_THREADS][3];
    double centroid_d[OPENMP_NUM_THREADS][3];
    double Si[OPENMP_NUM_THREADS][9];
    unsigned int n[OPENMP_NUM_THREADS];

    for (int i = 0; i < OPENMP_NUM_THREADS; i++) {
      sum[i] = centroid_m[i][0] = centroid_m[i][1] = centroid_m[i][2] = 0.0;
      centroid_d[i][0] = centroid_d[i][1] = centroid_d[i][2] = 0.0;
      Si[i][0] = Si[i][1] = Si[i][2] = Si[i][3] = Si[i][4] = 0.0;
      Si[i][5] = Si[i][6] = Si[i][7] = Si[i][8] = 0.0;
      n[i] = 0;
    }

#pragma omp parallel
    {
      int thread_num = omp_get_thread_num();

      Scan::getPtPairsParallel(pairs, PreviousScan, CurrentScan,
			       thread_num, step,
			       rnd, max_dist_match2,
			       sum, centroid_m, centroid_d, pairing_mode);

      n[thread_num] = (unsigned int)pairs[thread_num].size();

      if ((my_sc_ICPminimizer->getAlgorithmID() == 1) ||
          (my_sc_ICPminimizer->getAlgorithmID() == 2)) {
        for (unsigned int i = 0; i < n[thread_num]; i++) {

          double pp[3] = {pairs[thread_num][i].p1.x - centroid_m[thread_num][0],
			  pairs[thread_num][i].p1.y - centroid_m[thread_num][1],
			  pairs[thread_num][i].p1.z - centroid_m[thread_num][2]};
          double qq[3] = {pairs[thread_num][i].p2.x - centroid_d[thread_num][0],
			  pairs[thread_num][i].p2.y - centroid_d[thread_num][1],
			  pairs[thread_num][i].p2.z - centroid_d[thread_num][2]};
          // formula (6)
          Si[thread_num][0] += pp[0] * qq[0];
          Si[thread_num][1] += pp[0] * qq[1];
          Si[thread_num][2] += pp[0] * qq[2];
          Si[thread_num][3] += pp[1] * qq[0];
          Si[thread_num][4] += pp[1] * qq[1];
          Si[thread_num][5] += pp[1] * qq[2];
          Si[thread_num][6] += pp[2] * qq[0];
          Si[thread_num][7] += pp[2] * qq[1];
          Si[thread_num][8] += pp[2] * qq[2];
        }
      }
    } // end parallel

    // do we have enough point pairs?
    unsigned int pairssize = 0;
    for (int i = 0; i < OPENMP_NUM_THREADS; i++) {
      pairssize += n[i];
    }
    //add the number of point pair
    nr_pointPair = pairssize;

    if (pairssize > 3) {
      if ((my_sc_ICPminimizer->getAlgorithmID() == 1) ||
          (my_sc_ICPminimizer->getAlgorithmID() == 2) ) {
        ret = my_sc_ICPminimizer->Align_Parallel(OPENMP_NUM_THREADS,
						n, sum,
						centroid_m, centroid_d,
						Si, alignxf);
      } else if (my_sc_ICPminimizer->getAlgorithmID() == 6) {
        ret = my_sc_ICPminimizer->Align_Parallel(OPENMP_NUM_THREADS,
						n, sum,
						centroid_m, centroid_d,
						pairs,
						alignxf);
      } else {
        cout << "This parallel minimization algorithm is not implemented !!!"
	     << endl;
        exit(-1);
      }
    } else {
      //break;
    }
#else

    double centroid_m[3] = {0.0, 0.0, 0.0};
    double centroid_d[3] = {0.0, 0.0, 0.0};
    vector<sc_PtPair> pairs;

    Scan::getPtPairs(&pairs, PreviousScan, CurrentScan, 0, rnd,
		     max_dist_match2, ret, centroid_m, centroid_d, pairing_mode);

    //set the number of point paira
    nr_pointPair = pairs.size();

    // do we have enough point pairs?
    if (pairs.size() > 3) {
      if (my_sc_ICPminimizer->getAlgorithmID() == 3 ||
	  my_sc_ICPminimizer->getAlgorithmID() == 8 ) {
        memcpy(alignxf, CurrentScan->get_transMat(), sizeof(alignxf));
      }
      ret = my_sc_ICPminimizer->Align(pairs, alignxf, centroid_m, centroid_d);
    } else {
      break;
    }

#endif

//#define PLANAR
#ifdef  PLANAR
    double t_rPosTheta[3], t_rPos[3];
    Matrix4ToEuler(alignxf, t_rPosTheta, t_rPos);
    t_rPos[1] = 0.0;
    t_rPosTheta[0] = 0.0;
    t_rPosTheta[2] = 0.0;
    EulerToMatrix4(t_rPos, t_rPosTheta, alignxf);
#endif //PLANAR

    if ((iter == 0 && anim != -2) || ((anim > 0) && (iter % anim == 0))) {
      // transform the current scan
      CurrentScan->transform(alignxf, Scan::ICP, 0);
    } else {
      // transform the current scan
      CurrentScan->transform(alignxf, Scan::ICP, -1);
    }

    if (((fabs(ret - prev_ret) < epsilonICP) &&
	 (fabs(ret - prev_prev_ret) < epsilonICP)) ||
         (iter == max_num_iterations - 1) ) {
      double id[16];
      M4identity(id);
      if(anim == -2) {
	// write end pose
        CurrentScan->transform(id, Scan::ICP, -1);
      } else {
	// write end pose
        CurrentScan->transform(id, Scan::ICP, 0);
      }
      break;
    }
  }

  long endtime = GetCurrentTimeInMilliSec() - time;
  cout << "TIME  " << endtime << "   ITER " << iter <<  endl;
  return iter;
}


/**
 * Computes the point to point error between two scans
 *
 *
 */
double sc_ICP::Point_Point_Error(Scan* PreviousScan,
				Scan* CurrentScan,
				double max_dist_match,
				unsigned int *np,
				double scale_max)
{
  double scale = log(scale_max) / (max_dist_match*max_dist_match);
  double error = 0;
  unsigned int nr_ppairs = 0;

#ifdef _OPENMP
  omp_set_num_threads(OPENMP_NUM_THREADS);

  int max = (int)CurrentScan->size<DataXYZ>("xyz reduced");
  int step = ceil(max / (double)OPENMP_NUM_THREADS);

  vector<sc_PtPair> pairs[OPENMP_NUM_THREADS];
  double sum[OPENMP_NUM_THREADS];
  double centroid_m[OPENMP_NUM_THREADS][3];
  double centroid_d[OPENMP_NUM_THREADS][3];

  for (int i = 0; i < OPENMP_NUM_THREADS; i++) {
    sum[i] = centroid_m[i][0] = centroid_m[i][1] = centroid_m[i][2] = 0.0;
    centroid_d[i][0] = centroid_d[i][1] = centroid_d[i][2] = 0.0;
  }

#pragma omp parallel
  {
    int thread_num = omp_get_thread_num();
    Scan::getPtPairsParallel(pairs, PreviousScan, CurrentScan,
			     thread_num, step,
			     rnd, sqr(max_dist_match),
			     sum, centroid_m, centroid_d, CLOSEST_POINT);

  }

  for (unsigned int thread_num = 0;
       thread_num < OPENMP_NUM_THREADS;
       thread_num++) {
    for (unsigned int i = 0;
	 i < (unsigned int)pairs[thread_num].size();
	 i++) {
      double dist = sqr(pairs[thread_num][i].p1.x - pairs[thread_num][i].p2.x)
	+ sqr(pairs[thread_num][i].p1.y - pairs[thread_num][i].p2.y)
	+ sqr(pairs[thread_num][i].p1.z - pairs[thread_num][i].p2.z);
      error -= 0.39894228 * exp(dist*scale);
    }
    nr_ppairs += (unsigned int)pairs[thread_num].size();
  }
#else

  double centroid_m[3] = {0.0, 0.0, 0.0};
  double centroid_d[3] = {0.0, 0.0, 0.0};
  vector<sc_PtPair> pairs;

  Scan::getPtPairs(&pairs, PreviousScan, CurrentScan, 0,
		   rnd, sqr(max_dist_match),
		   error, centroid_m, centroid_d,
		   CLOSEST_POINT);

  // getPtPairs computes error as sum of squared distances
  error = 0;

  for (unsigned int i = 0; i < pairs.size(); i++) {
    double dist = sqr(pairs[i].p1.x - pairs[i].p2.x)
      + sqr(pairs[i].p1.y - pairs[i].p2.y)
      + sqr(pairs[i].p1.z - pairs[i].p2.z);
    error -= 0.39894228 * exp(dist*scale);
  }
  nr_ppairs = pairs.size();
#endif

  if (np) *np = nr_ppairs;
  return error/nr_ppairs;
}

/**
 * This function matches the scans only with ICP
 *
 * @param allScans Contains all necessary scans.
 */
void sc_ICP::doICP(vector <Scan *> allScans, PairingMode pairing_mode)
{
  double id[16];
  M4identity(id);

  vector < Scan* > meta_scans;
  Scan* my_MetaScan = 0;


  for(unsigned int i = 0; i < allScans.size(); i++) {
    cout << i << "*" << endl;

    Scan *CurrentScan = allScans[i];
    Scan *PreviousScan = 0;

    if (i > 0) {
      PreviousScan = allScans[i-1];
      if (eP) {                             // extrapolate odometry
        CurrentScan->mergeCoordinatesWithRoboterPosition(PreviousScan);
      }
    }

    if (i > 0) {
    //  if (meta) {
    //    match(my_MetaScan, CurrentScan, pairing_mode);
    //  } else
	//  if (cad_matching) {
	//  match(allScans[0], CurrentScan, pairing_mode);
	// } else {
	  match(PreviousScan, CurrentScan, pairing_mode);
	// }
    }

    // push processed scan
    //if ( meta && my_MetaScan) {
    //  delete my_MetaScan;
    //}
    //if ( meta && i != allScans.size()-1 ) {
    //  meta_scans.push_back(CurrentScan);

      // only keep last n scans as metascans
    //  if(max_num_metascans > 0) {
    //    while(meta_scans.size() > max_num_metascans) {
    //      meta_scans.erase(meta_scans.begin());
    //    }
    //    std::cout << meta_scans.size() << " scans in metascan" << std::endl;
    //  }

    //  my_MetaScan = new MetaScan(meta_scans, nns_method);
    //}
  }
}

