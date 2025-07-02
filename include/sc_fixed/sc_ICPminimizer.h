/** @file
 *  @brief Implementation of the virtual functor
 *         for ICP error function minimization
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_ICPMINIMIZER__
#define __SC_ICPMINIMIZER__

#ifdef _MSC_VER
#if !defined _OPENMP && defined OPENMP
#define _OPENMP
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sc_fixed/sc_PtPair.h"
#include "sc_fixed/sc_fixed_math.h"

#include <iostream>
#include <stdlib.h>

class sc_ICPminimizer {

public:
  /**
   * Constructor
   */
  sc_ICPminimizer(bool quiet = false) { this->quiet = quiet; };
  /**
   * Destructor
   */
  virtual ~sc_ICPminimizer() {};

  /**
   * aligning the point pairs
   */
  // a detailed discussion of the minimization techniques used for this
  // function is given in:
  // Andreas Nuechter, Jan Elseberg, Peter Schneider, and Dietrich Paulus.
  // Study of Parameterizations for the Rigid Body Transformations of The
  // Scan Registration Problem, Journal Computer Vision and Image
  // Understanding (CVIU), Elsevier Science, Volume 114, Issue 8,
  // pp. 963-980, August 2010.
  virtual double Align(const std::vector<sc_PtPair>& Pairs,
		       double *alignxf,
		       const double centroid_m[3],
		       const double centroid_d[3]) = 0;

  virtual void Align(const std::vector<std::array<f_float, 3>>& source,
                     const std::vector<std::array<f_float, 3>>& matchedTarget,
                     f_float *alignxf,
                     const std::array<f_float, 3> centerSource,
                     const std::array<f_float, 3> centerTarget) = 0;

  /**
   * aligning the point pairs parallel algorithms
   */
  virtual double Align_Parallel(const int openmp_num_threads,
				const unsigned int n[OPENMP_NUM_THREADS],
				const double sum[OPENMP_NUM_THREADS],
				const double centroid_m[OPENMP_NUM_THREADS][3],
				const double centroid_d[OPENMP_NUM_THREADS][3],
				const double Si[OPENMP_NUM_THREADS][9],
				double *alignxf)
  {
    std::cout << "this function is not implemented!!!" << std::endl;
    exit(-1);
  }
  virtual double Align_Parallel(const int openmp_num_threads,
						  const unsigned int n[OPENMP_NUM_THREADS],
						  const double sum[OPENMP_NUM_THREADS],
						  const double centroid_m[OPENMP_NUM_THREADS][3],
						  const double centroid_d[OPENMP_NUM_THREADS][3],
						  const std::vector<sc_PtPair> pairs[OPENMP_NUM_THREADS],
						  double *alignxf)
  {
    std::cout << "this function is not implemented!!!" << std::endl;
    exit(-1);
  }

  virtual int getAlgorithmID() = 0;

protected:
  bool quiet; ///< determines the verbosity
};

#endif
