/** @file
 *  @brief Definition of the ICP error function minimization
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_ICPAPX_H__
#define __SC_ICPAPX_H__

#include "sc_fixed/sc_fixed_math.h"
#include "sc_fixed/sc_ICPminimizer.h"

/**
 * @brief Implementation of the ICP error function minimization via
 *        approximation using the small angle approximation
 */
class sc_ICPapx : public sc_ICPminimizer
{
public:
  /**
   * Constructor
   */
  sc_ICPapx(bool quiet = false) : sc_ICPminimizer(quiet) {};
  /**
   * Destructor
   */
  virtual ~sc_ICPapx() {};

  double Align(const std::vector<sc_PtPair>& Pairs,
			double *alignxf,
			const double centroid_m[3],
			const double centroid_d[3]);
  double Align_Parallel(const int openmp_num_threads,
				    const unsigned int n[OPENMP_NUM_THREADS],
				    const double sum[OPENMP_NUM_THREADS],
				    const double centroid_m[OPENMP_NUM_THREADS][3],
				    const double centroid_d[OPENMP_NUM_THREADS][3],
				    const std::vector<sc_PtPair> pairs[OPENMP_NUM_THREADS],
				    double *alignxf);

  static void computeRt(const double *x, const double *dx, double *alignxf);

  inline int getAlgorithmID() { return 6; };
};

#endif
