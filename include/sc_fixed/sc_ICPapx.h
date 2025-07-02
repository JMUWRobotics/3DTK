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

  void Align(const std::vector<std::array<f_float, 3>>& source,
             const std::vector<std::array<f_float, 3>>& matchedTarget,
             f_float *alignxf,
             const std::array<f_float, 3> centerSource,
             const std::array<f_float, 3> centerTarget) override;

  static void computeRt(const double *x, const double *dx, double *alignxf);

  inline int getAlgorithmID() { return 6; };
};

#endif
