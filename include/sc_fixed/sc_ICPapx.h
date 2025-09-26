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

  fixed_val Align(const std::vector<std::array<fixed_val, 3>>& matchedSource,
             const std::vector<std::array<fixed_val, 3>>& matchedTarget,
             fixed_val *alignxf,
             const std::array<fixed_val, 3> centerSource,
             const std::array<fixed_val, 3> centerTarget) override;

  inline int getAlgorithmID() { return 6; };
};

#endif
