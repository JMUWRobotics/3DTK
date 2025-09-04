/*
 * sc_ICPminimizer implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Implementation of the virtual functor for ICP error function minimization
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer. Institute of Computer Science, University of Wuerzburg, Germany.
 */

#ifndef __SC_ICPMINIMIZER__
#define __SC_ICPMINIMIZER__

#include <vector>
#include "sc_fixed/sc_fixed_math.h"

#include <iostream>
#include <stdlib.h>
#include <array>

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

  virtual f_float Align(const std::vector<std::array<f_float, 3>>& matchedSource,
                     const std::vector<std::array<f_float, 3>>& matchedTarget,
                     f_float *alignxf,
                     const std::array<f_float, 3> centerSource,
                     const std::array<f_float, 3> centerTarget) = 0;

  virtual int getAlgorithmID() = 0;

protected:
  bool quiet; ///< determines the verbosity
};

#endif
