/**
 * @file
 * @brief Definition of point pairs
 *
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_PTPAIR_H__
#define __SC_PTPAIR_H__

#include "sc_fixed/sc_Point.h"
#include "sc_fixed/sc_fixed_math.h"
#include <iostream>
#include <fstream>

/**
 * @brief Representing point pairs
 */
class sc_PtPair {
public:
  /**
   * Constructor, by two 'sc_point' pointers
   */
  inline sc_PtPair(double *_p1, double *_p2);

  inline sc_PtPair(double *_p1, double *_p2, double *_norm);

  inline sc_PtPair(sc_Point &p1, sc_Point &p2);

  inline sc_PtPair();

  inline friend std::ostream& operator<<(std::ostream& os, const sc_PtPair& pair);

  sc_Point p1,  ///< The two points forming the pair
        p2;  ///< The two points forming the pair
};

#include "sc_fixed/sc_PtPair.icc"
#endif
