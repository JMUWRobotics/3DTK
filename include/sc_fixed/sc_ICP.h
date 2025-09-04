/** @file
 *  @brief Representation of 3D scan matching with ICP
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_ICP_H__
#define __SC_ICP_H__

#include <vector>
#include <array>

#include "sc_fixed/sc_fixed_math.h"
#include "sc_fixed/sc_ICPminimizer.h"

/**
 * Manages the matching of 3D scans.
 * Important values, such as maximal matching distance,
 * maximal number of iterations, etc.
 * are specified in the constructor.
 */
class sc_ICP {

public:
  /**
   * Constructor
   */
  sc_ICP(sc_ICPminimizer *my_sc_ICPminimizer,
	   f_float max_dist_match = 25.0,
	   int max_num_iterations = 50,
	   bool quiet = false,
	   int epsilonICPexp = 3);

  /**
   * Destructor (empty, but needed, because virtual)
   */
  virtual ~sc_ICP() {};

  virtual int match(std::vector<std::array<f_float, 3>>& source, std::vector<std::array<f_float, 3>>& target, std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf, std::ofstream& frame);

  inline f_float get_max_dist_match2();
  inline void set_max_dist_match2(f_float max_dist_match2);
  inline void set_max_num_iterations(int max_num_iterations);
  inline int get_nr_pointPair();

protected:

  /**
   * suppress output to cout
   */
  bool quiet;

  /**
   * the maximal distance (^2 !!!) for matching
   */
  f_float max_dist_match2;

  /**
   * the maximal number of iterations
   */
  int max_num_iterations;

  /**
   * epsilon for stopping ICP algorithm ( convergence criterium )
   */
  int epsilonICPexp;
  f_float epsilonICP;

  /**
   * ptr to ICP error function minimizer functor
   */
  sc_ICPminimizer *my_sc_ICPminimizer;

  /**
   * Maximum number of points in all scans
   */
  unsigned int max_scn_size; //FIXME -> update with metascan

  /**
   * number of matched points in ICP
   */
  int nr_pointPair;
};

#include "sc_fixed/sc_ICP.icc"

#endif //__SC_ICP_H__
