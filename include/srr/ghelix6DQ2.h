/** @file
 *  @brief The global correction based upon 6D Lu Milios style SLAM in 6D but with a helix correction step
 *
 *  @author Peter Schneider. Institute of Computer Science, University of Koblenz , Germany.
 */

#ifndef __GHELIX_Q2L_H__
#define __GHELIX_Q2L_H__

#include <vector>
using std::vector;

#include "slam6d/graphSlam6D.h"
#include "graphSlam6DL.h"
#include "slam6d/icp6D.h"
#include "slam6d/graph.h"

#include "newmat/newmatio.h"


/*
 * @brief Representation of 3D scan matching based upon Lu/Milios in 6D but using helix correction step.
 *
 *            H. Pottman, et al., Simultaneous registration of multiple
 *            views of a 3D object,
 */
class ghelix6DQ2L : public graphSlam6DL {

public:
  /**
   * Constructor (default)
   */
  ghelix6DQ2L() {};
  ghelix6DQ2L(icp6Dminimizer *my_icp6Dminimizer,
	   double mdm = 25.0,
	   double max_dist_match = 25.0,
	   int max_num_iterations = 50,
	   bool quiet = false,
	   bool meta = false,
	   int rnd = 1,
	   bool eP = true,
	   int anim = -1,
	   double epsilonICP = 0.0000001,
	   bool use_cache = false,
	   double epsilonLUM = 0.5);

  virtual ~ghelix6DQ2L();

  double doGraphSlam6D(Graph gr, vector <LScan*> MetaScan, int nrIt);
  double doGraphSlam6D(Graph gr, vector <Scan*> MetaScan, int nrIt) {return 0.0;}
  static double LUM[4];

private:

  double genBBdForLinkedPair( int firstScanNum, int secondScanNum, vPtPair *ptpairs,
						NEWMAT::Matrix *B, NEWMAT::ColumnVector *Bd);
};

#endif
