/** @file
 *  @brief The global correction based upon APX transform
 *
 *  @author Andreas Nuechter. Institute of Computer Science, University of Osnabrueck, Germany.
 *  @author Jan Elseberg. Institute of Computer Science, University of Osnabrueck, Germany.
 */

#ifndef __GAPX6D_H__
#define __GAPX6D_H__

#include <vector>
using std::vector;

#include "slam6d/graphSlam6D.h"
#include "graphSlam6DL.h"
#include "slam6d/icp6D.h"
#include "slam6d/graph.h"

#include "newmat/newmatio.h"

/*
 * @brief
 *
 */
class gapx6DL : public graphSlam6DL {

public:
  /**
   * Constructor (default)
   */
  gapx6DL() {};
  gapx6DL(icp6Dminimizer *my_icp6Dminimizer,
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

  virtual ~gapx6DL();

  double doGraphSlam6D(Graph gr, vector <Scan*> MetaScan, int nrIt) {return 0.0;}
  double doGraphSlam6D(Graph gr, vector <LScan*> MetaScan, int nrIt);

  static double LUM[4];

private:
  double genBArotForLinkedPair(int firstScanNum, int secondScanNum, vPtPair *ptpairs,
						 double *centroids_m, double *centroids_d, NEWMAT::Matrix *B, NEWMAT::ColumnVector *A);
  double genBAtransForLinkedPair(int firstScanNum, int secondScanNum,
						   double *centroids_m, double *centroids_d,
						   NEWMAT::SymmetricMatrix *B, NEWMAT::ColumnVector *A, NEWMAT::ColumnVector &X);

};

#endif
