/** @file
 *  @brief The global correction based upon 6D Lu Milios style SLAM in 6D but with a helix correction step
 *
 *  @author Peter Schneider. Institute of Computer Science, University of Koblenz , Germany.
 */

#ifndef __SIMPLEREG_H__
#define __SIMPLEREG_H__

#include <vector>
using std::vector;

#include "slam6d/graphSlam6D.h"
#include "graphSlam6DL.h"
#include "slam6d/icp6D.h"
#include "slam6d/graph.h"

#include "newmat/newmatio.h"

//typedef vector <PtPair> vPtPair;  ///< just a typedef: vPtPair = vector of type PtPair


/*
 * @brief Representation of 3D scan matching based upon Lu/Milios in 6D but using helix correction step.
 */
class simpleRegZ : public graphSlam6DL {

public:
  /**
   * Constructor (default)
   */
  simpleRegZ() {};
  simpleRegZ(icp6Dminimizer *my_icp6Dminimizer,
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

  virtual ~simpleRegZ();

  double doGraphSlam6D(Graph gr, vector <LScan*> MetaScan, int nrIt);
  double doGraphSlam6D(Graph gr, vector <Scan*> MetaScan, int nrIt) {return 0.0;}
  static double LUM[4];

private:

  double genBBdForLinkedPair( int firstScanNum, int secondScanNum, vPtPair *ptpairs,
						NEWMAT::Matrix *B, NEWMAT::ColumnVector *Bd);
};


class simpleReg : public graphSlam6DL {

public:
  /**
   * Constructor (default)
   */
  simpleReg() {};
  simpleReg(icp6Dminimizer *my_icp6Dminimizer,
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

  virtual ~simpleReg();

  double doGraphSlam6D(Graph gr, vector <LScan*> MetaScan, int nrIt);
  double doGraphSlam6D(Graph gr, vector <Scan*> MetaScan, int nrIt) {return 0.0;}
  static double LUM[4];

private:

  double genBBdForLinkedPair( int firstScanNum, int secondScanNum, vPtPair *ptpairs,
      double *centroids_m, double *centroids_d,
						NEWMAT::Matrix *B, NEWMAT::ColumnVector *Bd, double weight = 1.0);
};

#endif
