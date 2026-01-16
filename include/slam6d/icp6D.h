/** @file
 *  @brief Representation of 3D scan matching with ICP
 *  @author Kai Lingemann. Inst. of CS, University of Osnabrueck, Germany.
 *  @author Andreas Nuechter. Inst. of CS, University of Osnabrueck, Germany.
 */

#ifndef __ICP6D_H__
#define __ICP6D_H__

#include <vector>

#include "newmat/newmat.h"
//using namespace NEWMAT;

#include "slam6d/scan.h"
#include "slam6d/icp6Dminimizer.h"
#include "slam6d/pairingMode.h"

/**
 * @brief Representation of 3D scan matching with ICP.
 *
 * Manages the matching of 3D scans.
 * Important values, such as maximal matching distance,
 * maximal number of iterations, etc.
 * are specified in the constructor.
 */
class icp6D {

public:
  icp6D(icp6Dminimizer *my_icp6Dminimizer,
	   double max_dist_match = 25.0,
	   int max_num_iterations = 50,
	   bool quiet = false,
	   bool meta = false,
	   int rnd = 1,
	   bool eP = true,
	   int anim = -1,
	   double epsilonICP = 0.0000001,
	   int nns_method = simpleKD,
	   bool cuda_enabled = false,
	   bool cad_matching = false,
     int max_num_metascans = -1);

  /**
   * Destructor (empty, but needed, because virtual)
   */
  virtual ~icp6D() {};

  void doICP(vector <Scan *> allScans,
		   PairingMode pairing_mode = CLOSEST_POINT);
  virtual int match(Scan* PreviousScan,
				Scan* CurrentScan,
				PairingMode pairing_mode = CLOSEST_POINT);
  void covarianceEuler(Scan *scan1, Scan *scan2, NEWMAT::Matrix *C);
  void covarianceQuat(Scan *scan1, Scan *scan2, NEWMAT::Matrix *C);
  double Point_Point_Error(Scan* PreviousScan,
					  Scan* CurrentScan,
					  double max_dist_match,
					  unsigned int *nrp=0,
            double scale_max = 0.000001);

  inline int  get_rnd();
  inline bool get_meta();
  inline int  get_anim();
  inline int get_nns_method();
  inline void set_anim(int anim);
  inline double get_max_dist_match2();
  inline void set_max_dist_match2(double max_dist_match2);
  inline void set_max_num_iterations(int max_num_iterations);
  inline void set_cad_matching (bool cad_matching, int cad_index = 0);
  inline bool get_cad_matching (void);
  inline void set_meta(bool meta);
  inline int get_nr_pointPair();

protected:

  /**
   * suppress output to cout
   */
  bool quiet;

  /**
   * take every rnd point for matching
   */
  int rnd;

  /**
   * extrapolate odometry
   */
  bool eP;

  /**
   * match against all scans (= meta scan), or against the last scan only
   */
  bool   meta;

  /**
   * specifies which NNS method should be used
   */
  int nns_method;

  /**
   * specifies if the ANN trees have to be built
   */
  bool cuda_enabled;

  /**
   * the maximal distance (^2 !!!) for matching
   */
  double max_dist_match2;

  /**
   * the maximal number of iterations
   */
  int    max_num_iterations;

  /**
   * write anim'th animation frame
   */
  int anim;

  /**
   * epsilon for stopping ICP algorithm ( convergence criterium )
   */
  double epsilonICP;

  /**
   * ptr to ICP error function minimizer functor
   */
  icp6Dminimizer *my_icp6Dminimizer;

  /**
   * Maximum number of points in all scans
   */
  unsigned int max_scn_size; //FIXME -> update with metascan

  /**
   * determines if CAD models are matched against one scan
   */
  bool cad_matching;

  /**
   * the index of the scan to match against in cad_matching
   */
  int cad_index;

  /**
   * number of matched points in ICP
   */
  int nr_pointPair;

  /**
   * Window size for ICP with metascans
   */
  int max_num_metascans;
};

#include "icp6D.icc"

#endif
