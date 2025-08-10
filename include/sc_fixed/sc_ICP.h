/** @file
 *  @brief Representation of 3D scan matching with ICP
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_ICP_H__
#define __SC_ICP_H__

#include <vector>
#include <array>

#include "sc_fixed/sc_fixed_math.h"
#include "newmat/newmat.h"			//ggf. TODO: ersetzen?
//using namespace NEWMAT;

#include "slam6d/scan.h"             //TODO: ersetzen
#include "sc_fixed/sc_ICPminimizer.h"
#include "slam6d/pairingMode.h"      //kann so bleiben

/**
 * @brief Representation of 3D scan matching with ICP.
 *
 * Manages the matching of 3D scans.
 * Important values, such as maximal matching distance,
 * maximal number of iterations, etc.
 * are specified in the constructor.
 */
class sc_ICP {

public:
  sc_ICP(sc_ICPminimizer *my_sc_ICPminimizer,
	   double max_dist_match = 25.0,
	   int max_num_iterations = 50,
	   bool quiet = false,
	   bool meta = false,
	   int rnd = 1,
	   bool eP = true,
	   int anim = -1,
	   f_float epsilonICP = 0.0000001,
	   int nns_method = BruteForce,
	   bool cuda_enabled = false,
	   bool cad_matching = false,
     int max_num_metascans = -1);

  /**
   * Destructor (empty, but needed, because virtual)
   */
  virtual ~sc_ICP() {};

  void doICP(std::vector<std::vector<std::array<f_float, 3>>> allScans);
  //virtual int match2(std::vector<std::array<f_float, 3>>& PreviousScan, std::vector<std::array<f_float, 3>>& CurrentScan);
  virtual int match(std::vector<std::array<f_float, 3>>& source, std::vector<std::array<f_float, 3>>& target, std::array<f_float, 16>& transMat, std::array<f_float, 16>& dalignxf);
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
  inline void set_cad_matching (bool cad_matching);
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
  sc_ICPminimizer *my_sc_ICPminimizer;

  /**
   * Maximum number of points in all scans
   */
  unsigned int max_scn_size; //FIXME -> update with metascan

  /**
   * determines if CAD models are matched against one scan
   */
  bool cad_matching;

  /**
   * number of matched points in ICP
   */
  int nr_pointPair;

  /**
   * Window size for ICP with metascans
   */
  int max_num_metascans;
};

#include "sc_fixed/sc_ICP.icc"

#endif //__SC_ICP_H__
