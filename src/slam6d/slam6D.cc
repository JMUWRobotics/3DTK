/*
 * slam6D implementation
 *
 * Copyright (C) by the 3DTK contributors
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief Main programm for matching 3D scans (6D SLAM)
 *
 * Main programm to match 3D scans with ICP and the globally
 * consistent matching approach.
 * Use -i from the command line to match with ICP,
 * and -I to match 3D Scans using the global algorithm.
 *
 * @author Andreas Nuechter. Jacobs University Bremen gGmbH, Germany
 * @author Kai Lingemann. Inst. of CS, University of Osnabrueck, Germany.
 * @author Jochen Sprickerhof. Inst. of CS, University of Osnabrueck, Germany.
 */

#include "slam6d/scan.h"
#include "slam6d/metaScan.h"
#include "slam6d/io_utils.h"

#include "slam6d/icp6Dapx.h"
#include "slam6d/icp6Dsvd.h"
#include "slam6d/icp6Dquat.h"
#include "slam6d/icp6Dortho.h"
#include "slam6d/icp6Dhelix.h"
#include "slam6d/icp6Ddual.h"
#include "slam6d/icp6Dlumeuler.h"
#include "slam6d/icp6Dlumquat.h"
#include "slam6d/icp6Dquatscale.h"
#include "slam6d/icp6Dnapx.h"
#include "slam6d/icp6D.h"
#include "slam6d/lum6Deuler.h"
#include "slam6d/lum6Dquat.h"
#include "slam6d/ghelix6DQ2.h"
#include "slam6d/elch6Deuler.h"
#include "slam6d/elch6Dquat.h"
#include "slam6d/elch6DunitQuat.h"
#include "slam6d/elch6Dslerp.h"
#include "slam6d/graphSlam6D.h"
#include "slam6d/gapx6D.h"
#include "slam6d/graph.h"
#include "slam6d/globals.icc"

#include <csignal>

#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

#ifdef WITH_METRICS
#include "slam6d/metrics.h"
#endif //WITH_METRICS


#ifdef _MSC_VER
#if !defined _OPENMP && defined OPENMP
#define _OPENMP
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define WANT_STREAM ///< define the WANT stream :)

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <fstream>
using std::ifstream;

#include <boost/filesystem.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

//  Handling Segmentation faults and CTRL-C
void sigSEGVhandler (int v)
{
  static bool segfault = false;
  if(!segfault) {
    segfault = true;
    cout << endl
         << "# **************************** #" << endl
         << "  Segmentation fault or Ctrl-C" << endl
         << "# **************************** #" << endl
         << endl;
    // save frames and close scans
    for(ScanVector::iterator it = Scan::allScans.begin();
     it != Scan::allScans.end();
     ++it) {
      (*it)->saveFrames(Scan::continue_processing);
    }
    cout << "Frames saved." << endl;
    Scan::closeDirectory();
  }
  exit(-1);
}

void validate(boost::any& v, const std::vector<std::string>& values,
              IOType*, int) {
  if (values.size() == 0)
    throw std::runtime_error("Invalid model specification");
  std::string arg = values.at(0);
  try {
    v = formatname_to_io_type(arg.c_str());
  } catch (...) { // runtime_error
    throw std::runtime_error("Format " + arg + " unknown.");
  }
}


/** A function that parses the command-line arguments and sets the respective flags.
 * @param argc the number of arguments
 * @param argv the arguments
 * @param dir the directory
 * @param red using point reduction?
 * @param rand use randomized point reduction?
 * @param mdm maximal distance match
 * @param mdml maximal distance match for SLAM
 * @param mni maximal number of iterations
 * @param start starting at scan number 'start'
 * @param end stopping at scan number 'end'
 * @param maxDist - maximal distance of points being loaded
 * @param minDist - minimal distance of points being loaded
 * @param quiet switches on/off the quiet mode
 * @param veryQuiet switches on/off the 'very quiet' mode
 * @param extrapolate_pose - i.e., extrapolating the odometry by the last transformation
 *        (vs. taking the pose file as <b>exact</b>)
 * @param meta match against all scans (= meta scan), or against the last scan only???
 * @param anim selects the rotation representation for the matching algorithm
 * @param mni_lum sets the maximal number of iterations for SLAM
 * @param net specifies the file that includes the net structure for SLAM
 * @param cldist specifies the maximal distance for closed loops
 * @param epsilonICP stop ICP iteration if difference is smaller than this value
 * @param epsilonSLAM stop SLAM iteration if average difference is smaller than this value
 * @param algo specfies the used algorithm for rotation computation
 * @param lum6DAlgo specifies the used algorithm for global SLAM correction
 * @param loopsize defines the minimal loop size
 * @param bucketSize defines the k-d treeleaf bucket size
 * @return 0, if the parsing was successful. 1 otherwise
 */

int parse_options(int argc, char **argv, string &dir, double &red, int &rand,
              double &mdm, double &mdml, double &mdmll,
              int &mni, int &start, int &end, int &maxDist, int &minDist, string &customFilter, bool &quiet, bool &veryQuiet,
              bool &extrapolate_pose, bool &meta, int &algo, int &loopSlam6DAlgo, int &lum6DAlgo, int &anim,
              int &mni_lum, string &net, double &cldist, int &clpairs, int &loopsize,
              double &epsilonICP, double &epsilonSLAM,  int &nns_method, bool &exportPts, double &distLoop,
              int &iterLoop, double &graphDist, int &octree, IOType &type,
              bool& scanserver, PairingMode &pairing_mode, bool &continue_processing, int &bucketSize,
              boost::filesystem::path &loopclosefile, int &max_num_metascans, int &cad_index)
{

po::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "output this help message");

  bool point_to_plane = false;
  bool normal_shoot = false;
  po::options_description input("Input options");
  input.add_options()
    ("format,f", po::value<IOType>(&type)->default_value(UOS, "uos"),
     "using shared library <arg> for input. (chose F from {uos, uos_map, "
     "uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, "
     "riegl_txt, riegl_rgb, riegl_bin, zahn, ply, las})")
    ("start,s", po::value<int>(&start)->default_value(0),
     "start at scan <arg> (i.e., neglects the first <arg> scans) "
     "[ATTENTION: counting naturally starts with 0]")
    ("end,e", po::value<int>(&end)->default_value(-1),
     "end after scan <arg>")
    ("algo,a", po::value<int>(&algo)->default_value(1),
     "selects the minimization method for the ICP matching algorithm\n"
     "1 = unit quaternion based method by Horn\n"
     "2 = singular value decomposition by Arun et al.\n"
     "3 = orthonormal matrices by Horn et al.\n"
     "4 = dual quaternion method by Walker et al.\n"
     "5 = helix approximation by Hofer & Potmann\n"
     "6 = small angle approximation\n"
     "7 = Lu & Milios style, i.e., uncertainty based, with Euler angles\n"
     "8 = Lu & Milios style, i.e., uncertainty based, with Quaternion\n"
     "9 = unit quaternion with scale method by Horn")
    ("nns_method,t", po::value<int>(&nns_method)->default_value(simpleKD),
    "selects the Nearest Neighbor Search Algorithm\n"
    "0 = simple k-d tree\n"
    "1 = ANNTree\n"
    "2 = BOCTree\n"
    "3 = BruteForce")
    ("loop6DAlgo,L", po::value<int>(&loopSlam6DAlgo)->default_value(0),
     "selects the method for closing the loop explicitly\n"
     "0 = no loop closing technique\n"
     "1 = euler angles\n"
     "2 = quaternions\n"
     "3 = unit quaternions\n"
     "4 = SLERP (recommended)")
    ("graphSlam6DAlgo,G", po::value<int>(&lum6DAlgo)->default_value(0),
    "selects the minimizazion method for the SLAM matching algorithm\n"
    "0 = no global relaxation technique\n"
    "1 = Lu & Milios extension using euler angles due to Borrmann et al.\n"
    "2 = Lu & Milios extension using using unit quaternions\n"
    "3 = HELIX approximation by Hofer and Pottmann\n"
    "4 = small angle approximation")
    ("net,n", po::value<string>(&net)->default_value("none"),
    "specifies the file that includes the net structure for SLAM")
    ("iter,i", po::value<int>(&mni)->default_value(50),
    "sets the maximal number of ICP iterations to <NR>")
    ("iterSLAM,I", po::value<int>(&mni_lum)->default_value(-1),
    "sets the maximal number of iterations for SLAM to <NR>"
    "(if not set, graphSLAM is not executed)")
    ("max,m", po::value<int>(&maxDist)->default_value(-1),
    "neglegt all data points with a distance larger than NR 'units'")
    ("customFilter,u", po::value<string>(&customFilter),
    "Apply a custom filter. Filter mode and data are specified as a "
    "semicolon-seperated string:\n"
    "{filterMode};{nrOfParams}[;param1][;param2][...]\n"
    "Multiple filters can be specified in a file (syntax in file is same as"
    "direct specification)\n"
    "FILE;{fileName}\n"
    "See filter implementation in src/slam6d/pointfilter.cc for more detail.")
    ("loopsize,l", po::value<int>(&loopsize)->default_value(20),
    "sets the size of a loop, i.e., a loop must exceed <NR> of scans")
    ("cldist,c", po::value<double>(&cldist)->default_value(500),
    "specifies the maximal distance for closed loops")
    ("clpairs,C", po::value<int>(&clpairs)->default_value(-1),
    "specifies the minimal number of points for an overlap. If not specified"
    "cldist is used instead")
    ("min,M", po::value<int>(&minDist)->default_value(-1),
    "neglegt all data points with a distance smaller than NR 'units'")
    ("dist,d", po::value<double>(&mdm)->default_value(25.0),
    "sets the maximal point-to-point distance for matching with ICP to <NR> 'units'"
    "(unit of scan data, e.g. cm)")
    ("distSLAM,D", po::value<double>(&mdml)->default_value(25.0),
    "sets the maximal point-to-point distance for matching with SLAM to <NR> 'units'"
    "(unit of scan data, e.g. cm)")
    ("reduce,r", po::value<double>(&red)->default_value(-1.0),
    "turns on octree based point reduction (voxel size=<NR>)")
    ("octree,O", po::value<int>(&octree)->default_value(0),
    "use randomized octree based point reduction (pts per voxel=<NR>)")
    ("random,R", po::value<int>(&rand)->default_value(-1),
    "turns on randomized reduction, using about every <NR>-th point only")
    ("quiet,q", po::bool_switch(&quiet)->default_value(false),
    "Quiet mode. Suppress (most) messages")
    ("veryquiet,Q", po::bool_switch(&veryQuiet)->default_value(false),
    "Very quiet mode. Suppress all messages, except in case of error.")
    ("trustpose,p", po::bool_switch(&extrapolate_pose)->default_value(false),
    "Trust the pose file, do not extrapolate the last transformation."
    "(just for testing purposes, or gps input.)")
    ("anim,A", po::value<int>(&anim)->default_value(-1),
    "if specified, use only every NR-th frame for animation")
    ("metascan,2", po::bool_switch(&meta)->default_value(false),
    "Match current scan against a meta scan of all previous scans (default match against the last scan only)")
    ("DlastSLAM,4", po::value<double>(&mdmll)->default_value(-1.0),
    "sets the maximal point-to-point distance for the final SLAM correction,"
    "if final SLAM is not required don't set it.")
    ("epsICP,5", po::value<double>(&epsilonICP)->default_value(0.00001),
    "stop ICP iteration if difference is smaller than NR")
    ("epsSLAM,6", po::value<double>(&epsilonSLAM)->default_value(0.5),
    "stop SLAM iteration if average difference is smaller than NR")
    ("normal_shoot-simple,7", po::bool_switch(&normal_shoot)->default_value(false),
    "use closest point along normal for point correspondences'")
    ("point-to-plane-simple,z", po::bool_switch(&point_to_plane)->default_value(false),
    "use point to plane distance for correspondences'")
    ("exportAllPoints,8", po::bool_switch(&exportPts)->default_value(false),
    "writes all registered reduced points to the file points.pts before"
    "slam6D terminated")
    ("distLoop,9", po::value<double>(&distLoop)->default_value(700),
    "sets the maximal point-to-point distance for the loop closing")
    ("iterLoop,1", po::value<int>(&iterLoop)->default_value(100),
    "sets the maximal number of iterations for the loop closing")
    ("graphDist,3", po::value<double>(&graphDist)->default_value(500),
    "specifies the maximal distance for closed loops for the final SLAM correction")
    ("scanserver,S", po::bool_switch(&scanserver)->default_value(false),
    "Use the scanserver as an input method and handling of scan data")
    ("continue,0", po::bool_switch(&continue_processing)->default_value(false),
    "continue using last frames entry as starting pose")
    ("bucketSize,b", po::value<int>(&bucketSize)->default_value(20),
    "specifies the bucket size for leafs of the k-d tree. During construction of the"
    "tree, any subtree of at most this size will be replaced by an array.")
    ("loopclosefile", po::value<boost::filesystem::path>(&loopclosefile),
    "filename to write scan poses")
    ("maxmeta", po::value<int>(&max_num_metascans)->default_value(-1),
    "maximum nr of previous scans to combine to a metascan in scan matching")
    ("cad_index", po::value<int>(&cad_index)->default_value(-1),
    "If provided, match all scans against the scan with index <NR>");

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-dir", po::value<std::string>(&dir), "input dir");

  // all options
  po::options_description all;
  all.add(generic).add(input).add(hidden);

  // options visible with --help
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(input);

  // positional argument
  po::positional_options_description pd;
  pd.add("input-dir", 1);

  // process options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).positional(pd).run(), vm);

  // display help
  if (vm.count("help")) {
    std::cout << cmdline_options;
    std::cout << std::endl
         << "Example usage:" << std::endl
         << "\t./bin/slam6d dat" << std::endl
         << "\t./bin/slam6d --max=500 -r 10.2 -i 20 dat" << std::endl
         << "\t./bin/slam6d -s 2 -e 10 dat" << std::endl;
    exit(0);
  }


  po::notify(vm);

  // CAD matching is not always possible:
  if (cad_index >= 0 && meta) {
    cerr << "Error: CAD matching is not possible in MetaScan mode." << endl;
    exit(1);
  }
  if (cad_index >= 0 && loopSlam6DAlgo != 0) {
    cerr << "Error: CAD matching is not possible with Graph SLAM." << endl;
    exit(1);
  }
  if (cad_index >= 0 && lum6DAlgo != 0) {
    cerr << "Error: CAD matching is not possible with Global Relaxiation." << endl;
    exit(1);
  }
  if (cad_index >= 0 && mni_lum > 0) {
    cerr << "Error: CAD matching is not possible with Loop Closing." << endl;
    exit(1);
  }

#ifndef _MSC_VER
  if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
  if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif

  extrapolate_pose = !extrapolate_pose;

  if(point_to_plane) pairing_mode = CLOSEST_PLANE_SIMPLE;
  if(normal_shoot) pairing_mode = CLOSEST_POINT_ALONG_NORMAL_SIMPLE;

  if(!meta) max_num_metascans = -1;

  return 0;
}

/**
 * This function is does all the matching stuff
 * it iterates over all scans using the algorithm objects to calculate new poses
 * objects could be NULL if algorithm should not be used
 *
 * @param cldist maximal distance for closing loops
 * @param loopsize minimal loop size
 * @param allScans Contains all laser scans
 * @param my_icp6D the ICP implementation
 * @param meta_icp math ICP against a metascan
 * @param nns_method Indicates the nearest neigbor search method to be used
 * @param my_loopSlam6D used loopoptimizer
 * @param my_graphSlam6D used global optimization
 * @param nrIt The number of iterations the global SLAM-algorithm will run
 * @param epsilonSLAM epsilon for global SLAM iteration
 * @param mdml maximal distance match for global SLAM
 * @param mdmll max distance match for global SLAM after all scans are matched
 */
void matchGraph6Dautomatic(double cldist,
                           int loopsize,
                           vector <Scan *> allScans,
                           icp6D *my_icp6D,
                           bool meta_icp,
                           int nns_method,
                           loopSlam6D *my_loopSlam6D,
                           graphSlam6D *my_graphSlam6D,
                           int nrIt,
                           double epsilonSLAM,
                           double mdml,
                           double mdmll,
                           double graphDist,
                           bool &eP,
                           IOType type,
                           int max_num_metascans = -1)
{
  double cldist2 = sqr(cldist);

  // list of scan for metascan
  vector < Scan* > metas;

  // graph for loop optimization
  graph_t g;

  int n = allScans.size();

  int loop_detection = 0;
  double dist, min_dist = -1;
  int first = 0, last = 0;

  for(int i = 1; i < n; i++) {
    cout << i << "/" << n << endl;

    add_edge(i-1, i, g);

    if(eP) {
      allScans[i]->mergeCoordinatesWithRoboterPosition(allScans[i-1]);
    }

    //Hack to get all icp transformations into the .frames Files
    if(i == n-1 && my_icp6D != NULL && my_icp6D->get_anim() == -2) {
      my_icp6D->set_anim(-1);
    }

    /*if(i == 85 || i == 321 || i == 533) {
      my_icp6D->set_anim(1);
      }*/

    if(my_icp6D != NULL){
      cout << "ICP" << endl;
      // Matching strongly linked scans with ICPs
      if(meta_icp) {
        metas.push_back(allScans[i - 1]);
        if(max_num_metascans > 0) {
          while(metas.size() > max_num_metascans) {
            metas.erase(metas.begin());
          }
        }
        MetaScan* meta_scan = new MetaScan(metas);
        my_icp6D->match(meta_scan, allScans[i]);
        delete meta_scan;
      } else {
        switch(type) {
        case UOS_MAP:
        case UOS_MAP_FRAMES:
          my_icp6D->match(allScans[0], allScans[i]);
          break;
        case RTS_MAP:
          //untested (and could not work)
          //if(i < 220-22 && i > 250-22) match(allScans[0], CurrentScan);
          my_icp6D->match(allScans[0], allScans[i]);
          break;
        default:
          my_icp6D->match(allScans[i - 1], allScans[i]);
          break;
        }
      }
    } else {
      double id[16];
      M4identity(id);
      allScans[i]->transform(id, Scan::ICP, 0);
    }

    /*if(i == 85 || i == 321 || i == 533) {
      my_icp6D->set_anim(-2);
      }*/

    if(loop_detection == 1) {
      loop_detection = 2;
    }

    for(int j = 0; j < i - loopsize; j++) {
      dist = Dist2(allScans[j]->get_rPos(), allScans[i]->get_rPos());
      if(dist < cldist2) {
        loop_detection = 1;
        if(min_dist < 0 || dist < min_dist) {
          min_dist = dist;
          first = j;
          last = i;
        }
      }
    }

    if(loop_detection == 2) {
      loop_detection = 0;
      min_dist = -1;

      if(my_loopSlam6D != NULL) {
        cout << "Loop close: " << first << " " << last << endl;
        my_loopSlam6D->close_loop(allScans, first, last, g);
        add_edge(first, last, g);
      }

      if(my_graphSlam6D != NULL && mdml > 0) {
        int j = 0;
        double ret;
        do {
          // recalculate graph
          Graph *gr = new Graph(i + 1, cldist2, loopsize);
          cout << "Global: " << j << endl;
          ret = my_graphSlam6D->doGraphSlam6D(*gr, allScans, 1);
          delete gr;
          j++;
        } while (j < nrIt && ret > epsilonSLAM);
      }
    }
  }

  if(loop_detection == 1 && my_loopSlam6D != NULL) {
    cout << "Loop close: " << first << " " << last << endl;
    my_loopSlam6D->close_loop(allScans, first, last, g);
    add_edge(first, last, g);
  }

  if(my_graphSlam6D != NULL && mdml > 0.0) {
    int j = 0;
    double ret;
    do {
      // recalculate graph
      Graph *gr = new Graph(n, cldist2, loopsize);
      cout << "Global: " << j << endl;
      ret = my_graphSlam6D->doGraphSlam6D(*gr, allScans, 1);
      delete gr;
      j++;
    } while (j < nrIt && ret > epsilonSLAM);
  }

  if(my_graphSlam6D != NULL && mdmll > 0.0) {
    my_graphSlam6D->set_mdmll(mdmll);
    int j = 0;
    double ret;
    do {
      // recalculate graph
      Graph *gr = new Graph(n, sqr(graphDist), loopsize);
      cout << "Global: " << j << endl;
      ret = my_graphSlam6D->doGraphSlam6D(*gr, allScans, 1);
      delete gr;
      j++;
    } while (j < nrIt && ret > epsilonSLAM);
  }
}

/**
 * Main program for 6D SLAM.
 * Usage: bin/slam6D 'dir',
 * with 'dir' the directory of a set of scans
 * ...
 */
int main(int argc, char **argv)
{
  signal (SIGSEGV, sigSEGVhandler);
  signal (SIGINT,  sigSEGVhandler);

  cout << "slam6D - "
       << "A highly efficient SLAM implementation based on scan matching"
       << endl
       << "         with 6 degrees of freedom" << endl
       << "(c) University of Wuerzburg, Germany, since 2013" << endl
       << "    Jacobs University Bremen gGmbH, Germany, 2009 - 2013" << endl
       << "    University of Osnabrueck, Germany, 2006 - 2009" << endl << endl;

  // parsing the command line parameters
  // init, default values if not specified
  string dir;
  double red   = -1.0, mdmll = -1.0, mdml = 25.0, mdm = 25.0;
  int    rand  = -1,   mni = 50;
  int    start = 0,   end = -1;
  bool   quiet      = false;
  bool   veryQuiet  = false;
  int    maxDist    = -1;
  int    minDist    = -1;
  bool   eP         = true;   // should we extrapolate the pose??
  bool   meta       = false;  // match against meta scan,
                              // or against LAST scan only?
  int    algo       = 1;
  int    mni_lum    = -1;
  double cldist     = 500;
  int    clpairs    = -1;
  int    loopsize   = 20;
  string net        = "none";
  string customFilter;
  bool customFilterActive = false;
  int    anim       = -1;
  double epsilonICP = 0.00001;
  double epsilonSLAM = 0.5;
  int    nns_method = simpleKD;
  bool   exportPts  = false;
  int    loopSlam6DAlgo  = 0;
  int    lum6DAlgo  = 0;
  double distLoop   = 700.0;
  int iterLoop      = 100;
  double graphDist  = cldist;
  int octree       = 0;       // employ randomized octree reduction?
  IOType type    = UOS;
  bool scanserver = false;
  PairingMode pairing_mode = CLOSEST_POINT;
  bool continue_processing = false;
  int bucketSize = 20;
  boost::filesystem::path loopclose("loopclose.pts");
  int max_num_metascans = -1;
  int cad_index = -1;

  parse_options(argc, argv, dir, red, rand, mdm, mdml, mdmll, mni, start, end,
            maxDist, minDist, customFilter, quiet, veryQuiet, eP, meta,
            algo, loopSlam6DAlgo, lum6DAlgo, anim,
            mni_lum, net, cldist, clpairs, loopsize, epsilonICP, epsilonSLAM,
            nns_method, exportPts, distLoop, iterLoop, graphDist, octree, type,
            scanserver, pairing_mode, continue_processing, bucketSize,
            loopclose, max_num_metascans, cad_index);

  /* writing frames in zip archives is not supported by BasicScan */
  if(!boost::filesystem::is_directory(dir)) {
    std::cerr << argv[0] << " does not support scans in archives" << std::endl;
    exit(-1);
  }

  cout << "slam6D will proceed with the following parameters:" << endl;
  //@@@ to do :-)
  // TODO: writer a proper TODO ^

  if (continue_processing) Scan::continueProcessing();
  Scan::setProcessingCommand(argc, argv);

  Scan::openDirectory(scanserver, dir, type, start, cad_index > end ? cad_index : end);

  if(Scan::allScans.size() == 0) {
    cerr << "No scans found. Did you use the correct format?" << endl;
    exit(-1);
  }
  // custom filter set? quick check, needs to contain at least one ';'
  // (proper checking will be done case specific in pointfilter.cc)
  size_t pos = customFilter.find_first_of(";");
  if (pos != std::string::npos) {
    customFilterActive = true;

    // check if customFilter is specified in file
    if (customFilter.find("FILE;") == 0) {
      std::string selection_file_name = customFilter.substr(5, customFilter.length());
      std::ifstream selectionfile;
      // open the input file
      selectionfile.open(selection_file_name, std::ios::in);

      if (!selectionfile.good()) {
        std::cerr << "Error loading custom filter file " << selection_file_name << "!" << std::endl;
        std::cerr << "Data will NOT be filtered.!" << std::endl;
        customFilterActive = false;
      }
      else {
        std::string line;
        std::string custFilt;
        while (std::getline(selectionfile, line)) {
          // allow comment or empty lines
          if (line.find("#") == 0) continue;
          if (line.length() < 1) continue;
          custFilt = custFilt.append(line);
          custFilt = custFilt.append("/");
        }
        if (custFilt.length() > 0) {
          // last '/'
          customFilter = custFilt.substr(0, custFilt.length() - 1);
        }
      }
      selectionfile.close();
    }
  }
  else {
    // give a warning if custom filter has been inproperly specified
    if (customFilter.length() > 0){
      std::cerr << "Custom filter: specifying string has not been set properly, data will NOT be filtered." << std::endl;
    }
  }

  cout << start << " " << end << " " << Scan::allScans.size() << endl;
  for(ScanVector::iterator it = Scan::allScans.begin();
      it != Scan::allScans.end();
      ++it) {
    Scan* scan = *it;
    if (customFilterActive) scan->setCustomFilter(customFilter);
    scan->setRangeFilter(maxDist, minDist);
    unsigned int types = 0;
    if ((pairing_mode == CLOSEST_POINT_ALONG_NORMAL_SIMPLE) ||
        (pairing_mode == CLOSEST_PLANE_SIMPLE)) {
      types = PointType::USE_NORMAL;
    }
     scan->setReductionParameter(red, octree, PointType(types));
     scan->setSearchTreeParameter(nns_method, bucketSize);
  }
  icp6Dminimizer *my_icp6Dminimizer = 0;
  switch (algo) {
  case 1 :
    my_icp6Dminimizer = new icp6D_QUAT(quiet);
    break;
  case 2 :
    my_icp6Dminimizer = new icp6D_SVD(quiet);
    break;
  case 3 :
    my_icp6Dminimizer = new icp6D_ORTHO(quiet);
    break;
  case 4 :
    my_icp6Dminimizer = new icp6D_DUAL(quiet);
    break;
  case 5 :
    my_icp6Dminimizer = new icp6D_HELIX(quiet);
    break;
  case 6 :
    my_icp6Dminimizer = new icp6D_APX(quiet);
    break;
  case 7 :
    my_icp6Dminimizer = new icp6D_LUMEULER(quiet);
    break;
  case 8 :
    my_icp6Dminimizer = new icp6D_LUMQUAT(quiet);
    break;
  case 9 :
    my_icp6Dminimizer = new icp6D_QUAT_SCALE(quiet);
    break;
  case 10 :
    my_icp6Dminimizer = new icp6D_NAPX(quiet);
    break;
  }
  // match the scans and print the time used
  long starttime = GetCurrentTimeInMilliSec();

#ifdef WITH_METRICS
  Timer t = ClientMetric::matching_time.start();
#endif //WITH_METRICS

  if (mni_lum == -1 && loopSlam6DAlgo == 0) {
    icp6D *my_icp = 0;
    my_icp = new icp6D(my_icp6Dminimizer, mdm, mni, quiet, meta, rand, eP,
                       anim, epsilonICP, nns_method,false,false,max_num_metascans);
    // UOS_CAD is not available as file type, keep it for legacy reasons
    if (type == UOS_CAD)
    {
      my_icp->set_cad_matching(true); // will match against index 0
    }
    // Match against given index
    if (cad_index != -1)
    {
      my_icp->set_cad_matching (true, cad_index);
      if (cad_index > end) {
        vector<Scan*> cad_scanvec;
        cad_scanvec.reserve(end + 2);
        for (int i = 0; i <= end; ++i) {
          cad_scanvec.push_back(Scan::allScans[i]);
        }
        cad_scanvec.push_back(Scan::allScans[cad_index]);
        my_icp->set_cad_matching(true, end + 1);
        my_icp->doICP(cad_scanvec, pairing_mode);
      } else {
        my_icp->doICP(Scan::allScans, pairing_mode);
      }
    } // The usual successive matching
    else
    {
      if (my_icp) my_icp->doICP(Scan::allScans, pairing_mode);
    }

    delete my_icp;
  } else if (clpairs > -1) {
    //!!!!!!!!!!!!!!!!!!!!!!!!
    icp6D *my_icp = 0;
    my_icp = new icp6D(my_icp6Dminimizer, mdm, mni, quiet, meta, rand, eP,
                       anim, epsilonICP, nns_method,false,false,max_num_metascans);
    my_icp->doICP(Scan::allScans, pairing_mode);
    graphSlam6D *my_graphSlam6D = new lum6DEuler(my_icp6Dminimizer,
                                                 mdm, mdml, mni, quiet, meta,
                                                 rand, eP, anim, epsilonICP,
                                                 nns_method, epsilonSLAM);
    my_graphSlam6D->matchGraph6Dautomatic(Scan::allScans, mni_lum,
                                          clpairs, loopsize);

    //!!!!!!!!!!!!!!!!!!!!!!!!
  } else {
    graphSlam6D *my_graphSlam6D = 0;
    switch (lum6DAlgo) {
    case 1 :
      my_graphSlam6D = new lum6DEuler(my_icp6Dminimizer, mdm, mdml, mni,
                                      quiet, meta, rand, eP,
                                      anim, epsilonICP, nns_method,
                                      epsilonSLAM);
      break;
    case 2 :
      my_graphSlam6D = new lum6DQuat(my_icp6Dminimizer, mdm, mdml, mni,
                                     quiet, meta, rand, eP,
                                     anim, epsilonICP, nns_method, epsilonSLAM);
      break;
    case 3 :
      my_graphSlam6D = new ghelix6DQ2(my_icp6Dminimizer, mdm, mdml, mni,
                                      quiet, meta, rand, eP,
                                      anim, epsilonICP, nns_method,
                                      epsilonSLAM);
      break;
    case 4 :
      my_graphSlam6D = new gapx6D(my_icp6Dminimizer, mdm, mdml, mni,
                                  quiet, meta, rand, eP,
                                  anim, epsilonICP, nns_method, epsilonSLAM);
      break;
    }
    // Construct Network
    if (net != "none") {
      icp6D *my_icp = 0;
      my_icp = new icp6D(my_icp6Dminimizer, mdm, mni, quiet, meta, rand, eP,
                         anim, epsilonICP, nns_method);
      my_icp->doICP(Scan::allScans, pairing_mode);

      Graph* structure;
      structure = new Graph(net);
      my_graphSlam6D->doGraphSlam6D(*structure, Scan::allScans, mni_lum);
      if(mdmll > 0.0) {
        my_graphSlam6D->set_mdmll(mdmll);
        my_graphSlam6D->doGraphSlam6D(*structure, Scan::allScans, mni_lum);
      }

    } else {
      icp6D *my_icp = 0;
      if(algo > 0) {
        my_icp = new icp6D(my_icp6Dminimizer, mdm, mni, quiet, meta, rand, eP,
                           anim, epsilonICP, nns_method);

        loopSlam6D *my_loopSlam6D = 0;
        switch(loopSlam6DAlgo) {
        case 1:
          my_loopSlam6D = new elch6Deuler(veryQuiet, my_icp6Dminimizer,
                                          distLoop, iterLoop,
                                          rand, eP, 10, epsilonICP,
                                          nns_method);
          break;
        case 2:
          my_loopSlam6D = new elch6Dquat(veryQuiet, my_icp6Dminimizer,
                                         distLoop, iterLoop,
                                         rand, eP, 10, epsilonICP,
                                         nns_method);
          break;
        case 3:
          my_loopSlam6D = new elch6DunitQuat(veryQuiet, my_icp6Dminimizer,
                                             distLoop, iterLoop,
                                             rand, eP, 10, epsilonICP,
                                             nns_method);
          break;
        case 4:
          my_loopSlam6D = new elch6Dslerp(veryQuiet, my_icp6Dminimizer,
                                          distLoop, iterLoop,
                                          rand, eP, 10, epsilonICP,
                                          nns_method);
          break;
        }
        matchGraph6Dautomatic(cldist, loopsize, Scan::allScans, my_icp, meta,
                              nns_method, my_loopSlam6D, my_graphSlam6D,
                              mni_lum, epsilonSLAM, mdml, mdmll, graphDist, eP,
                              type, max_num_metascans);
        delete my_icp;
        if(my_loopSlam6D) {
          delete my_loopSlam6D;
        }
      }
      if(my_graphSlam6D) {
        delete my_graphSlam6D;
      }
    }
  }

#ifdef WITH_METRICS
  ClientMetric::matching_time.end(t);
#endif //WITH_METRICS

  long endtime = GetCurrentTimeInMilliSec() - starttime;
  cout << "Matching done in " << endtime << " milliseconds!!!" << endl;

  if (exportPts) {
    cout << "Export all 3D Points to file \"points.pts\"" << endl;
    ofstream redptsout("points.pts");
    for(unsigned int i = 0; i < Scan::allScans.size(); i++) {
      DataXYZ xyz_r(Scan::allScans[i]->get("xyz reduced"));
      // DataNormal normal_r(Scan::allScans[i]->get("normal reduced"));
      for(unsigned int i = 0; i < xyz_r.size(); ++i) {
        // int r,g,b;
        // r = (int)(normal_r[i][0] * (127.5) + 127.5);
        // g = (int)(normal_r[i][1] * (127.5) + 127.5);
        // b = (int)(fabs(normal_r[i][2]) * (255.0));
        redptsout << xyz_r[i][0] << ' ' << xyz_r[i][1] << ' ' << xyz_r[i][2]
        // << ' ' << r << ' ' << g << ' ' << b
            << endl;
      }
      redptsout << std::flush;
    }
    redptsout.close();
    redptsout.clear();
  }

  const double* p;
  ofstream redptsout(loopclose.string());
  for(ScanVector::iterator it = Scan::allScans.begin();
      it != Scan::allScans.end();
      ++it)
  {
    Scan* scan = *it;
    p = scan->get_rPos();
    Point x(p[0], p[1], p[2]);
    redptsout << x << endl;
    scan->saveFrames(continue_processing);
  }
  redptsout.close();

  Scan::closeDirectory();
  delete my_icp6Dminimizer;

  cout << endl << endl;
  cout << "Normal program end." << endl
       << (red < 0 && rand < 0 ? "(-> HINT: For a significant speedup, please use the '-r' or '-R' parameter <-)\n"
                               : "")
       << endl;

  // print metric information
#ifdef WITH_METRICS
  ClientMetric::print(scanserver);
#endif //WITH_METRICS
}
