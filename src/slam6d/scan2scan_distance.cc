/*
 * scan2scan_distance implementation
 *
 * Copyright (C) by the 3DTK contributors
 * Copyright (C) Dorit Borrmann, Helge Lauterbach
 *
 * Released under the GPL version 3.
 *
 */


/**
 * @file
 * @brief Main program for computing scan2scan distance.
 *
 * Program to compute the distance between scans using different distance criteria
 * ATTENTION: Scans will be written to the dir/diff if no 'outdir' is given
 *
 * @author Dorit Borrmann. University of Wuerzburg, Germany.
 * @author Helge Lauterbach. University of Wuerzburg, Germany.
 */

#include <iostream>
#include <iomanip>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "slam6d/globals.icc"
#include "slam6d/scan.h"
#include "slam6d/searchTree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace po = boost::program_options;

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

void parseArgs(int argc, char** argv, std::string &dir, int &start, int &end, IOType &type,
  double &max_dist, std::string &outdir, std::string &custom_filter, PairingMode &pairing_mode) {
    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "produce help message");

    bool point_to_plane = false;
    bool normal_shoot = false;
    po::options_description input("Input options");
    input.add_options()
            ("start,s", po::value<int>(&start)->default_value(0),
             "reference scan <arg> (i.e., neglects the first <arg> scans) "
             "[ATTENTION: counting naturally starts with 0]")
            ("end,e", po::value<int>(&end)->default_value(-1),
             "compared scan <arg>")
            ("format,f", po::value<IOType>(&type)->default_value(UOS, "uos"),
             "using shared library <arg> for input. (chose F from {uos, uos_map, "
             "uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, "
             "riegl_txt, riegl_rgb, riegl_bin, zahn, ply, las})")
            ("dist,d", po::value<double>(&max_dist)->default_value(100), "maximum search distance")
    	    ("normal_shoot-simple,7", po::bool_switch(&normal_shoot)->default_value(false),
             "use closest point along normal for point correspondences'")
            ("point-to-plane-simple,z", po::bool_switch(&point_to_plane)->default_value(false),
             "use point to plane distance for correspondences'")
            ("output,o",po::value<std::string>(&outdir)->default_value("distance"),"output directory")
            ("customFilter,u", po::value<std::string>(&custom_filter),
             "Apply a custom filter. Filter mode and data are specified as a "
             "semicolon-seperated string:\n"
             "\"{filterMode};{nrOfParams}[;param1][;param2][...]\"\n"
             "Multiple filters can be specified in a file (syntax in file is same as "
             "direct specification)\n"
             "\"FILE;{fileName}\"\n"
             "See filter implementation in src/slam6d/pointfilter.cc for more detail.")
            ;

    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("input-dir", po::value<std::string>(&dir), "input dir");

    po::options_description all;
    all.add(generic).add(input).add(hidden);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input);

    po::positional_options_description pd;
    pd.add("input-dir", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), vm);

    po::notify(vm);

    if (vm.count("help")) {
        std::cout << cmdline_options;
        exit(1);
    }

#ifndef _MSC_VER
    if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
    if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif

    if(point_to_plane) pairing_mode = CLOSEST_PLANE_SIMPLE;
    if(normal_shoot) pairing_mode = CLOSEST_POINT_ALONG_NORMAL_SIMPLE;

}

int main(int argc, char* argv[])
{
  std::string dir;
  std::string outdir;
  IOType type = UOS;
  int start;
  int end;
  double max_dist;

  std::string customFilter;
  bool customFilterActive = false;

  PairingMode pairing_mode = CLOSEST_POINT;

  parseArgs(argc, argv, dir, start, end, type, max_dist, outdir, customFilter, pairing_mode);

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

  double max_dist2 = max_dist*max_dist;

  boost::filesystem::create_directory(outdir);
  std::cout << "Writing files to " << outdir << std::endl;
  Scan::openDirectory(false, dir, type, start, end);

  Scan* reference = Scan::allScans[0];
  Scan* scan      = Scan::allScans[Scan::allScans.size() - 1];

  reference->setRangeFilter(-1,-1);
  if (customFilterActive) reference->setCustomFilter(customFilter);
  unsigned int types = 0;
  if ((pairing_mode == CLOSEST_POINT_ALONG_NORMAL_SIMPLE) ||
      (pairing_mode == CLOSEST_PLANE_SIMPLE)) {
    types = PointType::USE_NORMAL;
  }
  reference->setReductionParameter(-1, 0); //, PointType(types);
  reference->setSearchTreeParameter(simpleKD,20);

  scan->setRangeFilter(-1,-1);
  if (customFilterActive) scan->setCustomFilter(customFilter);
  scan->setReductionParameter(-1, 0, PointType(types));
  scan->setSearchTreeParameter(simpleKD,20);

  std::cout << "loading scans finished " << types << std::endl;

  double max_dist_match2 = max_dist*max_dist;
  unsigned int nr_pointPair = 0;
#ifdef _OPENMP
    // Implementation according to the paper
    // "The Parallel Iterative Closest Point Algorithm"
    // by Langis / Greenspan / Godin, IEEE 3DIM 2001
    //
    // The same information are given in (ecrm2007.pdf)
    // Andreas Nüchter. Parallelization of Scan Matching
    // for Robotic 3D Mapping. In Proceedings of the 3rd
    // European Conference on Mobile Robots (ECMR '07),
    // Freiburg, Germany, September 2007
    omp_set_num_threads(OPENMP_NUM_THREADS);

    int max = (int)scan->size<DataXYZ>("xyz reduced");
    //std::cout << max << " points" << std::endl;
    int step = max / OPENMP_NUM_THREADS;
    //std::cout << "Step: " << step << " " << OPENMP_NUM_THREADS << std::endl;

    std::vector<PtPair> pairs[OPENMP_NUM_THREADS];
    double sum[OPENMP_NUM_THREADS];
    double centroid_m[OPENMP_NUM_THREADS][3];
    double centroid_d[OPENMP_NUM_THREADS][3];
    unsigned int n[OPENMP_NUM_THREADS];

    for (int i = 0; i < OPENMP_NUM_THREADS; i++) {
      sum[i] = centroid_m[i][0] = centroid_m[i][1] = centroid_m[i][2] = 0.0;
      centroid_d[i][0] = centroid_d[i][1] = centroid_d[i][2] = 0.0;
      n[i] = 0;
    }
    #pragma omp parallel
    {
      int thread_num = omp_get_thread_num();
      Scan::getPtPairsParallel(pairs, reference, scan,
			       thread_num, step,
			       -1, max_dist_match2,
			       sum, centroid_m, centroid_d, pairing_mode);

      n[thread_num] = (unsigned int)pairs[thread_num].size();

    } // end parallel

    // do we have enough point pairs?
    unsigned int pairssize = 0;
    for (int i = 0; i < OPENMP_NUM_THREADS; i++) {
      pairssize += n[i];
    }
    //add the number of point pair
    nr_pointPair = pairssize;
    std::cout << "Writing points to " << outdir << std::endl;
    std::ofstream scan_out(outdir + "/scan" + to_string(scan->getIdentifier(),3) + ".3d");
    for(int i=0; i < OPENMP_NUM_THREADS; i++) {
      for(size_t j = 0; j < pairs[i].size(); j++) {
        PtPair &pp = pairs[i][j];
        double p12[3] = {
          pp.p1.x - pp.p2.x,
          pp.p1.y - pp.p2.y,
          pp.p1.z - pp.p2.z };
        scan_out << pp.p2.x << " " << pp.p2.y << " " << pp.p2.z
            << " " << Len(p12) << std::endl;
      }
    }
    scan_out.flush();
    scan_out.close();

#else
    std::cout << "single core version" << std::endl;

    std::cout << "has not been implemented yet, aborting..." << std::endl;
    return 1;
/** unused code with erros during building without OpenMP
    double ret = 0.0;
    double centroid_m[3] = {0.0, 0.0, 0.0};
    double centroid_d[3] = {0.0, 0.0, 0.0};
    std::vector<PtPair> pairs;

    Scan::getPtPairs(&pairs, PreviousScan, CurrentScan, 0, rnd,
		     max_dist_match2, ret, centroid_m, centroid_d, CLOSEST_POINT);

    //set the number of point paira
    nr_pointPair = pairs.size();
**/
#endif
    std::cout << nr_pointPair << " pairs" << std::endl;

  return 0;
}
