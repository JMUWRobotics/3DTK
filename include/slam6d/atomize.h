#ifndef _ATOMIZE_H
#define _ATOMIZE_H

#include <string>

using std::string;

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "slam6d/scan.h"
#include "globals.icc"

namespace po = boost::program_options;

// Enable bold characters on the outstream
ostream& bold_on(ostream& os) { return os << "\e[1m"; }

// Disable bold characters on the coutstream
ostream& bold_off(ostream& os) { return os << "\e[0m"; }

// Convert from string to IOType
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

int parse_options(int argc, char **argv, std::string &cond_dir, std::string &orig_dir,
            double &red, int &rand, int &start, int &end, int &maxDist, int &minDist,
            int &octree, IOType &type, std::string& customFilter, int &split, bool& rm_scatter,
            bool& skip_empty, bool& trustpose)
{
po::options_description generic("Generic options");
  generic.add_options()
    ("help,h", "output this help message");

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
    ("split,S", po::value<int>(&split)->default_value(end),
    "Iterativley put <arg> scans together. If not used, condense will put all scans in one file.")
    ("customFilter,u", po::value<string>(&customFilter),
    "Apply a custom filter. Filter mode and data are specified as a semicolon-seperated string:"
    "{filterMode};{nrOfParams}[;param1][;param2][...]\n"
    "Multiple filters can be specified in a file (syntax in file is same as direct specification\n"
    "FILE;{fileName}\n"
    "See filter implementation in src/slam6d/pointfilter.cc for more detail.")
    ("reduce,r", po::value<double>(&red)->default_value(-1.0),
    "turns on octree based point reduction (voxel size=<NR>)")
    ("octree,O", po::value<int>(&octree)->default_value(1),
    "use randomized octree based point reduction (pts per voxel=<NR>)")
    ("rm_scatter,d", po::bool_switch(&rm_scatter)->default_value(false),
     "Note: -r and -O are needed. Removes any voxel that has less than the specified number of points in it.")
    ("min,M", po::value<int>(&minDist)->default_value(-1),
    "neglegt all data points with a distance smaller than NR 'units'")
    ("max,m", po::value<int>(&maxDist)->default_value(-1),
    "neglegt all data points with a distance larger than NR 'units'")
    ("skipEmpty", po::bool_switch(&skip_empty)->default_value(false),
     "Skip empty scans (use it if you have condensed with it!)")
    ("trustpose", po::bool_switch(&trustpose)->default_value(true),
     "Does not do anything. Just to support copy/paste of condense args");

  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("condensed-dir", po::value<std::string>(&cond_dir), "condensed dir")
    ("original-dir", po::value<std::string>(&orig_dir), "original dir");
  // all options
  po::options_description all;
  all.add(generic).add(input).add(hidden);

  // options visible with --help
  po::options_description cmdline_options;
  cmdline_options.add(generic).add(input);

  // positional arguments
  po::positional_options_description pd;
  pd.add("condensed-dir", 1).add("original-dir", 1);

  // process options
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
            options(all).positional(pd).run(), vm);

  // display help
  if (vm.count("help")) {
    std::cout << cmdline_options;
    std::cout << std::endl
         << "This programm is the opposing program to condense." << std::endl
         << "Once you've condensed S subsequent scans into multiple MetaScans and matched them," << std::endl
         << "each MetaScan has a .frames file. The purpose of this program is to split the" << std::endl
         << ".frames files back and apply the relative transformation onto the corresponding single scans." << std::endl
         << bold_on << "Important: Use the same filters and reduction params as with condense!" << bold_off << std::endl
         << "Example usage:" << std::endl
         << "\tbin/atomize /condensed/dir /original/dir [params used for condense]" << std::endl
         << "\tbin/atomize dat/test/cond dat/test -S 10 -s 200 -e 1000 --min 10 --max 100 -r 10 -O 1 " << std::endl;
    exit(0);
  }
  po::notify(vm);

#ifndef _MSC_VER
  if (cond_dir[cond_dir.length()-1] != '/') cond_dir = cond_dir + "/";
  if (orig_dir[orig_dir.length()-1] != '/') orig_dir = orig_dir + "/";
#else
  if (cond_dir[cond_dir.length()-1] != '\\') cond_dir = cond_dir + "\\";
  if (orig_dir[orig_dir.length()-1] != '\\') orig_dir = orig_dir + "\\";
#endif

  return 0;
}

// Does what it says. outMat must be double[16]
void readTransformFromPose(std::string dir, int fileCounter, double* outMat) {
  // Setup filename and filestream
  std::ifstream poseIn;
  std::string poseFileName = dir + "scan" + to_string(fileCounter, 3) + ".pose";
  cout << "Reading transform from " << poseFileName << endl;

  // Read original transform from .pose file
  poseIn.open(poseFileName);
  double pose[6];
  // read 6 plain doubles
  for (int i = 0; i < 6; ++i) poseIn >> pose[i];
  // convert angles from deg to rad
  for (int i = 3; i < 6; ++i) pose[i] = rad(pose[i]);

  // Convert [x y z rotX rotY rotZ] to homogeneous transform matrix
  double rPos[3] = { pose[0], pose[1], pose[2] };
  double rPosTheta[3] = { pose[3], pose[4], pose[5] };
  EulerToMatrix4(rPos, rPosTheta, outMat);
  poseIn.close();
  poseIn.clear();
}

// Does what it says. outMat must be double[16]
void readTransformFromFrames(std::string dir, int fileCounter, double* outMat) {
  // Setup filename and filestream
  std::ifstream frameIn;
  std::string frameFileName = dir + "scan" + to_string(fileCounter, 3) + ".frames";
  cout << "Reading transform from " << frameFileName << endl;
  frameIn.open( frameFileName.c_str() );

  // Read transform from .frames file
  int algoTypeDummy; // used as trash bin for last frames entry
  while( frameIn.good() ) {
    try {
      frameIn >> outMat >> algoTypeDummy;
    } catch (const std::exception &e) {
      break;
    }
  }
  frameIn.close();
  frameIn.clear();
}

void writeFrame(std::string dir, const char* identifier, double* trans)
{
  std::ofstream frameOut;
  std::string frameFileName = dir + "scan" + identifier + ".frames";
  frameOut.open( frameFileName, std::ofstream::out | std::ofstream::trunc );

  if (frameOut.good()) {
    cout << "Writing to " << frameFileName << endl;
    for (int i = 0; i < 16; ++i)
      frameOut << trans[i] << " ";
    frameOut << Scan::ICPINACTIVE << endl;
  } else {
    cout << "ERR: could not open " << frameFileName << endl;
  }

  frameOut.flush();
  frameOut.close();
  frameOut.clear();
}

#endif //_ATOMIZE_H