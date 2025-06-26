/*
 * ICP Fixpoint implementation
 *
 * Released under the GPL version 3.
 *
 */


/**
 * @file
 * @author Tom Fleischmann, Yannik Winzer, Jonas Wiesner. Institute of Computer Science, University of Wuerzburg, Germany.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
using std::string;

#include <vector>
#include <map>

//#include "slam6d/point.h" //wird scheinbar nicht gebraucht
//#include "slam6d/scan.h" //wird scheinbar nicht gebraucht
#define int64 opencv_int64
#define uint64 opencv_uint64
#include "scanio/writer.h"
#include "scanio/framesreader.h"
#undef int64
#undef uint64
#define int64 systemc_int64
#define uint64 systemc_uint64
//#include "slam6d/globals.icc" //wird scheinbar nicht gebraucht
#include "sc_fixed/sc_Point.h" //neu, wird aber scheinbar nicht gebraucht
#include "sc_fixed/sc_ICP.h" //neu
#include "sc_fixed/sc_ICPapx.h" //neu
//#undef int64
//#undef uint64

#include "sc_fixed/sc_fixed_converter.h" //include DataXYZ to vector array converter and otherwise


#ifdef _MSC_VER
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#else
#include <strings.h>
#endif

#include <boost/program_options.hpp>
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

int parse_options(int argc, char **argv, std::string &dir,
            int &start, int &end, bool &use_pose,
        IOType &type, double &scaleFac)
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
    ("scale,y", po::value<double>(&scaleFac)->default_value(0.01),
    "scale factor for point cloud in m (be aware of the different units for uos (cm) and xyz (m), (default: 0.01 means that input and output remain the same)")
    ("trustpose,p", po::bool_switch(&use_pose)->default_value(false),
    "Trust the pose file, do not use the transformation from the .frames files.");
  

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
         << "\t./bin/icpFixpoint -s 0 -e 1 /Your/directory" << std::endl;
    exit(0);
  }
  po::notify(vm);

#ifndef _MSC_VER
  if (dir[dir.length()-1] != '/') dir = dir + "/";
#else
  if (dir[dir.length()-1] != '\\') dir = dir + "\\";
#endif

  return 0;
}

/**
 * program for ICP
 * Usage: bin/icpFixpoint 'dir',
 * with 'dir' the directory of a set of scans
 * ...
 */
int sc_main(int argc, char **argv)
{
  // parsing the command line parameters
  // init, default values if not specified
  std::string dir;
  int    start = 0,   end = -1;
  bool   uP         = false;  // should we use the pose information instead of the frames??
  IOType iotype    = UOS;
  double scaleFac = 0.01;
  bool quiet = false;


  try {
    parse_options(argc, argv, dir, start, end,
     uP, iotype, scaleFac);
  } catch (std::exception& e) {
    std::cerr << "Error while parsing settings: " << e.what() << std::endl;
    exit(1);
  }

  // Get Scans
  Scan::openDirectory(false, dir, iotype, start, end);
  if(Scan::allScans.size() == 0) {
    std::cerr << "No scans found. Did you use the correct format?" << std::endl;
    exit(-1);
  }

  //unsigned int types = PointType::USE_NONE;
  

  readFramesAndTransform(dir, start, end, -1 , uP, false);


  std::cout << "Export all 3D Points to file \"points.pts\"" << std::endl;
  std::cout << "Export all 6DoF poses to file \"positions.txt\"" << std::endl;
  std::cout << "Export all 6DoF matrices to file \"poses.txt\"" << std::endl;
  FILE *redptsout = fopen("points.pts", "wb");
  std::ofstream posesout("positions.txt");
  std::ofstream matricesout("poses.txt");

  //ab hier ICP
  sc_ICPminimizer *minimizer = new sc_ICPapx(quiet);
  sc_ICP icp(minimizer, 500, 500, false, false, 1, false, -1, 0.00001, 1, false, false, 0);

  std::cout << "Minimizer and sc_ICP object created" << std::endl;

  std::cout << Scan::allScans.size() << " scans detected" << std::endl;
  
  for(unsigned int i = 1; i < Scan::allScans.size(); i++){
    std::cout << std::to_string(i) + " match iteration" << std::endl;
    Scan *prevScan = Scan::allScans[i-1];
    Scan *nextScan = Scan::allScans[i];
    if(!prevScan || !nextScan) {
      std::cerr << "Scan" << i << " oder " << i-1 << " ist null!" <<std::endl;
      continue;
    }
    DataPointer prevXYZ = prevScan->get("xyz");
    DataPointer nextXYZ = nextScan->get("xyz");
    if(!prevXYZ.valid() || !nextXYZ.valid()) {
      std::cerr << " Leere Daten bei Index " << i << std::endl;
    }

    DataXYZ prevDat(prevXYZ);
    DataXYZ nextDat(nextXYZ);
    
    std::vector<std::array<f_float, 3>> prevFixed = array2fixedArray(prevDat);
    std::vector<std::array<f_float, 3>> nextFixed = array2fixedArray(nextDat);
    //icp.match(...)
    std::cout << std::to_string(i) + " match iteration" << std::endl;
  }

  //ab hier wieder Ausgabe
  /*
  for(unsigned int i = 0; i < Scan::allScans.size(); i++) {
    Scan *source = Scan::allScans[i];
    //std::string red_string = red > 0 ? " reduced" : "";
    DataXYZ xyz  = source->get("xyz" + std::string(""));
    
    // write_uos(xyz, redptsout, scaleFac*100.0, false, false);
    //writeTrajectoryUOS(posesout, source->get_transMat(), false, scaleFac*100.0);
    //writeTrajectoryUOS(matricesout, source->get_transMat(), true, scaleFac*100.0);
  }
  */

  //fclose(redptsout);
  //posesout.close();
  //posesout.clear();
  //matricesout.close();
  //matricesout.clear();
}
