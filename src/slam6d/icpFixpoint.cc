/*
 * ICP-fixpoint main implementation
 *
 * Copyright (C) Tom Fleischmann, Jonas Wiesner, Yannik Winzer
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file
 * @brief main file for starting ICP-fixpoint matching
 * @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer. Institute of Computer Science, University of Wuerzburg, Germany.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
using std::string;

#include <vector>
#include <map>

#define int64 opencv_int64
#define uint64 opencv_uint64
#include "scanio/writer.h"
#include "scanio/framesreader.h"
#undef int64
#undef uint64
#define int64 systemc_int64
#define uint64 systemc_uint64
#include "sc_fixed/sc_ICP.h"
#include "sc_fixed/sc_ICPapx.h"
#include "sc_fixed/sc_fixed_math.h"
#include "sc_fixed/sc_fixed_converter.h"

#include <strings.h>

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

// parse console input
int parse_options(int argc, char **argv, std::string &dir, int &mni, int &mdm, int &epsilonICPexp, int &start, int &end, IOType &type, double &scaleFac)
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
    ("iter,i", po::value<int>(&mni)->default_value(50),
     "sets the maximal number of ICP iterations to <NR>")
    ("max,m", po::value<int>(&mdm)->default_value(25),
     "neglegt all data points with a distance larger than <NR> 'units'")
    ("start,s", po::value<int>(&start)->default_value(0),
     "start at scan <arg> (i.e., neglects the first <arg> scans) "
     "[ATTENTION: counting naturally starts with 0]")
    ("end,e", po::value<int>(&end)->default_value(-1),
     "end after scan <arg>")
    ("epsICPexp,5", po::value<int>(&epsilonICPexp)->default_value(3),
     "stop ICP iteration if difference is smaller than 1e-<NR>")
    ("scale,y", po::value<double>(&scaleFac)->default_value(0.01),
    "scale factor for point cloud in m (be aware of the different units for uos (cm) and xyz (m), (default: 0.01 means that input and output remain the same)");

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
         << "\t./bin/icpFixpoint -s 0 -e 1 /your/directory" << std::endl;
    exit(0);
  }
  po::notify(vm);

  if (dir[dir.length()-1] != '/') dir = dir + "/";

  return 0;
}

// füllt Ganzzahlen von vorne mit Nullen auf, um 3 Stellen zu erreichen
std::string format_number(int number) {
    std::stringstream stream;
    stream << std::setw(3) << std::setfill('0') << number;
    return stream.str();
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
  int start = 0, end = -1; // start and end of input data
  int mni = 50; // maximum number of iterations
  int mdm = 25; // maximum distance match
  IOType iotype = UOS;
  double scaleFac = 0.01;
  bool quiet = false;
  int epsilonICPexp = 3; // exponent for fixpoint-epsilon termination criterion

  try {
    parse_options(argc, argv, dir, mni, mdm, epsilonICPexp, start, end, iotype, scaleFac);
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

  readFramesAndTransform(dir, start, end, -1, true, false);

  // ab hier ICP

  sc_ICPminimizer *minimizer = new sc_ICPapx(quiet);
  sc_ICP icp(minimizer, mdm, mni, quiet, epsilonICPexp);

  std::cout << Scan::allScans.size() << " scans detected" << std::endl;

  // lege die transMat's und dalignxf's für alle Scans an
  std::vector<std::array<fixed_val, 16>> transMats;
  std::vector<std::array<fixed_val, 16>> dalignxfs;
  transMats.reserve(Scan::allScans.size());
  dalignxfs.reserve(Scan::allScans.size());

  // füge alle Matrizen jeweils ans Ende des aktuellen Vektors (push_back)
  for(unsigned int i = 0; i < Scan::allScans.size(); i++){
    transMats.push_back(array2fixedArray16(Scan::allScans[i]->get_transMat()));
    dalignxfs.push_back(array2fixedArray16(Scan::allScans[i]->getDAlign()));
  }

  // gib die Ergebnis-Transformationsmatrix für den 0. Scan in .frames-Datei aus
  std::ofstream frame(dir + "/scan000.frames");
  for (unsigned int i = 0; i < (3*(Scan::allScans.size()-1)); i++) {
    writeFrame(frame, transMats[0], 2);
  }
  frame.close();

  // match mit ICP
  Scan *initialScan = Scan::allScans[0];
  DataPointer prevXYZ = initialScan->get("xyz");
  if(!prevXYZ.valid()) {
    std::cerr << "Leere Daten bei Index 0" << std::endl;
  }
  DataXYZ prevDat(prevXYZ);
  std::vector<std::array<fixed_val, 3>> prevFixed = array2fixedArray(prevDat);

  for(unsigned int i = 1; i < Scan::allScans.size(); i++){
    Scan *currentScan = Scan::allScans[i];
    if(!currentScan) {
      std::cerr << "Scan" << i << " ist null!" << std::endl;
      continue;
    }

    DataPointer currentXYZ = currentScan->get("xyz");
    if(!currentXYZ.valid()) {
      std::cerr << "Leere Daten bei Index " << i << std::endl;
    }

    DataXYZ currentDat(currentXYZ);
    std::vector<std::array<fixed_val, 3>> currentFixed = array2fixedArray(currentDat);

    std::cout <<std::endl << "RUNNING..." << std::endl;

    // erstelle den Output-Stream für die .frames-Datei des aktuellen Scans
    std::ofstream frame(dir + "/scan" + format_number(i) + ".frames");

    // fülle Zeilen im Stream für vorherige Scans auf (mit ViewFactor 0)
    for (unsigned int j = 0; j < 3*(i - 1); j++) {
      writeFrame(frame, transMats[i], 0);
    }

    // matche, mit Veränderung der transMat und dalignxf des aktuellen Scans
    int iter = icp.match(prevFixed, currentFixed, transMats[i], dalignxfs[i], frame);

    // fülle Zeilen im Stream für nachfolgende Scans auf (mit ViewFactor 2)
    unsigned int nextScans = (Scan::allScans.size() - i);
    for (unsigned int j = 0; j < 3*(nextScans - 1); j++) {
      writeFrame(frame, transMats[i], 2);
    }

    // schließe den Output-Stream
    frame.close();

    // nutze die berechneten Rotationen für die Ausrichtung des nächsten Scans
    prevFixed.clear();
    std::copy(currentFixed.begin(), currentFixed.end(), std::back_inserter(prevFixed));

    // Ausgabe der Anzahl der nötigen Iterationen
    std::cout << "ITER " << iter << std::endl;
  }

  Scan::closeDirectory();
  std::cout << "Matching done!!!" << std::endl << "Normal program end." << std::endl;
  return 0;
}
