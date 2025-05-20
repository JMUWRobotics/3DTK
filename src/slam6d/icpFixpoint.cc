/*
 * exportPoints implementation
 *
 * Copyright (C) Jochen Sprickerhof
 *
 * Released under the GPL version 3.
 *
 */


/**
 * @file
 * @author Jochen Sprickerhof. Institute of Computer Science, University of Osnabrueck, Germany.
 */

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
using std::string;

#include <vector>
#include <map>

#include "slam6d/point.h"
#include "slam6d/scan.h"
#include "scanio/writer.h"
#include "scanio/framesreader.h"
#include "slam6d/globals.icc"

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

int parse_options(int argc, char **argv, std::string &dir, double &red, int &rand,
            int &start, int &end, int &maxDist, int &minDist, bool &use_pose,
            bool &use_xyz, bool &use_reflectance, bool &use_type, bool &use_color, int &octree, IOType &type, std::string& customFilter, double &scaleFac,
	    bool &hexfloat, bool &high_precision, int &frame, bool &use_normals)
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
    ("scale,y", po::value<double>(&scaleFac)->default_value(0.01),
    "scale factor for point cloud in m (be aware of the different units for uos (cm) and xyz (m), (default: 0.01 means that input and output remain the same)")
    ("min,M", po::value<int>(&minDist)->default_value(-1),
    "neglegt all data points with a distance smaller than NR 'units'")
    ("max,m", po::value<int>(&maxDist)->default_value(-1),
    "neglegt all data points with a distance larger than NR 'units'")
    ("color,c", po::bool_switch(&use_color)->default_value(false),
     "export in color as RGB")
    ("reflectance,R", po::bool_switch(&use_reflectance)->default_value(false),
     "export reflectance values")
    ("type,T", po::bool_switch(&use_type)->default_value(false),
     "export type/class values")
    ("normals,N", po::bool_switch(&use_normals)->default_value(false),
     "export point normals")
    ("trustpose,p", po::bool_switch(&use_pose)->default_value(false),
    "Trust the pose file, do not use the transformation from the .frames files.")
    ("xyz,x", po::bool_switch(&use_xyz)->default_value(false),
     "export in xyz format (right handed coordinate system in m)")
    ("hexfloat,0", po::bool_switch(&hexfloat)->default_value(false),
     "export points with hexadecimal digits")
    ("highprecision,H", po::bool_switch(&high_precision)->default_value(false),
     "export points with full double precision")
    ("frame,n", po::value<int>(&frame)->default_value(-1),
     "uses frame NR for export");

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
         << "\t./bin/exportPoints -s 0 -e 1 /Your/directory" << std::endl;
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
 * program for point export
 * Usage: bin/exportPoints 'dir',
 * with 'dir' the directory of a set of scans
 * ...
 */
int main(int argc, char **argv)
{
  // parsing the command line parameters
  // init, default values if not specified
  std::string dir;
  double red   = -1.0;
  int    rand  = -1;
  int    start = 0,   end = -1;
  int    maxDist    = -1;
  int    minDist    = -1;
  bool   uP         = false;  // should we use the pose information instead of the frames??
  bool   use_xyz = false;
  bool   use_color = false;
  bool   use_reflectance = false;
  bool   use_type = false;
  bool   use_normals = false;
  int octree       = 0;  // employ randomized octree reduction?
  IOType iotype    = UOS;
  bool rangeFilterActive = false;
  bool customFilterActive = false;
  std::string customFilter;
  double scaleFac = 0.01;
  bool hexfloat = false;
  bool high_precision = false;
  int frame = -1;

  try {
    parse_options(argc, argv, dir, red, rand, start, end,
      maxDist, minDist, uP, use_xyz, use_reflectance, use_type, use_color, octree, iotype, customFilter, scaleFac,
      hexfloat, high_precision, frame, use_normals);
  } catch (std::exception& e) {
    std::cerr << "Error while parsing settings: " << e.what() << std::endl;
    exit(1);
  }

  if(!supportsNormals(iotype) && use_normals) {
    std::cerr << "WARNING File format does not support normals. Normals are not exported" << std::endl;
    use_normals = false;
  }

  rangeFilterActive = minDist > 0 || maxDist > 0;

  // custom filter set? quick check, needs to contain at least one ';'
  // (proper chsecking will be done case specific in pointfilter.cc)
  size_t pos = customFilter.find_first_of(";");
  if (pos != std::string::npos){
      customFilterActive = true;

      // check if customFilter is specified in file
      if (customFilter.find("FILE;") == 0){
          std::string selection_file_name = customFilter.substr(5, customFilter.length());
          std::ifstream selectionfile;
          // open the input file
          selectionfile.open(selection_file_name, std::ios::in);

          if (!selectionfile.good()){
              std::cerr << "Error loading custom filter file " << selection_file_name << "!" << std::endl;
              std::cerr << "Data will NOT be filtered!" << std::endl;
              customFilterActive = false;
          }
          else {
              std::string line;
              std::string custFilt;
              while (std::getline(selectionfile, line)){
                  if (line.find("#") == 0) continue;
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

  // Get Scans
  Scan::openDirectory(false, dir, iotype, start, end);
  if(Scan::allScans.size() == 0) {
    std::cerr << "No scans found. Did you use the correct format?" << std::endl;
    exit(-1);
  }

  unsigned int types = PointType::USE_NONE;
  if(supportsReflectance(iotype)) types |= PointType::USE_REFLECTANCE;
  if(supportsColor(iotype)) types |= PointType::USE_COLOR;
  if(supportsType(iotype)) types |= PointType::USE_TYPE;

  // if specified, filter scans
  for (size_t i = 0; i < Scan::allScans.size(); i++)  {
     if(rangeFilterActive) Scan::allScans[i]->setRangeFilter(maxDist, minDist);
     if(customFilterActive) Scan::allScans[i]->setCustomFilter(customFilter);
  }

//
  int end_reduction = (int)Scan::allScans.size();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int iterator = 0; iterator < end_reduction; iterator++) {
    if (red > 0) {
      PointType pointtype(types);
      std::cout << "Reducing Scan No. " << iterator << std::endl;
      Scan::allScans[iterator]->setReductionParameter(red, octree, pointtype);
      Scan::allScans[iterator]->calcReducedPoints();
    } else {
      std::cout << "Copying Scan No. " << iterator << std::endl;
    }
    // reduction filter for current scan!
  }

  readFramesAndTransform(dir, start, end, frame, uP, red > -1);

 std::cout << "Export all 3D Points to file \"points.pts\"" << std::endl;
 std::cout << "Export all 6DoF poses to file \"positions.txt\"" << std::endl;
 std::cout << "Export all 6DoF matrices to file \"poses.txt\"" << std::endl;
 FILE *redptsout = fopen("points.pts", "wb");
 std::ofstream posesout("positions.txt");
 std::ofstream matricesout("poses.txt");

  for(unsigned int i = 0; i < Scan::allScans.size(); i++) {
    Scan *source = Scan::allScans[i];
    std::string red_string = red > 0 ? " reduced" : "";

    DataXYZ xyz  = source->get("xyz" + red_string);

    if(use_reflectance) {
      DataReflectance xyz_reflectance =
          (((DataReflectance)source->get("reflectance" + red_string)).size() == 0) ?

          source->create("reflectance" + red_string, sizeof(float)*xyz.size()) :
          source->get("reflectance" + red_string);

      if (!(types & PointType::USE_REFLECTANCE)) {
        for(unsigned int i = 0; i < xyz.size(); i++) xyz_reflectance[i] = 255;
      }
      if(use_xyz) {
        write_xyzr(xyz, xyz_reflectance, redptsout, scaleFac, hexfloat, high_precision);
      } else {
        write_uosr(xyz, xyz_reflectance, redptsout, scaleFac*100.0 , hexfloat, high_precision);
      }

    } else if(use_type) {
      DataType xyz_type =
          (((DataType)source->get("type" + red_string)).size() == 0) ?

          source->create("type" + red_string, sizeof(int)*xyz.size()) :
          source->get("type" + red_string);

      if (!(types & PointType::USE_TYPE)) {
        for(unsigned int i = 0; i < xyz.size(); i++) xyz_type[i] = 0;
      }
      if(use_xyz) {
        write_xyzc(xyz, xyz_type, redptsout, scaleFac, hexfloat, high_precision);
      } else {
        write_uosc(xyz, xyz_type, redptsout, scaleFac*100.0 , hexfloat, high_precision);
      }

    } else if(use_color) {
      std::string data_string = red > 0 ? "color reduced" : "rgb";
      DataRGB xyz_color =
          (((DataRGB)source->get(data_string)).size() == 0) ?
          source->create(data_string, sizeof(unsigned char)*3*xyz.size()) :
          source->get(data_string);
      if (!(types & PointType::USE_COLOR)) {
          for(unsigned int i = 0; i < xyz.size(); i++) {
            xyz_color[i][0] = 0;
            xyz_color[i][1] = 0;
            xyz_color[i][2] = 0;
        }
      }
      if(use_xyz) {
        write_xyz_rgb(xyz, xyz_color, redptsout, scaleFac, hexfloat, high_precision);
      } else {
        write_uos_rgb(xyz, xyz_color, redptsout, scaleFac*100.0, hexfloat, high_precision);
      }

    } else if(use_normals) {
      std::string data_string = red > 0 ? "normal reduced" : "normal";
      DataNormal normals =
          (((DataNormal)source->get(data_string)).size() == 0) ?
          source->create(data_string, sizeof(double)*3*xyz.size()) :
          source->get(data_string);
      if(use_xyz) {
        write_xyz_normal(xyz, normals, redptsout, scaleFac, hexfloat, high_precision);
      } else {
        write_uos_normal(xyz, normals, redptsout, scaleFac*100.0, hexfloat, high_precision);
      }

    } else {
      if(use_xyz) {
        write_xyz(xyz, redptsout, scaleFac, hexfloat, high_precision);
      } else {
        write_uos(xyz, redptsout, scaleFac*100.0, hexfloat, high_precision);
      }

    }
    if(use_xyz) {
      writeTrajectoryXYZ(posesout, source->get_transMat(), false, scaleFac);
      writeTrajectoryXYZ(matricesout, source->get_transMat(), true, scaleFac);
    } else {
      writeTrajectoryUOS(posesout, source->get_transMat(), false, scaleFac*100.0);
      writeTrajectoryUOS(matricesout, source->get_transMat(), true, scaleFac*100.0);
    }

  }

  fclose(redptsout);
  posesout.close();
  posesout.clear();
  matricesout.close();
  matricesout.clear();
}
