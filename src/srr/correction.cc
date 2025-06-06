#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <queue>
using std::queue;
#include <vector>
using std::vector;
#include <sstream>
using std::stringstream;
using std::string;
using std::cout;
using std::endl;
using std::ifstream;

#include "slam6d/point_type.h"
#include "slam6d/globals.icc"
#include "slam6d/scan.h"
#include "slam6d/icp6Dquat.h"
#include "slam6d/icp6D.h"
#include "slam6d/graphSlam6D.h"
#include "slam6d/lum6Deuler.h"
#include "slam6d/io_types.h"
#include "simplereg.h"
#include "gapx6D.h"
#include "graphSlam6DL.h"
#include "lum6Deuler.h"
#include "ghelix6DQ2.h"
#include "linescan.h"
#include "continuousreg.h"
#include "lsegment.h"
#include "srr_program_options.h"

#include <signal.h>


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
      << endl
      << "Saving registration information in .frames files" << endl;
    for (unsigned int i = 0; i < LineScan::allLineScans.size(); i++) {
      if (i % 2) cout << ".";
      cout.flush();
      delete LineScan::allLineScans[i];
    }
    cout << endl;
  }
  exit(-1);
}

bool quickCheckCustomFilterStr(std::string &customFilter) {
  // custom filter set? quick check, needs to contain at least one ';'
  // (proper checking will be done case specific in pointfilter.cc)
  size_t pos = customFilter.find_first_of(";");
  if (pos != std::string::npos) {
    // check if customFilter is specified in file
    if (customFilter.find("FILE;") == 0) {
      std::string selection_file_name = customFilter.substr(5, customFilter.length());
      std::ifstream selectionfile;
      // open the input file
      selectionfile.open(selection_file_name, std::ios::in);

      if (!selectionfile.good()) {
        std::cerr << "Error loading custom filter file " << selection_file_name << "!" << std::endl;
        std::cerr << "Data will NOT be filtered.!" << std::endl;
        return false;
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
    return true;
  }
  else {
    // give a warning if custom filter has been inproperly specified
    if (customFilter.length() > 0){
      std::cerr << "Custom filter: specifying string has not been set properly, data will NOT be filtered." << std::endl;
    }
    return false;
  }
}


bool readScans(const std::string &dir, const IOType type,
    const int _start, const int _end,
    const double filterMinDist, const double filterMaxDist,
    std::vector<LineScan*> &linescans, const bool customFilterActive,std::string &customFilterStr) {

  bool scanserver = false;
  Scan::openDirectory(scanserver, dir, type, _start, _end);
  if(Scan::allScans.size() == 0) {
	  std::cerr << "No scans found. Did you use the correct format?" << std::endl;
    return false;
  }

  std::vector<Scan*>::iterator it = Scan::allScans.begin();
  for( ; it != Scan::allScans.end(); ++it) {
	  Scan* scan = *it;
	  if (customFilterActive) scan->setCustomFilter(customFilterStr);
	  scan->setRangeFilter(filterMaxDist, filterMinDist);
	  LineScan *ls;

	  //adds new LineScan to linescans
	  ls = new LineScan(scan);
	  if (!LineScan::error.empty()) {
	    std::cerr << LineScan::error << std::endl;
	    std::cerr << "Correction only up to this scan." << std::endl;
	    std::cout << "Read scans " << _start << " to " << scan->getIdentifier() << std::endl;
	    break;
	  }
  }
  Scan::closeDirectory();

  if (linescans.size() == 0) {
	  std::cerr << "Zero scans read" << std::endl;
    return false;
  }

  return true;

}

int main(int argc, char** argv) {

  signal (SIGSEGV, sigSEGVhandler);
  signal (SIGINT,  sigSEGVhandler);


  SRRSettings settings;
  parseArgs(argc, argv, settings);

  if(!boost::filesystem::is_directory(settings.dir)) {
    std::cerr << "Error: " << argv[0] << " does not support scans in archives." << std::endl;
    exit(-1);
  }

  string basedir = settings.dir;
  int start = settings.start;
  int end = settings.end;
  double voxelsize = settings.voxelsize;
  int octree = settings.octree;
  double mdm = settings.mdm;
  double mdml = settings.mdml;
  int iterations = settings.mni;
  int jterations = settings.mnj;
  int slamiterations = settings.slamiterations;

  int interval = settings.interval;
  int size = settings.size;

  int anim = settings.anim;

  double filterMinDist = settings.filterMinDist;
  double filterMaxDist = settings.filterMaxDist;
  int preRegInterval = settings.preRegInterval;

  // custom filter
  std::string customFilter = settings.customFilter;
  bool customFilterActive = quickCheckCustomFilterStr(customFilter);

  double epsICP = settings.epsICP;
  double odomWeightFactor = settings.odomWeightFactor;
  double gapOdomWeightFactor = settings.gapOdomWeightFactor;
  int    clpairs    = settings.clPairs;
  IOType type = settings.format;

  bool prereg = settings.prereg;
  int  preRegIter = settings.preRegIter;
  bool meta   = settings.meta;
  bool continue_processing = settings.continue_processing;

  /* With each call of LineScan::WriteFrames, the complete history of frames
   * kept in a string is written to the .frames files.
   * To reduce the size of the frames history kept, the idea is to append only
   * the new history since the last call to the .frames and then reset the
   * string.
   * However if not running in continue_processing mode, existing .frames files
   * need to be overwritten by the first call of WriteFrames.
   */
  bool appendFrames = continue_processing;

  if (mdml < 0) mdml = mdm;

  std::string segment_filename = settings.segment_filename;
  bool segment_mode = false;
  if(segment_filename != "") {
    std::cout << "running in segment mode" << std::endl;
    segment_mode = true;
  }

  cout << start << endl
	  << end << endl
	  << "voxelsize: " << voxelsize << endl
	  << "octree: " << octree << endl
	  << "icp max_dist_match: " << mdm << endl
	  << "slam max_dist_match: " << endl
	  << "filterMinDist: " << filterMinDist << endl
	  << "filterMaxDist: " << filterMaxDist << endl
	  << "anim " << anim << endl
          << "N: " << interval << endl;

  bool evenRepDistribution = false;
  if(settings.overlap != -1) {
    evenRepDistribution = true;
    size = interval/2 + interval*settings.overlap;
    std::cout << "relative overlap: " << settings.overlap << ", S: " << size << std::endl;
  }

  vector<LineScan*> &linescans = LineScan::allLineScans;

  if (continue_processing) Scan::continueProcessing();
  Scan::setProcessingCommand(argc, argv);

  std::vector<LSegment*> segments;

  if(!segment_mode) {
    std::cout << "Reading scans..." << std::endl;
    if(!readScans(basedir, type, start, end, filterMinDist, filterMaxDist, linescans, customFilterActive, customFilter)) {
      exit(1);
    }
    //segments.push_back(new LSegment(0,linescans.size()));
  } else {
    LSegmentParser segmentParser;
    std::cout << "reading segments ...";
    if(!segmentParser.loadCfg(segment_filename)) {
      std::cout << "failed" << std::endl;
      exit(1);
    }
    std::cout << "done." << std::endl;

    std::cout << "Reading scans..." << std::endl;
    for(std::pair<int, int> &segment : segmentParser.indices) {
      std::cout << "  ... segment from " << segment.first << " to " << segment.second << ". ";

      int firstLineScanIdx = linescans.size();
      if(!readScans(basedir, type, segment.first, segment.second,
                    filterMinDist, filterMaxDist, linescans, customFilterActive, customFilter)) {
        exit(1);
      }
      int lastLineScanIdx  = linescans.size()-1;
      segments.push_back(new LSegment(firstLineScanIdx,lastLineScanIdx));
      std::cout << "Done." << std::endl;
    }
  }


  icp6Dminimizer *my_icp6Dminimizer = new icp6D_QUAT(false);
  icp6D *icp = new icp6D(my_icp6Dminimizer,
					mdm,
					iterations,
					true, false, 1, true, anim,
					epsICP);

  ///////////////////////////////////////////////////////
  ///////////////   sequential registration
  if (prereg) {

    int parts = ceil((double) linescans.size() / preRegInterval);
    std::vector<std::shared_ptr<LScan>> lScans; // store structure of submaps
    for (unsigned int i = 0; i < parts; i++) {
      unsigned int firstIndex =  i    * preRegInterval;
      unsigned int lastIndex  = (i+1) * preRegInterval - 1;
      if(lastIndex >= linescans.size())
        lastIndex = linescans.size() - 1;
      unsigned int referenceIndex = firstIndex + (lastIndex - firstIndex)/2;
      std::shared_ptr<LScan> lscan(new LScan(firstIndex,lastIndex,referenceIndex));
      lScans.push_back(lscan);
    }

    for (int J=0; J < preRegIter; J++) {
      for (size_t i = 1; i < lScans.size(); i++) {

        // current linescan that we wish to
        // match and correct with previous
        // linescans
        if (meta) {
          LScan* metaLScan = new LScan(0,lScans[i-1]->getEnd(),lScans[i-1]->getRepresentative());
          preRegistration(linescans,
              metaLScan,
              lScans[i].get(),
              icp,
              voxelsize,
              octree);
          delete metaLScan;
        } else {
          preRegistration(linescans,
              lScans[i-1].get(),
              lScans[i].get(),
              icp,
              voxelsize,
              octree);
        }
      }

      if(preRegIter > 1) {
        std::cout << "pre-registration iteration " << J+1 << "/" << preRegIter << " done" << std::endl;
      }

      LineScan::writeFrames(basedir, "p", type, appendFrames);
      appendFrames = true;
    }

    lScans.clear();

    cout << "preregistration done!" << endl;
  }

  ////////// global semi rigid registration  // mdm was 5 for some reason?
  graphSlam6DL *luml = new lum6DEulerL(my_icp6Dminimizer, mdml, mdml);
  luml->set_quiet(true);

  ////////// global registration
  {
    graphSlam6D *my_graphSlam6D = new lum6DEuler(my_icp6Dminimizer, mdml, mdml);
    my_graphSlam6D->set_quiet(true);

    Graph *gr = 0;

    ////////////////////////// first frame
    double id[16];
    M4identity(id);
    for (unsigned int i = 0; i < linescans.size(); i++)  {
      linescans[i]->transform(id, Scan::ICP, 1 );
      if(!continue_processing) {
        linescans[i]->transform(id, Scan::ICP, 1 );
      }
    }

    double odomweight = odomWeightFactor * sqr(LineScan::allLineScans.size() );
    double gapodomweight = gapOdomWeightFactor * sqr(LineScan::allLineScans.size() );
    ///////////////////////

    for (int iters = 0; iters < jterations; iters++) {
      gr = 0;
      luml->setOdomWeight(odomweight,gapodomweight);
      if(!segment_mode)
        SSRR(linescans, my_graphSlam6D, luml, gr, 0,
            slamiterations,
            clpairs,
            mdml,
            interval,  // scan interval
            size,      // scan size
            voxelsize,
            octree);
      else  {
        for(LSegment* segment : segments) {
          segment->createLScans(interval,size, evenRepDistribution);
        }
        SSRR(linescans, segments, my_graphSlam6D, luml, gr, 0,
            slamiterations, clpairs, mdml, voxelsize, octree);
      }

      cout << "Matching with SLAM and distributing error... " << iters << endl;

      std::string suffix = prereg ? "p." + std::to_string(iters) : std::to_string(iters);
      LineScan::writeFrames(basedir, "p", type, appendFrames);
      appendFrames = true;
    }
    LineScan::setWriteFrames(true);
    for (unsigned int i = 0; i < linescans.size(); i++) {
      delete linescans[i];
    }
    return 0;
  }
  ///////////

}
