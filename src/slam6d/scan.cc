/*
 * scan implementation
 *
 * Copyright (C) by the 3DTK contributors
 *
 * Released under the GPL version 3.
 *
 */

/**
 * @file scan.cc
 * @brief the implementation for all scans (basic/managed)
 * @author Andreas Nuechter. Jacobs University Bremen gGmbH, Germany.
 * @author Kai Lingemann. Inst. of CS. University of Osnabrueck, Germany.
 * @author Dorit Borrmann. Jacobs University Bremen gGmbH, Germany.
 * @author Jan Elseberg. Jacobs University Bremen gGmbH, Germany.
 * @author Thomas Escher. Inst. of CS. University of Osnabrueck, Germany.
 */

#include "slam6d/scan.h"

#include "slam6d/basicScan.h"
#include "slam6d/managedScan.h"
#include "slam6d/metaScan.h"
#include "slam6d/searchTree.h"
#include "slam6d/kd.h"
#include "slam6d/Boctree.h"
#include "slam6d/globals.icc"

#include "slam6d/normals.h"

#ifdef WITH_METRICS
#include "slam6d/metrics.h"
#endif

#ifdef _MSC_VER
#define _NO_PARALLEL_READ
#endif

#ifdef __APPLE__
#define _NO_PARALLEL_READ
#endif

std::vector<Scan*> Scan::allScans;
unsigned int Scan::maxScanNr = 0;
bool Scan::scanserver = false;
bool Scan::continue_processing = false;
std::string Scan::processing_command;


void Scan::openDirectory(bool scanserver,
                         const std::string& path,
                         IOType type,
                         int start,
                         int end
#ifdef WITH_MMAP_SCAN
                         , boost::filesystem::path cache
#endif
                         )
{
  Scan::scanserver = scanserver;
  if (scanserver)
    ManagedScan::openDirectory(path, type, start, end);
  else
    BasicScan::openDirectory(path, type, start, end
#ifdef WITH_MMAP_SCAN
                             , cache
#endif
                            );
}

void Scan::openDirectory(dataset_settings& dss
#ifdef WITH_MMAP_SCAN
  , boost::filesystem::path cache
#endif
 )
{
  // custom filter set? quick check, needs to contain at least one ';'
// (proper checking will be done case specific in pointfilter.cc)
  size_t pos = dss.custom_filter.find_first_of(";");
  bool customFilterActive = false;
  std::string custom_filter = dss.custom_filter;
  if (pos != std::string::npos) {
    customFilterActive = true;

    // check if customFilter is specified in file
    if (dss.custom_filter.find("FILE;") == 0) {
      std::string selection_file_name = custom_filter.substr(5, custom_filter.length());
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
          custom_filter = custFilt.substr(0, custFilt.length() - 1);
        }
      }
      selectionfile.close();
    }
  }
  else {
    // give a warning if custom filter has been inproperly specified
    if (custom_filter.length() > 0) {
      std::cerr << "Custom filter: specifying string has not been set properly, data will NOT be filtered." << std::endl;
    }
  }

  std::vector<Scan*> valid_scans = Scan::allScans;
  int scan_nr = Scan::allScans.size();

  if (dss.use_scanserver)
    ManagedScan::openDirectory(dss.data_source, dss.format, dss.scan_numbers.min, dss.scan_numbers.max);
  else {
    BasicScan::openDirectory(dss
#ifdef WITH_MMAP_SCAN
      , cache
#endif
    );
  }

  for (; scan_nr < Scan::allScans.size(); ++scan_nr) {
    Scan* scan = Scan::allScans[scan_nr];
    scan->setRangeFilter(dss.distance_filter.max, dss.distance_filter.min);
    if (customFilterActive) scan->setCustomFilter(custom_filter);
    if (dss.sphere_radius > 0.0) scan->setRangeMutation(dss.sphere_radius);
    if (dss.octree_reduction_voxel > 0) {
      // scanserver differentiates between reduced for slam and
      // reduced for show, can handle both at the same time
      if (dss.use_scanserver) {
        if ((scan_nr - 1) % dss.skip_files != 0) delete scan;
        else {
          valid_scans.push_back(scan);
          dynamic_cast<ManagedScan*>(scan)->setShowReductionParameter(dss.octree_reduction_voxel, dss.octree_reduction_randomized_bucket);
        }
      }
      else {
        scan->setReductionParameter(dss.octree_reduction_voxel, dss.octree_reduction_randomized_bucket);
      }
    }
  }
  if (dss.use_scanserver && Scan::allScans.size() > valid_scans.size()) Scan::allScans = valid_scans;
}

void Scan::closeDirectory()
{
  if (scanserver)
    ManagedScan::closeDirectory();
  else
    BasicScan::closeDirectory();
}

void Scan::continueProcessing(bool continue_processing)
{
  Scan::continue_processing = continue_processing;
}

void Scan::setProcessingCommand(int argc, char** argv)
{
  std::string cmd;
  for(int arg_idx = 0; arg_idx < argc; arg_idx++) {
    cmd += argv[arg_idx];
    if(arg_idx < argc - 1) cmd += " ";
  }
  Scan::processing_command = cmd;
}

Scan::Scan()
{
  scanNr = maxScanNr;
  maxScanNr++;

  // pose and transformations
  for(size_t i = 0; i < 3; ++i) rPos[i] = 0;
  for(size_t i = 0; i < 3; ++i) rPosTheta[i] = 0;
  for(size_t i = 0; i < 4; ++i) rQuat[i] = 0;
  M4identity(transMat);
  M4identity(transMatOrg);
  M4identity(dalignxf);

  // trees and reduction methods
  nns_method = -1;
  kd = 0;

  // reduction on-demand
  reduction_voxelSize = 0.0;
  reduction_nrpts = 0;
  reduction_pointtype = PointType();

  // flags
  m_has_reduced = false;

  // octtree
  octtree_reduction_voxelSize = 0.0;
  octtree_voxelSize = 0.0;
  octtree_pointtype = PointType();
  octtree_loadOct = false;
  octtree_saveOct = false;
  octtree_autoOct = false;
}

Scan::~Scan()
{
  if (kd) delete kd;
}

void Scan::setReductionParameter(double voxelSize,
                                 int nrpts,
                                 PointType pointtype)
{
  reduction_voxelSize = voxelSize;
  reduction_nrpts = nrpts;
  reduction_pointtype = pointtype;
}

void Scan::setUpsamplingParameter(double voxelSize, double factor, PointType pointtype) {
    upsampling_voxelSize = voxelSize;
    upsampling_factor = factor;
    upsampling_pointtype = pointtype;
}

void Scan::setSearchTreeParameter(int nns_method, int bucketSize)
{
  searchtree_nnstype = nns_method;
  searchtree_bucketsize = bucketSize;
}

void Scan::setOcttreeParameter(double reduction_voxelSize,
                               double voxelSize,
                               PointType pointtype,
                               bool loadOct,
                               bool saveOct,
                               bool autoOct)
{
  octtree_reduction_voxelSize = reduction_voxelSize;
  octtree_voxelSize = voxelSize;
  octtree_pointtype = pointtype;
  octtree_loadOct = loadOct;
  octtree_saveOct = saveOct;
  octtree_autoOct = autoOct;
}

void Scan::clear(IODataType types)
{
  if(types & DATA_XYZ) clear("xyz");
  if(types & DATA_RGB) clear("rgb");
  if(types & DATA_REFLECTANCE) clear("reflectance");
  if(types & DATA_TEMPERATURE) clear("temperature");
  if(types & DATA_AMPLITUDE) clear("amplitude");
  if(types & DATA_TYPE) clear("type");
  if(types & DATA_DEVIATION) clear("deviation");
}

SearchTree* Scan::getSearchTree()
{
  // if the search tree hasn't been created yet, calculate everything
  if(kd == 0) {
    createSearchTree();
  }
  return kd;
}

void Scan::toGlobal() {
  calcReducedPoints();
  transform(transMatOrg, INVALID);
}

/**
 * Computes a search tree depending on the type.
 */
void Scan::createSearchTree()
{
  // multiple threads will call this function at the same time because they
  // all work on one pair of Scans, just let the first one
  // (who sees a null pointer)
  // do the creation
  boost::lock_guard<boost::mutex> lock(m_mutex_create_tree);
  if(kd != 0) return;

  // make sure the original points are created before starting the measurement
  DataXYZ xyz_orig(get("xyz reduced original"));

#ifdef WITH_METRICS
  Timer tc = ClientMetric::create_tree_time.start();
#endif //WITH_METRICS

  createSearchTreePrivate();

#ifdef WITH_METRICS
  ClientMetric::create_tree_time.end(tc);
#endif //WITH_METRICS
}

void Scan::calcReducedOnDemand()
{
  // multiple threads will call this function at the same time
  // because they all work on one pair of Scans,
  // just let the first one (who sees count as zero) do the reduction
  boost::lock_guard<boost::mutex> lock(m_mutex_reduction);
  if(m_has_reduced) return;

#ifdef WITH_METRICS
  Timer t = ClientMetric::on_demand_reduction_time.start();
#endif //WITH_METRICS

  calcReducedOnDemandPrivate();

  m_has_reduced = true;

#ifdef WITH_METRICS
  ClientMetric::on_demand_reduction_time.end(t);
#endif //WITH_METRICS
}

void Scan::calcNormalsOnDemand()
{
  // multiple threads will call this function at the same time
  // because they all work on one pair of Scans,
  // just let the first one (who sees count as zero) do the reduction
  boost::lock_guard<boost::mutex> lock(m_mutex_normals);
  if(m_has_normals) return;
  calcNormalsOnDemandPrivate();
  m_has_normals = true;
}

void Scan::copyReducedToOriginal()
{
#ifdef WITH_METRICS
  Timer t = ClientMetric::copy_original_time.start();
#endif //WITH_METRICS

  DataXYZ xyz_reduced(get("xyz reduced"));
  // check if we can create a large enough array. The maximum size_t on 32 bit
  // is around 4.2 billion which is too little for scans with more than 179
  // million points
  if (sizeof(size_t) == 4 && xyz_reduced.size() > ((size_t)(-1))/sizeof(double)/3) {
      throw std::runtime_error("Insufficient size of size_t datatype");
  }
  size_t size = xyz_reduced.size();
  DataXYZ xyz_reduced_orig(create("xyz reduced original",
                                  sizeof(double)*3*size));
  for(size_t i = 0; i < size; ++i) {
    for(size_t j = 0; j < 3; ++j) {
      xyz_reduced_orig[i][j] = xyz_reduced[i][j];
    }
  }

#ifdef WITH_METRICS
  ClientMetric::copy_original_time.end(t);
#endif //WITH_METRICS
}

void Scan::copyOriginalToReduced()
{
#ifdef WITH_METRICS
  Timer t = ClientMetric::copy_original_time.start();
#endif //WITH_METRICS

  DataXYZ xyz_reduced_orig(get("xyz reduced original"));
  // check if we can create a large enough array. The maximum size_t on 32 bit
  // is around 4.2 billion which is too little for scans with more than 179
  // million points
  if (sizeof(size_t) == 4 && xyz_reduced_orig.size() > ((size_t)(-1))/sizeof(double)/3) {
      throw std::runtime_error("Insufficient size of size_t datatype");
  }
  size_t size = xyz_reduced_orig.size();
  DataXYZ xyz_reduced(create("xyz reduced", sizeof(double)*3*size));
  for(size_t i = 0; i < size; ++i) {
    for(size_t j = 0; j < 3; ++j) {
      xyz_reduced[i][j] = xyz_reduced_orig[i][j];
    }
  }

#ifdef WITH_METRICS
  ClientMetric::copy_original_time.end(t);
#endif //WITH_METRICS
}


/**
 * Computes normals for all points<YY
 */
void Scan::calcNormals()
{
  std::cout << "calcNormals " << scanNr << std::endl;
  DataXYZ xyz(get("xyz"));
  // check if we can create a large enough array. The maximum size_t on 32 bit
  // is around 4.2 billion which is too little for scans with more than 179
  // million points
  if (sizeof(size_t) == 4 && xyz.size() > ((size_t)(-1))/sizeof(double)/3) {
      throw std::runtime_error("Insufficient size of size_t datatype");
  }
  DataNormal xyz_normals(create("normal", sizeof(double)*3*xyz.size()));
  if(xyz.size() == 0) throw
      std::runtime_error("Could not calculate reduced points, XYZ data is empty");

  std::vector<Point> points;
  points.reserve(xyz.size());
  std::vector<Point> normals;
  normals.reserve(xyz.size());
  for(size_t j = 0; j < xyz.size(); j++) {
    points.push_back(Point(xyz[j][0], xyz[j][1], xyz[j][2]));
  }
  const int K_NEIGHBOURS = 10;  //@FIXME
  calculateNormalsApxKNN(normals, points, K_NEIGHBOURS, get_rPos(), 1.0);
  // calculateNormalsKNN(normals, points, K_NEIGHBOURS, get_rPos());
  for (size_t i = 0; i < normals.size(); ++i) {
    xyz_normals[i][0] = normals[i].x;
    xyz_normals[i][1] = normals[i].y;
    xyz_normals[i][2] = normals[i].z;
  }
}

/**
 * Computes an octtree of the current scan, then getting the
 * reduced points as the centers of the octree voxels.
 */
void Scan::calcReducedPoints(bool rm_scatter)
{
#ifdef WITH_METRICS
  Timer t = ClientMetric::scan_load_time.start();
#endif //WITH_METRICS

  // get xyz to start the scan load, separated here for time measurement
  DataXYZ xyz(get("xyz"));
  DataXYZ xyz_normals(DataPointer(0, 0));
  if (reduction_pointtype.hasNormal()) {
    DataXYZ my_xyz_normals(get("normal"));
    xyz_normals =  my_xyz_normals;
  }
  DataReflectance reflectance(DataPointer(0, 0));
  if (reduction_pointtype.hasReflectance()) {
    DataReflectance my_reflectance(get("reflectance"));
    reflectance = my_reflectance;
  }
  DataType type(DataPointer(0, 0));
  if (reduction_pointtype.hasType()) {
    DataType my_type(get("type"));
    type = my_type;
  }
  DataRGB rgb(DataPointer(0, 0));
  if (reduction_pointtype.hasColor()) {
    DataRGB my_rgb(get("rgb"));
    rgb = my_rgb;
  }
  //Return if empty
  if(xyz.size() < 1) {
    DataXYZ xyz_reduced(create("xyz reduced", sizeof(double)*3*xyz.size()));
    if (reduction_pointtype.hasNormal()) {
      DataNormal normal_reduced(create("normal reduced", sizeof(double)*3*xyz.size()));
    }
    if (reduction_pointtype.hasReflectance()) {
      DataReflectance reflectance_reduced(create("reflectance reduced", sizeof(float)*reflectance.size()));
    }
    if (reduction_pointtype.hasType()) {
      DataType type_reduced(create("type reduced", sizeof(int)*type.size()));
    }
    if (reduction_pointtype.hasColor()) {
      /** @author: Fabian Arzberger
       * FIXME: this is inconsistent with the other field identifiers.
       * Why is color always refered to with get("rgb"), but when reduced it becomes get("color reduced") ?
       */
      DataRGB rgb_reduced(create("color reduced", sizeof(unsigned char)*3*xyz.size()));
    }
    return;
  }

#ifdef WITH_METRICS
    ClientMetric::scan_load_time.end(t);
    Timer tl = ClientMetric::calc_reduced_points_time.start();
#endif //WITH_METRICS

  if(reduction_voxelSize <= 0.0) {
    // copy the points
    // check if we can create a large enough array. The maximum size_t on 32 bit
    // is around 4.2 billion which is too little for scans with more than 179
    // million points
    if (sizeof(size_t) == 4 && xyz.size() > ((size_t)(-1))/sizeof(double)/3) {
        throw std::runtime_error("Insufficient size of size_t datatype");
    }
    DataXYZ xyz_reduced(create("xyz reduced", sizeof(double)*3*xyz.size()));
    for(size_t i = 0; i < xyz.size(); ++i) {
      for(size_t j = 0; j < 3; ++j) {
        xyz_reduced[i][j] = xyz[i][j];
      }
    }
    if (reduction_pointtype.hasReflectance()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 1.07
      // billion points
      if (sizeof(size_t) == 4 && reflectance.size() > ((size_t)(-1))/sizeof(float)) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataReflectance reflectance_reduced(create("reflectance reduced",
                                        sizeof(float)*reflectance.size()));
      for(size_t i = 0; i < xyz.size(); ++i) {
           reflectance_reduced[i] = reflectance[i];
        }
    }
    if (reduction_pointtype.hasType()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion
      if (sizeof(size_t) == 4 && type.size() > ((size_t)(-1))/sizeof(int)) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataType type_reduced(create("type reduced",
                                        sizeof(int)*type.size()));
      for(size_t i = 0; i < xyz.size(); ++i) {
           type_reduced[i] = type[i];
        }
    }
    if (reduction_pointtype.hasColor()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 1.4
      // billion points
      if (sizeof(size_t) == 4 && xyz.size() > ((size_t)(-1))/sizeof(unsigned char)/3) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      /** @author: Fabian Arzberger
       * FIXME: this is inconsistent with the other field identifiers.
       * Why is color always refered to with get("rgb"), but when reduced it becomes get("color reduced") ?
       */
      DataRGB rgb_reduced(create("color reduced", sizeof(unsigned char)*3*xyz.size()));
      for(size_t i = 0; i < xyz.size(); ++i) {
          for(size_t j = 0; j < 3; ++j) {
            rgb_reduced[i][j] = rgb[i][j];
          }
        }
    }
    if (reduction_pointtype.hasNormal()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 179
      // million points
      if (sizeof(size_t) == 4 && xyz.size() > ((size_t)(-1))/sizeof(double)/3) {
          throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataNormal normal_reduced(create("normal reduced",
                                       sizeof(double)*3*xyz.size()));
        for(size_t i = 0; i < xyz.size(); ++i) {
          for(size_t j = 0; j < 3; ++j) {
            normal_reduced[i][j] = xyz_normals[i][j];
          }
        }
    }

  } else {

    double **xyz_in = new double*[xyz.size()];
    for (size_t i = 0; i < xyz.size(); ++i) {
      xyz_in[i] = new double[reduction_pointtype.getPointDim()];
      size_t j = 0;
      for (; j < 3; ++j)
        xyz_in[i][j] = xyz[i][j];
      if (reduction_pointtype.hasReflectance())
        xyz_in[i][j++] = reflectance[i];
      if (reduction_pointtype.hasType())
        xyz_in[i][j++] = type[i];
      if (reduction_pointtype.hasColor())
        memcpy(&xyz_in[i][j++], &rgb[i][0], 3);
      if (reduction_pointtype.hasNormal())
        for (size_t l = 0; l < 3; ++l)
          xyz_in[i][j++] = xyz_normals[i][l];
    }

    // start reduction
    // build octree-tree from CurrentScan
    // put full data into the octtree
    BOctTree<double> *oct = new BOctTree<double>(xyz_in,
                                                 xyz.size(),
                                                 reduction_voxelSize,
                                                 reduction_pointtype);

    std::vector<double*> center;
    center.clear();
    if (reduction_nrpts != 0) {
      if (reduction_nrpts == -1) {
	   oct->GetOctTreeAvg(center);
      } else if (reduction_nrpts == 1) {
        oct->GetOctTreeRandom(center);
      } else {
     // if rm_scatter is true, all voxels with less than
     // reduction_nrpts will be removed
        oct->GetOctTreeRandom(center, reduction_nrpts, rm_scatter);
      }
    } else {
        oct->GetOctTreeCenter(center);
    }

    // storing it as reduced scan
    // check if we can create a large enough array. The maximum size_t on 32 bit
    // is around 4.2 billion which is too little for scans with more than 179
    // million points
    if (sizeof(size_t) == 4 && center.size() > ((size_t)(-1))/sizeof(double)/3) {
        throw std::runtime_error("Insufficient size of size_t datatype");
    }
    size_t size = center.size();
    DataXYZ xyz_reduced(create("xyz reduced", sizeof(double)*3*size));
    DataReflectance reflectance_reduced(DataPointer(0, 0));
    DataType type_reduced(DataPointer(0, 0));
    DataRGB rgb_reduced(DataPointer(0, 0));
    DataNormal normal_reduced(DataPointer(0, 0));
    if (reduction_pointtype.hasReflectance()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 1.07
      // billion points
      if (sizeof(size_t) == 4 && size > ((size_t)(-1))/sizeof(float)) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataReflectance my_reflectance_reduced(create("reflectance reduced",
                                                    sizeof(float)*size));
      reflectance_reduced = my_reflectance_reduced;
    }
    if (reduction_pointtype.hasType()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion
      if (sizeof(size_t) == 4 && size > ((size_t)(-1))/sizeof(int)) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataType my_type_reduced(create("type reduced",
                                                    sizeof(int)*size));
      type_reduced = my_type_reduced;
    }
    if (reduction_pointtype.hasColor()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 1.4
      // billion points
      if (sizeof(size_t) == 4 && size > ((size_t)(-1))/sizeof(unsigned char)/3) {
              throw std::runtime_error("Insufficient size of size_t datatype");
      }
      /** @author: Fabian Arzberger
       * FIXME: this is inconsistent with the other field identifiers.
       * Why is color always refered to with get("rgb"), but when reduced it becomes get("color reduced") ?
       */
      DataRGB my_rgb_reduced(create("color reduced",
                                          sizeof(unsigned char)*3*size));
      rgb_reduced = my_rgb_reduced;
    }
    if (reduction_pointtype.hasNormal()) {
      // check if we can create a large enough array. The maximum size_t on 32 bit
      // is around 4.2 billion which is too little for scans with more than 179
      // million points
      if (sizeof(size_t) == 4 && center.size() > ((size_t)(-1))/sizeof(double)/3) {
          throw std::runtime_error("Insufficient size of size_t datatype");
      }
      DataNormal my_normal_reduced(create("normal reduced",
                                          sizeof(double)*3*size));
      normal_reduced = my_normal_reduced;
    }
    for(size_t i = 0; i < size; ++i) {
      size_t j = 0;
      for (; j < 3; ++j)
        xyz_reduced[i][j] = center[i][j];
      if (reduction_pointtype.hasReflectance())
        reflectance_reduced[i] = center[i][j++];
      if (reduction_pointtype.hasType())
        type_reduced[i] = center[i][j++];
      if (reduction_pointtype.hasColor())
        memcpy(&rgb_reduced[i][0], &center[i][j++], 3);
      if (reduction_pointtype.hasNormal())
        for (size_t l = 0; l < 3; ++l)
          normal_reduced[i][l] = center[i][j++];
    }
    delete oct;
    for(size_t i = 0; i < xyz.size(); i++) {
      delete[] xyz_in[i];
    }
    delete[] xyz_in;
  }

#ifdef WITH_METRICS
    ClientMetric::calc_reduced_points_time.end(tl);
#endif //WITH_METRICS
}


/**
 * Computes the covariance matrix of points in leaf.
 * Attention: leaf must contain the mean of all points at the first position.
 * @param leaf A vector containing points of one leaf.
 * @return The covariance matrix (symmetric).
 */
NEWMAT::SymmetricMatrix Scan::calcCovarianceMatrix(std::vector<double*>& leaf) {

    unsigned int dim = upsampling_pointtype.getPointDim();

    NEWMAT::SymmetricMatrix C(dim);
    C = 0;

    // For all points in the leaf except the mean
    // (because of that we start at pointI=1)
    for(size_t pointI = 1; pointI < leaf.size(); pointI++) {

        // For all values in lower triangular covariance matrix C
        for(size_t k = 1; k <= 0.5 * dim * (dim + 1); k++) {

            // 2D indices for lower triangular matrix
            size_t i = std::ceil(std::sqrt(2.0 * k + 0.25) - 0.5);
            size_t j = k - (i - 1) * i / 2;

            // Summing up the covariance terms for current point and element
            C(i, j) += ((leaf.at(pointI)[i - 1] - leaf.front()[i - 1]) *
                   (leaf.at(pointI)[j - 1] - leaf.front()[j - 1]));

        }
    }

    // Multiply by 1 / (n - 1) for sample variance
    // --> additional mean does not belongs to the original sampling
    // --> 1 / (n - 2)
    C *= 1.0 / (leaf.size() - 2);

    return C;

}


/**
 * Computes an octree of the current scan, then upsample the
 * current scan such that the calculated points fit to the gaussian
 * distribution of one leaf containing the points of a certain area
 * of the scan.
 */
void Scan::calcUpsampledPoints() {

    unsigned int dim = upsampling_pointtype.getPointDim();
    size_t writeIndex = 0;

    DataXYZ xyz(get("xyz"));
    DataXYZ xyz_upsampled(create("xyz upsampled", sizeof(double)*upsampling_factor*dim*xyz.size()));

    // Copy all points from xyz to xyz_in (for octree) and xyz_upsampled
    double **xyz_in = new double*[xyz.size()];
    for (size_t i = 0; i < xyz.size(); ++i) {
        xyz_in[i] = new double[upsampling_pointtype.getPointDim()];
        for (size_t j = 0; j < dim; ++j) {
            xyz_in[i][j] = xyz[i][j];
            xyz_upsampled[i][j] = xyz[i][j];
        }
        ++writeIndex;
    }

    // start upsampling
    // build octree-tree from CurrentScan
    // put full data into the octree
    BOctTree<double> *oct = new BOctTree<double>(xyz_in,
                                                 xyz.size(),
                                                 upsampling_voxelSize,
                                                 upsampling_pointtype);

    // Determine all leafs and its related points
    std::vector<std::vector<double*>> leafPoints;
    oct->AllLeafPointsWithAvg(leafPoints);


    // Calculate upsampled points for every leaf considering its gaussian distribution
    for(std::vector<double*> leaf : leafPoints) {

        // Check if real leaf size > 3 has to be check for 4 because leaf contains mean at position 0
        NEWMAT::SymmetricMatrix C = (leaf.size() > 4) ? calcCovarianceMatrix(leaf) : NEWMAT::SymmetricMatrix(NEWMAT::IdentityMatrix(dim));
        NEWMAT::LowerTriangularMatrix L = Cholesky(C);

        // To create normal distributed point (for Box-Muller transform)
        std::default_random_engine generator;
        std::normal_distribution<double> distribution(0,1);

        // Create column vector and copy mean values to it
        NEWMAT::ColumnVector currentMean(dim);
        for(size_t i = 0; i < dim; i++)
            currentMean(i+1) = leaf.front()[i];

        // Create normal distributed point x by using Box-Muller transorm to calculate a point fitting to the gaussian distribution of the data set
        for(size_t i = 0; i < (leaf.size() - 1) * (upsampling_factor - 1); i++) {

            NEWMAT::ColumnVector x(dim);
            for(size_t j = 1; j <= dim; j++)
                x(j) = distribution(generator);

            NEWMAT::ColumnVector y = L * x + currentMean;

            for(size_t j = 0; j < dim; j++)
                xyz_upsampled[writeIndex][j] = y(j+1);
            ++writeIndex;

        }

    }

    // Delete memory on heap
    delete oct;
    for(size_t i = 0; i < xyz.size(); i++) {
        delete[] xyz_in[i];
    }
    delete[] xyz_in;

}


/**
 * Merges the scan's intrinsic coordinates with the robot position.
 * @param prevScan The scan that's transformation is extrapolated,
 * i.e., odometry extrapolation
 *
 * For additional information see the following paper (jfr2007.pdf):
 *
 * Andreas Nüchter, Kai Lingemann, Joachim Hertzberg, and Hartmut Surmann,
 * 6D SLAM - 3D Mapping Outdoor Environments Journal of Field Robotics (JFR),
 * Special Issue on Quantitative Performance Evaluation of Robotic and Intelligent
 * Systems, Wiley & Son, ISSN 1556-4959, Volume 24, Issue 8-9, pages 699 - 722,
 * August/September, 2007
 *
 */
void Scan::mergeCoordinatesWithRoboterPosition(Scan* prevScan)
{
  double tempMat[16], deltaMat[16];
  M4inv(prevScan->get_transMatOrg(), tempMat);
  MMult(prevScan->get_transMat(), tempMat, deltaMat);
  // apply delta transformation of the previous scan
  transform(deltaMat, INVALID);
}

/**
 * The method transforms all points with the given transformation matrix.
 */
void Scan::transformAll(const double alignxf[16])
{
  DataXYZ xyz(get("xyz"));
  size_t i=0 ;
  //  #pragma omp parallel for
  for(; i < xyz.size(); ++i) {
    transform3(alignxf, xyz[i]);
  }
  // TODO: test for ManagedScan compability,
  // may need a touch("xyz") to mark saving the new values
}

//! Internal function of transform which alters the reduced points
void Scan::transformReduced(const double alignxf[16])
{
#ifdef WITH_METRICS
  Timer t = ClientMetric::transform_time.start();
#endif //WITH_METRICS

  DataXYZ xyz_reduced(get("xyz reduced"));
  size_t i=0;
  // #pragma omp parallel for
  for( ; i < xyz_reduced.size(); ++i) {
    transform3(alignxf, xyz_reduced[i]);
  }

  if (reduction_pointtype.hasNormal()) {
    DataNormal normal_reduced(get("normal reduced"));
    for (size_t i = 0; i < normal_reduced.size(); ++i) {
      transform3normal(alignxf, normal_reduced[i]);
    }
  }


#ifdef WITH_METRICS
  ClientMetric::transform_time.end(t);
#endif //WITH_METRICS
}

//! Internal function of transform which handles the matrices
void Scan::transformMatrix(const double alignxf[16])
{
  double tempxf[16];

  // apply alignxf to transMat and update pose vectors
  MMult(alignxf, transMat, tempxf);
  memcpy(transMat, tempxf, sizeof(transMat));
  Matrix4ToEuler(transMat, rPosTheta, rPos);
  Matrix4ToQuat(transMat, rQuat);

#ifdef DEBUG
  std::cerr << "(" << rPos[0] << ", " << rPos[1] << ", " << rPos[2] << ", "
       << rPosTheta[0] << ", " << rPosTheta[1] << ", " << rPosTheta[2] << ")" << std::endl;

  std::cerr << transMat << std::endl;
#endif

  // apply alignxf to dalignxf
  MMult(alignxf, dalignxf, tempxf);
  memcpy(dalignxf, tempxf, sizeof(transMat));
}

/**
 * Transforms the scan by a given transformation and writes a new frame.
 * The idea is to write for every transformation in all files,
 * such that the show program is able to determine,
 * which scans have to be drawn in which color. Hidden scans
 * (or later processed scans) are written with INVALID.
 *
 * @param alignxf Transformation matrix
 * @param colour Specifies which colour should the written to the frames file
 * @param islum Is the transformtion part of LUM, i.e., all scans
 *              are transformed?
 *              In this case only LUM transformation is stored, otherwise all
 *              scans are processed
 *        -1  no transformation is stored
 *         0  ICP transformation
 *         1  LUM transformation, all scans except last scan
 *         2  LUM transformation, last scan only
 */
void Scan::transform(const double alignxf[16], const AlgoType type, int islum)
{
  MetaScan* meta = dynamic_cast<MetaScan*>(this);

  if(meta) {
    for(size_t i = 0; i < meta->size(); ++i) {
      meta->getScan(i)->transform(alignxf, type, -1);
    }
  }

#ifdef TRANSFORM_ALL_POINTS
  transformAll(alignxf);
#endif //TRANSFORM_ALL_POINTS

#ifdef DEBUG
  std::cerr << alignxf << std::endl;
  std::cerr << "(" << rPos[0] << ", " << rPos[1] << ", " << rPos[2] << ", "
       << rPosTheta[0] << ", " << rPosTheta[1] << ", " << rPosTheta[2]
       << ") ---> ";
#endif

  // transform points
  transformReduced(alignxf);

  // update matrices
  transformMatrix(alignxf);

  // store transformation in frames
  if(type != INVALID) {
#ifdef WITH_METRICS
    Timer t = ClientMetric::add_frames_time.start();
#endif //WITH_METRICS
    bool in_meta;
    MetaScan* meta = dynamic_cast<MetaScan*>(this);
    int found = 0;
    size_t scans_size = allScans.size();

    switch (islum) {
    case -1:
      // write no tranformation
      break;
    case 0:
      for(size_t i = 0; i < scans_size; ++i) {
        Scan* scan = allScans[i];
        in_meta = false;
        if(meta) {
          for(size_t j = 0; j < meta->size(); ++j) {
            if(meta->getScan(j) == scan) {
              found = i;
              in_meta = true;
            }
          }
        }

        if(scan == this || in_meta) {
          found = i;
          scan->addFrame(type);
        } else {
          if(found == 0) {
            scan->addFrame(ICPINACTIVE);
          } else {
            scan->addFrame(INVALID);
          }
        }
      }
      break;
    case 1:
      addFrame(type);
      break;
    case 2:
      for(size_t i = 0; i < scans_size; ++i) {
        Scan* scan = allScans[i];
        if(scan == this) {
          found = i;
          addFrame(type);
          allScans[0]->addFrame(type);
          continue;
        }
        if (found != 0) {
          scan->addFrame(INVALID);
        }
      }
      break;
    default:
      std::cerr << "invalid point transformation mode" << std::endl;
    }

#ifdef WITH_METRICS
    ClientMetric::add_frames_time.end(t);
#endif //WITH_METRICS
  }
}

/**
 * Transforms the scan by a given transformation and writes a new frame.
 * The idea is to write for every transformation in all files, such that
 * the show program is able to determine, whcih scans have to be drawn
 * in which color. Hidden scans (or later processed scans) are written
 * with INVALID.
 *
 * @param alignQuat Quaternion for the rotation
 * @param alignt    Translation vector
 * @param colour Specifies which colour should the written to the frames file
 * @param islum Is the transformtion part of LUM, i.e., all scans are
 *              transformed?
 *              In this case only LUM transformation is stored, otherwise
 *              all scans are processed
 *        -1  no transformation is stored
 *         0  ICP transformation
 *         1  LUM transformation, all scans except last scan
 *         2  LUM transformation, last scan only
 */
void Scan::transform(const double alignQuat[4], const double alignt[3],
                     const AlgoType type, int islum)
{
  double alignxf[16];
  QuatToMatrix4(alignQuat, alignt, alignxf);
  transform(alignxf, type, islum);
}

/**
 * Transforms the scan, so that the given Matrix
 * prepresent the next pose.
 *
 * @param alignxf Transformation matrix to which this scan will be set to
 * @param islum Is the transformation part of LUM?
 */
void Scan::transformToMatrix(double alignxf[16], const AlgoType type, int islum)
{
  double tinv[16];
  M4inv(transMat, tinv);
  transform(tinv, INVALID);
  transform(alignxf, type, islum);
}

/**
 * Transforms the scan, so that the given Euler angles
 * prepresent the next pose.
 *
 * @param rP Translation to which this scan will be set to
 * @param rPT Orientation as Euler angle to which this scan will be set
 * @param islum Is the transformation part of LUM?
 */
void Scan::transformToEuler(double rP[3],
                            double rPT[3],
                            const AlgoType type,
                            int islum)
{
#ifdef WITH_METRICS
  // called in openmp context in lum6Deuler.cc:422
  ClientMetric::transform_time.set_threadsafety(true);
  ClientMetric::add_frames_time.set_threadsafety(true);
#endif //WITH_METRICS

  double tinv[16];
  double alignxf[16];
  M4inv(transMat, tinv);
  transform(tinv, INVALID);
  EulerToMatrix4(rP, rPT, alignxf);
  transform(alignxf, type, islum);

#ifdef WITH_METRICS
  ClientMetric::transform_time.set_threadsafety(false);
  ClientMetric::add_frames_time.set_threadsafety(false);
#endif //WITH_METRICS
}

/**
 * Transforms the scan, so that the given Euler angles
 * prepresent the next pose.
 *
 * @param rP Translation to which this scan will be set to
 * @param rPQ Orientation as Quaternion to which this scan will be set
 * @param islum Is the transformation part of LUM?
 */
void Scan::transformToQuat(double rP[3],
                           double rPQ[4],
                           const AlgoType type,
                           int islum)
{
  double tinv[16];
  double alignxf[16];
  M4inv(transMat, tinv);
  transform(tinv, INVALID);
  QuatToMatrix4(rPQ, rP, alignxf);
  transform(alignxf, type, islum);
}

/**
 * Calculates Source\Target
 * Calculates a set of corresponding point pairs and returns them. It
 * computes the k-d trees and deletes them after the pairs have been
 * found. This slow function should be used only for testing
 *
 * @param pairs The resulting point pairs (vector will be filled)
 * @param Target The scan to whiche the points are matched
 * @param thread_num number of the thread (for parallelization)
 * @param rnd randomized point selection
 * @param max_dist_match2 maximal allowed distance for matching
 */

void Scan::getNoPairsSimple(std::vector <double*> &diff,
                            Scan* Source, Scan* Target,
                            int thread_num,
                            double max_dist_match2)
{
  DataXYZ xyz_reduced(Source->get("xyz reduced"));
  KDtree* kd = new KDtree(
                 PointerArray<double>(Target->get("xyz reduced")).get(),
                 Target->size<DataXYZ>("xyz reduced"));

  std::cout << "Max: " << max_dist_match2 << std::endl;
  for (size_t i = 0; i < xyz_reduced.size(); i++) {

    double p[3];
    p[0] = xyz_reduced[i][0];
    p[1] = xyz_reduced[i][1];
    p[2] = xyz_reduced[i][2];


    double *closest = kd->FindClosest(p, max_dist_match2, thread_num);
    if (!closest) {
         diff.push_back(xyz_reduced[i]);
         //diff.push_back(closest);
    }
  }

  delete kd;
}

/**
 * Calculates a set of corresponding point pairs and returns them. It
 * computes the k-d trees and deletes them after the pairs have been
 * found. This slow function should be used only for testing
 *
 * @param pairs The resulting point pairs (vector will be filled)
 * @param Source The scan whose points are matched to Targets' points
 * @param Target The scan to whiche the points are matched
 * @param thread_num number of the thread (for parallelization)
 * @param rnd randomized point selection
 * @param max_dist_match2 maximal allowed distance for matching
 */
void Scan::getPtPairsSimple(std::vector <PtPair> *pairs,
                            Scan* Source, Scan* Target,
                            int thread_num,
                            int rnd, double max_dist_match2,
                            double *centroid_m, double *centroid_d)
{
  KDtree* kd = new KDtree(
                 PointerArray<double>(Source->get("xyz reduced")).get(),
                 Source->size<DataXYZ>("xyz reduced"));
  DataXYZ xyz_reduced(Target->get("xyz reduced"));

  for (size_t i = 0; i < xyz_reduced.size(); i++) {
    // take about 1/rnd-th of the numbers only
    if (rnd > 1 && rand(rnd) != 0) continue;

    double p[3];
    p[0] = xyz_reduced[i][0];
    p[1] = xyz_reduced[i][1];
    p[2] = xyz_reduced[i][2];

    double *closest = kd->FindClosest(p, max_dist_match2, thread_num);
    if (closest) {
      centroid_m[0] += closest[0];
      centroid_m[1] += closest[1];
      centroid_m[2] += closest[2];
      centroid_d[0] += p[0];
      centroid_d[1] += p[1];
      centroid_d[2] += p[2];
      PtPair myPair(closest, p);
      pairs->push_back(myPair);
    }
  }
  centroid_m[0] /= pairs[thread_num].size();
  centroid_m[1] /= pairs[thread_num].size();
  centroid_m[2] /= pairs[thread_num].size();
  centroid_d[0] /= pairs[thread_num].size();
  centroid_d[1] /= pairs[thread_num].size();
  centroid_d[2] /= pairs[thread_num].size();

  delete kd;
}


/**
 * Calculates a set of corresponding point pairs and returns them.
 * The function uses the k-d trees stored the the scan class, thus
 * the function createTrees and deletTrees have to be called before
 * resp. afterwards.
 * Here we implement the so called "fast corresponding points"; k-d
 * trees are not recomputed, instead the apply the inverse transformation
 * to to the given point set.
 *
 * @param pairs The resulting point pairs (vector will be filled)
 * @param Source The scan whose points are matched to Targets' points
 * @param Target The scan to whiche the points are matched
 * @param thread_num number of the thread (for parallelization)
 * @param rnd randomized point selection
 * @param max_dist_match2 maximal allowed distance for matching
 * @return a set of corresponding point pairs
 */
void Scan::getPtPairs(std::vector <PtPair> *pairs,
                      Scan* Source, Scan* Target,
                      int thread_num,
                      int rnd, double max_dist_match2, double &sum,
                      double *centroid_m, double *centroid_d,
                      PairingMode pairing_mode)
{
  // initialize centroids
  for(size_t i = 0; i < 3; ++i) {
    centroid_m[i] = 0;
    centroid_d[i] = 0;
  }

  // get point pairs
  DataXYZ xyz_reduced(Target->get("xyz reduced"));
  DataNormal normal_reduced(DataPointer(0, 0));
  if ((pairing_mode == CLOSEST_POINT_ALONG_NORMAL_SIMPLE) || (pairing_mode == CLOSEST_PLANE_SIMPLE)) {
    DataNormal my_normals(Target->get("normal reduced"));
    normal_reduced =  my_normals;
  }
  Source->getSearchTree()->getPtPairs(pairs, Source->dalignxf,
                                      xyz_reduced,
                                      normal_reduced,
                                      0,
                                      xyz_reduced.size(),
                                      thread_num,
                                      rnd,
                                      max_dist_match2,
                                      sum,
                                      centroid_m, centroid_d,
                                      pairing_mode);

  // normalize centroids
  size_t size = pairs->size();
  if(size != 0) {
    for(size_t i = 0; i < 3; ++i) {
      centroid_m[i] /= size;
      centroid_d[i] /= size;
    }
  }
}


/**
 * Calculates a set of corresponding point pairs and returns them.
 * The function uses the k-d trees stored the the scan class, thus
 * the function createTrees and delteTrees have to be called before
 * resp. afterwards.
 *
 * @param pairs The resulting point pairs (vector will be filled)
 * @param Source The scan whose points are matched to Targets' points
 * @param Target The scan to whiche the points are matched
 * @param thread_num The number of the thread that is computing ptPairs
 *                   in parallel
 * @param step The number of steps for parallelization
 * @param rnd randomized point selection
 * @param max_dist_match2 maximal allowed distance for matching
 * @param sum The sum of distances of the points
 *
 * These intermediate values are for the parallel ICP algorithm
 * introduced in the paper
 * "The Parallel Iterative Closest Point Algorithm"
 *  by Langis / Greenspan / Godin, IEEE 3DIM 2001
 *
 */
void Scan::getPtPairsParallel(std::vector <PtPair> *pairs,
                              Scan* Source, Scan* Target,
                              int thread_num, int step,
                              int rnd, double max_dist_match2,
                              double *sum,
                              double centroid_m[OPENMP_NUM_THREADS][3],
                              double centroid_d[OPENMP_NUM_THREADS][3],
                              PairingMode pairing_mode)
{
  // initialize centroids
  for(size_t i = 0; i < 3; ++i) {
    centroid_m[thread_num][i] = 0;
    centroid_d[thread_num][i] = 0;
  }

  // get point pairs
  SearchTree* search = Source->getSearchTree();
  // differentiate between a meta scan (which has no reduced points)
  // and a normal scan
  // if Source is also a meta scan it already has a special meta-kd-tree
  MetaScan* meta = dynamic_cast<MetaScan*>(Target);
  if(meta) {
    for(size_t i = 0; i < meta->size(); ++i) {
      // determine step for each scan individually
      DataXYZ xyz_reduced(meta->getScan(i)->get("xyz reduced"));
      DataNormal normal_reduced(DataPointer(0, 0));
      if ((pairing_mode == CLOSEST_POINT_ALONG_NORMAL_SIMPLE) || (pairing_mode == CLOSEST_PLANE_SIMPLE)) {
        DataNormal my_normals(meta->getScan(i)->get("normal reduced"));
        normal_reduced =  my_normals;
      }

      size_t max = xyz_reduced.size();
      size_t step = ceil(max / (double)OPENMP_NUM_THREADS);
      size_t endindex = thread_num == (OPENMP_NUM_THREADS - 1) ? max : step * thread_num + step;
      // call ptpairs for each scan and accumulate ptpairs, centroids and sum
      search->getPtPairs(&pairs[thread_num], Source->dalignxf,
                         xyz_reduced, normal_reduced,
                         step * thread_num, endindex,
                         thread_num,
                         rnd, max_dist_match2, sum[thread_num],
                         centroid_m[thread_num], centroid_d[thread_num],
                         pairing_mode);
    }
  } else {
    DataXYZ xyz_reduced(Target->get("xyz reduced"));
    DataNormal normal_reduced(DataPointer(0, 0));
    if ((pairing_mode == CLOSEST_POINT_ALONG_NORMAL_SIMPLE) || (pairing_mode == CLOSEST_PLANE_SIMPLE)) {
      DataNormal my_normals(Target->get("normal reduced"));
      normal_reduced =  my_normals;
    }
    size_t endindex = thread_num == (OPENMP_NUM_THREADS - 1) ?  xyz_reduced.size() : step * thread_num + step;
    search->getPtPairs(&pairs[thread_num], Source->dalignxf,
                       xyz_reduced, normal_reduced,
                       thread_num * step, endindex,
                       thread_num,
                       rnd, max_dist_match2, sum[thread_num],
                       centroid_m[thread_num], centroid_d[thread_num],
                       pairing_mode);
  }

  // normalize centroids
  size_t size = pairs[thread_num].size();
  if(size != 0) {
    for(size_t i = 0; i < 3; ++i) {
      centroid_m[thread_num][i] /= size;
      centroid_d[thread_num][i] /= size;
    }
  }
}

size_t Scan::getMaxCountReduced(ScanVector& scans)
{
  size_t max = 0;
  for(std::vector<Scan*>::iterator it = scans.begin();
      it != scans.end();
      ++it) {
    size_t count = (*it)->size<DataXYZ>("xyz reduced");
    if(count > max)
      max = count;
  }
  return max;
}
