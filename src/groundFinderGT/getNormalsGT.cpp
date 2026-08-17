/**
 * @file
 * @brief Program to compute ground truth normal vectors for every pose of a trajectory.
 *
    Uses e.g. RIEGL scan as reference and computes for every pose of the trajectory the nearest ground point & its k
 nearest neighbors to get the related ground truth normal vector from the RIEGL scan. Specifically designed for
 Spherical Mobile Mapping System.
 *
 * @author Tim Schubert. University of Wuerzburg, Germany.
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "slam6d/frame.h"
#include "slam6d/globals.icc"
#include "slam6d/icp6D.h"
#include "slam6d/icp6Dsvd.h"
#include "slam6d/kd.h"
#include "slam6d/normals.h"
#include "slam6d/scan.h"
#include "slam6d/searchTree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace po = boost::program_options;

void parseArgs(int argc, char **argv, std::string &dir, int &start, int &end, IOType &type, int &knn,
	       double &knn_radius, double &max_dist, double &max_dist_match, int &max_iters, std::string &outfile,
	       std::string &custom_filter, bool &save_matched, double &red, double &red_input, int &octree,
	       bool &use_frames, bool &skip_matching)
{

	po::options_description generic("Generic options");
	generic.add_options()("help,h", "produce help message");

	po::options_description input("Input options");
	input.add_options()("start,s", po::value<int>(&start)->default_value(0), "start scan number (Pandar)")(
	    "end,e", po::value<int>(&end)->default_value(-1),
	    "end scan number for Pandar scans (-1 means process up to last Pandar scan)")(
	    "format,f", po::value<IOType>(&type)->default_value(UOS, "uos"), "scan format")(
	    "knn,k", po::value<int>(&knn)->default_value(200), "number of nearest neighbors for normal calculation")(
	    "knn-radius, R", po::value<double>(&knn_radius)->default_value(5.0),
	    "search for knn in this radius [cm]")("max-dist,d", po::value<double>(&max_dist)->default_value(5.0),
						  "maximum search distance along gravity direction (meters)")(
	    "max-dist-match,m", po::value<double>(&max_dist_match)->default_value(20.0),
	    "maximum distance for ICP matching (cm)")(
	    "max-iterations,i", po::value<int>(&max_iters)->default_value(100), "maximum ICP iterations")(
	    "output,o", po::value<std::string>(&outfile)->default_value("normals_gt.csv"),
	    "output CSV file")("save-matched", po::bool_switch(&save_matched),
			       "save ICP-matched Pandar scans for debugging (in matched_scans/ subdirectory)")(
	    "reduce,r", po::value<double>(&red)->default_value(-1.0),
	    "turns on octree based point reduction for RIEGL scan (voxel size= <arg>)")(
	    "reduce_input", po::value<double>(&red_input)->default_value(-1.0))(
	    "octree,O", po::value<int>(&octree)->default_value(0),
	    "Use randomized octree based point reduction for RIEGL scan (pts per voxel=<arg>")(
	    "continue", po::bool_switch(&use_frames)->default_value(false),
	    "Use pose specified in .frames file instead of .pose file.")(
	    "skip-matching", po::bool_switch(&skip_matching)->default_value(false),
	    "Skip the ICP matching part and only compute normals.")(
	    "customFilter,u", po::value<std::string>(&custom_filter),
	    "Apply a custom filter. Filter mode and data are specified as a "
	    "semicolon-seperated string:\n"
	    "\"{filterMode};{nrOfParams}[;param1][;param2][...]\"\n"
	    "Multiple filters can be specified in a file (syntax in file is same as "
	    "direct specification)\n"
	    "\"FILE;{fileName}\"\n"
	    "See filter implementation in src/slam6d/pointfilter.cc for more detail.");

	po::options_description hidden("Hidden options");
	hidden.add_options()("input-dir", po::value<std::string>(&dir), "input dir");

	po::options_description all;
	all.add(generic).add(input).add(hidden);

	po::options_description cmdline_options;
	cmdline_options.add(generic).add(input);

	po::positional_options_description pd;
	pd.add("input-dir", 1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), vm);

	if (vm.count("help")) {
		std::cout << cmdline_options << std::endl;
		exit(0);
	}

	po::notify(vm);

#ifndef _MSC_VER
	if (dir[dir.length() - 1] != '/')
		dir = dir + "/";
#else
	if (dir[dir.length() - 1] != '\\')
		dir = dir + "\\";
#endif
}

void validate(boost::any &v, const std::vector<std::string> &values, IOType *, int)
{
	if (values.size() == 0)
		throw std::runtime_error("Invalid model specification");
	std::string arg = values.at(0);
	try {
		v = formatname_to_io_type(arg.c_str());
	} catch (...) {
		throw std::runtime_error("Format " + arg + " unknown.");
	}
}

// Returns the angle in degrees between the normal and world-up (+Y)
// Assumes normal is already sign-corrected to point "upward" (normal[1] >= 0).
double angleFromVertical(const double normal[3])
{
	double len = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
	if (len < 1e-12)
		return 90.0; // invalid due to zero length

	double cos_angle = normal[1] / len; // dot(normal, up) / (|normal|*|up|)
	if (cos_angle > 1.0)
		cos_angle = 1.0;
	if (cos_angle < -1.0)
		cos_angle = -1.0;

	return std::acos(cos_angle) * 180.0 / M_PI;
}

void computeNormalFromKNN(const std::vector<Point> &neighbors, double normal[3])
{
	std::vector<Point> temp_points;
	for (size_t i = 0; i < neighbors.size(); i++) {
		temp_points.push_back(neighbors[i]);
	}

	double eigen[3];
	calculateNormal(temp_points, normal, eigen);

	// Ensure normal points UPWARDS (ground plane normal in 3DTK Y-up coordinate system)
	if (normal[1] < 0) {
		normal[0] = -normal[0];
		normal[1] = -normal[1];
		normal[2] = -normal[2];
	}
}

void clearScanBuffers(Scan *scan)
{
	if (!scan)
		return;
	scan->clear("xyz");
	scan->clear("xyz original");
	scan->clear("xyz reduced");
	scan->clear("xyz reduced original");
	scan->clear("xyz show");
	scan->clear("xyz show reduced");
	scan->clear("reflectance");
	scan->clear("amplitude");
	scan->clear("type");
	scan->clear("deviation");
	scan->clear("rgb");
	scan->clear("normals");
}

// Writes one transformation to the frames file
inline void writeTransformToFrameFile(std::ofstream &file, const double *mat, unsigned int type)
{
	for (int i = 0; i < 16; ++i) {
		file << mat[i] << " ";
	}
	file << type << "\n" << std::flush;
}

// saves only the final transformation but NOT an animation for a single scan
inline void saveLastFrame(Scan *scan, bool append)
{
	const double *last_matrix;
	Scan::AlgoType last_atype;
	scan->getFrame(scan->getFrameCount() - 1, last_matrix, last_atype);
	std::string filename = scan->getPath() + "scan" + scan->getIdentifier() + ".frames";
	std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;
	std::ofstream file(filename.c_str(), mode);
	writeTransformToFrameFile(file, last_matrix, last_atype);
	file.close();
}

// saves only the final transformation but NOT an animation for a vector of scans
inline void saveLastFrame(ScanVector &scans, bool append)
{
	for (size_t i = 0; i < scans.size(); ++i) {
		saveLastFrame(scans.at(i), append);
	}
}

//  Handling Segmentation faults and CTRL-C
void sigSEGVhandler(int v)
{
	static bool segfault = false;
	if (!segfault) {
		segfault = true;
		cout << endl
		     << "# **************************** #" << endl
		     << "  Segmentation fault or Ctrl-C" << endl
		     << "# **************************** #" << endl
		     << endl
		     << "Saving Frames... ";
		saveLastFrame(Scan::allScans, true);
		cout << " done!" << endl;
	}
	exit(-1);
}

int main(int argc, char *argv[])
{
	std::string dir;
	int start, end;
	IOType type;
	int knn;
	double knn_radius;
	double max_dist;
	double max_dist_match;
	int max_iters;
	std::string outfile;
	std::string custom_filter;
	bool save_matched = false;
	double red;
	double red_input;
	int octree;
	bool use_frames = false;
	bool skip_matching = false;

	parseArgs(argc, argv, dir, start, end, type, knn, knn_radius, max_dist, max_dist_match, max_iters, outfile,
		  custom_filter, save_matched, red, red_input, octree, use_frames, skip_matching);

	std::cout << "Loading scans from: " << dir << std::endl;

	if (use_frames)
		Scan::continueProcessing();

	Scan::openDirectory(false, dir, type, start, -1);

	if (Scan::allScans.empty()) {
		std::cerr << "Error: No scans found in directory " << dir << std::endl;
		return 1;
	}

	Scan *rieglScan = Scan::allScans[Scan::allScans.size() - 1]; // riegl is last in directory

	// If end scan index was specified, prune Pandar scans with IDs > end
	if (end != -1) {
		std::vector<Scan *> filtered_scans;
		for (size_t i = 0; i < Scan::allScans.size() - 1; ++i) {
			Scan *scan = Scan::allScans[i];
			if (std::atoi(scan->getIdentifier()) <= end) {
				filtered_scans.push_back(scan);
			} else {
				delete scan;
			}
		}
		filtered_scans.push_back(rieglScan);
		Scan::allScans = filtered_scans;
	}

	if (Scan::allScans.size() < 2) {
		std::cerr << "Error: Need at least 2 scans (Pandar scan(s) + RIEGL reference scan)!" << std::endl;
		return 1;
	}

	std::cout << "Using scan " << rieglScan->getIdentifier() << " as RIEGL reference scan" << std::endl;

	// Config RIEGL scan
	rieglScan->setRangeFilter(-1, -1);
	if (custom_filter.length() > 0) {
		rieglScan->setCustomFilter(custom_filter);
	}
	rieglScan->setReductionParameter(red, octree);
	rieglScan->setSearchTreeParameter(simpleKD, 20);

	// Config Pandar scans
	for (size_t i = 0; i < Scan::allScans.size() - 1; i++) {
		Scan *scan = Scan::allScans[i];
		scan->setRangeFilter(-1, -1); // no filter by default
		if (custom_filter.length() > 0) {
			scan->setCustomFilter(custom_filter);
		}
		scan->setReductionParameter(red_input, octree);
		scan->setSearchTreeParameter(simpleKD, 20);
	}

	if (!skip_matching) {
		// 1. Setup ICP minimizer
		icp6Dminimizer *my_icp6Dminimizer = new icp6D_SVD(true);

		// 2. Initial alignment: Match RIEGL to first Pandar scan to bring RIEGL into world frame
		icp6D *icp_init = new icp6D(my_icp6Dminimizer, max_dist_match, max_iters, false, false, 0, false, -2,
					    0.00001, simpleKD, false, false, 1);

		Scan *firstPandarScan = Scan::allScans[0];

		firstPandarScan->createSearchTree();

		std::cout << "Matching RIEGL to first Pandar scan..." << std::endl;
		int initial_iters = icp_init->match(firstPandarScan, rieglScan, CLOSEST_POINT);

		saveLastFrame(Scan::allScans, false);

		clearScanBuffers(firstPandarScan);

		// Clean up first scan search tree
		delete icp_init;

		// 3. Setup CAD ICP instance for remaining Pandar scans with cad matching
		icp6D *icp_cad = new icp6D(my_icp6Dminimizer, max_dist_match, max_iters, false, false, 0, false, -2,
					   0.00001, simpleKD, false, true, 1);

		std::cout << "\nRunning sequential CAD ICP..." << std::endl;

		// Build search tree for transformed RIEGL scan
		rieglScan->createSearchTree();

		for (size_t i = 1; i < Scan::allScans.size() - 1; i++) {
			Scan *pandarScan = Scan::allScans[i];

			pandarScan->mergeCoordinatesWithRoboterPosition(Scan::allScans[i - 1]);

			//: Align pandarScan against rieglScan
			icp_cad->match(rieglScan, pandarScan, CLOSEST_POINT);

			std::cout << "scan: " << i << endl;

			saveLastFrame(pandarScan, false);
			clearScanBuffers(pandarScan); // overload for only 1 scan
		}
		// saveLastFrame(Scan::allScans, false); // overload for all scans in directory

		delete icp_cad;
		delete my_icp6Dminimizer;
	} else {
		std::cout << "Skipping ICP matching. Using existing poses from .pose or .frames files." << std::endl;
	}

	/*
	NORMAL CALCULATION
	*/

	for (size_t i = 0; i < Scan::allScans.size() - 1; i++) {
		Scan::allScans[i]->clear("xyz");
		Scan::allScans[i]->clear("xyz reduced");
	}

	// Build KDtree with RIEGL now in world frame
	DataXYZ riegl_xyz_world = rieglScan->get("xyz reduced");
	KDtree *rieglTree = new KDtree(PointerArray<double>(riegl_xyz_world).get(), riegl_xyz_world.size());

	std::cout << "RIEGL scan has " << riegl_xyz_world.size() << " points (now in world frame)" << std::endl;
	std::cout << "Processing remaining Pandar scans..." << std::endl;

	boost::filesystem::path output_path(outfile);
	boost::filesystem::path output_dir = output_path.parent_path();
	if (!output_dir.empty() && !boost::filesystem::exists(output_dir)) {
		boost::filesystem::create_directories(output_dir);
		std::cout << "Created new output directory: " << output_dir << std::endl;
	}

	std::ofstream out(outfile.c_str());
	if (!out.good()) {
		std::cerr << "Error: Cannot open output file: " << outfile << std::endl;
		return 1;
	}

	out << "scan_id,timestamp,pos_x,pos_y,pos_z,GT_nx,GT_ny,GT_nz,num_neighbors" << std::endl;

	double direction[3] = {0.0, -1.0, 0.0}; // Search downward along -Y in world frame
	double max_dist2 = max_dist * max_dist;
	int valid_normals = 0;

	// Process Pandar scans 1 to end
	for (size_t i = 1; i < Scan::allScans.size() - 1; i++) {
		Scan *pandarScan = Scan::allScans[i];
		std::cout << "\nProcessing scan " << pandarScan->getIdentifier() << " (" << i << "/"
			  << (Scan::allScans.size() - 2) << ")" << std::endl;

		// Get ICP-matched pose
		const double *icpTransMat = pandarScan->get_transMat();
		double original_pose[3] = {icpTransMat[12], icpTransMat[13], icpTransMat[14]};
		std::cout << "  Pandar pose (after ICP): (" << original_pose[0] << ", " << original_pose[1] << ", "
			  << original_pose[2] << ")" << std::endl;

		// Use transformed Pandar pose to query ground point in RIEGL scan
		double query_point[3] = {original_pose[0], original_pose[1], original_pose[2]};

		// Find ground point in RIEGL scan
		double *closest = rieglTree->FindClosestAlongDir(query_point, direction, max_dist2, 0);

		if (closest == nullptr) {
			std::cout << "Warning: No ground point found" << std::endl;
			std::cout << "  Original pose: (" << original_pose[0] << ", " << original_pose[1] << ", "
				  << original_pose[2] << ")" << std::endl;
			std::cout << "  ICP-corrected pose: (" << query_point[0] << ", " << query_point[1] << ", "
				  << query_point[2] << ")" << std::endl;
			continue;
		}

		double sqRad2 = knn_radius * knn_radius;

		std::vector<Point> knn_points = rieglTree->kNearestRangeSearch(closest, knn, sqRad2, 0);

	std:
		cout << "  Found " << knn_points.size() << " neighbors for normal calculation" << std::endl;

		// if less than 20 neighbors found, no reliable normal is calculated so skip
		if (knn_points.size() < 20) {
			std::cout << "Warning: Not enough neighbors found" << std::endl;
			continue;
		}

		double normal[3];
		computeNormalFromKNN(knn_points, normal);

		// Ensure normal points upward in world left-handed frame (Y-up)
		if (normal[1] < 0) {
			normal[0] = -normal[0];
			normal[1] = -normal[1];
			normal[2] = -normal[2];
		}

		// Check tilt angle of normal vector. If too steep normal does not belong to ground
		const double max_tilt_deg = 45.0;
		double tilt_deg = angleFromVertical(normal);
		if (tilt_deg > max_tilt_deg) {
			std::cout << "  Warning: rejected normal at scan " << pandarScan->getIdentifier() << " - tilt "
				  << tilt_deg << " deg exceeds " << max_tilt_deg
				  << " deg threshold (nx,ny,nz = " << normal[0] << ", " << normal[1] << ", "
				  << normal[2] << ")" << std::endl;
			continue; // skip this scan's row entirely
		}

		// Transform normal from 3DTK left-handed to ROS right-handed coordinate system
		// 3DTK left-handed (world): x, y, z
		// ROS right-handed: X=z, Y=-x, Z=y
		double normal_world_rh[3];
		normal_world_rh[0] = normal[2];	 // X_rh = z_lh
		normal_world_rh[1] = -normal[0]; // Y_rh = -x_lh
		normal_world_rh[2] = normal[1];	 // Z_rh = y_lh

		// Read timestamp from .pose file
		double timestamp = 0.0;
		std::stringstream pose_filename;
		pose_filename << dir << "scan" << std::setfill('0') << std::setw(3) << pandarScan->getIdentifier()
			      << ".pose";
		std::ifstream pose_file(pose_filename.str().c_str());
		if (pose_file.good()) {
			std::string line;
			// Line 1: xyz
			std::getline(pose_file, line);
			// Line 2: rpy
			std::getline(pose_file, line);
			// Line 3: timestamp
			if (std::getline(pose_file, line)) {
				std::istringstream iss(line);
				iss >> timestamp;
			}
			pose_file.close();
		}

		if (timestamp == 0.0) {
			std::cout << "Warning: Could not read timestamp, will be 0.00" << std::endl;
		}

		out << std::fixed << std::setprecision(6) << pandarScan->getIdentifier() << "," << std::setprecision(9)
		    << timestamp << "," << original_pose[0] << "," << original_pose[1] << "," << original_pose[2] << ","
		    << normal_world_rh[0] << "," << normal_world_rh[1] << "," << normal_world_rh[2] << ","
		    << knn_points.size() << std::endl;

		valid_normals++;
	}

	out.close();

	std::cout << "\nProcessed " << (Scan::allScans.size() - 1) << " Pandar scans" << std::endl;
	std::cout << "Found " << valid_normals << " valid ground truth normals" << std::endl;
	std::cout << "Normals written to: " << outfile << std::endl;

	delete rieglTree;

	Scan::closeDirectory();

	return 0;
}