/**
 * subgraphicp implementation
 *
 * The program puts S (--size) subsequent scans into a Graph.
 * First, GraphSLAM is performed on these so called "Sub-Graphs" of size S.
 * Then, another Graph is constructed where each Sub-Graph is connected to other Sub-Graphs.
 * GraphSLAM or plain ICP is then performed on the subsequent Sub-Graphs.
 *
 * @file
 * @brief An experimental program which performs scan registration.
 *
 * Motivation:
 * This program plays with the inheritance of the MetaScan class and its base Scan class.
 * The program is not intended to compete with srr (SemiRigidRegistration / correction).
 * While this program potentially runs faster, srr is a much more powerfull tool.
 * My personal use-case for this program is as a fast and robust pre-registration step to srr.
 * The Scan class defines the
 *
 *  virtual DataPointer get(const std::string &identifier)
 *
 * and
 *
 *  virtual DataPointer create(const std::string& identifier, size_t size)
 *
 * functions which MetaScan has to implement.
 * I've made an attempt in include/slam6d/metaScan.h
 * There might be better ways to do that.
 * However, it allows us to write programs like this one,
 * without breaking MetaScan's previous functionality.
 *
 * Example usage of extended MetaScan, essentially you can now do:
 *
 *  // Create a vector of some scans (previously possible):
 *
 *  ScanVector* someScans = new ScanVector();
 *  someScans->push_back( Scan::allScans[0] );
 *  someScans->push_back( Scan::allScans[1] );
 *  someScans->push_back( Scan::allScans[2] ); // push any number of Scan*
 *
 *  // Create a vector of some MetaScans (previously possible but useless):
 *
 *  ScanVector* someMetaScans = new ScanVector();
 *  someMetaScans->push_back( new MetaScan(someScans) );
 *  someMetaScans->push_back( new MetaScan(someMoreScans) ); // push any number of MetaScan*
 *
 *  // Perform GraphSLAM or ICP on MetaScans (previously impossible):
 *
 *  // [...] (--> code to initialize graphSlam6D* gSlam and Graph* gr, see e.g. this file)
 *  gSlam->doGraphSlam6D(*gr, *someMetaScans, iter);
 *
 * @author Fabian Arzberger, Institute of CS, University of Wuerzburg, Germany
 *
 */

// System libs
#include <iostream>
#include <signal.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

// 3dtk internal
#include "slam6d/Boctree.h"
#include "slam6d/globals.icc"
#include "slam6d/icp6D.h"
#include "slam6d/icp6Dquat.h"
#include "slam6d/lum6DeulerS.h"
#include "slam6d/metaScan.h"
#include "slam6d/point_type.h"
#include "slam6d/subgraph_options.h"

using std::cout;
using std::endl;

// Writes one transformation (mat) to the file
inline void writeTransformToFrameFile(std::ofstream &file, const double *mat, unsigned int type)
{
	for (int i = 0; i < 16; ++i) {
		file << mat[i] << " ";
	}
	file << type << "\n" << std::flush;
}

// saves only the final transformation but NOT an animation
inline void saveLastFrame(ScanVector &scans, bool append)
{
	// Save frames
	for (size_t i = 0; i < scans.size(); ++i) {
		Scan *scan = scans.at(i);
		const double *last_matrix;
		Scan::AlgoType last_atype;
		scan->getFrame(scan->getFrameCount() - 1, last_matrix, last_atype);
		std::string filename = scan->getPath() + "scan" + scan->getIdentifier() + ".frames";
		std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;
		std::ofstream file(filename.c_str(), mode);
		writeTransformToFrameFile(file, last_matrix, last_atype);
		file.close();
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

int main(int argc, char **argv)
{
	// Signal handling (CTRL+C or segfault)
	signal(SIGSEGV, sigSEGVhandler);
	signal(SIGINT, sigSEGVhandler);

	// Parse program options
	std::string scandir;
	int start, end, maxdist, mindist, octree, iter, iter_lum, clpairs, maxmeta;
	size_t size;
	IOType type;
	bool quiet, scanserver, use_frames, meta, icp_only;
	double red, mdm, mdml, epsICP, epsSLAM;
	parse_options(argc, argv, scandir, start, end, type, quiet, maxdist, mindist, octree, red, scanserver,
		      use_frames, mdm, mdml, iter, iter_lum, epsICP, epsSLAM, size, clpairs, meta, maxmeta, icp_only);

	// Read the scan directory
	if (use_frames)
		Scan::continueProcessing();

	Scan::setProcessingCommand(argc, argv);
	Scan::openDirectory(scanserver, scandir, type, start, end);

	// Create the subgraphs
	vector<ScanVector *> subgraphs;
	ScanVector *subgraph = new ScanVector();
	subgraphs.reserve(ceil(Scan::allScans.size() / size));
	for (size_t i = 0; i < Scan::allScans.size(); ++i) {
		Scan *scan = Scan::allScans[i];
		scan->setRangeFilter(maxdist, mindist);
		scan->setReductionParameter(red, octree);
		scan->setSearchTreeParameter(nns_type::BOCTree);
		scan->calcReducedPoints();
		subgraph->push_back(scan);
		if (subgraph->size() >= size) {
			cout << "Reduced points for segment " << subgraph->front()->getIdentifier() << "-"
			     << subgraph->back()->getIdentifier() << endl;
			subgraphs.push_back(subgraph);
			subgraph = new ScanVector();
		}
		double id[16];
		M4identity(id);
		scan->transform(id, Scan::ICP, 1);
		scan->transform(id, Scan::ICP, 1);
	}
	saveLastFrame(Scan::allScans, false);

	// Setup ICP minimizer and graph SLAM objects
	icp6Dminimizer *icp_minimizer_functor = new icp6D_QUAT();
	graphSlam6D *gSlam = new lum6DEuler(icp_minimizer_functor, mdml, mdml, iter_lum, true, false, 1, true, -1,
					    epsICP, BOCTree, epsSLAM);
	gSlam->set_quiet(true);

	// Optimize subgraphs
	for (size_t i = 0; i < subgraphs.size(); ++i) {
		ScanVector *sg = subgraphs[i];
		cout << "Optimizing graph between " << sg->front()->getIdentifier() << "-"
		     << sg->back()->getIdentifier() << endl;
		Graph *gr = gSlam->computeGraph6Dautomatic(*sg, clpairs);
		gSlam->doGraphSlam6D(*gr, *sg, iter_lum);
	}
	saveLastFrame(Scan::allScans, true);

	// Create bundled scan packages , i.e. one Scan object per subgraph
	cout << "Creating metascans... ";
	ScanVector *metascans = new ScanVector();
	metascans->reserve(subgraphs.size());
	for (size_t i = 0; i < subgraphs.size(); ++i) {
		// Get each subgraph
		ScanVector *subgraph = subgraphs[i];
		MetaScan *scan = new MetaScan(*subgraph);
		scan->setRangeFilter(maxdist, mindist);
		scan->setReductionParameter(red, octree);
		scan->setSearchTreeParameter(nns_type::BOCTree, 20);
		scan->calcReducedPoints();
		scan->createSearchTree();
		metascans->push_back(scan);
	}

	if (icp_only) {
		// Do ICP on the metascans
		icp6D *icp = new icp6D(icp_minimizer_functor, mdm, iter, false, meta, 1, true, 1, epsICP,
				       nns_type::BOCTree, false, false, maxmeta);
		icp->doICP(*metascans, PairingMode::CLOSEST_POINT);
	} else {
		// Do GraphSLAM on the metascans
		delete icp_minimizer_functor;
		icp_minimizer_functor = new icp6D_QUAT();
		gSlam = new lum6DEuler(icp_minimizer_functor, mdm, mdm, iter, true, false, 1, true, -1, epsICP, BOCTree,
				       epsSLAM);
		gSlam->set_quiet(true);
		Graph *gr = gSlam->computeGraph6Dautomatic(*metascans, clpairs);
		gSlam->doGraphSlam6D(*gr, *metascans, iter);
	}
	// Metascans store the transformation on the metascan objects,
	// we have to apply them manually on the individual scans
	double id[16];
	M4identity(id);
	for (size_t i = 0; i < metascans->size(); ++i) {
		cout << "Writing result for " << metascans->at(i)->getIdentifier() << endl;
		metascans->at(i)->transform(id, Scan::ICP, 0); // write end pose
	}

	// Save frames (TODO: export animation as well)
	saveLastFrame(Scan::allScans, true);
}