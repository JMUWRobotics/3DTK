/** @file
 *  @brief Point-to-plane based post registration correction algorithm.
 *  Requires globaly registered scan archive and extracted planes.
 *
 *  Takes the scans and registers them against the predefined extracted planes.
 *
 *  You can use bin/slam6d or bin/correction for registration.
 *  Use bin/planes to extract planes from a scan.
 *  To extract planes from the globaly registered scan, you have
 *  to condense the scan archive into one globaly consistent .3d scan-file.
 *  Use bin/condense for that.
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 */

#include <csignal>
#include <string>
#include <map>
#include <omp.h>
#include <boost/program_options.hpp>

#include "newmat/newmat.h"
#include "newmat/newmatio.h"

#include "slam6d/scan.h"
#include "slam6d/io_types.h"

#include "io/ioplanes.h"

#include "model/planescan.h"
#include "match/planematcher.h"

#include "opt/optimizer.h"
#include "opt/adadelta6d.h"
#include "opt/newtons6d.h"
#include "opt/svd.h"

#include "util.h"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// Boost needs to convert from string to IOType. Use this validation function.
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

Dimensions formatname_to_dimensions(const char* str)
{
    if (strcasecmp(str, "all") == 0) return ALL;
    else if (strcasecmp(str, "rolling") == 0) return ROLLING;
    else if (strcasecmp(str, "descending") == 0) return DESCENDING;
    else if (strcasecmp(str, "rotating") == 0) return ROTATING;
    else if (strcasecmp(str, "udlock") == 0) return UDLOCK;
    else if (strcasecmp(str, "jasper") == 0) return UDLOCK;
    else if (strcasecmp(str, "jrbm") == 0) return UDLOCK;
    else if (strcasecmp(str, "sliding") == 0) return SLIDING;
    else if (strcasecmp(str, "xytrans") == 0) return UDLOCK;
    else if (strcasecmp(str, "radler") == 0) return RADLER;
    else throw std::runtime_error(string("Invalid Dimension type."));
}

// Same for Dimensions enum
void validate(boost::any& v, const std::vector<std::string>& values,
              Dimensions*, int) {
  if (values.size() == 0)
    throw std::runtime_error("Invalid model specification");
  std::string arg = values.at(0);
  try {
    v = formatname_to_dimensions(arg.c_str());
  } catch (...) { // runtime_error
    throw std::runtime_error("Format " + arg + " unknown.");
  }
}

/*
 * Use boost to set and parse command line options.
 */
int parse_options(  int argc,
                    char** argv,
                    string &scandir,
                    string &planedir,
                    IOType &type,
                    bool &quiet,
                    int &maxDist,
                    int &minDist,
                    int &octree,
                    double &red,
                    bool &scanserver,
                    double &eps_dist,
                    double &eps_ppd,
                    int& start,
                    int& end,
                    int& n_iter,
                    int& min_scansize,
                    double& eps_crit,
                    double& alpha,
                    double& rPos_alpha_scale,
                    int& anim,
                    int& nthreads,
                    Dimensions& dims,
                    bool& autoAlpha,
                    int& updateCor,
                    bool& continuous,
                    bool& use_min_cor,
                    bool& use_frames,
                    int& k_nearest,
                    bool& use_clustering,
                    double& d_growth,
                    double& d_growth_max_adapt,
                    int& n_min_clusterpoints,
                    double& eps_norm_similarity,
                    string& clusters_outputpath,
                    bool& use_normal_cor,
                    bool& color,
                    double& eigratio)
{
    // Add program option descriptions and reference variables accordingly
    po::options_description generic("Generic options");
    po::options_description input("Program options");
    po::options_description hidden("Hidden options");

    // Add program options in a curried way, using the overloaded ()-operator
    generic.add_options()
        ("help,h", "Display a very helpful message");
    input.add_options()
        ("format,f", po::value<IOType>(&type)->default_value(UOS, "uos"),
            "chose input format from {uos, uosr, uos_rgb, ous_normal, old, xyz}")
        ("start,s", po::value<int>(&start)->default_value(0),
            "Skip the first <arg> scans.")
        ("end,e", po::value<int>(&end)->default_value(-1),
            "Stop at scan index <arg>")
        ("planes,p", po::value<string>(&planedir)->default_value("./planes/"),
            "Directory of the globaly consistent planes (either convex hull or normal representation)")
        ("quiet,q", po::bool_switch(&quiet)->default_value(false),
            "Supress any output (except segfault) of the program.")
        ("anim", po::value(&anim)->default_value(-2),
            "If set, store animation in .frames files. Can be visualized with bin/show")
        ("max,m", po::value<int>(&maxDist)->default_value(-1),
            "neglegt all data points with a distance larger than <arg> 'units'")
        ("min,M", po::value<int>(&minDist)->default_value(-1),
            "neglegt all data points with a distance smaller than <arg> 'units'")
        ("minsize,n", po::value<int>(&min_scansize)->default_value(0),
            "Only uses optimization for scans with more pts than <arg>.")
        ("reduce,r", po::value<double>(&red)->default_value(-1.0),
            "turns on octree based point reduction (voxel size= <arg>)")
        ("octree,O", po::value<int>(&octree)->default_value(0),
            "Use randomized octree based point reduction (pts per voxel=<arg>)")
        ("scanserver,S", po::value<bool>(&scanserver)->default_value(false),
            "Use the scanserver as an input method and handling of scan data")
        ("eps_dist,d", po::value<double>(&eps_dist)->default_value(25),
            "Maximum Hesse distance threshold for correspondence.")
        ("eps_ppd,P", po::value<double>(&eps_ppd)->default_value(50),
            "Maximum projected polygon distance threshold for correspondence.")
        ("iter,i", po::value<int>(&n_iter)->default_value(-1),
            "Number <arg> of maximum gradient descent iterations before updating correspondences.")
        ("eps_crit", po::value<double>(&eps_crit)->default_value(-1),
            "Quickstop gradient descent if convergence falls below <arg> threshold.")
        ("alpha,a", po::value<double>(&alpha)->default_value(0.1),
            "Gradient descent factor / convergence rate. Use carefully.")
        ("pos_scale_alpha,c", po::value<double>(&rPos_alpha_scale)->default_value(1),
            "Scales the <alpha> value applied to the robots position about <arg>.")
        ("nthreads,t", po::value<int>(&nthreads)->default_value(1),
            "Use <arg> threads for gradient descent. Each thread minimizes one sub-scan.")
        ("dimensions,D", po::value<Dimensions>(&dims)->default_value(ALL, "all"),
            "Specifies which dimensions should be optimized. Chose from {all, rolling, descending, rotating}.")
        ("auto", po::bool_switch(&autoAlpha)->default_value(false),
            "Auto-detect optimum alpha to initialize AdaDelta. Convergence will be guaranteed to be fast. However, use with care as sometimes smaller steps on the gradient are prefered.")
        ("update_correspondences,k", po::value<int>(&updateCor)->default_value(10),
            "Update correspondences <arg> times. Total number of iterations will therefor be i*k")
        ("sequential,S", po::bool_switch(&continuous)->default_value(false),
            "Do not parallelize iterations. Instead, apply the last iteration pose change onto the next scan.")
        ("use_min_cor", po::bool_switch(&use_min_cor)->default_value(false),
            "Use minimum distance to plane for correspondence in case of ambiguity.")
        ("use_normal_cor", po::bool_switch(&use_normal_cor)->default_value(false),
            "Use direct normals for better point-2-plane correspondences. Uses k nearest neighbours for normal calculation. Or you can read normals using -f uos_normal")
        ("use_clustering", po::bool_switch(&use_clustering)->default_value(false),
            "Use clustering for better point-2-plane correspondences. Utilizes local normal calculation using k-nearest neighbours. Or you can read cluster using -f uos_rgb")
        ("k_nearest,K", po::value(&k_nearest)->default_value(20),
            "Use <arg> nearest neighbours in local normal calculation.")
        ("d_growth,g", po::value(&d_growth)->default_value(50),
            "Threshold for cluster region growing. Stop region growing if there is no new point nearer than <arg>.")
        ("d_max,G", po::value(&d_growth_max_adapt)->default_value(-1.0),
            "If set, use <arg> as maximum adaptation for region growing parameter (-g)")
        ("min_cluster_size,z", po::value(&n_min_clusterpoints)->default_value(500),
            "Do not consider a cluster unless it has more then <arg> points in it.")
        ("eps_similarity,y", po::value(&eps_norm_similarity)->default_value(10),
            "Two normals are considered similar if their angle is smaller then <arg> in degrees.")
        ("eigenratio,R", po::value(&eigratio)->default_value(0.05),
            "Maximum eigenvalue ratio { e1/(e1+e2+e3) }. Everything above is considered a 'bad quality plane'")
        ("write_clusters,w", po::value(&clusters_outputpath)->default_value(""),
            (string("When used, the program writes the results of the local normal clustering ")+
            string("to the specified path in <arg>. You can visualize the clusters with nice colors as follows:\n")+
            string("bin/show -c -f uos_rgb /your/dir/ -s 0")).c_str())
        ("color,C", po::bool_switch(&color)->default_value(false),
            "Exports clusters with RGB values.")
        ("continue", po::bool_switch(&use_frames)->default_value(false),
            "Use pose specified in .frames file instead of .pose file.");
    hidden.add_options()
        ("input-dir", po::value<string>(&scandir), "input dir");

    // All options together
    po::options_description alloptions;
    alloptions.add(generic).add(input).add(hidden);

    // Only commandline, visible with --help
    po::options_description cmdoptions;
    cmdoptions.add(generic).add(input);

    // positional argument for input directory
    po::positional_options_description pos;
    pos.add("input-dir", 1); // max 1 pos arg

    // Map and store option inputs to variables
    po::variables_map vars;
    po::store( po::command_line_parser(argc, argv).
                options(alloptions).
                positional(pos).
                run(),
                vars);

    // help display msg
    if ( vars.count("help") )
    {
        cout << cmdoptions;
        cout << endl << "Example usage:" << endl
            << "\t bin/preg6d dat/ -f uos —planes dat/planes/ —reduce 100 —octree 10 -d 100 -i 50 -k 50 -a 0.01 -c 10 -t 20 —dimensions rolling —anim —max 1500"
            << endl;
        exit(0);
    }
    po::notify(vars);

    // Add trailing directory slash if there is none. Works differently when compiling under Windows
#ifndef _MSC_VER
    if (scandir[ scandir.length()-1 ] != '/') scandir = scandir + "/";
    if (planedir[ planedir.length()-1 ] != '/') planedir = planedir + "/";
#else
    if (scandir[ scandir.length()-1]  != '\\') scandir = scandir + "\\";
    if (planedir[ planedir.length()-1]  != '\\') planedir = planedir + "\\";
#endif

    // return with success exit code
    return 0;
}

int main(int argc, char **argv)
{
    // Parse program options first:
    string scandir;
    int startScan = 0;
    int endScan = -1;
    string planedir;
    IOType type = UOS;
    bool quiet = false;
    int maxDist = -1, minDist = -1;
    int octree = 0;
    double red = -1.0;
    bool scanserver = false;
    double eps_dist = 0;
    double eps_ppd = 0;
    int n_iter = -1;
    int min_scansize = 1;
    double eps_crit = -1;
    double alpha = 0.001;
    double rPos_alpha_scale = 100;
    int anim = 1;
    int nthreads = 1;
    Dimensions dims;
    bool autoAlpha = false;
    int updateCor = 10;
    bool continuous = false;
    bool use_min_cor = false;
    bool use_frames = false;
    int k_nearest = 20;
    bool use_clustering = false;
    double d_growth = 5.0;
    double d_growth_max_adapt = -1.0;
    int n_min_clusterpoints = 100;
    double eps_norm_similarity = 5;
    string cluster_output_path = "";
    bool use_normal_cor = false;
    bool color = false;
    double eigratio = 0.05;

    parse_options(argc, argv, scandir, planedir, type,
        quiet, maxDist, minDist, octree, red, scanserver,
        eps_dist, eps_ppd, startScan, endScan, n_iter, min_scansize, eps_crit, alpha,
        rPos_alpha_scale, anim, nthreads, dims, autoAlpha, updateCor,
        continuous, use_min_cor, use_frames, k_nearest, use_clustering,
        d_growth, d_growth_max_adapt ,n_min_clusterpoints, eps_norm_similarity, cluster_output_path,
        use_normal_cor, color, eigratio);

    Optimizer::setUseP2P(false);
    autoAlpha = false;
    AdaDelta::setAuto(false);
    // Open and setup all planes
    PlaneIO::read( planedir );
    std::cout << "Read " << PlaneIO::allPlanes.size() << " planes:" << std::endl;
    {
        Planes::iterator plane = PlaneIO::allPlanes.begin();
        Planes::iterator end = PlaneIO::allPlanes.end();
        for (; plane != end; ++plane)
        {
            std::cout << "Plane normal ["
                << (*plane)->n[0] << ", "
                << (*plane)->n[1] << ", "
                << (*plane)->n[2] << "]" << std::endl;
        }
    }

#ifdef _OPENMP
    omp_set_num_threads(OPENMP_NUM_THREADS);
#endif

    // Output of given command line arguments for optimization
    if (min_scansize > 0) {
        std::cout << "Consider a scan as 'empty' when it has less/eq. than " << min_scansize << " pts." << std::endl;
        PlaneScan::setMinScansize( min_scansize );
    }
    if (n_iter != -1) {
        std::cout << "Using maximum " << n_iter << " iterations." << std::endl;
        Optimizer::setMaxIter(n_iter);
    }
    //if (eps_crit != -1) {
        std::cout << "Quickstop if convergence < " << eps_crit << std::endl;
        Optimizer::setEpsConvergence(eps_crit);
    //}

    std::cout << "Scaling alpha applied to position about " << rPos_alpha_scale << std::endl;
    AdaDelta::setRPosAlphaScale( rPos_alpha_scale );
    Optimizer::setUpdateCor( updateCor ); // k iteration value
    if (quiet) {
        std::cout << "Quiet mode." << std::endl;
        Optimizer::setQuiet(true);
    }
    if (anim != -2) {
        std::cout << "Storing animation in .frames." << std::endl;
        Optimizer::setAnim(anim);
    }
    std::cout << "Updating correspondences " << updateCor << " times." << std::endl;
    if (use_normal_cor && use_clustering)
    {
        cout << "Both normals and clustering correspondences are activated." << endl;
        cout << "Please deactivate one of them. You cannot use both of them at the same time." << endl;
        return 1;
    }
    AdaDelta::setAlphaInit(alpha);
    PlaneScan::setUseNormalCor( use_normal_cor );
    PlaneScan::setUseClustering( use_clustering );
    PlaneScan::setKNearest( k_nearest );
    PlaneScan::setGrowthThresh( d_growth );
    PlaneScan::setMaxEigRatio( eigratio );
    Optimizer::setEpsKernel(eps_dist);
    // Growthrate adaption for clustering
    if (d_growth_max_adapt != -1) {
        double mind2 = __DBL_MAX__, maxd2 = 0.0;
        double tmp_mind2, tmp_maxd2;
        for (uint i = 0; i < Scan::allScans.size(); ++i ) {
            Scan *scan = Scan::allScans.at(i);
            DataXYZ xyz(scan->get("xyz"));
            minMaxD2(xyz, tmp_mind2, tmp_maxd2);
            mind2 = min(mind2, tmp_mind2);
            maxd2 = max(maxd2, tmp_maxd2);
        }
        PlaneScan::setGrowthMaxAdapt(d_growth_max_adapt);
        PlaneScan::adaptGrowth = {mind2, maxd2,
                                sqr(d_growth), sqr(d_growth_max_adapt),
                                QUADRATIC };
        cout << "Region growing max. adaptation: " << d_growth_max_adapt << endl;
    }
    PlaneScan::setEpsSimilarity( eps_norm_similarity );
    PlaneScan::setMinClusterPoints( n_min_clusterpoints );

    // You can read normals to save time during clustering.
    if (type == IOType::UOS_NORMAL) {
        PlaneScan::setReadNormals(true);
        cout << "Using normals." << endl;
    }
    // You can read entire clusters
    else if (type == IOType::UOS_RGB) {
        PlaneScan::setReadClusters(true);
        cout << "Reading clusters." << endl;
    }
    else if (use_clustering) {
        std::cout << "Using " << k_nearest << " nearest neighbours for normal calculation." << endl;
    }

    // Open and pre-process all scans
    if (use_frames) Scan::continueProcessing(true);
    Scan::setProcessingCommand(argc, argv);
    Scan::openDirectory(scanserver, scandir, type, startScan, endScan);

    // Converting the Scan objects into PlaneScan objects
    PlaneScan::setEpsDist( eps_dist );
    PlaneScan::setEpsPPD( eps_ppd );
    PlaneScan::setReduce( red != -1 );
    PlaneScan::setPlanes( PlaneIO::allPlanes );
    PlaneScan::setUseCorrespondenceMin( use_min_cor );
    PlaneScan::setClusterOutputPath( cluster_output_path );
    PlaneScan::setExportClusterOutputRGB( color );

    PlaneScan::setWeightAlpha(1);
    PlaneScan::setWeightOverlap(1);
    PlaneScan::setWeightHesse(0.8);
    PlaneScan::setWeightPPD(0.3);
    PlaneScan::setWeightEigen(0);

    bool use_point_2_point = false;
    Optimizer::setUseP2P( use_point_2_point );

    if ( use_normal_cor ) {
        PlaneScan::_match_type = 3;
    } else {
        PlaneScan::_match_type = use_min_cor ? 2 : 1;
    }
    //PlaneScan::_match_type = 4;
    PlaneScan *ps;
    {
        #pragma omp parallel for schedule(dynamic)
        for( unsigned int k = 0; k < Scan::allScans.size() ; ++k)
        {
            Scan *scan = Scan::allScans.at(k);
            cout << "Reading scan" << scan->getIdentifier() << endl;
            scan->setRangeFilter( maxDist, minDist );
            scan->setReductionParameter( red, octree );
            // Finds point-2-plane correspondences
            ps = new PlaneScan( scan );
            if (!continuous) ps->labelPoints(PlaneScan::_match_type, false);
        }
    }

    // Iterate all the scans. Transformations are buffered between iterations
    PlaneScan *prevScan;
    Optimizer *iter;
    if (continuous)
    {
        for ( unsigned int k = 0; k < PlaneScan::allPlaneScans.size(); ++k)
        {
            PlaneScan *pScan = PlaneScan::allPlaneScans.at(k);
            if (k != 0)
            {
                cout << "Merge scan" << pScan->identifier << " to scan" << prevScan->identifier << endl;
                bool* dims_array = Optimizer::convertDimensionsToBoolArray(dims);
                pScan->mergeCoordinatesWithRoboterPosition(prevScan, dims_array);
            }
            prevScan = pScan;

            cout << "Optimizing scan" << pScan->identifier << endl;
            iter = new PlaneSVD(pScan, dims);
            //iter.relabel();
            pScan->labelPoints(PlaneScan::_match_type, quiet);
            iter->operator()(); // do the iterations
            delete iter;
        }
    }
    else
    {
        // #pragma omp parallel for schedule(dynamic)
        for ( unsigned int k = 0; k < PlaneScan::allPlaneScans.size(); ++k)
        {
            // Finds point-2-plane correspondences
            PlaneScan *pScan = PlaneScan::allPlaneScans.at(k);
            iter = new PlaneSVD(pScan, dims);
            iter->operator()(); // do the iterations
            delete iter;
        }
    }

    // cout << "Interpolating 6D path..." << endl;
    // int length = PlaneScan::allPlaneScans.size();
    // double *transmatold = PlaneScan::allPlaneScans.at(length-1)->orig_transMat;
    // double *transmatnew = PlaneScan::allPlaneScans.at(length-1)->transMat;
    // // compute change between old and corrected pose
    // double transMatDIFF[16];
    // double transMatOldInv[16];
    // M4inv(transmatold,transMatOldInv);
    // MMult(transmatnew,transMatOldInv,transMatDIFF);
    // double qDiff[4];     // change in rotation as quternion
    // double transDiff[3]; // change in translation (linear)
    // Matrix4ToQuat(transMatDIFF,qDiff,transDiff);
    // double *transMat;
    // double qStart[4]{1,0,0,0};
    // // interpolate rotation with slerp and linear iterpolate translation
    // for (int i = 0; i < length; ++i)
    // {
    //     PlaneScan *pScan = PlaneScan::allPlaneScans.at(i);
    //     transMat = pScan->transMat;
    //     double t = (static_cast<double>(i) / static_cast<double>(length));
    //     double qDiff_i[4];
    //     double transDiff_i[3];
    //     slerp(qStart,qDiff,t,qDiff_i);
    //     for(unsigned int j=0; j<3; j++) {
    //         transDiff_i[j] = transDiff[j]*t;
    //     }
    //     double transMatDiff_i[16];
    //     QuatToMatrix4(qDiff_i,transDiff_i,transMatDiff_i);
    //     pScan->transform(transMatDiff_i);
    //     pScan->addFrame(Scan::ICPINACTIVE);
    // }

    cout << "Writing transformations to .frames..." << endl;
    // In order for the animation to work, all .frames files must have same length
    //PlaneScan::fillFrames();
    // Write .frames files to disk from the buffered transformations
    PlaneScan::writeFrames();

    return 0;
}
