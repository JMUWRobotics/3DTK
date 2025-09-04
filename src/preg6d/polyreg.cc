/**
 * @brief
 * Flexible program that performs SLAM given laserscans.
 * Based on finding and matching polygons to each other.
 * Combines point-based and plane-based correspondences.
 * You can either include planar clusters (use bin/cluster for export)
 * or let this program do the job.
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 *
 */

#include <csignal>
#include <string>
#include <omp.h>
#include <boost/program_options.hpp>

#include "newmat/newmatio.h"

#include "slam6d/scan.h"
#include "slam6d/io_types.h"

#include "model/planescan.h"
#include "model/detectplanes.h"

#include "opt/optimizer.h"
#include "opt/adadelta6d.h"
#include "opt/newtons6d.h"
#include "opt/gaussnewton.h"
#include "opt/svd.h"

#include "match/planematcher.h"

#include "io/ioplanes.h"

#include "util.h"

using namespace std;
namespace po = boost::program_options;

//  Handling Segmentation faults and CTRL-C with csignal lib
void sigSEGVhandler (int v)
{
  static bool segfault = false;
  if(!segfault) {
    segfault = true;
    cout << endl
         << "# **************************** #" << endl
         << "  Segmentation fault or Ctrl-C" << endl
         << "# **************************** #" << endl
         << endl;
    // save frames and close scans
    PlaneScan::writeFrames();
    cout << "Frames saved." << endl;
    if (PlaneScan::allPlanes.size() != 0) {
        PlaneIO::write( PlaneScan::allPlanes, PlaneScan::allPlaneScans.at(0)->basedir );
        cout << "Planes saved." << endl;
    }
    Scan::closeDirectory();
  }
  exit(-1);
}

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

void validate(boost::any& v, const std::vector<std::string>& values,
              plane_alg*, int) {
  if (values.size() == 0)
    throw std::runtime_error("Invalid model specification");
  std::string arg = values.at(0);
  try {
    v = formatname_to_plane_alg(arg.c_str());
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
    // evil floor hack
    else if (strcasecmp(str, "udlock") == 0) return UDLOCK;
    else if (strcasecmp(str, "jasper") == 0) return UDLOCK;
    else if (strcasecmp(str, "jrbm") == 0) return UDLOCK;
    else if (strcasecmp(str, "sliding") == 0) return SLIDING;
    else if (strcasecmp(str, "xytrans") == 0) return SLIDING;
    else if (strcasecmp(str, "radler") == 0) return RADLER;
    else if (strcasecmp(str, "dorit") == 0) return RADLER;
    else if (strcasecmp(str, "xytransyaw") == 0) return RADLER;
    else if (strcasecmp(str, "norot") == 0) return NOROT;
    else if (strcasecmp(str, "trans") == 0) return NOROT;
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
                    double& eps_crit,
                    double& alpha,
                    double& rPos_alpha_scale,
                    int& anim,
                    Dimensions& dims,
                    bool& autoAlpha,
                    int& updateCor,
                    bool& use_frames,
                    int& k_nearest,
                    double& d_growth,
                    double& d_growth_max_adapt,
                    int& n_max_clusters,
                    int& n_min_clusterpoints,
                    int& n_min_clusterpoints_far,
                    double& eps_norm_similarity,
                    int& eps_thickness,
                    string& clusters_outputpath,
                    int& min_scansize,
                    bool& color,
                    double& eigratio,
                    double& w_alpha,
                    double& w_overlap,
                    double& w_hesse,
                    double& w_eigen,
                    double& w_ppd,
                    bool& trustpose,
                    double& percentile,
                    int& jter,
                    int& opt_type,
                    int& match_type,
                    plane_alg& alg,
                    bool& p2p
                    )
{
    // Add program option descriptions and reference variables accordingly
    po::options_description generic("Generic options");
    po::options_description input("Program options");
    po::options_description reduction("Reduction options");
    po::options_description optimization("Optimization options");
    po::options_description clustering("Clustering options");
    po::options_description hidden("Hidden options");

    // Add program options in a curried way, using the overloaded ()-operator
    generic.add_options()
        ("help,h", "Display a very helpful message")
        ("quiet,q", po::bool_switch(&quiet)->default_value(false),
            "Supress any output (except segfault) of the program.");
    input.add_options()
        ("format,f", po::value<IOType>(&type)->default_value(UOS, "uos"),
            "chose input format from {uos, uosr, uos_rgb, ous_normal, old, xyz}")
        ("start,s", po::value<int>(&start)->default_value(0),
            "Skip the first <arg> scans.")
        ("end,e", po::value<int>(&end)->default_value(-1),
            "Stop at scan index <arg>")
        ("scanserver,S", po::value<bool>(&scanserver)->default_value(false),
            "Use the scanserver as an input method and handling of scan data")
        ("use_frames", po::bool_switch(&use_frames)->default_value(false),
            "Use pose specified in .frames file instead of .pose file.");
    reduction.add_options()
        ("max,m", po::value<int>(&maxDist)->default_value(-1),
            "neglegt all data points with a distance larger than <arg> 'units'")
        ("min,M", po::value<int>(&minDist)->default_value(-1),
            "neglegt all data points with a distance smaller than <arg> 'units'")
        ("reduce,r", po::value<double>(&red)->default_value(-1.0),
            "turns on octree based point reduction (voxel size= <arg>)")
        ("octree,O", po::value<int>(&octree)->default_value(0),
            "Use randomized octree based point reduction (pts per voxel=<arg>)");
    optimization.add_options()
        ("optimizer", po::value<int>(&opt_type)->default_value(4),
            "Sets the type of optimization. Chose from\n"
            "1 = gradient descent with adadelta\n"
            "2 = relaxed newtons method with adadelta fallback\n"
            "3 = newton-gauss\n"
            "4 = SVD"
            )
        ("matcher", po::value<int>(&match_type)->default_value(4),
            "Sets the type of point-2-plane matching. Chose from\n"
            "1 = single euklid dist correspondence\n"
            "2 = min. euklid dist correspondence\n"
            "3 = min. dist normal similarity\n"
            "4 = cluster similarity scores\n"
            "5 = min-energy convex-hull matching"
            )
        ("eps_dist,d", po::value<double>(&eps_dist)->default_value(50),
            "Maximum distance threshold for correspondence in normal direction.")
        ("eps_ppd,P", po::value<double>(&eps_ppd)->default_value(50),
            "Maximum distance threshold for correspondence, orthogonal to normal direction.")
        ("eps_similarity,y", po::value(&eps_norm_similarity)->default_value(15),
            "Two normals are considered similar if their angle is smaller then <arg> in degrees. Also used for clustering if needed.")
        ("iter,i", po::value<int>(&n_iter)->default_value(200),
            "Number <arg> of maximum gradient descent iterations before updating correspondences.")
        ("update_correspondences,k", po::value<int>(&updateCor)->default_value(25),
            "Update correspondences <arg> times to mimic ICP-like behavior. Number of iterations will be i*k")
        ("jter,j", po::value(&jter)->default_value(-1),
            "Corrects the global plane model after <arg> iterations.")
        ("eps_slam,E", po::value<double>(&eps_crit)->default_value(0.001),
            "Quickstop gradient descent if convergence falls below <arg> threshold.")
        ("alpha,a", po::value<double>(&alpha)->default_value(0.01),
            "Gradient descent factor / convergence rate. Use carefully.")
        ("pos_scale_alpha,C", po::value<double>(&rPos_alpha_scale)->default_value(10),
            "Scales the <alpha> value applied to the robots position about <arg>.")
        ("dimensions,D", po::value<Dimensions>(&dims)->default_value(ALL, "all"),
            "Specifies which dimensions should be optimized. Chose from {all, rolling, descending, rotating, udlock, xytrans, xytransyaw, norot}.")
        // ("auto", po::bool_switch(&autoAlpha)->default_value(false),
        //     "Auto-detect optimum alpha to initialize AdaDelta. Convergence will be guaranteed to be fast. However, use with care as sometimes smaller steps on the gradient are prefered.")
        ("anim", po::value(&anim)->default_value(-2),
            "If set, store animation in .frames files. Stores every <arg>th iteration. Can be visualized with bin/show.")
        ("minsize,n", po::value<int>(&min_scansize)->default_value(0),
            "Only uses optimization for scans with more pts than <arg>.")
        ("w_alpha,1", po::value<double>(&w_alpha)->default_value(0.9),
            "Score function weight for angle between planes. Scores are usually ok, you dont have to change them.")
        ("w_overlap,2", po::value<double>(&w_overlap)->default_value(1),
            "Score function weight for overlapping 3D region.")
        ("w_hesse,3", po::value<double>(&w_hesse)->default_value(0.8),
            "Score function weight for point to plane distance.")
        ("w_eigen,4", po::value<double>(&w_eigen)->default_value(0.1),
            "Score function weight for eigenvalue similarity.")
        ("w_ppd,5", po::value<double>(&w_ppd)->default_value(0.3),
            "Score function weight for projected polygon distance (PPD)")
        ("trustpose", po::bool_switch(&trustpose)->default_value(false),
            "Do not extrapolate previous pose changes onto subsequent scans.")
        ("use_point2point", po::value<bool>(&p2p)->default_value(false),
            "Uses point-to-point correspondences for translation. Set false to only use point-to-plane.");
    clustering.add_options()
        ("plane_alg", po::value<plane_alg>(&alg)->default_value(RG),
            "If region growing (rg) is selected, use above params.\n"
            "Else, the program ignores the above and uses params defined in \"bin/hough.cfg\".\n"
            "Choose from {rg, rht, sht, pht, ppht, apht, ran}")
        ("k_nearest,K", po::value(&k_nearest)->default_value(20),
            "Use <arg> nearest neighbours in local normal calculation.")
        ("d_growth,g", po::value(&d_growth)->default_value(50),
            "Threshold for cluster region growing. Stop region growing if there is no new point nearer than <arg> cm.")
        ("d_max,G", po::value(&d_growth_max_adapt)->default_value(-1.0),
            "If set, use <arg> cm as maximum adaptation for region growing parameter (-g)")
        ("min_cluster_size,z", po::value(&n_min_clusterpoints)->default_value(100),
            "Do not consider a cluster unless it has more then <arg> points in it. The value can be adapted for far away parts of the scan.")
        ("min_cluster_size_far,Z", po::value(&n_min_clusterpoints_far)->default_value(-1.0),
            "If set, use <arg> as minimum adaption for far away parts of the scan.")
        ("max_clusters,x", po::value<int>(&n_max_clusters)->default_value(-1.0),
            "Upper bound for number of clusters.")
        ("thickness,T", po::value(&eps_thickness)->default_value(25),
            "If a point is further away than <arg> from the cluster in normal direction, dont include it in the cluster.")
        ("color,c", po::bool_switch(&color)->default_value(false),
            "Exports clusters with color. Use for visualization only!")
        ("eigenratio,R", po::value(&eigratio)->default_value(0.05),
            "Maximum eigenvalue ratio { e1/(e1+e2+e3) }. Everything above is considered a 'bad quality plane'")
        ("percentile,p", po::value<double>(&percentile)->default_value(-1.0),
            "E.g. <arg> = 0.9 -> Use 90 percent of points in clusters. Filter smaller clusters.")
        ("write_clusters,w", po::value(&clusters_outputpath)->default_value(""),
            (string("When used, the program writes the results of the local normal clustering ")+
            string("to the specified path in <arg>. You can visualize the clusters with nice colors as follows:\n")+
            string("bin/show -c -f uos_rgb /your/dir/ -s 0")).c_str());
    hidden.add_options()
        ("input-dir", po::value<string>(&scandir), "input dir");

    // All options together
    po::options_description alloptions;
    alloptions.add(generic).add(input).add(reduction).add(clustering).add(optimization).add(hidden);

    // Only commandline, visible with --help
    po::options_description cmdoptions;
    cmdoptions.add(generic).add(input).add(reduction).add(clustering).add(optimization);

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
            << "bin/polyreg dat/example/clusters -f uos_rgb -r 10 -O 1 -d 50 -P 10 -y 44 --use_point2point=false"
            << endl;
        exit(0);
    }
    po::notify(vars);

    // Add trailing directory slash if there is none. Works differently when compiling under Windows
#ifndef _MSC_VER
    if (scandir[ scandir.length()-1 ] != '/') scandir = scandir + "/";
#else
    if (scandir[ scandir.length()-1]  != '\\') scandir = scandir + "\\";
#endif

    // return with success exit code
    return 0;
}

inline void setOptimizer(Optimizer*& opt_iter,
                  int opt_type,
                  PlaneScan*& pScan,
                  Dimensions dim)
{
    switch(opt_type) {
        case 1:
        opt_iter = new AdaDelta( pScan, dim );
        break;
        case 2:
        opt_iter = new Newtons( pScan, dim );
        break;
        case 3:
        opt_iter = new GaussNewton( pScan, dim);
        break;
        case 4:
        opt_iter = new PlaneSVD( pScan, dim);
        break;
        default:
        opt_iter = 0;
        break;
    }
}

int main(int argc, char **argv)
{
    signal (SIGSEGV, sigSEGVhandler);
    signal (SIGINT,  sigSEGVhandler);

    // Parse program options first:
    string scandir;
    int startScan;
    int endScan;
    IOType type;
    bool quiet;
    int maxDist, minDist;
    int octree;
    double red;
    bool scanserver;
    double eps_dist;
    double eps_ppd;
    int n_iter;
    double eps_crit;
    double alpha;
    double rPos_alpha_scale;
    int anim;
    Dimensions dims;
    bool autoAlpha;
    int updateCor;
    bool use_frames;
    int k_nearest;
    double d_growth;
    double d_growth_max_adapt;
    int n_min_clusterpoints;
    int n_min_clusterpoints_far;
    int n_max_clusters;
    double percentile;
    double eps_norm_similarity;
    int eps_thickness;
    string cluster_output_path;
    int min_scansize;
    bool color;
    double eigratio;
    double w_alpha;
    double w_overlap;
    double w_hesse;
    double w_eigen;
    double w_ppd;
    bool trustpose; // extrapolate odometry by pose?
    int jter; // update global planes after this number of iterations
    bool use_p2p; // use point-2-point?
    int opt_type;
    int match_type;
    plane_alg alg;

    parse_options(
        argc, argv, scandir, type,
        quiet, maxDist, minDist, octree, red, scanserver,
        eps_dist, eps_ppd, startScan, endScan, n_iter, eps_crit, alpha,
        rPos_alpha_scale, anim, dims, autoAlpha, updateCor,
        use_frames, k_nearest, d_growth, d_growth_max_adapt, n_max_clusters ,n_min_clusterpoints,
        n_min_clusterpoints_far, eps_norm_similarity, eps_thickness, cluster_output_path,
        min_scansize, color, eigratio, w_alpha, w_overlap, w_hesse, w_eigen, w_ppd, trustpose,
        percentile, jter, opt_type, match_type, alg, use_p2p
    );

    // Write the command to a file so you dont forget it
    ofstream call_file(scandir + "lastcall.txt");
    for (int i = 0; i < argc; ++i)
        call_file << argv[i] << " ";
    call_file.flush();
    call_file.clear();
    call_file.close();


    // Set all the programm options:
    // ========================================================
    Optimizer::setUseP2P(use_p2p);
    PlaneScan::setUseClustering( true );
    PlaneScan::setPlaneAlgo( alg );
    PlaneScan::setExportClusterOutputRGB( true );
    PlaneScan::setUseRegionGrowingBruteForce( false );
    autoAlpha=false; // TODO: That smells. Auto learnrate isnt even used…

    // Set point type for reduction
    unsigned int ptype = 0;

    // You can read normals to save time during clustering.
    if (type == IOType::UOS_NORMAL) {
        ptype = PointType::USE_NORMAL;
        PlaneScan::setReadNormals(true);
        cout << "Using normals." << endl;
    }

    // You can read entire clusters
    else if (type == IOType::UOS_RGB || type == IOType::UOSC) {
        PlaneScan::setReadClusters(true);
#ifdef _OPENMP
        omp_set_num_threads(OPENMP_NUM_THREADS);
#endif
        cout << "Reading clusters." << endl;
        if (type == IOType::UOS_RGB)
            ptype = PointType::USE_COLOR;
        else
            ptype = PointType::USE_TYPE;
    }

    // None of the above. No Normals nor Clusters.
    else if (type == IOType::UOSR) {
        ptype = PointType::USE_REFLECTANCE;
    }
    else if (type == IOType::UOS || type == IOType::XYZ) {
        ptype = PointType::USE_NONE;
    }

    // Output of given command line arguments for optimization
    // set optimization options
    if (min_scansize > 0) {
        std::cout << "Consider a scan as 'empty' when it has less/eq. than " << min_scansize << " pts." << std::endl;
        PlaneScan::setMinScansize( min_scansize );
    }
    if (!autoAlpha) {
        AdaDelta::setAuto( false );
        std::cout << "Using alpha = " << alpha << std::endl;
        AdaDelta::setAlphaInit( alpha );
    } else {
        std::cout << "Auto detecting initial alpha." << std::endl;
        AdaDelta::setAuto( true );
    }
    std::cout << "Scaling alpha applied to position about " << rPos_alpha_scale << std::endl;
    AdaDelta::setRPosAlphaScale( rPos_alpha_scale );

    if (n_iter != -1) {
        std::cout << "Using maximum " << n_iter << " iterations." << std::endl;
        Optimizer::setMaxIter(n_iter);
    }
    //if (eps_crit != -1) {
        std::cout << "Quickstop if convergence < " << eps_crit << std::endl;
        Optimizer::setEpsConvergence(eps_crit);
    //}
    Optimizer::setUpdateCor( updateCor ); // k iteration value
    if (quiet) {
        std::cout << "Quiet mode." << std::endl;
        Optimizer::setQuiet(true);
    }
    if (anim != -2) {
        std::cout << "Storing animation in .frames." << std::endl;
        std::cout << "Using every " << anim << "th step for animation." << std::endl;
        Optimizer::setAnim( anim );
    }
    std::cout << "Updating correspondences " << updateCor << " times." << std::endl;
    Optimizer::setEpsKernel( eps_dist );
    PlaneScan::setEpsDist( eps_dist );
    PlaneScan::setEpsPPD( eps_ppd );
    PlaneScan::setReduce( red != -1 );

    // set clustering options
    PlaneScan::setKNearest( k_nearest );
    PlaneScan::setGrowthThresh( d_growth );
    PlaneScan::setMaxPlaneThickness( eps_thickness );
    PlaneScan::setEpsSimilarity( eps_norm_similarity );
    PlaneScan::setMaxClusters( n_max_clusters );
    PlaneScan::setMinClusterPoints( n_min_clusterpoints );
    PlaneScan::setMinClusterPointsFar( n_min_clusterpoints_far );
    PlaneScan::setPercentileFilter( percentile );
    PlaneScan::setExportClusterOutputRGB( color );
    PlaneScan::setWeightAlpha( w_alpha );
    PlaneScan::setWeightOverlap( w_overlap );
    PlaneScan::setWeightHesse( w_hesse );
    PlaneScan::setWeightEigen( w_eigen );
    PlaneScan::setWeightPPD( w_ppd );
    PlaneScan::setMaxEigRatio( eigratio );

    unsigned long time_start = GetCurrentTimeInMilliSec();
    // Growthrate adaption for clustering
    if (use_frames) Scan::continueProcessing(true);
    Scan::setProcessingCommand(argc, argv);

    double mind2 = __DBL_MAX__, maxd2 = 0.0;
    // Open directory
    Scan::openDirectory(scanserver, scandir, type, startScan, endScan);
    for (uint i = 0; i < Scan::allScans.size(); ++i ) {

        Scan *scan = Scan::allScans.at(i);
        scan->setRangeFilter( maxDist, minDist );
        scan->setReductionParameter( red, octree, PointType(ptype) );

        // Find the most near and most distant point in all scans.
        string red_string = (red != -1) ? "xyz reduced show" : "xyz";

        // If adaptive parameters are set, find minimum and maximum points in all scans.
        if (d_growth_max_adapt != -1 || n_min_clusterpoints_far != -1) {
            DataXYZ xyz(scan->get(red_string));
            double tmp_mind2, tmp_maxd2;
            minMaxD2(xyz, tmp_mind2, tmp_maxd2);
            mind2 = min(mind2, tmp_mind2);
            maxd2 = max(maxd2, tmp_maxd2);
        }
    }

    if (d_growth_max_adapt != -1 || n_min_clusterpoints_far != -1) {
        cout << "Dist of nearest point: " << sqrt(mind2) << endl;
        cout << "Dist of furthest point: " << sqrt(maxd2) << endl;
    }

    // Adapt growth rate
    PlaneScan::setGrowthMaxAdapt(d_growth_max_adapt);
    if (d_growth_max_adapt != -1) {
        PlaneScan::adaptGrowth = {
            mind2, maxd2, // map these squared distances ...
            sqr(d_growth), sqr(d_growth_max_adapt), // ... to these growth radii ...
            LINEAR // ... in a linear fashion.
        };
        cout << "Region growing max. adaptation: " << d_growth_max_adapt << endl;
    }

    // Adapt minimum number of cluster points theshold
    PlaneScan::setMinClusterPointsFar(n_min_clusterpoints_far);
    if (n_min_clusterpoints_far != -1) {
        PlaneScan::adaptNrPts = {
            mind2, maxd2, // map these squared distances ...
            (double)n_min_clusterpoints, (double)n_min_clusterpoints_far, // ... to these thresholds ...
            LINEAR // ... in a linear fashion.
        };
        cout << "Minimum cluster-size threshold for far points: " << n_min_clusterpoints_far << endl;
    }

    // Setup cluster output directory if necessary
    if ( cluster_output_path != "" )
        PlaneScan::setClusterOutputPath( cluster_output_path );
    else if ( !(type == IOType::UOS_RGB || type == IOType::UOSC) )
    {
        std::string save_dir = scandir + "clusters/";
        if ( !existsDir( save_dir.c_str() ) )
        {
            boost::filesystem::create_directory(save_dir);
            cout << "Creating \"" << save_dir << "\"." << endl;
        }
        else cout << save_dir << " exists allready." << endl;
        PlaneScan::setExportClusterOutputRGB(true);
        PlaneScan::setClusterOutputPath( save_dir );
    }

    // Process Scan directory and apply polyreg:
    // ========================================================
    cout << "Building Planescan objects... ";
    cout.flush();
    PlaneScan *pScan;
    for (uint i = 0; i < Scan::allScans.size(); i++)
    {
        cout << i << " ";
        cout.flush();
        Scan *scan = Scan::allScans.at(i);
        pScan = new PlaneScan( scan );
    }
    cout << "done." << endl;

    Optimizer *opt_iter;
    PlaneScan *prevScan = PlaneScan::allPlaneScans.at(0);
    for (uint i = 0; i < PlaneScan::allPlaneScans.size(); i++)
    {
        pScan = PlaneScan::allPlaneScans[i];
        if (i != 0) {
            if (!trustpose) {
                cout << "Merge scan" << pScan->identifier << " to scan" << prevScan->identifier << endl;
                bool* dimensions = Optimizer::convertDimensionsToBoolArray( dims );
                pScan->mergeCoordinatesWithRoboterPosition(prevScan, dimensions);
            }
        } else {
            pScan->labelPoints(match_type, true);
            PointClusters clusters = pScan->getGlobalClusters();
            for (uint j = 0; j < clusters.size(); ++j) {
                NormalPlane *plane = new NormalPlane( clusters[j] );
                if ( !eigenValueOK(plane->eigen, eigratio) )
                {
                    delete plane;
                    break;
                } else
                    PlaneScan::allPlanes.push_back( plane );
            }

            prevScan = pScan;
            pScan->labelPoints(match_type, quiet);
            // TODO :: pScan->labelPoints(5, quiet); ?
            setOptimizer( opt_iter, opt_type, pScan, dims );
            opt_iter->updateScan();
            // opt_iter->operator()(); // function operator 'cause it looks cool, starts iteration
            // delete opt_iter;

            continue;
        }

        if (pScan->nrpts < (size_t)min_scansize) continue;

        pScan->labelPoints(match_type, quiet); //quiet = false
        cout << "Scan" << pScan->identifier << " has " << pScan->global_matches.size() << " matches"
             << " and " << pScan->correspondences.size() << " correspondences." << endl;

        prevScan = pScan;

        setOptimizer( opt_iter, opt_type, pScan, dims );
        opt_iter->operator()(); // function operator 'cause it looks cool, starts iteration
        //delete opt_iter;
        cout << "Global model update..." << endl;

        // Merge corresponding Cluster-2-Plane matches
        cout << "Merge planes..." << endl;
        for (const auto& match : pScan->global_matches)
        {
            NormalPlane *plane = match.first;
            PointCluster cluster = *(match.second);
            NormalPlane *new_plane = new NormalPlane( cluster );
            plane->mergeWith( new_plane ); // auto deletes local plane after merging
        }

        // Add new planes (Clusters that did not get matched to planes)
        for (const auto& mismatch : pScan->global_mismatches)
        {
            PointCluster cluster = *mismatch;
            NormalPlane* new_plane = new NormalPlane( cluster );
            if ( eigenValueOK( new_plane->eigen, eigratio ) )
            {
                PlaneScan::allPlanes.push_back( new_plane );
            }
        }

    } // end outer for i

    unsigned long time_end = GetCurrentTimeInMilliSec();
    cout << "ELAPSED TIME: " << (time_end - time_start) / 1000.0 << "s." << endl;

    // Write .frames files to disk from the buffered transformations
    cout << "Writing transformations to .frames..." << endl;
    PlaneScan::writeFrames();

    cout << "Writing global planes to hulls/planes.list ..." << endl;
    PlaneIO::write( PlaneScan::allPlanes, PlaneScan::allPlaneScans.at(0)->basedir );

    return 0;
}