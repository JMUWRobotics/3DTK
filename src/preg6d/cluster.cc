/**
 * @file
 * @brief Program to cluter 3d points based on their planarity.
 * Exports clusters to be read later by preg6d.
 *
 * @author Fabian Arzberger, JMU, Germany
 *
 * Released under GPL version 3.
 */

#include "cluster.h"

int main(int argc, char **argv)
{
    string scandir;
    IOType type;
    bool quiet;
    int maxDist;
    int minDist;
    int octree;
    double red;
    int start;
    int end;
    int k_nearest;
    int n_threads;
    double d_growth;
    double d_growth_max_adapt;
    int n_max_clusters;
    int n_min_clusterpoints;
    int n_min_clusterpoints_far;
    double eps_alpha_similarity;
    int eps_thickness;
    bool scanserver;
    bool color;
    double eigratio;
    bool bruteforce;
    double percentile;
    plane_alg alg;

    parse_options(argc, argv, scanserver, scandir, type, quiet, maxDist,
                  minDist, octree, red, start, end, k_nearest, n_threads,
                  d_growth, d_growth_max_adapt, n_max_clusters ,n_min_clusterpoints, n_min_clusterpoints_far,
                  eps_alpha_similarity, eps_thickness, color, eigratio, bruteforce, percentile, alg);

    // We use a PlaneScan object to do clustering

    unsigned int ptype = 0;
    if (type == IOType::UOS_NORMAL || alg != RG) {
        PlaneScan::setReadNormals(true);
        ptype = PointType::USE_NORMAL;
#ifdef _OPENMP
        omp_set_num_threads(n_threads);
#endif
        cout << "Using normals." << endl;
        cout << "Using " << n_threads << " threads." << endl;
    }
    else {
        std::cout << "Using " << k_nearest << " nearest neighbours for normal calculation." << endl;
        std::cout << "Multithreading requires precalculated normals. Use bin/calc_normals to export into uos_normal format." << endl;
    }
    // None of the above
    if (type == IOType::UOSR) {
        ptype = PointType::USE_REFLECTANCE;
    }
    else if (type == IOType::UOS || type == IOType::XYZ) {
        ptype = PointType::USE_NONE;
    }
    // Set static POs
    // This is a cluster program. So that one should be obvious.
    PlaneScan::setReadClusters(false);
    PlaneScan::setUseClustering(true);
    PlaneScan::setPlaneAlgo( alg );
    PlaneScan::setKNearest(k_nearest);
    PlaneScan::setReduce( red != -1);
    PlaneScan::setGrowthThresh(d_growth);
    PlaneScan::setEpsSimilarity(eps_alpha_similarity);
    PlaneScan::setMaxClusters(n_max_clusters);
    PlaneScan::setMinClusterPoints(n_min_clusterpoints);
    PlaneScan::setUseRegionGrowingBruteForce( bruteforce );
    PlaneScan::setMaxEigRatio( eigratio );
    PlaneScan::setMaxPlaneThickness( eps_thickness );
    PlaneScan::setPercentileFilter( percentile );

    double mind2 = __DBL_MAX__, maxd2 = 0.0;
    // Open directory
    cout << "Reading Scans from disk:" << endl;
    Scan::openDirectory(scanserver, scandir, type, start, end);
    for (uint i = 0; i < Scan::allScans.size(); ++i ) {

        Scan *scan = Scan::allScans.at(i);
        cout << scan->getIdentifier() << " ";
        cout.flush();
        scan->setRangeFilter( maxDist, minDist );
        scan->setReductionParameter( red, octree, PointType(ptype) );
        //scan->calcReducedPoints();

        // Find the most near and most distant point in all scans.
        //TODO: Check if we should use "xyz show" here.
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
    cout << endl;

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
            QUADRATIC // ... in a quadratic fashion.
        };
        cout << "Region growing max. adaptation: " << d_growth_max_adapt << endl;
    }

    // Adapt minimum number of cluster points theshold
    PlaneScan::setMinClusterPointsFar(n_min_clusterpoints_far);
    if (n_min_clusterpoints_far != -1) {
        PlaneScan::adaptNrPts = {
            mind2, maxd2,
            (double)n_min_clusterpoints, (double)n_min_clusterpoints_far,
            QUADRATIC
        };
        cout << "Minimum cluster-size threshold for far points: " << n_min_clusterpoints_far << endl;
    }

    //setup save directory
    std::string save_dir = scandir + "clusters/";
    if ( !existsDir( save_dir.c_str() ) )
    {
        boost::filesystem::create_directory(save_dir);
        cout << "Creating \"" << save_dir << "\"." << endl;
    }
    else cout << save_dir << " exists allready." << endl;
    PlaneScan::setClusterOutputPath( save_dir );
    if (color) {
        cout << "WARN: Exporting clusters as RGB. Visualize them with 'bin/show -f uos_rgb -c'" << endl;
    }
    PlaneScan::setExportClusterOutputRGB( color );

    // Actually setting up the planescan objects now
    PlaneScan *ps;
    {
        #pragma omp parallel for num_threads(n_threads)
        for( unsigned int k = 0; k < Scan::allScans.size() ; ++k)
        {
            Scan *scan = Scan::allScans.at(k);
            // Finds point-2-plane correspondences
            ps = new PlaneScan( scan );
        }
    }
}