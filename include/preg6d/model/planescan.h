/** @file
 *  @brief this file contains the PlaneScan Class that represents a Scan
 *  and holds the information about Point-to-Plane correspondences.
 *  The class has many utilities for labeling the points based on various models.
 *  It can hold matches and correspondences between points and planes.
 *
 *  The source files for this class are split.
 *  You find the basic class-related stuff, as well as I/O and some other
 *  utilities such as program options, in "planescan.cc".
 *  Everything that has to do with clustering is defined in "clustering.cc"
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 */

#ifndef __PLANESCAN_H_
#define __PLANESCAN_H_

#include <deque>
#include <map>
#include <chrono>

#include "io/ioplanes.h"

#include "model/normalplane.h"
#include "model/detectplanes.h"

#include "match/matcher.h"
#include "match/simplematchers.h"
#include "match/clustermatcher.h"
#include "match/planematcher.h"

#include "tree/kdplanes.h"
#include "tree/bkdplanes.h"

#include "ANN/ANN.h"
#include "slam6d/kdIndexed.h"
#include "slam6d/scan.h"
#include "slam6d/globals.icc"
#include "slam6d/bkd.h"

/**
 * The difference between Cluster and PointCluster is the following:
 * Cluster stores indices of points, thus is very efficient.
 * It's the internal way of how clustering is processed.
 * PointCluster stores points directly, can be used for I/O between classes.
 *
 * The difference between correspondence and match is the following:
 * Correspondence is established between a Point and a global plane.
 * Match is established between a PointCluster and a global plane.
 */
typedef std::vector<size_t> Cluster;
typedef std::vector<double*> PointCluster;

typedef std::vector< std::vector<size_t> > Clusters;
typedef std::vector< std::vector<double*> > PointClusters;

typedef std::pair<double*, NormalPlane*> Correspondence;
typedef std::vector<Correspondence> Correspondences;

typedef std::pair<double*, double*> PointPair;
typedef std::vector<PointPair> PointPairs;

typedef std::pair<NormalPlane*, PointCluster* > Match;

typedef std::vector<BkdTree*> Forest;
typedef std::vector<double*> Cache;

/*
 * Some inline util functions
 */

double inline distSqr(double *p1, double *p2) {
    return sqr(p1[0]-p2[0]) + sqr(p1[1]-p2[1]) + sqr(p1[2]-p2[2]);
}

double inline distEuklid(double *p1, double *p2) {
    return sqrt( distSqr(p1, p2) );
}

// Returns the smallest angle between n1 and n2 in range [0, pi/2]. Looks at all 4 quadrants. Returns in rad
double inline angleBetweenNormals(double* n1, double* n2) {
    double dot = Dot(n1,n2); // normalized vectors may not be perfect
    dot = min( max(dot, -1.0), 1.0); // clamp to prevent NaN in acos
    double alpha = acos(dot);
    // If between 270 and 360 degrees, just take the conjugate angle.
    if (alpha > M_PI + M_PI_2) alpha = M_2_PI - alpha;
    // If between 180 and 270, flip normal.
    else if (alpha > M_PI) alpha = alpha - M_PI;
    // If between 90 and 180, flip normal.
    else if (alpha > M_PI_2) alpha = M_PI - alpha;
    return fabs(alpha);
}

// PlaneScan class
class PlaneScan
{
public:

    // The class itself stores all instances and all planes in a static manner
    static std::vector<PlaneScan*> allPlaneScans;
    static Planes allPlanes;

    ~PlaneScan();
    PlaneScan() {};
    PlaneScan(Scan*);
    PlaneScan(Scan*, Planes&);

    Matcher* matcher;
    static int _match_type;

    // for each point, store on which plane(s) it lies
    // good for efficient iteration when pairs are needed
    std::vector< Correspondence > correspondences;
    std::vector< PointPair > point_pairs;

    // Used if clustering is available
    std::map< NormalPlane*, PointCluster* > local_matches; // needed for centroid gradient
    std::map< NormalPlane*, PointCluster* > global_matches; // needed for hessian gradient

    std::set< PointCluster* > global_mismatches; // needed for new planes

    // Scan data (indexed)
    double **points; // raw point data
    double **normals; // raw normal data
    size_t nrpts;
    double timestamp;
    double orig_transMat[16];
    // Pose Data
    double rPos[3];
    double rPosTheta[3];
    double transMat[16];
    // Disk data
    std::string basedir;
    std::string identifier;

    // We use Quadratic Interpolation Maps for adapting clustering parameters
    static InterpolMap adaptGrowth; // Adaptor for growth rate
    static InterpolMap adaptNrPts; // Adaptor for number of points

    // Stores planar clusters (stores indices of pts)
    Clusters clusters;

    // Get point
    const double* operator()(int);

    /**
     * @brief Converts the Clusters (which store indices of the points) to
     * globaly consistent PointClusters (which store the points in global ref. frame).
     * @return A vector of vector of Points, representing global Point-Clusters.
     */
    vector<vector<double*> > getGlobalClusters();
    vector<vector<double*> > getLocalClusters(); // same but local

    vector<PointPlane*> getGlobalPointPlanes();
    vector<PointPlane*> getLocalPointPlanes();

    /**
     * Label points to planes.
     * This is initially done by the constructor, so you dont have to explicitly call the function
     * unless you want to update correspondences after transforming the scan.
     */
    void labelPoints(int mtype, bool quiet);

    bool isOkay();
    bool isEmpty();
    void reset();
    size_t index();

    // Static functions

    // writes all the frame buffers to .frames files
    static bool writeFrames();
    static double totalError();
    static double dist2Plane(const double*, const double*, const NormalPlane*);
    static double dist2Plane(const double*, const NormalPlane*);
    static double projDist2Plane(const double*, const double*, NormalPlane*);
    static double projDist2Plane(const double*, NormalPlane*);

    // Program option setters (also static)

    static void setUseCorrespondenceMin(bool);
    static void setPlanes(Planes&);
    static void setEpsDist(double);
    static void setEpsPPD(double);
    static void setReduce(bool);
    static void setUseClustering(bool);
    static void setEpsSimilarity(double);
    static void setMaxEigRatio(double);
    static void setClusterOutputPath(string);
    static void setExportClusterOutputRGB(bool);
    static void setReadNormals(bool);
    static void setReadClusters(bool);
    static void setUseNormalCor(bool);
    static void setMinScansize(size_t);
    static void setPlaneAlgo(plane_alg);

    // Region growing related program options:
    static void setPercentileFilter(double);
    static void setKNearest(int);
    static void setGrowthThresh(double);
    static void setGrowthMaxAdapt(double);
    static void setMaxPlaneThickness(double);
    static void setMinClusterPoints(int);
    static void setMinClusterPointsFar(int);
    static void setMaxClusters(int);

    // dont use this, please. I beg you. its shit. I use it only for evaluation purpose.
    static void setUseRegionGrowingBruteForce(bool);

    // Weight setters for scoring function
    static void setWeightAlpha(double);
    static void setWeightOverlap(double);
    static void setWeightHesse(double);
    static void setWeightEigen(double);
    static void setWeightPPD(double);

    double calcErr();
    double calcErrPtPairs();
    void addFrame(Scan::AlgoType); // adds the current transMat to the .frames file
    string getFrameData();

    // Use this to quickly transform the scan.
    void transform(const double*);

    // Use this to quickly apply the intrinsic coordinates with the robot position.
    // Therefore, apply the complete pose change from the previous scan, i.e. odometry extrapolation.
    void mergeCoordinatesWithRoboterPosition(PlaneScan*, bool* dims = 0);

    //static helper functions
    static bool isInPlane(double*, double* trans, NormalPlane*);
    static bool similar(double*, double*);

    static bool use_correspondence_min;
    static bool use_clustering;
    static bool use_normal_cor;
    static double _eps_dist;

private:
    static unsigned int idx_global;

    bool empty = false;
    unsigned int idx;

    /**
     * The stringstream sout buffers the frame file. It will be written to disk at
     * once at the end of the program. This reduces file operations and saves time.
     */
    stringstream sout;
    // Cache for mean cluster normals
    vector<double*> ncm_cache;
    // Cache for cluster centroids
    vector<double*> centroid_cache;

    /**
     * Dynamic searchtree with all points in the scan.
     * Finds k nearest neighbours for region growing, then deletes those
     * points from the tree so you dont visit these points again.
     */
    KDtreeIndexed* tree;
    // Normal calculation from all scan data
    void calcNormals(DataXYZ&);

    // Instance helper functions for clustering

    void calcLocalClusters(Scan *scan);

    // Read clusters from RGB scans
    int readClusters(DataRGB&);
    // Read cluster from UOSC (type) scans
    int readClusters(DataType&);

    void writeClustersToScanfileRGB(string path);
    void writeClustersToScanfileType(string path);

    // Programm options

    static double _eps_ppd;
    static bool reduce;
    static string cluster_outpath;
    static bool export_cluster_rgb;
    static bool read_normals;
    static bool read_clusters;
    static size_t min_scansize;
    static plane_alg plane_algo;


    /* ------------------------------------------------
     *            REGION GROWING RELATED
     * the below fields are only used if rg is selected
     * ------------------------------------------------ */

    /*
     * Multiple dynamic searchtrees, each tree for each cluster.
     * We iterativley put the pts in the cluster in a Bkd tree, i.e. multiple balanced kd-trees.
     * We get closest point for point-2-cluster distance very fast this way.
     */
    vector<BkdTree*> bkd_forest;
    // Only used when there is a cluster without normal information
    double* calcNormalOnDemand(Cluster&);
    // Updates the centroid cache for a cluster c with point j
    void update_centroid_cache(int c, size_t j);
    // Get cached centroid of cluster at index c
    double * centroidClusterMeanCached(int c);
    // Updates the normal cache for a cluster c with point j
    void update_ncm_cache(int c, size_t j);
    // Get normal vector of cluster at index c (cached, requires use of the above)
    double * normalClusterMeanCached(int c);
    // Get normal vector of the param. cluster (slow, not cached)
    double * normalClusterMean(Cluster&);
    /**
     * @brief Performs the region growing step.
     * @param threadNum: Thread number identifier
     */
    void regionGrowing(int threadNum = 0);
    /**
     * @brief Performs filtering on clusters.
     * Removes small clusters, and clusters of bad quality.
     * @param hulls: The clusters in NormalPlane-ptr form.
     */
    void filter(vector<NormalPlane*> &hulls);
    /**
     * @brief Performs filtering on clusters based on histogram percentiles.
     * Only keeps the largest clusters that contain a given percentage of points in them.
     * @param hulls: The clusters in NormalPlane-ptr form.
     * @param histogram: The histogram of the clusters, in descending order, ordered by size.
     * @param percentile: Between 0 and 1. The amount of points you want to keep.
     */
    void filterPercentile(vector<NormalPlane*> &hulls, vector<size_t> &histogram, double percentile);
    /**
     * @brief BF! Calculates distance from a point to a cluster in O(N). WARN: very slow.
     * @param cluster: vector of point indices, corresponding to a cluster
     * @param p: query point
     * @return Distance from <p> to <cluster>.
     */
    double clusterPointDist(Cluster&, double*);
    /**
     * @brief Calculates distance from a point to a cluster in O(log(N))
     * @param cluster: Pointer to a BkdTree of points, corresponding to a cluster
     * @param p: query point
     * @param threadNum: Thread number index
     * @return Squared distance from <p> to <cluster>.
     */
    double clusterPointDist2(BkdTree*, double*, int threadNum = 0);
    // Calculates distance from a cluster to another cluster, using bkd trees. O(n*log(k))
    double clusterDistBkd(Cluster&, BkdTree*, int threadNum = 0);
    // Calculates distance from a cluster to another cluster, using brute force. O(n*k)
    double clusterDistBF(Cluster&, Cluster&);
    // Calculates local clusters, if <use_clustering> is set.
    void initRG();
    void detectRG();

    // Region growing program options
    static int k_neighbours;
    static double d_growth;
    static double d_growth_max_adapt;
    static double d_thickness;
    static int n_min_clusterpoints;
    static int n_min_clusterpoints_far;
    static int n_max_clusters;
    static double eps_similarity;
    static double eigratio;
    static double score_w_alpha;
    static double score_w_overlap;
    static double score_w_hesse;
    static double score_w_ppd;
    static double score_w_eigen;
    static bool use_bruteforce;
    static double percentile_filter;
};

typedef std::vector<PlaneScan*> PlaneScans;

#endif //__PLANESCAN_H_
