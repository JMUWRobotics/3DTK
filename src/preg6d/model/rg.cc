#include "shapes/convexplane.h"
#include "model/planescan.h"
#include "util.h"

typedef vector<double*>* HullPtr;

int PlaneScan::readClusters(DataRGB &scan_color)
{
    size_t nrclusterpts = scan_color.size();
    if (nrclusterpts == 0) return -1;

    // Init colors with size 1 for the first point
    vector<rgb> colors(0);
    rgb color = {
        (double)scan_color[0][0],   // I know. These are some
        (double)scan_color[0][1],   // dirty "unsigned char" to "double" conversions.
        (double)scan_color[0][2]    // Your compiler may complain about that.
    };
    colors.push_back(color);
    // cout << color.r << ", " << color.g << ", " << color.b << endl;
    // Init cluster with size 1, corresponding to color of first point
    clusters = Clusters(1);
    clusters[0].push_back(0);
    for (size_t i = 1; i < nrclusterpts; ++i)
    {
        bool rgb_equals = false;
        rgb color1 = {
            (double)scan_color[i][0],   // However, i ensure you that
            (double)scan_color[i][1],   // they will always work in that usecase,
            (double)scan_color[i][2]    // so you don't actually have to worry about it.
        };
        for (size_t j = 0; j < colors.size(); ++j)
        {
            rgb color2 = colors[j];
            rgb_equals = rgbEquals(color1, color2);
            if ( rgb_equals ) {
                clusters[j].push_back(i);
                break;
            }
        }
        if ( !rgb_equals ) {
            //cout << color1.r << ", " << color1.g << ", " << color1.b << endl;
            clusters.push_back( Cluster(0) );
            colors.push_back(color1);
            clusters[ clusters.size() - 1].push_back(i);
        }
    }
    //cout << clusters.size() << " clusters found." << endl;
    return 0;
}

int PlaneScan::readClusters(DataType &scan_type)
{
    size_t nrclusterpts = scan_type.size();
    if (nrclusterpts == 0) return -1;

    // Init possible indices. They do not have to be ordered!
    vector<size_t> indices(0);
    size_t first_index = scan_type[0];
    indices.push_back( first_index );

    // Init cluster with size 1, corresponding to color of first point
    clusters = Clusters(1);
    clusters[0].push_back(0);

    for (size_t i = 1; i < nrclusterpts; ++i) {
        bool index_equals = false;
        size_t index = scan_type[i];
        for (size_t j = 0; j < indices.size(); ++j) {
            index_equals = (index == indices[j]);
            if (index_equals) {
                clusters[j].push_back(i);
                break;
            }
        }
        if (!index_equals) {
            clusters.push_back( Cluster(0) );
            indices.push_back(index);
            clusters[ clusters.size() -1].push_back(i);
        }
    }
    return 0;
}

/**
 * @brief Calculates the distance between two clusters.
 * It does so by looking at each point of the first cluster, calculate its
 * closest point in the second cluster and store the min distance. Thus, the
 * complexity of the algorithm is O(nlog(k)) when n is the number of points in
 * the first cluster, and k is the number of points in the second cluster.
 * @param c1: The first cluster with n pts.
 * @param c2: The second cluster with k pts.
 */
double PlaneScan::clusterDistBkd(Cluster &c1, BkdTree *c2, int threadNum)
{
    double min_dist = DBL_MAX;
    for (size_t i : c1) {
        double *p = c2->FindClosest( points[i], DBL_MAX, threadNum );
        min_dist = min( min_dist, Dist2( p, points[i] ));
    }
    return min_dist;
}

// Uses brute force to calculate cluster to cluster distance.
// Thus, its complexity is O(n*k) (which is shit)
double PlaneScan::clusterDistBF(Cluster &c1, Cluster &c2)
{
    double min_dist = DBL_MAX;
    double dist;
    for (size_t i : c1) {
        for (size_t j : c2) {
            dist = Dist2(points[i], points[j]);
            if (i != j && dist < min_dist )
                min_dist = dist;
        }
    }
    return min_dist;
}

void PlaneScan::update_centroid_cache(int idxc, size_t idxp)
{
    for(int i = 0; i < 3; i++)
        centroid_cache[idxc][i] += points[idxp][i];
}

double * PlaneScan::centroidClusterMeanCached(int idxc)
{
    double *centroid = new double[3]; // Code smell, have to free mem manually
    for (int i = 0; i < 3; i++)
        centroid[i] = 0;
    const size_t csize = clusters[idxc].size();
    if (csize == 0) return centroid;
    for (int i = 0; i < 3; i++)
        centroid[i] = centroid_cache[idxc][i] / csize;
    return centroid;
}

// Adds the normal for the point at index <idxp> to the normal cluster mean cache for cluster at index <idxc>
void PlaneScan::update_ncm_cache(int idxc, size_t idxp)
{
    for(int i = 0; i < 3; i++)
        ncm_cache[idxc][i] += normals[idxp][i];
}

// Cached approach for calculating clusters. To be used ONLY and ONLY for CALCULATION of clusters.
// Does not work at all after some clusters are merged together or deleted.
double * PlaneScan::normalClusterMeanCached(int idxc)
{
    double *nmean = new double[3]; // Code smell... have to free mem manually
    nmean[0] = 0; nmean[1] = 0; nmean[2] = 0;
    const size_t csize = clusters[idxc].size();
    if (csize == 0) return nmean;
    nmean[0] = ncm_cache[idxc][0] / csize;
    nmean[1] = ncm_cache[idxc][1] / csize;
    nmean[2] = ncm_cache[idxc][2] / csize;
    Normalize3(nmean);
    return nmean;
}

// Calculate normal mean for a cluster, given only the points but no normal information.
double* PlaneScan::calcNormalOnDemand(Cluster &indices)
{
    vector<Point> temp;
    for (size_t i : indices)
    {
        double tmp[3];
        transform3(transMat, points[i], tmp);
        temp.push_back( Point(tmp) );
    }
    double *norm = new double[3];
    double eigen[3];
    calculateNormal(temp, norm, eigen);
    return norm;
}

// Calculates the mean normal of a cluster
double* PlaneScan::normalClusterMean(Cluster &indices)
{
    double *nmean = new double[3];
    nmean[0] = 0; nmean[1] = 0; nmean[2] = 0;
    if (indices.size() == 0) return nmean;

    // When clusters are read, there is no normal information.
    if (read_clusters) return calcNormalOnDemand(indices);

    for (size_t i = 0; i < indices.size(); ++i) {
        nmean[0] += normals[indices[i]][0];
        nmean[1] += normals[indices[i]][1];
        nmean[2] += normals[indices[i]][2];
    }
    nmean[0] /= indices.size();
    nmean[1] /= indices.size();
    nmean[2] /= indices.size();
    return nmean;
}

// Calculates the distance from point "p" to the cluster "cluster",
// i.e. the distance from the point "p" to a point in the cluster, where the distance is minimal
double PlaneScan::clusterPointDist(Cluster &cluster, double *p)
{
    double dist = DBL_MAX, d;
    for (size_t i = 0; i < cluster.size(); ++i) {
        d = distSqr(p, points[cluster[i]]);
        if ( d < dist) dist = d;
    }
    return dist;
}

double PlaneScan::clusterPointDist2(BkdTree *cluster, double *p, int threadNum)
{
    double *closest = cluster->FindClosest(p, DBL_MAX, threadNum);
    double d = Dist2(closest, p);
    return d;
}

void PlaneScan::initRG()
{

    // Setup clusters, initialize with one clusterlabelUseNorm
    clusters = Clusters();
    // Setup cluster caches
    ncm_cache = vector<double*>();
    centroid_cache = vector<double*>();
    // Setup dynamic cluster search-trees which support insertion of points
    bkd_forest = vector<BkdTree*>();
    // Setup static searchtree that is able to remove points
    tree = new KDtreeIndexed(points, nrpts);

    // Init first cluster with first point
    // clusters[0].push_back(0);
    // bkd_forest.push_back( new BkdTree() );
    // bkd_forest[0]->insert(points[0]);
    // // Also init caches
    // update_ncm_cache(0, 0);
    // update_centroid_cache(0, 0);
}

/**
 * Each point has the following states:
 * 1. It is in the queue.
 * 2. Point was visited, but not labled.
 * 3. Point was visited and labled.
 */
void PlaneScan::regionGrowing(int threadNum)
{
    // Use double ended queue to ensure correct iteration order
    deque<size_t> idx_queue;

    vector<bool> labled(nrpts);
    vector<bool> visited(nrpts);
    vector<bool> inqueue(nrpts);
    vector<bool> searched(nrpts);

    for (size_t i = 0; i < nrpts; ++i) {
        labled[i] = false;
        inqueue[i] = false;
        visited[i] = false;
        searched[i] = false;
    }
    // Only allow one cluster to grow at a time
    int current_growing_cluster = 0;
    double d2 = sqr(d_growth); // squared region growing param..

    /**
     * This outer for-loop ensures we always hit points where we can grow a region.
     * The selection of a new starting point for region growing happens if a cluster
     * is completed, i.e. it can no longer grow. Then, we pick the next best point.
     */
    for (size_t k = 0; k < nrpts; ++k) {

        // We skip points that have already been labled
        if (labled[k]) continue;

        // Visited points that have not been labeld are usually edge points (where normals suck)
        // Skip them, too.
        if (visited[k]) continue;

        // Put back the edge points, as they could potentially belong to another cluster
        for (size_t i = 0; i < nrpts; ++i) visited[i] = false;
        visited[k] = true;
        labled[k] = true;

        // Remove point from search tree, so we dont "grow into it" again
        tree->Remove( points[k], threadNum );

        // Open new cluster for new point
        clusters.push_back( Cluster() );
        // Each cluster gets a BkdTree (dynamicaly scalable kd-tree)
        bkd_forest.push_back( new BkdTree() );

        // Open new Normal Cache for cluster
        double *ncm_tail = new double[3];
        for(int i = 0; i < 3; i++)
            ncm_tail[i] = 0;
        ncm_cache.push_back( ncm_tail );

        // Open new Centroid Cache for cluster
        double *centroid_tail = new double[3];
        for(int i = 0; i < 3; i++)
            centroid_tail[i] = 0;
        centroid_cache.push_back( centroid_tail );

        // Update cluster and caches with current point
        current_growing_cluster = clusters.size() - 1;
        update_ncm_cache( current_growing_cluster, k );
        update_centroid_cache( current_growing_cluster, k);
        clusters[ current_growing_cluster ].push_back(k);
        if (!use_bruteforce) {
            bkd_forest[ current_growing_cluster ]->insert(points[k], threadNum);
        }

        // Grow from the new point using its k-nearest neighbors
        vector<size_t> knns = tree->kNearestNeighbors( points[k], k_neighbours, threadNum);
        searched[k] = true;
        for (size_t const i : knns) {
            if (!inqueue[i]) {
                idx_queue.push_front( i );
                inqueue[i] = true;
            }
        }

        /*
         * Queue points to grow a single cluster.
         * Once the cluster becomes full and can no longer grow, break this loop.
         */
        size_t grown = 0;
        while( !idx_queue.empty() ) {

            // Look at current point (front of the queue):
            int j = idx_queue.front(); idx_queue.pop_front();
            inqueue[j] = false;

            // Case: Point with index j is last element in the queue.
            // Then, we have to start a new growing process as follows
            if ( idx_queue.empty() ) {

                // We did not cluster any point: Finished cluster, break loop
                if (grown == 0) {
                    idx_queue.clear();
                    for (size_t i = 0; i < nrpts; ++i)
                        inqueue[i] = false;
                    break;

                // Some points got clustered in the last run. Continue growing
                } else {
                    PointCluster pts;
                    // Grow (almost) all points from the current cluster...
                    for (size_t i = 0; i < clusters[ current_growing_cluster ].size(); ++i) {
                        size_t growing_idx = clusters[ current_growing_cluster ][i];

                        // skip the ones already grown.
                        if (searched[growing_idx]) continue;

                        // Grow remaining (newly aquired) clusterpoints...
                        double *p = points[ growing_idx ];
                        vector<size_t> knns = tree->kNearestNeighbors( p, k_neighbours, threadNum );
                        // mark them as grown.
                        searched[growing_idx] = true;

                        // For all the new points in the range search...
                        for (size_t const index : knns) {
                            // if they aren't already in the queue and are, in fact, new...
                            if (!inqueue[index] && !visited[index]) {
                                // add them to the iteration queue.
                                idx_queue.push_front( index );
                                // mark point as 'in queue'
                                inqueue[index] = true;
                            }
                        }
                    }
                }
                // Reset, count growing points for new region-grow now.
                grown = 0;
            }

            // The point is already labled, skip!
            if (labled[j]) continue;
            if (visited[j]) continue;
            visited[j] = true;

            // Here we check if the point with index j belongs to the currently growing cluster
            //---------------------------------------------------------------------------------

            // Get point p with normal n
            double *p = points[j];
            double *n = normals[j];

            // Get BkdTree of current cluster for distance calculation
            BkdTree *clust = bkd_forest[current_growing_cluster];

            // Calculate angle alpha between point normal and mean cluster normal
            double *ncm = normalClusterMeanCached( current_growing_cluster );
            double alpha = angleBetweenNormals(ncm, n);

            // Calculate distance from point to cluster (neglect direction)
            double cpd2;
            if (use_bruteforce)
                cpd2 = clusterPointDist( clusters[current_growing_cluster], p );
            else
                cpd2 = clusterPointDist2( clust, p, threadNum );

            // Check distance of the point in normal direction (should be small)
            //TODO: DEBUG THIS LOCAL/GLOBAL MESS. SHOULD BE LOCAL!
            // double hesse2cluster = 0;
            // if (clusters[current_growing_cluster].size() >= (size_t)n_min_clusterpoints) {
            //     double* x = centroidClusterMeanCached( current_growing_cluster );
            //     // take cached cluster normal
            //     double v[3];
            //     for (int i = 0; i < 3; i++)
            //         v[i] = p[i] - x[i];
            //     hesse2cluster = Dot(v, ncm);
            //     delete[] x;
            //     x = NULL;
            // }

            // Dont forget to free this mem. code smells...
            delete[] ncm;
            ncm = NULL;

            // Adapt threshold for inclusion based on the distance of the point from the origin
            if (d_growth_max_adapt != -1.0) {
                double o[3] = {0, 0, 0};
                double r2_Orig = Dist2(p, o);
                d2 = adaptGrowth( r2_Orig );
            }
            //cout << "CPD: " << sqrt(cpd2) << ", Hesse: " << hesse2cluster << endl;
            // Check if we have found a good cluster for current point with index j
            if (deg(alpha) <= eps_similarity && cpd2 <= d2) { // && fabs(hesse2cluster) <= d_thickness) {

                // Update the clusters
                clusters[current_growing_cluster].push_back(j);
                if (!use_bruteforce){
                    bkd_forest[current_growing_cluster]->insert(points[j], threadNum);
                }

                // Remove point from static tree so we dont find it again with next knn search
                tree->Remove( points[j], threadNum );

                // update cluster caches
                update_ncm_cache(current_growing_cluster, j);
                update_centroid_cache(current_growing_cluster, j);

                // Mark point as labeled
                labled[j] = true;
                grown++;
            }
        }
    }
}

void PlaneScan::filter(Planes &hulls)
{
    // delete clusters that have too few points or bad eigenvals
    cout << "Deleting low quality clusters." << endl;

    Cache::iterator ncm_it = ncm_cache.begin();
    Forest::iterator bkd_it = bkd_forest.begin();
    Planes::iterator hulls_it = hulls.begin();
    for (Clusters::iterator it = clusters.begin();  // start at the beginning
         it != clusters.end();                      // stop at the end
         )                                          // no increment here
    {
        // Quality check
        //double *ncm = (*hulls_it)->n;
        double *eig = (*hulls_it)->eigen;
        double *ncm_cached = new double[3];
        for (int i = 0; i < 3; i++)
            ncm_cached[i] = (*ncm_it)[i] / it->size();
        Normalize3(ncm_cached);

        // "Too few points for a cluster" can have a different meaning according to scanned distance,
        // for some sensors. Ususally, far measurements are less dense, and cover larger areas.
        int n_thresh = n_min_clusterpoints;
        if ( n_min_clusterpoints_far != -1
            && (size_t)n_min_clusterpoints_far < it->size()
            && it->size() < (size_t)n_min_clusterpoints )
        {
            double origin[3] = {0, 0, 0};
            double centroid[3] = {0, 0, 0};
            for (size_t i = 0; i < it->size(); ++i) {
                centroid[0] += points[ it->at(i) ][0];
                centroid[1] += points[ it->at(i) ][1];
                centroid[2] += points[ it->at(i) ][2];
            }
            centroid[0] /= it->size();
            centroid[1] /= it->size();
            centroid[2] /= it->size();
            double dist2ToOrigin = Dist2(origin, centroid);
            n_thresh = adaptNrPts( dist2ToOrigin );
        }

        // TODO: ADD AREA FILTER HERE maybe
        // We put a threshold on cluster size because we want big, planar areas, not non-feature gibberish
        if (it->size() < (size_t)n_thresh
            // we check eigenvalue of the plane. if the plane is not "thin" enough, delete it
            || !eigenValueOK(eig, eigratio) )
        {
            it = clusters.erase(it);
            ncm_it = ncm_cache.erase( ncm_it );
            bkd_it = bkd_forest.erase(bkd_it);
            hulls_it = hulls.erase(hulls_it);
        }
        else
        {
            ++it; // increment as usual only if you dont delete the elem
            ++ncm_it;
            ++bkd_it;
            ++hulls_it;
        }
    }
}

void PlaneScan::filterPercentile(vector<NormalPlane*> &hulls, vector<size_t> &histogram, double percentile)
{
    // delete clusters that have too few points or bad eigenvals
    cout << "Apply histogram percentile filter. Keeping " << percentile_filter*100 << "% of pts" << endl;

    int hist_idx = 0;
    size_t nr_all_pts = 0;
    size_t nr_current_pts = 0;
    // Calculate Maximum number of points for each cluster
    for (uint i = 0; i < histogram.size(); ++i)
        nr_all_pts += histogram.at(i);

    Cache::iterator ncm_it = ncm_cache.begin();
    Forest::iterator bkd_it = bkd_forest.begin();
    Planes::iterator hulls_it = hulls.begin();
    for (Clusters::iterator it = clusters.begin();  // start at the beginning
         it != clusters.end();                      // stop at the end
         )                                          // no increment here
    {
        // Check condition, if true -> delete cluster
        nr_current_pts += histogram.at(hist_idx);
        hist_idx++;
        if ( nr_current_pts > percentile * nr_all_pts
          || (hist_idx > n_max_clusters && n_max_clusters != -1) )
        {
            it = clusters.erase(it);
            ncm_it = ncm_cache.erase( ncm_it );
            bkd_it = bkd_forest.erase(bkd_it);
            hulls_it = hulls.erase(hulls_it);
        }
        else
        {
            ++it; // increment as usual only if you dont delete the elem
            ++ncm_it;
            ++bkd_it;
            ++hulls_it;
        }
    }
}

/*
 * Calculates local clusters in the scan, using normal information of the points.
 * Normals get calculated automatically if not read beforehand.
 */
void PlaneScan::detectRG()
{
    cout << "Calculating clusters for scan" << identifier << endl;
    initRG();

#ifdef _OPENMP
    int thread_num = omp_get_thread_num();
#else
    int thread_num = 0;
#endif

    // Start region growing
    Timer<> clock; // millisecond clock
    clock.tick();
    regionGrowing( thread_num );
    clock.tock();
    cout << "Duration was " << clock.duration().count() / 1000.0 << "s" << endl;

    // Calculate and write histogram
    sortClustersBySize(clusters);
    std::vector<size_t> histogram = clusterHist(clusters);
    if (cluster_outpath[ cluster_outpath.size() - 1] != '/')
            cluster_outpath += "/";
    string hist_path = cluster_outpath + "hist" + identifier + ".dat";
    std::cout << "Writing histogram to " << hist_path << std::endl;
    writeHist( histogram, hist_path );

    // Build convex hulls to check overlap
    vector<NormalPlane*> hulls(0);
    for (size_t i = 0; i < clusters.size(); ++i) {
        vector<double*> hull;
        for (size_t j = 0; j < clusters[i].size(); ++j)
            hull.push_back(points[clusters[i][j]]);

        hulls.push_back(new NormalPlane(hull));
    }

    // Filter the clusters
    if (percentile_filter == -1.0)
        filter(hulls);
    else
        filterPercentile(hulls, histogram, percentile_filter);

    // Write clusters to scanfile
    cout << clusters.size() << " clusters found for scan" << identifier << endl;
    if ( cluster_outpath != "")
    {
        if (cluster_outpath[ cluster_outpath.size() - 1] != '/')
            cluster_outpath += "/";
        if (export_cluster_rgb)
            writeClustersToScanfileRGB(cluster_outpath);

        else
            writeClustersToScanfileType(cluster_outpath);
    }
}

/**
 * Saves the result of Clustering to a scanfile.
 * Can be visualized with "show -f uos_rgb -c"
 */
void PlaneScan::writeClustersToScanfileRGB(string path)
{
    // int sign = 1;
    ofstream ofile((path+"scan"+identifier+".3d").c_str(), ios_base::out);
    hsv color_hsv = {0, 1, 1};
    for (uint i = 0; i < clusters.size(); ++i)
    {
        rgb color_rgb = hsv2rgb( color_hsv );
        for(size_t j = 0; j < clusters[i].size(); ++j)
        {
            ofile << points[clusters[i][j]][0] << " "
                << points[clusters[i][j]][1] << " "
                << points[clusters[i][j]][2] << " "
                << (int) (color_rgb.r*255) << " "
                << (int) (color_rgb.g*255) << " "
                << (int) (color_rgb.b*255) << endl;
        }
        //color_hsv.h += sign*180.0;
        color_hsv.h += 360.0 / (double)clusters.size();
        //sign *= -1;
    }
    // Write corresponding pose.
    ofstream ofilepose((path+"scan"+identifier+".pose").c_str(), ios_base::out);
    ofilepose << rPos[0] << " " << rPos[1] << " " << rPos[2] << " ";
    ofilepose << deg(rPosTheta[0]) << " " << deg(rPosTheta[1]) << " " << deg(rPosTheta[2]) << endl;
}

// Use this function to export clusters.
void PlaneScan::writeClustersToScanfileType(string path)
{
    ofstream ofile((path+"scan"+identifier+".3d").c_str(), ios_base::out);
    for (uint i = 0; i < clusters.size(); ++i)
    {
        for(size_t j = 0; j < clusters[i].size(); ++j)
        {
            ofile << points[clusters[i][j]][0] << " "
                << points[clusters[i][j]][1] << " "
                << points[clusters[i][j]][2] << " "
                << i << endl; // cluster index as type
        }
    }
    // Write corresponding pose.
    ofstream ofilepose((path+"scan"+identifier+".pose").c_str(), ios_base::out);
    ofilepose << rPos[0] << " " << rPos[1] << " " << rPos[2] << " ";
    ofilepose << deg(rPosTheta[0]) << " " << deg(rPosTheta[1]) << " " << deg(rPosTheta[2]) << endl;
}

vector<vector<double*> > PlaneScan::getGlobalClusters()
{
    vector<vector<double*> > result;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        vector<double*> cluster;
        for (size_t j = 0; j < clusters[i].size(); j++)
        {
            double *p = new double[3];
            transform3(transMat, points[clusters[i][j]], p);
            cluster.push_back( p );
        }
        result.push_back(cluster);
    }
    return result;
}

vector<vector<double*> > PlaneScan::getLocalClusters()
{
    vector<vector<double*> > result;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        vector<double*> cluster;
        for (size_t j = 0; j < clusters[i].size(); j++)
            cluster.push_back( points[clusters[i][j]] );
        result.push_back(cluster);
    }
    return result;
}


// Outdated. We dont need to merge anything.

// void PlaneScan::merge(Planes &hulls, int threadNum)
// {
//     // Cleanup, merge clusters with similar normal if they are near enough
//     cout << "Merge similar clusters" << endl;

//     // We need all cluster iterators, since we have indexed clusters:
//     // Bkd-tree iterators
//     Forest::iterator it_bkd = bkd_forest.begin();
//     Forest::iterator it_bkd2 = bkd_forest.begin();
//     // Normal mean cache iterators
//     Cache::iterator it_ncm = ncm_cache.begin();
//     Cache::iterator it_ncm2 = ncm_cache.begin();
//     // Convex hull iterators
//     Hulls::iterator it_hull = hulls.begin();
//     Hulls::iterator it_hull2 = hulls.begin();

//     // Iterate the clusters
//     for (Clusters::iterator it = clusters.begin();
//          it != clusters.end();
//          // we increment / delete manualy.
//          )
//     {
//         // Normal vector of first cluster
//         double *ncm1 = (*it_hull)->n; // faster than using ncm_cache (we dont need division)

//         // Reference iterator is next of current iterator
//         it_ncm2 = next(it_ncm);
//         it_hull2 = next(it_hull);
//         it_bkd2 = next(it_bkd);
//         for (Clusters::iterator it2 = next(it);
//              it2 != clusters.end();
//             // conditional increment below...
//             )
//         {
//             // Normal vector of second cluster
//             double *ncm2 = (*it_hull2)->n;

//             NormalPlane *hull1 = (*it_hull);
//             NormalPlane *hull2 = (*it_hull2);

//                 // If planes are not the same...
//             if (it != it2
//                 // ...have similar normals...
//                 && similar(ncm1, ncm2)
//                 // ...have small inter-point distance to each other...
//                 && clusterDistBkd((*it), (*it_bkd2), threadNum) <= sqr(d_growth)
//                 // ...and overlap...
//                 && hull1->overlaps(hull2) )
//             {
//             // Case "merge" (delete current, and set next elem)
//                 merge(it, it2);
//                 // Combine normal mean caches
//                 merge(it_ncm, it_ncm2);
//                 // Combine hulls
//                 merge(it_hull, it_hull2);
//                 // Merge bkd trees
//                 size_t npts_new = it->size();
//                 double **data_new = new double*[npts_new];
//                 for (size_t i = 0; i < npts_new; ++i)
//                     data_new[i] = points[it->at(i)];
//                 (*it_bkd) = new BkdTree(data_new, npts_new);

//                 // Erase ref elem and go on with the next one.
//                 it2 = clusters.erase(it2);
//                 it_ncm2 = ncm_cache.erase(it_ncm2);
//                 it_hull2 = hulls.erase(it_hull2);
//                 it_bkd2 = bkd_forest.erase(it_bkd2);
//             }
//             // do not merge (just iterate further)
//             else
//             {
//                 ++it2;
//                 ++it_ncm2;
//                 ++it_hull2;
//                 ++it_bkd2;
//             }
//         }
//         // prepare next outer loop
//         ++it;
//         ++it_ncm;
//         ++it_hull;
//         ++it_bkd;
//     }
// }


// void PlaneScan::merge(Clusters::iterator &it, Clusters::iterator &it2)
// {
//     it->insert( // insert
//         it->end(), // at the end of the current cluster
//         make_move_iterator(it2->begin()), // everything from begining
//         make_move_iterator(it2->end()) // to end of the ref. cluster.
//     );
// }

// void PlaneScan::merge(Cache::iterator &it_ncm, Cache::iterator &it_ncm2)
// {
//     (*it_ncm)[0] += (*it_ncm2)[0];
//     (*it_ncm)[1] += (*it_ncm2)[1];
//     (*it_ncm)[2] += (*it_ncm2)[2];
// }

// void PlaneScan::merge(Hulls::iterator &it_hull, Hulls::iterator &it_hull2)
// {
//     (*it_hull)->mergeWith((*it_hull2));
// }