#include "match/clustermatcher.h"

ClusterMatcher::ClusterMatcher(PlaneScan* p,
    double w_ov,
    double w_al,
    double w_he,
    double w_pp,
    double w_ei,
    double eps_d,
    double eps_p,
    double eps_s,
    double eigr
)
: Matcher(p)
{
    this->w_overlap = w_ov;
    this->w_alpha = w_al;
    this->w_hesse = w_he;
    this->w_ppd = w_pp;
    this->w_eigen = w_ei;
    this->eps_dist = eps_d;
    this->eps_ppd = eps_p;
    this->eps_sim = eps_s;
    this->eigratio = eigr;
}

void ClusterMatcher::match()
{
    /*
     * Setting up memory for plane-to-plane simliarity
     */
    const int nplanes = PlaneScan::allPlanes.size();
    const int nclusters = ps->clusters.size();
    int cccount[nclusters][nplanes]; // cluster correspondence count
    double overlap[nclusters][nplanes]; // overlap is percentage of cluster correspondence count
    double distancecount[nclusters][nplanes][2]; // minimum distance measures between polygons
    double eigendiffs[nclusters][nplanes]; // differences between eigenvalues
    double* hull_normals[nclusters];
    double* hull_eigs[nclusters];
    // Initialize with zero
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int c = 0; c < nclusters; ++c) {

        // Calculate the global Normal vector of the cluster
        // use eigenvalues to check quality of the normal
        vector<Point> temp;
        for (size_t i : ps->clusters[c]) {
            double tmp[3];
            transform3(ps->transMat, ps->points[i], tmp);
            temp.push_back( Point(tmp) );
        }
        double* tmp_norm = new double[3];
        double* tmp_eig = new double[3];
        calculateNormal(temp, tmp_norm, tmp_eig);
        hull_normals[c] = tmp_norm;
        hull_eigs[c] = tmp_eig;

        // Initialize above arrays to some dummy vals
        for (int p = 0; p < nplanes; ++p) {
            cccount[c][p] = 0;
            eigendiffs[c][p] = eigendiff(hull_eigs[c], PlaneScan::allPlanes[p]->eigen);
            distancecount[c][p][0] = -1; // hesse distance
            distancecount[c][p][1] = -1; // ppd (-1 is dummy for uninitialized)
        }
    }

    // Iterate the clusters...
    for (int k = 0; k < nclusters; ++k) {
        // ...check if it belongs to planes...
        for (int j = 0; j < nplanes; ++j) {
            // ...for each point in the cluster...
#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic)
#endif
            for (size_t q = 0; q < ps->clusters[k].size(); ++q) {
                // ...and if so...
                double *p = ps->points[ps->clusters[k][q]]; // Point p
                double hesse = fabs( PlaneScan::dist2Plane(p, ps->transMat, PlaneScan::allPlanes[j]) );
                double ppd = fabs( PlaneScan::projDist2Plane(p, ps->transMat, PlaneScan::allPlanes[j]) );
                // We seek to store min distance, then find minimum of those later for optimal score function
                if (distancecount[k][j][0] == -1 || hesse < distancecount[k][j][0])
                    #pragma omp critical
                    distancecount[k][j][0] = hesse;
                if (distancecount[k][j][1] == -1 || ppd < distancecount[k][j][1])
                    #pragma omp critical
                    distancecount[k][j][1] = ppd;
                // If the point fullfills the overlap criteria...
                if ( fabs(hesse) < eps_dist && fabs( ppd ) < eps_ppd ) {
                    // ...increment the k-th cluster count for the j-th plane.
                    #pragma omp critical
                    cccount[k][j]++;
                }
            }
        } // Uncomment following line if you need that. Shows cccount as a table
    } // debugCCCount(cccount, nplanes, nclusters);

    // Find max values for score function:
    double hesse_max = 0, ppd_max = 0;
    double eigdiff_max = 0;
    for (int c=0;c<nclusters;c++) {
        for (int p=0;p<nplanes;p++) {
            overlap[c][p] = (double)cccount[c][p] / ps->clusters[c].size();
            hesse_max = max( hesse_max, distancecount[c][p][0] );
            ppd_max = max( ppd_max, distancecount[c][p][1] );
            eigdiff_max = max( eigdiff_max, eigendiffs[c][p]);
        }
    }

    // Look at each cluster
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int c = 0; c < nclusters; ++c)
    {
        double *ncm = hull_normals[c];
        double *eigen = hull_eigs[c];

        // Skip bad clusters (just to be sure)
        if (!eigenValueOK(eigen, eigratio)) continue;

        int p_best_score = -1;
        // Search for a plane that has a similar normal and a big overlapping region
        double best_score = 0;
        for (int p = 0; p < nplanes; ++p)
        {
            // Setup a score function for every parameter
            // i.e. alpha, hesse, ppd and overlap
            double *n = PlaneScan::allPlanes[p]->n;
            double alpha = angleBetweenNormals(ncm, n);
            int count = cccount[c][p];
            double hesse = distancecount[c][p][0];
            double ppd = distancecount[c][p][1];
            double eigdiff = eigendiffs[c][p];

            // Score will have value between 0 (worst) and 1 (best)
            double score_alpha = 1 - (alpha / M_PI_2);
            double score_overlap = overlap[c][p];
            double score_hesse = hesse_max == 0 ? 1 : 1 - hesse / hesse_max;
            double score_ppd = ppd_max == 0 ? 1 : 1 - ppd / ppd_max;
            double score_eig = eigdiff_max == 0 ? 1 : 1 - eigdiff / eigdiff_max;

            double score =  w_alpha * score_alpha
                          + w_overlap * score_overlap
                          + w_hesse * score_hesse
                          + w_eigen * score_eig
                          + w_ppd * score_ppd;

            // Save best score if there is an overlap.
            if (count != 0 && score > best_score) {
                best_score = score;
                p_best_score = p;
            }
        }

        // Find corresponding parameters for the best plane
        double *n;
        double alpha = M_2_PI, hesse=DBL_MAX, ppd=DBL_MAX;
        int count = 0;
        if (p_best_score != -1)
        {
            n = PlaneScan::allPlanes[p_best_score]->n;
            alpha = angleBetweenNormals(ncm, n);
            hesse = distancecount[c][p_best_score][0];
            ppd = distancecount[c][p_best_score][1];
            count = cccount[c][p_best_score];
        }

        // Convert Local Cluster (which stores indices) to Global PointCluster (which stores global points).
        PointCluster *pcluster = new PointCluster();
        PointCluster *pcluster_local = new PointCluster();
        for (size_t i = 0; i < ps->clusters[c].size(); i++) {
            double *p_global = new double[3];
            transform3( ps->transMat, ps->points[ps->clusters[c][i]], p_global);
            pcluster->push_back( p_global );
            pcluster_local->push_back( ps->points[ps->clusters[c][i]]  );
        }
        // Check thresholds
        #pragma omp critical
        {
            if (count != 0
                && deg(alpha) <= eps_sim
                && fabs(hesse) <= eps_dist
                && fabs(ppd) <= eps_ppd)
            {
                // Insert cluster-2-plane match
                ps->global_matches.insert(
                    Match(
                        PlaneScan::allPlanes[p_best_score],
                        pcluster
                    )
                );
                ps->local_matches.insert(
                    Match(
                        PlaneScan::allPlanes[p_best_score],
                        pcluster_local
                    )
                );
                for (size_t i = 0; i < ps->clusters[c].size(); i++)
                {

                    ps->correspondences.push_back(
                        Correspondence(
                            ps->points[ps->clusters[c][i]],
                            PlaneScan::allPlanes[p_best_score]
                        )
                    );

                    double p_trans[3];
                    transform3(ps->transMat, ps->points[ps->clusters[c][i]], p_trans);
                    double *nearest = PlaneScan::allPlanes[p_best_score]
                                            ->search_tree
                                            ->FindClosest(
                                                p_trans,
                                                sqr(eps_dist));

                    if (nearest) {
                    ps->point_pairs.push_back(
                            PointPair(
                                ps->points[ps->clusters[c][i]], nearest
                            )
                        );
                    }

                }

            }
            else if ( eigenValueOK(eigen, eigratio) ) // insert "mismatches", i.e. clusters that did not get matched
                ps->global_mismatches.insert( pcluster );
        } // end pragma critical
    } // end for all custers

    // Cleanup
    for (int i = 0; i < nclusters; ++i) {
        delete hull_normals[i];
        delete hull_eigs[i];
    }

} // end func match()