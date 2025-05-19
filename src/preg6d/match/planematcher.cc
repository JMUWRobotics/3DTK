#include "match/planematcher.h"

#define INV_M_PI_2 0.63661977236

PlaneMatcher::PlaneMatcher(PlaneScan* ps, double eps_hesse, double eps_ppd, double eps_sim)
: Matcher(ps)
{
    this->eps_hesse = eps_hesse;
    this->eps_ppd = eps_ppd;
    this->eps_sim = eps_sim;

}

bool PlaneMatcher::sanityPass(EnergyPlanePair &epp) {
    // TODO: Think about that later
    return true;
}

void PlaneMatcher::match()
{
#ifdef _OPENMP
    int thread_num = omp_get_thread_num();
#else
    int thread_num = 0;
#endif

    // Setup iterators
    PointClusters global_clusters = ps->getGlobalClusters();
    PointClusters local_clusters = ps->getLocalClusters();
    Clusters cluster_indices = ps->clusters;
    PointClusters::iterator local_cluster = local_clusters.begin();
    Clusters::iterator indices = cluster_indices.begin();

    // List of candidates

    // Get global clusters
    int j = 1; // the j variable is only for debug output.
    for ( PointClusters::iterator cluster = global_clusters.begin();
          cluster != global_clusters.end();
          // the ++j is only for debug output.
          ++j// no increment here, we want to increment both global and local clusters at the bottom
        ) {

        // Calculate centroid and normal for each global cluster
        NormalPlane *data_plane = new NormalPlane( *cluster );
        PlanePairs candidates;

        // Find maximum distant plane
        for (NormalPlane* model_plane : PlaneScan::allPlanes)
        {

            double hesse = data_plane->hesseDist2Plane(model_plane);
            double ppd =  data_plane->projDist2Plane(model_plane);
            double alpha = deg(angleBetweenNormals(model_plane->n, data_plane->n));

            // Smart skip
            if (eps_sim < alpha
             || eps_hesse < hesse
             || eps_ppd < ppd)
                continue;

            // Else: calc energy for pair and store

            double hesse_min, hesse_max, ppd_min, ppd_max;
            data_plane->getMinMaxHesseTo(model_plane, hesse_min, hesse_max);
            double delta_hesse = hesse_max - hesse_min;
            data_plane->getMinMaxProjDistTo(model_plane, ppd_min, ppd_max);
            double delta_ppd1 = ppd_max - ppd_min;
            model_plane->getMinMaxProjDistTo(data_plane, ppd_min, ppd_max);
            double delta_ppd2 = ppd_max - ppd_min;
            double delta_ppd = max( delta_ppd1, delta_ppd2 );

            // This could actually be wrong.
            double energy_alpha = alpha ;
            double energy_hesse = delta_hesse;//hesse / (1e-9 + hesse + delta_hesse);
            double energy_ppd = delta_ppd;

            PlanePair pl_pair(data_plane, model_plane);
            // Store
            candidates.insert( end(candidates),
                EnergyPlanePair(
                    energy_alpha,
                    energy_hesse,
                    energy_ppd,
                    pl_pair
                )
            );
        }

        //cout << "Candidates found: " << candidates.size() << endl;
        candidates.sort(
            [=] (EnergyPlanePair &pair1, EnergyPlanePair &pair2) -> bool {
                return pair1.total_energy < pair2.total_energy;
            }
        );

        // // Use this to debug correspondences
        // int i = 0;
        // for (EnergyPlanePair &epp : candidates) {
        //     NormalPlane* tmp_plane = epp.plane_pair.second;

        //     cout << "------------Combi("<< j << ", " << ++i << ")--------------" << endl;
        //     cout << "dAlpha=" << epp.delta_alpha << endl;
        //     cout << "dHesse=" << epp.delta_hesse << endl;
        //     cout << "dPPD=" << epp.delta_ppd << endl;
        //     cout << "Hesse=" << data_plane->hesseDist2Plane(tmp_plane)<< endl;
        //     cout << "PPD=" << data_plane->projDist2Plane(tmp_plane) << endl;
        //     cout << "TOTAL ENERGY:" << epp.total_energy << endl;
        //     cout << "-----------------------------------" << endl;
        // }

        PointCluster* glob_clust = new PointCluster();
        for (double * p : (*cluster))
            glob_clust->push_back(p);

        // Setup correspondences
#pragma omp critical
{
        if ( candidates.size() > 0 )
        {
            NormalPlane* best = candidates.front().plane_pair.second;

            PointCluster* loc_clust = new PointCluster();
            for (double * p : (*local_cluster))
                loc_clust->push_back(p);

            // Plane-2-Plane Matching

            // Local Data -> Global Model Plane-2-Plane Match
            ps->local_matches.insert(
                Match(
                    best, loc_clust
                )
            );

            // Global Data -> Global Model Plane-2-Plane Match
            ps->global_matches.insert(
                Match(
                    best, glob_clust
                )
            );

            // Point-based Matching
            for (size_t index : (*indices))
            {
                double *p = ps->points[index];
                double p_trans[3];
                transform3(ps->transMat, p, p_trans);

                double *nearest = best->search_tree
                                ->FindClosest(
                                    p_trans,
                                    sqr(eps_hesse),
                                    thread_num
                                );

                if ( PlaneScan::dist2Plane(p, ps->transMat, best) <= eps_hesse
                  && PlaneScan::projDist2Plane(p, ps->transMat, best) <= eps_ppd )
                {
                    ps->correspondences.push_back(
                        Correspondence(
                            p, best
                        )
                    );

                    if(nearest) {
                        // Local Data -> Global Model Point-2-Point Match
                        ps->point_pairs.push_back(
                            PointPair(
                                p, nearest
                            )
                        );
                    }
                }
            }
            delete data_plane;

        // Case: No best Match found.
        } else {

            ps->global_mismatches.insert(
                glob_clust
            );

        }
} // end #pragma omp_critical
        candidates.clear();
        ++indices;
        ++cluster;
        ++local_cluster;
    }
}
