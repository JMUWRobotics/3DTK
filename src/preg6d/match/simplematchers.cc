#include "match/simplematchers.h"

EuklidMatcher::EuklidMatcher(PlaneScan* ps, bool use_cm)
: Matcher(ps)
{
    use_correspondence_min = use_cm;
}

void EuklidMatcher::match()
{
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < ps->nrpts; ++i)
    {
        Planes *planeList = new Planes(); // can be also empty or bigger than one
        for ( const auto& plane : PlaneScan::allPlanes ) {
            if ( PlaneScan::isInPlane(ps->points[i], ps->transMat, plane) )
                planeList->push_back( plane );
        }
        if (use_correspondence_min)
        {
            // find the plane with the minimum distance amongst all candidates
            if (!planeList->empty())
            {
                NormalPlane* min = planeList->at(0);

                double ppd_min = __DBL_MAX__, hesse_min = __DBL_MAX__;
                for ( size_t j = 0; j < planeList->size(); ++j)
                {
                    double hesse = fabs(PlaneScan::dist2Plane( ps->points[i], ps->transMat, planeList->at(j)));
                    double ppd = fabs(PlaneScan::projDist2Plane( ps->points[i], ps->transMat, planeList->at(j)));
                    if ( ppd < ppd_min )
                    {
                        min = planeList->at(j);
                        ppd_min = ppd;
                        hesse_min = hesse;
                    }
                    // ppd == 0 in this case
                    else if ( ppd == ppd_min && hesse < hesse_min)
                    {
                        min = planeList->at(j);
                        hesse_min = hesse;
                    }
                }

                // Store a correspondence (unordered)
                #pragma omp critical
                ps->correspondences.push_back( Correspondence( ps->points[i], min ) );

                // min plane might not have a search tree in preg6D.
                // This is if the plane model was read globaly from convex hull (e.g., export of bin/planes).
                // Works just fine with cluster-based planes, since points for the search tree are available.
                if (PlaneScan::use_clustering)
                {
                    double p_trans[3];
                    transform3(ps->transMat, ps->points[i], p_trans);
                    double *nearest = min->search_tree->FindClosest(p_trans, sqr(PlaneScan::_eps_dist));

                    if (nearest) {
                    #pragma omp critical
                    ps->point_pairs.push_back(
                            PointPair(
                                ps->points[i], nearest
                            )
                        );
                    }
                }
            }
        }
        else
        {
            if (planeList->size() == 1) {
                #pragma omp critical
                ps->correspondences.push_back( Correspondence( ps->points[i], planeList->at(0)) );

                if (PlaneScan::use_clustering)
                {
                    double p_trans[3];
                    transform3(ps->transMat, ps->points[i], p_trans);
                    double *nearest = planeList->at(0)->search_tree->FindClosest(p_trans, sqr(PlaneScan::_eps_dist));
                    if (nearest) {
                    #pragma omp critical
                    ps->point_pairs.push_back(
                            PointPair(
                                ps->points[i], nearest
                            )
                        );
                    }
                }
            }

        }
        delete planeList;

    }
}

NormalMatcher::NormalMatcher(PlaneScan* p, double eps)
: Matcher(p)
{
    this->eps_sim = eps;
}

void NormalMatcher::match()
{
    // For all points ...
    for (size_t i = 0; i < ps->nrpts; ++i)
    {
        double min_alpha = DBL_MAX;
        size_t min_j;

        for (size_t j = 0; j < PlaneScan::allPlanes.size(); ++j)
        {
            // ... that fall into a distance band of a plane ...
            if ( PlaneScan::isInPlane( ps->points[i], ps->transMat, PlaneScan::allPlanes[j]) )
            {  //TODO: NORMALS MIGHT NOT BE AVAILABLE!
                // ... look at the angle between point normal and plane normal ...
                double alpha = angleBetweenNormals(ps->normals[i], PlaneScan::allPlanes[j]->n);

                if (alpha < min_alpha)
                {
                    // ... and save the minimum.
                    min_j = j;
                    min_alpha = alpha;
                }
            }
        }
        // Establish correspondence between point and plane with min. angle
        if ( min_alpha != DBL_MAX && deg(min_alpha) < eps_sim)
            ps->correspondences.push_back( Correspondence(ps->points[i], PlaneScan::allPlanes[min_j]) );
    }
}