#ifndef __CLUSTERMATCHER_H_
#define __CLUSTERMATCHER_H_

#include <memory>

#include "model/planescan.h"
#include "match/matcher.h"

class ClusterMatcher : public Matcher
{
public:
    double w_overlap;
    double w_alpha;
    double w_hesse;
    double w_ppd;
    double w_eigen;
    double eps_dist;
    double eps_ppd;
    double eps_sim;
    double eigratio;
    ClusterMatcher(PlaneScan* p, double w_overlap = 0.8,
                                 double w_alpha = 1.0,
                                 double w_hesse =1.0,
                                 double w_ppd = 0.3,
                                 double w_eigen = 0.1,
                                 double eps_dist = 25,
                                 double eps_ppd = 25,
                                 double eps_sim = 5,
                                 double eigratio = 0.05);
    void match();
};

#endif //__CLUSTERMATCHER_H_