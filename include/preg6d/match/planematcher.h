#ifndef __PLANEMATCHER_H_
#define __PLANEMATCHER_H_

#include <queue>
#include <list>
#include <utility>

#include "model/planescan.h"
#include "match/matcher.h"
#include "tree/bkdplanes.h"

typedef std::pair<NormalPlane*, NormalPlane*> PlanePair;

struct EnergyPlanePair {
    EnergyPlanePair(double da, double dh, double dp, PlanePair &p) :
        delta_alpha(da), delta_hesse(dh), delta_ppd(dp), plane_pair(p)
        {
            // TODO: NORMALIZE SINGLE ENERGIES?
            total_energy = delta_alpha+delta_hesse+delta_ppd;
        }
    double delta_alpha;
    double delta_hesse;
    double delta_ppd;
    PlanePair plane_pair;

    double total_energy;
};

typedef std::list<EnergyPlanePair> PlanePairs;

class PlaneMatcher : public Matcher
{
public:

    PlaneMatcher(PlaneScan* ps, double eps_hesse, double eps_ppd, double eps_sim);
    void match();

private:
    double eps_hesse;
    double eps_ppd;
    double eps_sim; // similarity crit for normals in deg

    bool sanityPass(EnergyPlanePair &epp);
};

#endif // __PLANEMATCHER_H_