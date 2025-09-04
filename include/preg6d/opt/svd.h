#ifndef __OPT_SVD_H
#define __OPT_SVD_H

#include "opt/optimizer.h"

class PlaneSVD : public Optimizer
{
public:

    PlaneSVD(PlaneScan* p, Dimensions d);
    void operator()(void);

    PlaneSVD();
    PlaneSVD(PlaneScan*);
    ~PlaneSVD();

    void iterate();

protected:

    /**
     * @brief Set up correspondences from current scan to global model.
     * Then use correspondences to perform SVD to get alignment.
     * Used inside the iterate()-function of this class.
     * @param it: (Optional) Print this iteration along with RMSE.
     */
    void step(int it=-1);
    void rematch();
    void calc_centroids(double* cm, double* cd, double* dtr = 0);
};

#endif