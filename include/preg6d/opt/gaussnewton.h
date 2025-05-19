#ifndef __GAUSS_NEWTON_H_
#define __GAUSS_NEWTON_H_

#include <utility>
#include <cassert>

#include "opt/optimizer.h"

class GaussNewton : public Optimizer
{
public:
    GaussNewton(PlaneScan* p, Dimensions d); // overrides default constructor
    void operator()(void); // overrides default iteration mechanism

    void iterate();

    GaussNewton();
    GaussNewton(PlaneScan*);

    ~GaussNewton(){}
private:
    void updateGradient();
    Matrix J;
    DiagonalMatrix W;
    ColumnVector R;
    void rematch();

};

#endif