/** @file newtons6d.h
 *
 * @brief
 * the distance D composes the error function:
 * Error = sum of all D^2
 * Jacobian = sum of all 2 * D * D'
 * Hessian = sum of all 2 * (D'*D'^t + D*D'')
 * (hint: use product rule differentiation, D depends on six transform params)
 * We minimize the Error as follows:
 * X_k+1 = X_k - Hessian^-1 * Jacobian
 *
 * @author Fabian Arzberger, JMU, Germany.
 *
 * Released under the GPL version 3.
 *
 */

#ifndef __NEWTONS_H_
#define __NEWTONS_H_

#include <cassert>

// Use (void) to silent unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

#include "opt/optimizer.h"

class Newtons : public Optimizer
{
public:
    Newtons(PlaneScan* p, Dimensions d);
    void operator()(void);

    Newtons();
    Newtons(PlaneScan*);

    void iterate();

protected:
    double err_min;
    ColumnVector X_min; // stores state where error func. was minimal

    ColumnVector J; // Jacobian of hesse distance
    Matrix H; // Hessian of hesse distance

    ColumnVector J_c; // Jacobian of centroid distance
    Matrix H_c; // Hessian of centroid distance

    ColumnVector JJ; // combined jacobian
    Matrix HH; // combined hessian

    void updateHesseGradient();
    void updateCentroidGradient();

    void rematch();

    void overlapGradients();
};

#endif //__NEWTONS_H_
