/** @file
 *  @brief This file implements an interface for optimizers.
 *  The optimizer class exists to minimize point distances.
 *  See adadelta6d.h or newtons6d.h for examples on how to derive the class properly.
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 */

#ifndef __OPTIMIZER_H_
#define __OPTIMIZER_H_

#include <math.h>

#define WANT_STREAM // want to stream matrices on std output :)
#include <newmat/newmat.h>
#include <newmat/newmatio.h>
using namespace NEWMAT;

#include "model/planescan.h"

// Number of past states to compare with for stop condition
#define COMP_N_FRAMES 5

enum Dimensions {
    ALL,            // performs 6D optimization
    ROLLING,        // only optimizes rotX, as well as y and z (forward rolling motion, e.g. L.U.N.A sphere)        DESCENDING,     // only optimizes rotY and y (e.g. DAEDALUS descent phase)
    DESCENDING,     // only optimizes rotY and Y
    ROTATING,       // only optimizes rotY (e.g. terrestial stationary RIEGL laserscanner)
    UDLOCK,         // no correction in upwards / downwards direction. Usefull for Jaspers Sphere.
    SLIDING,        // only optimizes translation on XY plane (no up/down and no rotation)
    RADLER,         // only optimizes XY plane translation and yaw.
    NOROT           // only optimizes translation
};

class Optimizer
{
public:

    // Override constructors for initialatition
    Optimizer() {};
    Optimizer(PlaneScan*, Dimensions) {};

    /**
     * These interfaces need to be implemented by deriving classes:
     */
    virtual ~Optimizer() {};
    // Override operator() for iteration
    virtual void operator()(void) {};

    // converts dimensions enum to bool array
    static bool* convertDimensionsToBoolArray(Dimensions);

    // Update scan to state X and write .frames entry
    void updateScan();

    double getRSE();

    // Program options
    static void setMaxIter(int);
    static void setUpdateCor(int);
    static void setEpsConvergence(double);
    static void setEpsKernel(double);
    static void setAnim(int);
    static void setQuiet(bool);
    static void setUseP2P(bool);

protected:

    PlaneScan* ps;
    ColumnVector X; // state vector
    ColumnVector dX; // state change

    deque<ColumnVector> Xk; // previous states

    Dimensions dim;
    double lastErr;
    int frames_count;

    static double _eps_convergence;
    static double _eps_kernel;
    static int _max_iter; // run this many iterations before updating correspondences
    static int _update_cor; // update correspondences this often
    static int _anim;
    static bool _quiet;
    static bool _use_p2p;

    void lock(); // lock dX according to enum DIMENSIONS
    bool stop_condition(ColumnVector state, double eps = _eps_convergence);

    static double weightLoss(double);

};

#endif // __OPTIMIZER_H_