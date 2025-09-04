#ifndef __ADADELTA_H_
#define __ADADELTA_H_

#include <utility>

#include "opt/optimizer.h"

#define EPS_CONVERGENCE (-1)
#define MAX_ITER (-1)
#define ALPHA 0.01
#define RPOS_ALPHA_SCALE 10
#define MAX_FRAMES 20

class AdaDelta : public Optimizer
{
public:

    // Inheritance stuff:

    AdaDelta(PlaneScan* p, Dimensions d); // overrides default constructor
    void operator()(void); // overrides default iteration mechanism

    // Constructors:
    AdaDelta();
    AdaDelta(PlaneScan*);
    AdaDelta(Matrix&, Matrix&);
    ~AdaDelta();


    // Run <arg> iterations of gradient descent
    void iterate(int);
    void operator()(int);

    // Run until state difference between iterations is smaller then <arg>
    void iterate(double);
    void operator()(double);

    // Combination of both top approaches
    /*
     * The overloaded operator() allows you to initialize
     * and run an AdaDelta instance in the following manner:
     *
     * AdaDelta::setMaxIter(10000);
     * AdaDelta::setEpsConvergence(0.0001);
     * AdaDelta iter( my_planescan, Dimensions::ALL );
     * iter();
     *
     * Use -E and -i to set these.
     * Multiple criteria are supported.
     * Running without criteria runs 1000 iterations.
     */
    void iterate();
    void iterate(int, double);
    void iterateAuto();

    ColumnVector getResult();
    void relabel();

    // Instance setters
    void setPlaneScan(PlaneScan*);
    void setState(ColumnVector&);
    void setJacobian(Matrix&);

    // Static setters for programm options:
    static void setAuto(bool);
    static void setAlphaInit(double);
    static void setRPosAlphaScale(double);

private:

    ColumnVector X0; // initial state. keep in order to reset if needed.
    ColumnVector X1; // State vector of the very last iteration.
    ColumnVector Xmin; // always holds the state where RMSE is minimum.
    ColumnVector alpha; // gradient descent factor / convergence rate.
    ColumnVector Ga; // decaying squared accumulated gradient, used for adaptive alpha
    ColumnVector Xa; // decaying squared accumulated update, used for adaptive alpha
    Matrix J; // Jacobian (Gradient)

    double a; // weight that determines the descent agression

    static double _alpha; // static alpha used for initilazition
    static double _rPos_a_scale; // scales alpha applied to translation
    static bool _auto; // auto detect optimum alpha, used for speedup

    // Computes only one step of gradient descent. Updates Transformation and writes animation.
    void step();

    // Resets the PlaneScan object to its initial state, i.e. the orientation in the .pose file
    void reset();

    // Update jacobian from state vector X
    void updateJacobian();

    // Update alpha from accumulated gradients and outputs
    void adaptAlpha();

};

#endif //__ADADELTA_H_