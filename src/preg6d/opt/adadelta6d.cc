#include "opt/adadelta6d.h"


#define ALPHA_DECAY 1.0

double AdaDelta::_alpha;
double AdaDelta::_rPos_a_scale;
bool AdaDelta::_auto;

#define P_DECAY 0.95
#define E_SMALL_NUM 1e-6
/*
 * Static setters, program options
 */

void AdaDelta::setAlphaInit(double val)
{
    _alpha = val;
}

void AdaDelta::setRPosAlphaScale(double val)
{
    _rPos_a_scale = val;
}

void AdaDelta::setAuto(bool val)
{
    _auto = val;
}

// Base constructor
AdaDelta::AdaDelta()
{
    frames_count = 0;
    this->alpha = ColumnVector(6);
    /*
     * Typically, a change in orientation has much more impact on
     * the total error than a change in position. While translating
     * the scan in x, y or z, the error grows linearily for all points.
     * However, when rotating about x, y or z, points that are further
     * away from the robot are moved even more dramatic, leading to
     * a higher sensibility on the error function. For this reason,
     * the alpha applied on orientation must be much smaller than the
     * alpha applied on the position.
     */
    this->a = _alpha;
    this->alpha << a
                << a
                << a
                << a * _rPos_a_scale
                << a * _rPos_a_scale
                << a * _rPos_a_scale;
    this->J = ColumnVector(6);
    this->Ga = ColumnVector(6);
    this->Xa = ColumnVector(6);
    dX = ColumnVector(6);
    // Do not buffer these from RAM ! C++ optimization is shit
    dX << 0 << 0 << 0 << 0 << 0 << 0;
    Ga << 0 << 0 << 0 << 0 << 0 << 0;
    Xa << 0 << 0 << 0 << 0 << 0 << 0;
}

// PlaneScan* constructor. Calls base constructor.
AdaDelta::AdaDelta(PlaneScan* pScan) : AdaDelta()
{
    // Setup current state vector
    X0 = ColumnVector(6);
    X0 << pScan->rPosTheta[0]  // left-handed roll (rotX)
        << pScan->rPosTheta[1] // left-handed pitch (rotY)
        << pScan->rPosTheta[2] // left-handed yaw (rotZ)
        << pScan->rPos[0] //x
        << pScan->rPos[1] //y
        << pScan->rPos[2];//z
    this->X = X0;
    this->Xmin = X;
    this->ps = pScan;
    updateJacobian(); // init once before start iteration
}

AdaDelta::AdaDelta(Matrix &x0, Matrix &j) : AdaDelta()
{
    X = x0;
    J = j;
}

// Overrides constructor of base Optimizer-class
AdaDelta::AdaDelta(PlaneScan* p, Dimensions d) : AdaDelta(p)
{
    dim = d;
}

AdaDelta::~AdaDelta() { }

void AdaDelta::setPlaneScan(PlaneScan* p) {
    ps = p;
}

void AdaDelta::setState(ColumnVector& state)
{
    // actually, i dont know if they are all necessary.
    // however, it works, so i wont change it.
    X0 = state;
    X = state;
    X1 = state;
    updateScan();
}

void AdaDelta::setJacobian(Matrix &jac)
{
    J = jac;
}

ColumnVector AdaDelta::getResult()
{
    return X;
}

void AdaDelta::relabel()
{
    // After having found the best transformation, label points again and repeat
    ps->correspondences.clear();
    ps->global_matches.clear();
    ps->global_mismatches.clear();
    ps->point_pairs.clear();
    ps->labelPoints( ps->_match_type , _quiet);
}

void AdaDelta::adaptAlpha()
{
    // E_SMALL_NUM is a magic number that works quite well, used to avoid division with 0
    alpha(1) = a * sqrt(E_SMALL_NUM + Xa(1) ) / sqrt( E_SMALL_NUM + Ga(1) );
    alpha(2) = a * sqrt(E_SMALL_NUM + Xa(2) ) / sqrt( E_SMALL_NUM + Ga(2) );
    alpha(3) = a * sqrt(E_SMALL_NUM + Xa(3) ) / sqrt( E_SMALL_NUM + Ga(3) );
    alpha(4) = a * _rPos_a_scale * sqrt(E_SMALL_NUM + Xa(4) ) / sqrt( E_SMALL_NUM + Ga(4) );
    alpha(5) = a * _rPos_a_scale * sqrt(E_SMALL_NUM + Xa(5) ) / sqrt( E_SMALL_NUM + Ga(5) );
    alpha(6) = a * _rPos_a_scale * sqrt(E_SMALL_NUM + Xa(6) ) / sqrt( E_SMALL_NUM + Ga(6) );
}

void AdaDelta::step()
{
    X1 = X;
    lastErr = getRSE();
    /*
    * Gradient descent step, using AdaDelta:
    * https://arxiv.org/pdf/1212.5701.pdf
    */
    updateJacobian(); // set J(X_k) to J(X_k+1) for next iteration
    Ga = P_DECAY*Ga + (1-P_DECAY)*SP(J, J); // acumulate averaged squared gradient
    adaptAlpha(); // use Ga and Xa to update alpha
    dX = SP(alpha, J); // compute update step
    lock();
    // cout << "a=" << alpha << endl;
    //cout << dX << endl;
    Xa = P_DECAY*Xa + (1-P_DECAY)*SP(dX, dX); // accumulate averaged squared update
    X = X - dX; // apply update
    updateScan(); // write new transformation to planescan object
}

// Update jacobian matrix from state X
void AdaDelta::updateJacobian()
{
    if (ps->isEmpty()) return;

    // Sum over correspondences
    Matrix sum = ColumnVector(6);
    sum << 0 << 0 << 0 << 0 << 0 << 0; // just to be sure.

    double phi =  X(1) ; // rot x
    double theta = X(2) ; // rot y
    double psi =  X(3); // rot z
    double tx = X(4); // x
    double ty = X(5); // y
    double tz = X(6); // z

    double sph = sin(phi);
    double st = sin(theta);
    double sps = sin(psi);
    double cph = cos(phi);
    double ct = cos(theta);
    double cps = cos(psi);

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < ps->correspondences.size(); ++i)
    {
        Correspondence cor = ps->correspondences[i];
        //cout << "Hesse gradient! " << ps->correspondences.size() << endl;
    // Auto for-each iterator is cool, but OpenMP doesnt support it.
    // for (const auto& cor : ps->correspondences)
    // {
        double x = cor.first[0];
        double y = cor.first[1];
        double z = cor.first[2];

        double nx = cor.second->n[0];
        double ny = cor.second->n[1];
        double nz = cor.second->n[2];
        // a = {ax, ay, az} is the center point of the convex hull.
        double ax = cor.second->x[0];
        double ay = cor.second->x[1];
        double az = cor.second->x[2];

        // Coordinates of the transformed point:
        double Tpx = x*ct*cps - y*ct*sps + z*st + tx;
        double Tpy = x*(cph*sps+cps*sph*st) + y*(cph*cps-sph*st*sps) - z*ct*sph + ty;
        double Tpz = x*(sph*sps-cph*cps*st) + y*(cps*sph+cph*st*sps) + z*cph*ct + tz;

        double D = nx*(Tpx - ax)  // as you can see, this represents the distance
                 + ny*(Tpy - ay)  // of the point to an ever extending, infinite
                 + nz*(Tpz - az); // plane. Called "hesse distance"

        // The first order gradient of that distance is:

        double dEdPhi =
              ny*(x*(cps*cph*st-sps*sph) + y*(-sph*cps-cph*st*sps) - z*ct*cph)
            + nz*(x*(cph*sps+sph*cps*st) + y*(cps*cph-sph*st*sps) - z*ct*sph);

        double dEdTheta =
              nx*(-x*st*cps+y*st*sps+z*ct)
            + ny*(x*(ct*cps*sph) + y*(-ct*sph*sps) + z*st*sph)
            + nz*(-x*ct*cph*cps + y*ct*cph*sps - z*st*cph);

        double dEdPsi =
              nx*(-x*sps*ct - y*cps*ct)
            + ny*(x*(cps*cph-sps*sph*st) + y*(-sps*cph-cps*sph*st))
            + nz*(x*(cps*sph+sps*cph*st) + y*(-sps*sph+cps*cph*st));

        // centroid gradient

        Matrix jacobian = ColumnVector(6);

        if (_use_p2p) {
            jacobian << dEdPhi*2*D
                    << dEdTheta*2*D
                    << dEdPsi*2*D
                    << 0//nx*2*D
                    << 0//ny*2*D
                    << 0;//nz*2*D;
        } else {

            jacobian << dEdPhi*2*D
                    << dEdTheta*2*D
                    << dEdPsi*2*D
                    << nx*2*D
                    << ny*2*D
                    << nz*2*D;
            // cout << "using plane-based correspondences! J=" << endl;
            // cout << jacobian << endl;
        }

        double w = weightLoss(D);
        #pragma omp critical
        sum += w * jacobian;
    }

    // Centroid gradient ( USE THIS LATER FOR ICP-LIKE GAUSS-NEWTON )
    if (_use_p2p) {
        double mx, my, mz; // global plane centroid (model)
        double cx, cy, cz; // local plane centroid (data)
        for (const auto& match : ps->point_pairs)
        {
            double* c = match.first;
            double* m = match.second;

            if (!m || !c) {
                cout << "invalid pointer! ";
                if (!m)
                    cout <<"Model set:"<< &m << endl;
                if (!c)
                    cout << "Data set:"<<&c << endl;
                continue;
            }

            // Get global centroid of plane model set:
            mx = m[0];
            my = m[1];
            mz = m[2];

            // Calculate local centroid point of data set:
            cx = c[0];
            cy = c[1];
            cz = c[2];
            double Tcx, Tcy, Tcz; // transformed global plane centroid (data)
            Tcx = cx*ct*cps - cy*ct*sps + cz*st + tx;
            Tcy = cx*(cph*sps+cps*sph*st) + cy*(cph*cps-sph*st*sps) - cz*ct*sph + ty;
            Tcz = cx*(sph*sps-cph*cps*st) + cy*(cps*sph+cph*st*sps) + cz*cph*ct + tz;

            double Dcx = Tcx - mx;
            double Dcy = Tcy - my;
            double Dcz = Tcz - mz;

            // First oder derivatives

            // double dDcx_dPhi = 0;
            // double dDcx_dTheta = -cx*st*cps + cy*st*sps + cz*ct;
            // double dDcx_dPsi = -cx*ct*sps - cy*ct*cps;

            // double dDcy_dPhi = cx*(-sph*sps+cps*cph*st) + cy*(-sph*cps-cph*st*sps) - cz*ct*cph;
            // double dDcy_dTheta = cx*(cps*sph*ct) + cy*(-sph*ct*sps) + cz*st*sph;
            // double dDcy_dPsi = cx*(cph*cps-sps*sph*st) + cy*(-cph*sps-sph*st*cps);

            // double dDcz_dPhi = cx*(cph*sps+sph*cps*st) + cy*(cps*cph-sph*st*sps) - cz*sph*ct;
            // double dDcz_dTheta = cx*(-cph*cps*ct) + cy*(cph*ct*sps) - cz*cph*st;
            // double dDcz_dPsi = cx*(sph*cps+cph*sps*st) + cy*(-sps*sph+cph*st*cps);

            ColumnVector dDx(6);
            ColumnVector dDy(6);
            ColumnVector dDz(6);


                dDx << 0//dDcx_dPhi
                    << 0//dDcx_dTheta
                    << 0//dDcx_dPsi
                    << 1
                    << 0
                    << 0;

                dDy << 0//dDcy_dPhi
                    << 0//dDcy_dTheta
                    << 0//dDcy_dPsi
                    << 0
                    << 1
                    << 0;

                dDz << 0//dDcz_dPhi
                    << 0//dDcz_dTheta
                    << 0//dDcz_dPsi
                    << 0
                    << 0
                    << 1;

            double D = sqrt(sqr(Dcx)+sqr(Dcy)+sqr(Dcz));
            double w = weightLoss( D );
            ColumnVector dE(6);
            dE << 2 * (Dcx*dDx + Dcy*dDy + Dcz*dDz);
            sum += w*dE;
        }
    }

    // double norm = 1.0 / ps->correspondences.size();
    // this->J = norm * sum;
    this->J = sum;
   // cout << "sum=" << J << endl;
}

void AdaDelta::iterate(int n)
{
    if (ps->isEmpty()) return;
    int n_iter = 0;
    for (int k = 0; k < _update_cor; ++k)
    {
        this->a *= ALPHA_DECAY;
        if (k != 0 && _update_cor != 1) relabel();
        for (int i = 0; i < n; ++i)
        {
            if ( !_quiet )
            {
                cout << "i = " << i << " in scan " << ps->identifier;
                cout << ", RSME = "
                << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
                << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
                << std::setw(10) << std::setprecision(7)
                <<  getRSE() << endl;
            }
            step();

        }
        if (n_iter == 1) break;
        else n_iter = 0;
    }
}

void AdaDelta::operator()(int n)
{
    this->iterate(n);
}

// Iterates until convergence. Use with care!
void AdaDelta::iterate(double eps)
{
    if (ps->isEmpty()) return;
    int n_iter = 0;
    for (int k = 0; k < _update_cor ;++k)
    {
        this->a *= ALPHA_DECAY;
        if (k != 0 && _update_cor != 1) relabel();
        do
        {
            if ( !_quiet )
            {
                cout << "Scan " << ps->identifier;
                cout << ", RSME = "
                << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
                << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
                << std::setw(10) << std::setprecision(7)
                << getRSE() << endl;
            }
            step();

            n_iter++;
        } while ( !stop_condition(X, eps) );
        // After having found the best transformation, label points again and repeat
        if (n_iter == 1) break;
        else n_iter = 0;
    }
}

void AdaDelta::reset()
{
    X = X0;
    Xmin = X;
    frames_count = 0;
    updateScan();
    Ga << 0 << 0 << 0 << 0 << 0 << 0;
    Xa << 0 << 0 << 0 << 0 << 0 << 0;
}

// Iterates until convergence. Use with care!
void AdaDelta::operator()(double eps)
{
    this->iterate(eps);
}

// At maximum n iterations, quickstop if difference is too small
void AdaDelta::iterate(int n, double eps)
{
    if (ps->isEmpty()) return;
    int n_iter = 0;
    for (int k = 0; k < _update_cor; ++k)
    {
        this->a *= ALPHA_DECAY;
        if (k != 0 && _update_cor != 1) relabel();
        for (int i = 0; i < n; ++i)
        {
            if ( !_quiet )
            {
                cout << "i = " << i << " in scan " << ps->identifier;
                cout << ", RSME = "
                << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
                << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
                << std::setw(10) << std::setprecision(7)
                <<  getRSE() << endl;
            }
            step();

            n_iter++;
            //If eps criterium is met, stop iteration.
            if ( stop_condition(X, eps) ) break;
        }
        if (n_iter == 1) break;
        else n_iter = 0;
    }
}

void AdaDelta::iterate()
{
    if ( _max_iter != -1 && _eps_convergence != -1 )
        this->iterate(_max_iter, _eps_convergence);
    else if ( _eps_convergence == -1 && _max_iter != -1)
        this->iterate(_max_iter);
    else if ( _max_iter == -1 && _eps_convergence != -1)
        this->iterate(_eps_convergence);
    else
        this->iterate(1000);
    // this->ps->addFrame(Scan::ICPINACTIVE);
    size_t scans_size = PlaneScan::allPlaneScans.size();
    int found = 0;
    for (size_t i = 0; i < scans_size; ++i) {
        PlaneScan* scan = PlaneScan::allPlaneScans[i];
        if(scan == ps) {
            found = i;
            scan->addFrame( Scan::ICP );
        } else {
            if (found == 0) {
                scan->addFrame( Scan::ICPINACTIVE);
            } else {
                scan->addFrame( Scan::INVALID);
            }
        }
    }
}

// Checks program options and iterates in the specified way
void AdaDelta::operator()(void)
{
    // FIXME: Add better auto-mode
    // if (_auto) { this->iterateAuto(); }
    iterate();
    relabel();
}

// Finds the optimal alpha to initialize AdaDelta
// Experimental, usually this is very agressive.
void AdaDelta::iterateAuto()
{
    a = 1;
    double initErr, err;
    for (int i = 0; i < 10; ++i)
    {
        initErr = getRSE();
        for (int k = 0; k < 10; ++k) step();
        err = getRSE();
        reset();
        if ( err < initErr ) {
            a *= 0.1;
            if (!_quiet) cout << "Optimal alpha = " << a << endl;
            iterate(); // start the actual iteration, using specified criteria
            return;
        } else {
            a *= 0.5;
        }
    }
}

