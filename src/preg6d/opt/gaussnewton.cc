#include "opt/gaussnewton.h"

// TODO:
/*
 * Optimize in the following manner:
 *
 * X_k+1 = X_k - (J^t * W * J)^-1 * J^t * W * r
 * J should be a vector of N Jacobians (ColumnVector)
 * J^t is then a RowVector of the N Jacobians
 * W is the weight matrix containing N weights
 * r is a ColumnVector of N residuals (point2plane err)
 *
 * weights are:
 * w_i = 1 / r_i * p'(r_i)
 * where p' is the first order derivative of loss p
 *
 * For loss function table, see:
 * P. Babin et al.
 * https://arxiv.org/pdf/1810.01474.pdf
 *
 */

GaussNewton::GaussNewton() {
    J = 0;
    W = 0;
    R = 0;
    dX = ColumnVector(6);
    dX = 0;
    X = ColumnVector(6);
    X = 0;
}

GaussNewton::GaussNewton(PlaneScan* p) : GaussNewton() {
    this->ps = p;
    X << ps->rPosTheta[0]  // left-handed roll (rotX)
      << ps->rPosTheta[1] // left-handed pitch (rotY)
      << ps->rPosTheta[2] // left-handed yaw (rotZ)
      << ps->rPos[0] //x
      << ps->rPos[1] //y
      << ps->rPos[2];//z
    updateScan();
    updateGradient();
    this->dim = Dimensions::ALL;
}

GaussNewton::GaussNewton(PlaneScan *p, Dimensions d): GaussNewton(p)
{
    this->dim = d;
}

void GaussNewton::rematch()
{
    // Update Correspondences
    ps->correspondences.clear();
    ps->global_matches.clear();
    ps->global_mismatches.clear();
    ps->point_pairs.clear();
    ps->labelPoints( ps->_match_type, true );
}

void GaussNewton::iterate()
{
    if (ps->isEmpty()) return;
    int iter = 0;
    do {
        rematch();
        lastErr = ps->calcErr();
        // Print current error
        cout << "i = " << iter++ << " in scan " << ps->identifier << " with " << ps->correspondences.size() << " pts";
            cout << ", RMSE = "
            << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
            << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
            << std::setw(10) << std::setprecision(7)
            <<  lastErr  << endl;

        updateGradient();


        try {
            //cout << "Without residual we get: " << endl;
            LinearEquationSolver H = (J.t() * J);
            RowVector S = H.i() * J.t();
            //cout << S.t() << endl;
            dX = S.t() * lastErr;
        } catch(NEWMAT::SingularException &e) {
            break;
        }
        //cout << "dX= " << endl << dX << endl;
        lock();

        X = X - dX;

        updateScan();

    } while ( !stop_condition(X)
        && iter < _max_iter);
}

void GaussNewton::operator()(void)
{
    iterate();
}

void GaussNewton::updateGradient()
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
        jacobian << dEdPhi*2*D
                 << dEdTheta*2*D
                 << dEdPsi*2*D
                 << nx*2*D
                 << ny*2*D
                 << nz*2*D;

        // double w = fabs(D) <= _eps_kernel ? sqr(1-sqr(D*inv_kernel)) : 0; // Tukey loss
        //double w = 1.0 / (1 + sqr(D * inv_kernel)) ; // Cauchy loss
        double w = weightLoss(D);
        #pragma omp critical
        sum += w * jacobian;
    }

    // Centroid gradient ( USE THIS LATER FOR ICP-LIKE GAUSS-NEWTON )

    // double mx, my, mz; // global plane centroid (model)
    // double cx, cy, cz; // local plane centroid (data)
    // for (const auto& match : ps->point_pairs)
    // {
    //     double* c = match.first;
    //     double* m = match.second;

    //     if (!m || !c) {
    //         cout << "invalid pointer! ";
    //         if (!m)
    //             cout <<"Model set:"<< &m << endl;
    //         if (!c)
    //             cout << "Data set:"<<&c << endl;
    //         continue;
    //     }

    //     // Get global centroid of plane model set:
    //     mx = m[0];
    //     my = m[1];
    //     mz = m[2];

    //     // Calculate local centroid point of data set:
    //     cx = c[0];
    //     cy = c[1];
    //     cz = c[2];
    //     double Tcx, Tcy, Tcz; // transformed global plane centroid (data)
    //     Tcx = cx*ct*cps - cy*ct*sps + cz*st + tx;
    //     Tcy = cx*(cph*sps+cps*sph*st) + cy*(cph*cps-sph*st*sps) - cz*ct*sph + ty;
    //     Tcz = cx*(sph*sps-cph*cps*st) + cy*(cps*sph+cph*st*sps) + cz*cph*ct + tz;

    //     double Dcx = Tcx - mx;
    //     double Dcy = Tcy - my;
    //     double Dcz = Tcz - mz;

    //     // First oder derivatives

    //     double dDcx_dPhi = 0;
    //     double dDcx_dTheta = -cx*st*cps + cy*st*sps + cz*ct;
    //     double dDcx_dPsi = -cx*ct*sps - cy*ct*cps;

    //     double dDcy_dPhi = cx*(-sph*sps+cps*cph*st) + cy*(-sph*cps-cph*st*sps) - cz*ct*cph;
    //     double dDcy_dTheta = cx*(cps*sph*ct) + cy*(-sph*ct*sps) + cz*st*sph;
    //     double dDcy_dPsi = cx*(cph*cps-sps*sph*st) + cy*(-cph*sps-sph*st*cps);

    //     double dDcz_dPhi = cx*(cph*sps+sph*cps*st) + cy*(cps*cph-sph*st*sps) - cz*sph*ct;
    //     double dDcz_dTheta = cx*(-cph*cps*ct) + cy*(cph*ct*sps) - cz*cph*st;
    //     double dDcz_dPsi = cx*(sph*cps+cph*sps*st) + cy*(-sps*sph+cph*st*cps);

    //     ColumnVector dDx(6);
    //     dDx << 0//dDcx_dPhi
    //         << 0//dDcx_dTheta
    //         << 0//dDcx_dPsi
    //         << 1
    //         << 0
    //         << 0;

    //     ColumnVector dDy(6);
    //     dDy << 0//dDcy_dPhi
    //         << 0//dDcy_dTheta
    //         << 0//dDcy_dPsi
    //         << 0
    //         << 1
    //         << 0;

    //     ColumnVector dDz(6);
    //     dDz << 0//dDcz_dPhi
    //         << 0//dDcz_dTheta
    //         << 0//dDcz_dPsi
    //         << 0
    //         << 0
    //         << 1;

    //     double D = sqrt( sqr(Dcx) + sqr(Dcy) + sqr(Dcz));
    //     double w = weightLoss(D);
    //     // Jacobian of centroid error function:
    //     ColumnVector dE(6);
    //     dE << 2 * (Dcx*dDx + Dcy*dDy + Dcz*dDz);

    //     sum += w*dE;
    // }

    // double norm = 1.0 / ps->correspondences.size();
    // this->J = norm * sum;
    this->J = sum;
}

