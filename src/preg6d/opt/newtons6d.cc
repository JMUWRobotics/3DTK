#include "opt/newtons6d.h"
#include "opt/adadelta6d.h"

// All Newtons constructors call this:
Newtons::Newtons()
{
    frames_count = 0;
    X = ColumnVector(6);
    J = ColumnVector(6);
    H = Matrix(6,6);
    J_c = ColumnVector(6);
    H_c = Matrix(6,6);
    JJ = ColumnVector(6);
    HH = Matrix(6,6);
    err_min = __DBL_MAX__;
}

Newtons::Newtons(PlaneScan* pScan) : Newtons()
{
    this->ps = pScan;
    X << pScan->rPosTheta[0]  // left-handed roll (rotX)
      << pScan->rPosTheta[1] // left-handed pitch (rotY)
      << pScan->rPosTheta[2] // left-handed yaw (rotZ)
      << pScan->rPos[0] //x
      << pScan->rPos[1] //y
      << pScan->rPos[2];//z
    updateScan();
    updateHesseGradient();
    updateCentroidGradient();
    overlapGradients();
}

// OVERRIDE from Optimizer base class
Newtons::Newtons(PlaneScan* pScan, Dimensions d)
: Newtons(pScan) // call PlaneScan constructor, which calls base
{
    dim = d;
}

// OVERRIDE from Optimizer base class
void Newtons::operator()(void)
{
    iterate();
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
    rematch();
}

void Newtons::rematch()
{
    // Update Correspondences
    ps->correspondences.clear();
    ps->global_matches.clear();
    ps->global_mismatches.clear();
    ps->point_pairs.clear();
    ps->labelPoints( ps->_match_type, true );
}

void Newtons::iterate()
{
    if (ps->isEmpty()) return;
    int iter = 0;
    int break_iter = 0;
    size_t lastNrCor, nrCor_max = 0;
    ColumnVector jacobian(6);
    Matrix hessian(6,6);
    // use c and tau for Armijo conditions.
    // See:
    // https://en.wikipedia.org/wiki/Backtracking_line_search#Algorithm
    double c = 0.5, tau = 0.5;
    do {
        jacobian = JJ;
        hessian = HH;

        // update matches
        rematch();
        lastNrCor = ps->correspondences.size();

        // If 3 subsequent iterations ONLY increased error, break.
        // if (ps->calcErr() > lastErr) {
        //     if (break_iter++ >= 3) {
        //         X = X_min;
        //         updateScan();
        //         break;
        //     } else {
        //         break_iter = 0;
        //     }
        // }

        lastErr = ps->calcErr();
        if (lastNrCor > nrCor_max) {
            X_min = X;
            nrCor_max = lastNrCor;
            if (lastErr < err_min) {
                X_min = X;
                err_min = lastErr;
            }
        }

        // Print current error
        cout << "i = " << iter++ << " in scan " << ps->identifier << " with " << ps->correspondences.size() << " pts";
            cout << ", RMSE = "
            << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
            << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
            << std::setw(10) << std::setprecision(7)
            <<  lastErr  << endl;

        // Check if H is positive definite
        bool is_positive_definite = true;
        {
            SymmetricMatrix H_s(6);
            H_s << hessian;
            int nr = H_s.Nrows();
            LowerTriangularMatrix T(nr);
            Real* s = H_s.Store(); Real* t = T.Store(); Real* ti = t;
            for (int i=0; i<nr; i++)
            {
                Real* tj = t; Real sum; int k;
                for (int j=0; j<i; j++)
                {
                    Real* tk = ti; sum = 0.0; k = j;
                    while (k--) { sum += *tj++ * *tk++; }
                    *tk = (*s++ - sum) / *tj++;
                }
                sum = 0.0; k = i;
                while (k--) { sum += sqr(*ti++); }
                Real d = *s++ - sum;
                if (d<=0.0) {
                    is_positive_definite = false;
                    break;
                };
                *ti++ = sqrt(d);
            }
        }
        double a = 1;

        // Compute step direction dX
        LinearEquationSolver He = hessian;
        dX = He.i() * jacobian;
        // Put constraints on the step
        lock();

        // If hessian is not positive definite, we use AdaDelta as fallback
        // if (!is_positive_definite) {

        //     bool exec_fallback = false;

        //     // Perform line search in the specified direction (Armijo conditions)
        //     int saved_anim = _anim;
        //     setAnim(-2);
        //     double m = DotProduct(jacobian, dX);
        //     double t = c*m;
        //     if (fabs(t) < 0.0000001 ) {
        //         updateScan();
        //         rematch();
        //         setAnim(saved_anim);
        //         break;
        //     }
        //     ColumnVector saved_state = X;
        //     X = X - a*dX;
        //     updateScan();
        //     double tmpErr = ps->calcErr();
        //     X = saved_state;
        //     updateScan();
        //     while ( lastErr - tmpErr < a*t ) {
        //         a *= tau;
        //         X = X - a*dX;
        //         updateScan();
        //         tmpErr = ps->calcErr();
        //         X = saved_state;
        //         updateScan();
        //         // Line search was unsuccessfull at this point.
        //         if (a < 1e-12) {
        //             exec_fallback = true;
        //             break;
        //         }
        //     }
        //     // Reset to the state before line search
        //     X = saved_state;
        //     updateScan();
        //     setAnim(saved_anim);
        //     exec_fallback= false;
        //     // Use gradient descent based on Jacobian only
        //     if (exec_fallback) {

        //         cout << "Hessian not positive definite! Reset and Skip." << endl;
        //         X = X_min;
        //         updateScan();
        //         rematch();
        //         int saved_anim = _anim;
        //         AdaDelta::setAnim(-2);
        //         AdaDelta ada_iter(ps,dim);
        //         ada_iter();
        //         AdaDelta::setAnim(saved_anim);
        //         lastErr = ps->calcErr();
        //         X = ada_iter.getResult();
        //         X_min = X;
        //         updateScan();
        //         rematch();
        //         break;
        //     } else {
        //         cout << "Using armijos method with a=" << a << endl;
        //     }
        // }

        // ColumnVector saved_state = X;

        // Apply newton step
        X = X - a*dX;

        // Update Scan
        updateScan();
        // if ( ps->calcErr() > lastErr ) {
        //     X=saved_state;
        //     updateScan();
        //     int saved_anim = _anim;
        //     AdaDelta::setAnim(-2);
        //     AdaDelta ada_iter(ps,dim);
        //     ada_iter();
        //     AdaDelta::setAnim(saved_anim);
        //     lastErr = ps->calcErr();
        //     X = ada_iter.getResult();
        //     X_min = X;
        //     updateScan();
        //     rematch();
        //     break;
        // }

        // Update Jacobian and Hessian
        updateHesseGradient();
        updateCentroidGradient();
        overlapGradients();

    } while (
        !stop_condition(X, _eps_convergence)
        && iter < _max_iter
    );

    // Set state to minimal error
   // X = X_min;
    updateScan();
}

void Newtons::updateHesseGradient()
{
    if (ps->isEmpty()) return;

    // Update Jacobian (Vector6):

    Matrix sum = ColumnVector(6);
    sum = 0;
    J = 0;

    Matrix sum2 = Matrix(6,6);
    sum2 = 0;
    H = 0;

    double inv_kernel = 1.0 / _eps_kernel;

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (size_t i = 0; i < ps->correspondences.size(); ++i)
    {
        Correspondence cor = ps->correspondences[i];

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

        // Coordinates of the transformed point:
        double Tpx = x*ct*cps - y*ct*sps + z*st + tx;
        double Tpy = x*(cph*sps+cps*sph*st) + y*(cph*cps-sph*st*sps) - z*ct*sph + ty;
        double Tpz = x*(sph*sps-cph*cps*st) + y*(cps*sph+cph*st*sps) + z*cph*ct + tz;

        double D = nx*(Tpx - ax)  // as you can see, this represents the distance
                 + ny*(Tpy - ay)  // of the point to an ever extending, infinite
                 + nz*(Tpz - az); // plane. Called "hesse distance"

        // The first order gradient of that distance is:

        double dDdPhi =
              ny*(x*(cps*cph*st-sps*sph) + y*(-sph*cps-cph*st*sps) - z*ct*cph)
            + nz*(x*(cph*sps+sph*cps*st) + y*(cps*cph-sph*st*sps) - z*ct*sph);

        double dDdTheta =
              nx*(-x*st*cps+y*st*sps+z*ct)
            + ny*(x*(ct*cps*sph) + y*(-ct*sph*sps) + z*st*sph)
            + nz*(-x*ct*cph*cps + y*ct*cph*sps - z*st*cph);

        double dDdPsi =
              nx*(-x*sps*ct - y*cps*ct)
            + ny*(x*(cps*cph-sps*sph*st) + y*(-sps*cph-cps*sph*st))
            + nz*(x*(cps*sph+sps*cph*st) + y*(-sps*sph+cps*cph*st));

        Matrix jacobian = ColumnVector(6);
        if (_use_p2p) {
            jacobian << dDdPhi*2*D
                    << dDdTheta*2*D
                    << dDdPsi*2*D
                    << 0//nx*2*D
                    << 0//ny*2*D
                    << 0;//nz*2*D;
        } else {
            jacobian << dDdPhi*2*D
                    << dDdTheta*2*D
                    << dDdPsi*2*D
                    << nx*2*D
                    << ny*2*D
                    << nz*2*D;
        }

        // The second order gradient of that distance is:

        // Update Hessian
        double H11, H12, H13, H14, H15, H16,
               H21, H22, H23, H24, H25, H26,
               H31, H32, H33, H34, H35, H36,
               H41, H42, H43, H44, H45, H46,
               H51, H52, H53, H54, H55, H56,
               H61, H62, H63, H64, H65, H66;

        // d2D/(dPhi*dPhi) =
        H11 = ny*(x*(-cps*sph*st-sps*cph) + y*(-cph*cps+sph*st*sps) + z*ct*sph)
            + nz*(x*(-sph*sps+cph*cps*st) + y*(-cps*sph-cph*st*sps) - z*ct*cph);

        // d2D/(dPhi*dTheta) =
        H12 = ny*(x*(cps*cph*ct) + y*(-cph*ct*sps) + z*st*cph)
            + nz*(x*(sph*cps*ct) + y*(-sph*ct*sps) + z*st*sph);

        // d2D/(dPhi*dPsi) =
        H13 = ny*(x*(-sps*cph*st-cps*sph) + y*(sph*sps-cph*st*cps))
            + nz*(x*(cph*cps-sph*sps*st) + y*(-sps*cph-sph*st*cps));

        // d2D/(dPhi*dtx)
        H14 = 0;

        // d2D/(dPhi*dty)
        H15 = 0;

        // d2D/(dPhi*dtz)
        H16 = 0;

        // d2D/(dTheta*dPhi)
        H21 = H12;

        // d2D/(dTheta*dTheta)
        H22 = nx*(-x*ct*cps+y*ct*sps-z*st)
            + ny*(x*(-st*cps*sph) + y*(st*sph*sps) + z*ct*sph)
            + nz*(x*st*cph*cps - y*st*cph*sps - z*ct*cph);

        // dD/(dTheta*dPsi)
        H23 = nx*(x*st*sps+y*st*cps)
            + ny*(x*(-ct*sps*sph) + y*(-ct*sph*cps))
            + nz*(x*ct*cph*sps + y*ct*cph*cps);

        // dD/(dTheta*dtx)
        H24 = 0;

        // dD/(dTheta*dty)
        H25 = 0;

        // dD/(dTheta*dtz)
        H26 = 0;

        // dD/(dPsi*dPhi)
        H31 = H13;

        // dD/(dPsi*dTheta)
        H32 = H23;

        // dD/(dPsi*dPsi)
        H33 = nx*(-x*cps*ct + y*sps*ct)
            + ny*(x*(-sps*cph-cps*sph*st) + y*(-cps*cph+sps*sph*st))
            + nz*(x*(-sps*sph+cps*cph*st) + y*(-cps*sph-sps*cph*st));

        // dD/(dPsi*dtx)
        H34 = 0;

        // dD/(dPsi*dty)
        H35 = 0;

        // dD/(dPsi*dtz)
        H36 = 0;

        // The rest of the matrix is symetrical (Schwarzs Law)
        H41 = H14;
        H42 = H24;
        H43 = H34;
        H44 = 0; // dD/(dtx*dtx)
        H45 = 0; // dD/(dtx*dty)
        H46 = 0; // dD/(dtx*dtz)

        H51 = H15;
        H52 = H25;
        H53 = H35;
        H54 = H45;
        H55 = 0; // dD/(dty*dty)
        H56 = 0; // dD/(dty*dtz)

        H61 = H16;
        H62 = H26;
        H63 = H36;
        H64 = H46;
        H65 = H56;
        H66 = 0; // dD/(dtz*dtz)

        Matrix D2 = Matrix(6,6);
        D2 << H11 << H12 << H13 << H14 << H15 << H16
           << H21 << H22 << H23 << H24 << H25 << H26
           << H31 << H32 << H33 << H34 << H35 << H36
           << H41 << H42 << H43 << H44 << H45 << H46
           << H51 << H52 << H53 << H54 << H55 << H56
           << H61 << H62 << H63 << H64 << H65 << H66;

        // Assert symetric matrix
        Matrix D2t(6,6);
        D2t << D2.t();
        assertm( (D2.IsEqual(D2t)) , "Second order of distance not symmetric." );

        // pure first order distance derivative
        ColumnVector D1(6);

        if (_use_p2p) {
            D1 << dDdPhi
            << dDdTheta
            << dDdPsi
            << 0//nx
            << 0//ny
            << 0;//nz;
        } else {
            D1 << dDdPhi
            << dDdTheta
            << dDdPsi
            << nx
            << ny
            << nz;
        }

        // Transpose first order distance derivative
        RowVector D1_t(6);
        D1_t << D1.t();

        // Setup hessian of error function
        Matrix hessian = Matrix(6,6);
        hessian << 2 * ( (D1*D1_t) + (D*D2) );

        // Tukey loss
        double w = fabs(D) <= _eps_kernel ? sqr(1-sqr(D*inv_kernel)) : 0;

#pragma omp critical
{
        sum += w*jacobian;
        sum2 += w*hessian;
}

    }
    J = sum;
    H = sum2;
}

void Newtons::updateCentroidGradient()
{
    //cout << "update centroid" << endl;
    if (ps->isEmpty() || !_use_p2p) return;

    Matrix sum = ColumnVector(6);
    sum = 0;
    J_c = 0;

    Matrix sum2 = Matrix(6,6);
    sum2 = 0;
    H_c = 0;

    double inv_kernel = 1.0 / _eps_kernel;

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

    double mx, my, mz; // global point (model)
    double cx, cy, cz; // local point (data)
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
        dDx << 0//dDcx_dPhi
            << 0//dDcx_dTheta
            << 0//dDcx_dPsi
            << 1
            << 0
            << 0;

        ColumnVector dDy(6);
        dDy << 0//dDcy_dPhi
            << 0//dDcy_dTheta
            << 0//dDcy_dPsi
            << 0
            << 1
            << 0;

        ColumnVector dDz(6);
        dDz << 0//dDcz_dPhi
            << 0//dDcz_dTheta
            << 0//dDcz_dPsi
            << 0
            << 0
            << 1;

        ColumnVector dDx_trans(6);
        ColumnVector dDy_trans(6);
        ColumnVector dDz_trans(6);

        dDx_trans << 0 << 0 << 0 << dDx(4) << dDx(5) << dDx(6);
        dDy_trans << 0 << 0 << 0 << dDy(4) << dDy(5) << dDy(6);
        dDz_trans << 0 << 0 << 0 << dDz(4) << dDz(5) << dDz(6);

        // Second order derivatives

        // dDcx_dPhi = 0;

        //dDcx_dTheta = -cx*st*cps + cy*st*sps + cz*ct
        double d2Dx22 = -cx*ct*cps + cy*ct*sps - cz*st;
        double d2Dx23 = cx*st*sps + cy*st*cps;
        //dDcx_dPsi = -cx*ct*sps - cy*ct*cps;
        double d2Dx32 = d2Dx23;
        double d2Dx33 = -cx*ct*cps + cy*ct*sps;

        // dDcx_dT = const. -> d2Dcx_dT^2 = 0 for all x,y,z

        Matrix d2Dx(6,6);
        d2Dx << 0 << 0      << 0      << 0 << 0 << 0
             << 0 << d2Dx22 << d2Dx23 << 0 << 0 << 0
             << 0 << d2Dx32 << d2Dx33 << 0 << 0 << 0
             << 0 <<   0    <<   0    << 0 << 0 << 0
             << 0 <<   0    <<   0    << 0 << 0 << 0
             << 0 <<   0    <<   0    << 0 << 0 << 0;

        // dDcy_dPhi = cx*(-sph*sps+cps*cph*st) + cy*(-sph*cps-cph*st*sps) - cz*ct*cph;
        double d2Dy11 = cx*(-cph*sps-cps*sph*st) + cy*(-cph*cps+sph*st*sps) + cz*ct*sph;
        double d2Dy12 = cx*(cps*cph*ct) + cy*(-cph*ct*sps) + cz*st*cph;
        double d2Dy13 = cx*(-sph*cps-sps*cph*st) + cy*(sph*sps-cph*st*cps);

        // dDcy_dTheta = cx*(cps*sph*ct) + cy*(-sph*ct*sps) + cz*st*sph;
        double d2Dy21 = d2Dy12;
        double d2Dy22 = cx*(-cps*sph*st) + cy*(sph*st*sps) + cz*ct*sph;
        double d2Dy23 = cx*(-sps*sph*ct) + cy*(-sph*ct*cps);

        // dDcy_dPsi = cx*(cph*cps-sps*sph*st) + cy*(-cph*sps-sph*st*cps);
        double d2Dy31 = d2Dy13;
        double d2Dy32 = d2Dy23;
        double d2Dy33 = cx*(-cph*sps-cps*sph*st) + cy*(-cph*cps+sph*st*sps);

        // d2Dcy_dT^2 = 0

        Matrix d2Dy(6, 6);
        d2Dy << d2Dy11 << d2Dy12 << d2Dy13 << 0 << 0 << 0
             << d2Dy21 << d2Dy22 << d2Dy23 << 0 << 0 << 0
             << d2Dy31 << d2Dy32 << d2Dy33 << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0;

        //dDcz_dPhi = cx*(cph*sps+sph*cps*st) + cy*(cps*cph-sph*st*sps) - cz*sph*ct;
        double d2Dz11 = cx*(-sph*sps+cph*cps*st) + cy*(-cps*sph-cph*st*sps) - cz*cph*ct;
        double d2Dz12 = cx*(sph*cps*ct) + cy*(-sph*ct*sps) + cz*sph*st;
        double d2Dz13 = cx*(cph*cps-sph*sps*st) + cy*(-sps*cph-sph*st*cps);

        //dDcz_dTheta = cx*(-cph*cps*ct) + cy*(cph*ct*sps) - cz*cph*st;
        double d2Dz21 = d2Dz12;
        double d2Dz22 = cx*(cph*cps*st) + cy*(-cph*st*sps) - cz*cph*ct;
        double d2Dz23 = cx*(cph*sps*ct) + cy*(cph*ct*cps);

        //dDcz_dPsi = cx*(sph*cps+cph*sps*st) + cy*(-sps*sph+cph*st*cps);
        double d2Dz31 = d2Dz13;
        double d2Dz32 = d2Dz23;
        double d2Dz33 = cx*(-sph*sps+cph*cps*st) + cy*(-cps*sph-cph*st*sps);

        Matrix d2Dz(6, 6);
        d2Dz << d2Dz11 << d2Dz12 << d2Dz13 << 0 << 0 << 0
             << d2Dz21 << d2Dz22 << d2Dz23 << 0 << 0 << 0
             << d2Dz31 << d2Dz32 << d2Dz33 << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0
             <<   0    <<    0   <<    0   << 0 << 0 << 0;


         // Jacobian of centroid error function:
        ColumnVector dE(6);
        dE << 2 * (Dcx*dDx + Dcy*dDy + Dcz*dDz);

        // Tukey loss
        double D = sqrt(sqr(Dcx)+sqr(Dcy)+sqr(Dcz));
        double w = fabs(D) <= _eps_kernel ? sqr(1-sqr(D*inv_kernel)) : 0;

        ColumnVector dE_trans(6);
        dE_trans << 0 << 0 << 0 << dE(4) << dE(5) << dE(6);
        sum += w*dE_trans;

        // Hessian of centroid error function:
        Matrix d2E(6, 6);
        // d2E << 2 * ( (dDx*dDx.t()) + Dcx*d2Dx
        //            + (dDy*dDy.t()) + Dcy*d2Dy
        //            + (dDz*dDz.t()) + Dcz*d2Dz );
        d2E << 2 * ( (dDx_trans*dDx_trans.t()) + Dcx*d2Dx
                   + (dDy_trans*dDy_trans.t()) + Dcy*d2Dy
                   + (dDz_trans*dDz_trans.t()) + Dcz*d2Dz );
        sum2 += d2E;
    }
    J_c = sum;
    H_c = sum2;
}

void Newtons::overlapGradients()
{
    if (_use_p2p) {
        JJ = J +  J_c;
        HH = H +  H_c;
    } else {
        JJ = J;
        HH = H;
    }
}




