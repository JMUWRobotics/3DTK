#include "opt/svd.h"

PlaneSVD::~PlaneSVD() {}

PlaneSVD::PlaneSVD()
{
    frames_count = 0;
    X = ColumnVector(6);
    dX = ColumnVector(6);
}

PlaneSVD::PlaneSVD( PlaneScan* pScan ) : PlaneSVD()
{
    this->ps = pScan;
    X << pScan->rPosTheta[0]  // left-handed roll (rotX)
      << pScan->rPosTheta[1] // left-handed pitch (rotY)
      << pScan->rPosTheta[2] // left-handed yaw (rotZ)
      << pScan->rPos[0] //x
      << pScan->rPos[1] //y
      << pScan->rPos[2];//z
    updateScan();
}

PlaneSVD::PlaneSVD( PlaneScan* pScan, Dimensions d ) : PlaneSVD( pScan )
{
    dim = d;
}

void PlaneSVD::operator()(void)
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
}

void PlaneSVD::iterate()
{
    if (ps->isEmpty()) return;
    int n_iter = 0;
    do {
        rematch();
        step( n_iter );
    } while( !stop_condition(X, _eps_convergence)
        &&  (++n_iter <= _max_iter) );

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

void PlaneSVD::rematch()
{
  // Update Correspondences
    ps->correspondences.clear();
    ps->global_matches.clear();
    ps->global_mismatches.clear();
    ps->point_pairs.clear();
    ps->labelPoints( ps->_match_type, true );
}

// TODO: I think there is a bug here.
// We need to look at src/slam6D/searchTree.cc::109
// Transform the data points to be in the model reference frame first.
void PlaneSVD::calc_centroids(double* cm, double* cd, double* dtr)
{
    for(int j = 0; j < 3; ++j) {
        cd[j] = 0;
        cm[j] = 0;
        if (dtr != 0)
            dtr[j] = 0;
    }

    if (ps->isEmpty()) return;

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

        double x = cor.first[0];
        double y = cor.first[1];
        double z = cor.first[2];

        double nx = cor.second->n[0];
        double ny = cor.second->n[1];
        double nz = cor.second->n[2];
        // double norm[3] = {nx, ny, nz};
        // transform3normal(ps->transMat, norm);
        // nx = norm[0]; ny = norm[1]; nz = norm[2];


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

        // Calculate the point projection on the plane, which is used as correspondence.
        double Ppx = Tpx - D * nx;
        double Ppy = Tpy - D * ny;
        double Ppz = Tpz - D * nz;

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            if (dtr != 0) {
                dtr[0] -= 2*D*nx;
                dtr[1] -= 2*D*ny;
                dtr[2] -= 2*D*nz;
            }

            // Model centroid (plane projections)
            cm[0] += Ppx;
            cm[1] += Ppy;
            cm[2] += Ppz;

            // Data to align (world frame)
            cd[0] += Tpx;
            cd[1] += Tpy;
            cd[2] += Tpz;
        }
    }

    for (int i = 0; i < 3; ++i) {
        cm[i] /= ps->correspondences.size();
        cd[i] /= ps->correspondences.size();
        if (dtr != 0)
            dtr[i] /= ps->correspondences.size();
    }

}

void PlaneSVD::step(int it)
{
    if (ps->isEmpty()) return;
    double error = 0;
    double sum = 0.0;

    // Calculate centroid of correspondences
    double cm[3], cd[3];
    double dtr[3]; // this is the needed translation
    calc_centroids(cm, cd, dtr);

    // cout << "Centroid cm: " << cm[0] << " " << cm[1] << " " << cm[2] << endl;
    // cout << "Centroid cd: " << cd[0] << " " << cd[1] << " " << cd[2] << endl;

    // Allocate memory for centered correspondences
    double **m = new double*[ps->correspondences.size()];
    double **d = new double*[ps->correspondences.size()];

    // Get current state (transformation)
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

    // Fill the centered correspondences
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (unsigned int i = 0; i < ps->correspondences.size(); ++i) {
        Correspondence cor = ps->correspondences[i];

        double x = cor.first[0];
        double y = cor.first[1];
        double z = cor.first[2];

        double nx = cor.second->n[0];
        double ny = cor.second->n[1];
        double nz = cor.second->n[2];
        // double norm[3] = {nx, ny, nz};
        // transform3normal(ps->transMat, norm);
        // nx = norm[0]; ny = norm[1]; nz = norm[2];

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

        // Calculate the point projection on the plane, which is used as correspondence.
        double Ppx = Tpx - D * nx;
        double Ppy = Tpy - D * ny;
        double Ppz = Tpz - D * nz;

        m[i] = new double[3];
        d[i] = new double[3];
        m[i][0] = Ppx - cm[0];
        m[i][1] = Ppy - cm[1];
        m[i][2] = Ppz - cm[2];
        d[i][0] = Tpx - cd[0];
        d[i][1] = Tpy - cd[1];
        d[i][2] = Tpz - cd[2];

        #pragma omp critical
        sum += sqr(Ppx - Tpx) + sqr(Ppy - Tpy) + sqr(Ppz - Tpz);
    }
    error = sqrt(sum / (double) ps->correspondences.size());
    if ( !_quiet )
    {
        if (it != -1) cout << it << ": ";
        cout << "SVD point-2-plane RMSE in scan" << ps->identifier;
        cout << " = "
        << std::resetiosflags(std::ios::adjustfield) << std::setiosflags(std::ios::internal)
        << std::resetiosflags(std::ios::floatfield) << std::setiosflags(std::ios::fixed)
        << std::setw(10) << std::setprecision(7)
        << error << endl;
    }

    // The following is a rip-off from src/slam6D/icp6Dsvd.cc:

    // Fill H matrix
    Matrix H(3,3), R(3,3);
    for(int j = 0; j < 3; j++){
        for(int k = 0; k < 3; k++){
            H(j+1, k+1) = 0.0;
        }
    }

    for(unsigned int i = 0; i < ps->correspondences.size(); i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                H(j+1, k+1) += d[i][j]*m[i][k];
            }
        }
    }
    Matrix U(3,3);
    DiagonalMatrix Lamda(3);
    Matrix V(3,3);

    // Make SVD
    SVD(H, Lamda, U, V);
    // Get rotation
    R = V*(U.t());

    // Check if R turned out to be a reflection (likely because of large errors or coplanar input)
    // (cf. Arun, Huang, and Blostein: "Least-Squares Fitting of Two 3-D Point Sets")
    if(R.Determinant() < 0) {
        V(1, 3) = -V(1, 3);
        V(2, 3) = -V(2, 3);
        V(3, 3) = -V(3, 3);
        // unless input is extremely noisy or colinear, R should now be a rotation
        R = V*(U.t());

        if(R.Determinant() < 0) {
        // if we still failed, at least give a warning
        cerr << "Warning - PlaneSVD::iterate failed to compute a matching transformation - the returned solution is a reflection!" << endl;
        }
    }

    // Calculate translation
    double translation[3];
    ColumnVector col_vec(3);
    for(int j = 0; j < 3; j++)
        col_vec(j+1) = cd[j];
    ColumnVector r_time_colVec = ColumnVector(R*col_vec);
    translation[0] = cm[0] - r_time_colVec(1);
    translation[1] = cm[1] - r_time_colVec(2);
    translation[2] = cm[2] - r_time_colVec(3);

    // TODO: check definition of normal direction!
    // for (size_t i = 0; i < 3; i++)
    //     translation[i] = dtr[i];

    // Create 4x4 homogeneous transformation matrix from SVD result
    double alignfx[16];
    alignfx[0] = R(1,1);
    alignfx[1] = R(2,1);
    alignfx[2] = 0;
    alignfx[2] = R(3,1);
    alignfx[3] = 0;
    alignfx[4] = R(1,2);
    alignfx[5] = R(2,2);
    alignfx[6] = R(3,2);
    alignfx[7] = 0;
    alignfx[8] = R(1,3);
    alignfx[9] = R(2,3);
    alignfx[10] = R(3,3);
    alignfx[11] = 0;
    alignfx[12] = translation[0];
    alignfx[13] = translation[1];
    alignfx[14] = translation[2];
    alignfx[15] = 1;

    // Convert to 6DoF state vector
    double rot[3], tra[3];
    Matrix4ToEuler(alignfx, rot, tra);

    // for(int i = 0; i < 3; ++i) {
    //     X(i+1) = rot[i];
    //     X(i+4) = tra[i];
    // }
    // updateScan();

    // Take state X and convert to 4x4, then apply alignfx and
    // convert the result back to 6DoF state.
    for(int i = 0; i < 3; ++i) {
        dX(i+1) = rot[i];
        dX(i+4) = tra[i];
    }
    lock();
    for(int i = 0; i < 3; ++i) {
        rot[i] = dX(i+1);
        tra[i] = dX(i+4);
    }

    // Apply transformation to PlaneScan object and add frame to the animation
    EulerToMatrix4(tra, rot, alignfx);
    this->ps->transform(alignfx);

    // Save resulting transformation in state vector X
    Matrix4ToEuler(ps->transMat, rot, tra);
    for(int i = 0; i < 3; ++i) {
        X(i+1) = rot[i];
        X(i+4) = tra[i];
    }

    if (_anim != -2 && _anim > 0 && ++frames_count % _anim == 0) {
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

    // Free memory of centered correspondences again
    for(unsigned int i = 0; i <  ps->correspondences.size(); i++) {
        delete[] m[i];
        delete[] d[i];
    }
    delete[] m;
    delete[] d;
}