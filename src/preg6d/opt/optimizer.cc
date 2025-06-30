#include "opt/optimizer.h"

int Optimizer::_max_iter;
int Optimizer::_update_cor;
double Optimizer::_eps_convergence;
double Optimizer::_eps_kernel;
int Optimizer::_anim;
bool Optimizer::_quiet;
bool Optimizer::_use_p2p = true;

void Optimizer::setMaxIter(int val)
{
    _max_iter = val;
}

void Optimizer::setEpsConvergence(double val)
{
    _eps_convergence = val;
}

void Optimizer::setUpdateCor(int val)
{
    _update_cor = val;
}

void Optimizer::setAnim(int val)
{
    _anim = val;
}

void Optimizer::setUseP2P(bool val)
{
    _use_p2p = val;
}

void Optimizer::setQuiet(bool val)
{
    _quiet = val;
}

void Optimizer::setEpsKernel(double val)
{
    _eps_kernel = val;
}

double Optimizer::getRSE()
{
    if (_use_p2p)
        return ps->calcErr() + ps->calcErrPtPairs();
    else
        return ps->calcErr();
}

bool* Optimizer::convertDimensionsToBoolArray(Dimensions d)
{
    bool* res = new bool[6];
    for(int i = 0; i < 6; ++i) res[i] = true;
    switch(d)
    {
        case ALL:
            break;
        case UDLOCK:
            res[4]=false;
            break;
        case RADLER:
            res[0]=false;
            res[2]=false;
            res[4]=false;
            break;
        // TODO: FINISH
        default:
            break;

    }
    return res;
}

/* weights are:
 * w_i = 1 / r_i * p'(r_i)
 * where p' is the first order derivative of loss p
 *
 * For loss function table, see:
 * P. Babin et al.
 * https://arxiv.org/pdf/1810.01474.pdf
 */
double Optimizer::weightLoss(double e) {

    return 1; // L2

    //return 1/fabs(e); // L1

    // return fabs(e) <= _eps_kernel ? 1 : _eps_kernel / fabs(e); // Huber

    // return 1 / (1 + sqr(e/_eps_kernel)); // Cauchy

    // return fabs(e) <= _eps_kernel ? sqr(1-sqr(e/_eps_kernel)) : 0; // Tukey
}

void Optimizer::lock()
{
    //
    double drPosTheta[3] = {dX(1), dX(2), dX(3)};
    double drPos[3] = {dX(4), dX(5), dX(6)};
    switch(dim)
    {
        case ALL:
            break;
        case ROLLING:
            drPosTheta[1] = 0;
            drPosTheta[2] = 0;
            break;
        case DESCENDING:
            drPosTheta[0] = 0;
            drPosTheta[2] = 0;
            drPos[0] = 0;
            drPos[2] = 0;
            break;
        case ROTATING:
            drPosTheta[0] = 0;
            drPosTheta[2] = 0;
            drPos[0] = 0;
            drPos[1] = 0;
            drPos[2] = 0;
            break;
        //AKA. JASPER-MODE.
        case UDLOCK:
            drPos[1] = 0;
            break;
        case SLIDING:
            drPos[1] = 0;
            drPosTheta[0] = 0;
            drPosTheta[1] = 0;
            drPosTheta[2] = 0;
            break;
        case RADLER:
            drPosTheta[0] = 0;
            drPosTheta[2] = 0;
            drPos[1] = 0;
            break;
        case NOROT:
            drPosTheta[0] = 0;
            drPosTheta[1] = 0;
            drPosTheta[2] = 0;
            break;
    }
    dX(1) = drPosTheta[0];
    dX(2) = drPosTheta[1];
    dX(3) = drPosTheta[2];
    dX(4) = drPos[0];
    dX(5) = drPos[1];
    dX(6) = drPos[2];
}

void Optimizer::updateScan()
{
    ps->rPosTheta[0] = X(1); // roll
    ps->rPosTheta[1] = X(2); // pitch
    ps->rPosTheta[2] = X(3); // yaw
    ps->rPos[0] = X(4);
    ps->rPos[1] = X(5);
    ps->rPos[2] = X(6);
    EulerToMatrix4(ps->rPos, ps->rPosTheta, ps->transMat);
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
}

bool Optimizer::stop_condition(ColumnVector state, double eps)
{
    // Default false -> do not stop
    bool result = false;

    //cout << "States: " << endl;
    // Compare current state with other states
    for (ColumnVector& s : Xk) {

        //cout << "State" << s << " ";
        ColumnVector ds = state - s;

        // Loop or too small change detected
        if (ds.SumAbsoluteValue() <= eps) {
            cout << endl << "Loop or too small change detected for state " << endl << state;
            result = true;
            break;
        }
    }

    // Add state to deque
    Xk.push_back(state);
    if (Xk.size() > COMP_N_FRAMES)
        Xk.pop_front();

    return result;
}