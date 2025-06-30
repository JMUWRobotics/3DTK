/** @file
 *  @brief Representation of a plane with normal and reference point.
 *  Optionally stores the convex hull of the plane.
 *  This class also holds some static members to check if a point is
 *  inside the convex hull.
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 */

#include "model/normalplane.h"
#include "model/planescan.h"
#include "math.h"
#include "slam6d/normals.h"
#include <cassert>

// Use (void) to silent unused warnings.
#define assertm(exp, msg) assert(((void)msg, exp))

// Constructors

NormalPlane::NormalPlane(double *n, double *x)
{
    for(int i = 0; i < 3; i++)
    {
        this->n[i] = n[i];
        this->x[i] = x[i];
    }
    rho = Dot(n,x);
    if(rho < 0)
    {
        rho = -rho;
        for(int i = 0; i < 3; i++) {
            n[i] = -n[i];
        }
    }
    if(fabs(n[0]) < fabs(n[1])) {
        if(fabs(n[1]) < fabs(n[2])) {
        direction = 'z';
        } else {
        direction = 'y';
        }
    } else if (fabs(n[2]) < fabs(n[0])){
        direction = 'x';
    } else {
        direction = 'z';
    }
}

NormalPlane::NormalPlane(vector<double*> &ps)
{
    vector<Point> pts;
    search_tree = new BkdTree();
    search_tree->insert(ps);
    for (double* p : ps) {
        all_pts.push_back( p );
        pts.push_back( Point( p ) );
    }

    calculateConvexHullFromCluster(pts);

    // convert to ordered point array:
    _hull_parr = new Point[hull.size()];
    for ( size_t i = 0; i < hull.size(); ++i)
    {
        _hull_parr[i].x = hull[i][0];
        _hull_parr[i].y = hull[i][1];
        _hull_parr[i].z = hull[i][2];
    }
}

NormalPlane::NormalPlane(double *n, double *x, vector<double*> &ps)
    : NormalPlane(n, x)
{
    addConvexHull(ps);
}

// Calculate the Convex Hull from a cluster
// NormalPlane::NormalPlane(vector<Point> &cluster)
// {
//     std::vector<double*> tmp_pts;
//     for (const Point& p : cluster) {
//         double *d = new double[3]{p.x, p.y, p.z};
//         all_pts.push_back( d );
//         tmp_pts.push_back( d );
//     }
//     search_tree = new BkdTree();
//     search_tree->insert(tmp_pts);

//     calculateConvexHullFromCluster(cluster);

//     // convert to ordered point array:
//     _hull_parr = new Point[hull.size()];
//     for ( size_t i = 0; i < hull.size(); ++i)
//     {
//         _hull_parr[i].x = hull[i][0];
//         _hull_parr[i].y = hull[i][1];
//         _hull_parr[i].z = hull[i][2];
//     }
// }

// Operators

double& NormalPlane::operator()(int idx)
{
    return n[idx];
}

double& NormalPlane::operator[](int idx)
{
    return x[idx];
}

// Accessors / Mutators

void NormalPlane::addConvexHull(std::vector<double*> &ps)
{
    hull = ps;
    // Convert to Point array.
    _hull_parr = new Point[ps.size()];
    for ( size_t i = 0; i < ps.size(); ++i)
    {
        _hull_parr[i].x = ps[i][0];
        _hull_parr[i].y = ps[i][1];
        _hull_parr[i].z = ps[i][2];
    }
}

void NormalPlane::calculateConvexHullFromCluster(vector<Point> &cluster)
{
    // Getting plane parameters
    double center[3];
    for(int i = 0; i < 3; i++) {
        center[i] = 0.0;
    }
    //cout << "calculate centroid of the cluster" << endl;
    for(size_t i = 0; i < cluster.size(); ++i)
    {
        center[0] += cluster[i].x;
        center[1] += cluster[i].y;
        center[2] += cluster[i].z;
    }
    for(int i = 0; i < 3; i++) {
        center[i] /= cluster.size();
        x[i] = center[i];
    }
    //cout << "done." << endl;
    //cout << "Calc normal" << endl;
    calculateNormal(cluster, n, eigen);
    rho = Dot(x, n);

    // TODO: ADD PLANE UNCERTAINTY HERE!
    // See: https://www.ipb.uni-bonn.de/pdfs/foerstner17efficient.pdf
    //

    if(rho < 0)
    {
        rho = -rho;
        for(int i = 0; i < 3; i++) {
            n[i] = -n[i];
        }
    }
    if(fabs(n[0]) < fabs(n[1])) {
        if(fabs(n[1]) < fabs(n[2])) {
            direction = 'z';
        } else {
            direction = 'y';
        }
    } else if (fabs(n[2]) < fabs(n[0])){
        direction = 'x';
    } else {
        direction = 'z';
    }
    //cout << "Done" << endl;

    //cout << "Make convex hull 2d" << endl;
    // Calculating 2d convex hull
    std::list<double *> point_list;
    for (std::vector<Point>::iterator it = cluster.begin(); it != cluster.end(); it++) {
        Point p = (*it);
        double * point = new double[2];
        switch(direction) {
            case 'x': point[0] = p.y;
                      point[1] = p.z;
                    break;
            case 'y': point[0] = p.x;
                      point[1] = p.z;
                    break;
            case 'z': point[0] = p.x;
                      point[1] = p.y;
                    break;
            default: throw std::runtime_error("default branch taken");
        }

        point_list.push_back(point);
    }
    //cout << "List of 2d Points created. Convexing now..." << endl;
    vector<double*> hull2d;
    if (point_list.size() > 0)
        ConvexPlane::JarvisMarchConvexHull(point_list, hull2d);
    //cout << "Done." << endl;

    // TODO: calc area

    //cout << "Embedd 2d convex hull in 3d space" << endl;
    // Convert back to 3d
    for(std::vector<double*>::iterator it = hull2d.begin();
        it != hull2d.end();
        it++) {

        double *point = new double[3];

        switch(direction) {
            case 'x':
                point[0] = (rho - (*it)[0] * n[1] - (*it)[1] * n[2]) / n[0];
                point[1] = (*it)[0];
                point[2] = (*it)[1];
                break;
            case 'y':
                point[0] = (*it)[0];
                point[1] = (rho - (*it)[0] * n[0] - (*it)[1] * n[2]) / n[1];
                point[2] = (*it)[1];
                break;
            case 'z':
                point[0] = (*it)[0];
                point[1] = (*it)[1];
                point[2] = (rho - (*it)[0] * n[0] - (*it)[1] * n[1]) / n[2];
                break;
            default: throw std::runtime_error("default branch taken");
        }
        hull.push_back(point);
    }
    //cout << "Done." << endl;
    return ;
}

void NormalPlane::getBBox(double *mins, double *maxs)
{
    size_t hsize = hull.size();
    assertm( (hsize > 0), "Convex hull has no pts.");

    for (int i = 0; i < 3; ++i) {
        // Init with first pt
        mins[i] = hull[0][i];
        maxs[i] = hull[0][i];
    }
    for (size_t i = 1; i < hsize; ++i) {
        for(int j = 0; j < 3; ++j) {
            mins[j] = min(mins[j], hull[i][j]);
            maxs[j] = max(maxs[j], hull[i][j]);
        }
    }
}

double NormalPlane::projDist2Plane(NormalPlane *other)
{
    // Attention: ppd is order sensitive, as the direction of
    // projection matters. Therefore, take max of those.
    assertm( (this->hull.size() != 0
           && other->hull.size() != 0),
           "There exists no hull." );
    double res;
    double res1 = __DBL_MAX__, res2 = __DBL_MAX__;
    for (double* p : this->hull) {
        double tmp = PlaneScan::projDist2Plane(p, other);
        res1 = min(res1, tmp);
    }
    for (double* p : other->hull) {
        double tmp = PlaneScan::projDist2Plane(p, this);
        res2 = min(res2, tmp);
    }
    res = max(res1, res2);
    return res;
}

// Check convex hull overlap by
// querying every point of the convex hull of this object into "isInPlane",
// checking the convex hull of the other plane.
double NormalPlane::overlap(NormalPlane *other)
{
    assertm( (this->all_pts.size() != 0
          && other->all_pts.size() != 0),
            "Plane doesnt contain any pts." );
    double ratio1, ratio2;
    size_t hit = 0, all = 0;
    double Mid[16];
    M4identity(Mid);
    for (double *p : this->all_pts) {
        if ( PlaneScan::isInPlane(p, Mid, other) )
            hit++;
        all++;
    }
    ratio1 = all == 0 ? 0 : (double)hit/all;
    hit = all = 0;
    for (double *p : other->all_pts) {
        if ( PlaneScan::isInPlane(p, Mid, this) )
            hit++;
        all++;
    }
    ratio2 = all == 0 ? 0 : (double)hit/all;
    return max(ratio1, ratio2);
}

void NormalPlane::getMinMaxHesseTo(NormalPlane* other, double &mind, double &maxd)
{
    assertm( (this->hull.size() != 0),
           "There exists no hull." );
    mind = __DBL_MAX__;
    maxd = 0;
    for (double* p : this->hull) {
        double tmp = fabs(PlaneScan::dist2Plane(p, other));
        mind = min(tmp, mind);
        maxd = max(tmp, maxd);
    }
}

void NormalPlane::getMinMaxProjDistTo(NormalPlane* other, double &mind, double &maxd)
{
    assertm( (this->hull.size() != 0),
           "There exists no hull." );
    mind = __DBL_MAX__;
    maxd = 0;
    for (double* p : this->hull) {
        double tmp = PlaneScan::projDist2Plane(p, other);
        mind = min(tmp, mind);
        maxd = max(tmp, maxd);
    }
}

double NormalPlane::hesseDist2Plane(NormalPlane *other)
{
    // Attention: hesse is order sensitive, as the direction of
    // projection matters. Therefore, take max of those.
    assertm( (this->hull.size() != 0
           && other->hull.size() != 0),
           "There exists no hull." );
    double res;
    double res1 = __DBL_MAX__, res2 = __DBL_MAX__;
    for (double* p : this->hull) {
        double tmp = fabs(PlaneScan::dist2Plane(p, other));
        res1 = min(res1, tmp);
    }
    for (double* p : other->hull) {
        double tmp = fabs(PlaneScan::dist2Plane(p, this));
        res2 = min(res2, tmp);
    }
    res = max(res1, res2);
    return res;
}


void NormalPlane::mergeWith(NormalPlane *other)
{
    assertm( (this != other), "Trying to merge plane with itself.");
    vector<Point> pts;
    pts.insert( end(pts), begin(hull), end(hull));
    pts.insert( end(pts), begin(other->hull), end(other->hull));

    this->hull.clear();
    delete[] _hull_parr;

    this->all_pts.insert( end(all_pts), begin(other->all_pts), end(other->all_pts) );
    this->search_tree->insert(other->all_pts);

    calculateConvexHullFromCluster(pts);

    // convert to ordered point array:
    _hull_parr = new Point[hull.size()];
    for ( size_t i = 0; i < hull.size(); ++i)
    {
        _hull_parr[i].x = hull[i][0];
        _hull_parr[i].y = hull[i][1];
        _hull_parr[i].z = hull[i][2];
    }
    //cout << "Deleting other..." << endl;

    // delete other plane after merging with this
    delete other;
}

void NormalPlane::writePlane(string path)
{
    std::ofstream out;
  out.open(path.c_str());

  for(std::vector<double*>::iterator it = hull.begin();
      it != hull.end();
      it++) {
        //cout << "Export point [" << (*it)[0] << ", " << (*it)[1] << ", " << (*it)[2] << "]" << endl;
    switch(direction) {
      case 'x':
        out << (*it)[0] << " ";
        out << (*it)[1] << " ";
        out << (*it)[2] << std::endl;
        break;
      case 'y':
        out << (*it)[0] << " ";
        out << (*it)[1] << " ";
        out << (*it)[2] << std::endl;
        break;
      case 'z':
        out << (*it)[0] << " ";
        out << (*it)[1] << " ";
        out << (*it)[2] << std::endl;
        break;
      default: throw std::runtime_error("default branch taken");
    }

  }
  out.flush();
  out.close();
}

Point* NormalPlane::hullAsPointArr()
{
    return _hull_parr;
}

// Static functions

int NormalPlane::isLeft(Point P0, Point P1, Point P2 )
{
return ( (P1.x - P0.x) * (P2.y - P0.y)
                - (P2.x -  P0.x) * (P1.y - P0.y) );
}

int NormalPlane::cn_PnPoly(Point P, Point* V, int n)
{
    int cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
    if (((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
        || ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
            if (P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn&1);    // 0 if even (out), and 1 if  odd (in)
}

int NormalPlane::wn_PnPoly( Point P, Point* V, int n )
{
    int    wn = 0;    // the  winding number counter

    // loop through all edges of the polygon
    for (int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
        if (V[i].y <= P.y) {          // start y <= P.y
            if (V[i+1].y  > P.y)      // an upward crossing
                if (isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
                    ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if (V[i+1].y  <= P.y)     // a downward crossing
                if (isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
                    --wn;            // have  a valid down intersect
        }
    }
    return wn;
}

void NormalPlane::convert3Dto2D(Point p, char dir, Point& out)
{
    switch(dir) {
        case 'x': out.x = p.y;
                out.y = p.z;
                break;
        case 'y': out.x = p.x;
                out.y = p.z;
                break;
        case 'z': out.x = p.x;
                out.y = p.y;
                break;
        default: throw std::runtime_error("default branch taken");
    }
}

void NormalPlane::convert3Dto2D(Point* ps, int n, char dir, Point* out)
{
    for (int i = 0; i < n; ++i)
        convert3Dto2D(ps[i], dir, out[i]);
}

double NormalPlane::nearestLineSegment(const Point& p, const Point* polygon, int n, Point& p1, Point& p2)
{
    Point pch2_1;
    Point p_pch1;
    Point pch2_p;
    Point s;
    double D_min = __DBL_MAX__;
    double D = D_min;
    int j_min = 0, i_min = 1;
    int j = j_min;
    for (int i = i_min; i < n; )
    {
        p_pch1 = p - polygon[j];
        pch2_1 = polygon[i] - polygon[j];
        pch2_p = polygon[i] - p;

        // Projection factor onto the line segment. if between 0 and 1,
        // point projection is on the line segment. Else its not in the linesegment.
        double r = dot( pch2_1, p_pch1 ) / dot( pch2_1, pch2_1 );
        // Clamping between 0 and 1. If 0, calc distance to polygon[j].
        // If 1, calc distance to polygon[i]
        double t = min( max(r, 0.0), 1.0);
        // Parametrization of a line, starting at polygon[j] ending at polygon[i]
        s = polygon[j] + t*pch2_1 - p;
        // The shortest distance, either to the line or one of the points.
        D = sqrt(dot( s, s )) ;
        if (fabs(D) < D_min) {
            D_min = D;
            j_min = j;
            i_min = i;
        }
        j++;
        i++;
    }
    p1 = polygon[j_min];
    p2 = polygon[i_min];
    return D_min;
}

NormalPlane::~NormalPlane()
{
    // for(double* p : all_pts)
    //     delete p;
    // for(double* p : hull)
    //     delete p;
    if (_hull_parr) {
        delete[] _hull_parr;
        _hull_parr = 0;
    }
    if (search_tree) {
        delete search_tree;
        search_tree = 0;
    }
    hull.clear();
    all_pts.clear();
}
