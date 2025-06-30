/** @file
 *  @brief Representation of a plane in multiple formats.
 *
 * This class adapts to:
 *      - Planes that are defined by normal and reference point vector.
 *      - Planes that are defined by their convex hull.
 *      - Planes that are defined by the bestfit through a cluster of points.
 * See contructor briefings for reference how to do so.
 *
 * The class is self-backwards compatible, which means that you can, e.g., construct
 * a plane by a cluster of points while normal, reference point, and convex hull
 * representation get also automaticaly infered.
 *
 * @author Fabian Arzberger, JMU, Germany.
 *
 * Released under the GPL version 3.
 */

#ifndef __NORMALPLANE_H_
#define __NORMALPLANE_H_

#include <vector>
#include "slam6d/globals.icc"
#include "slam6d/scan.h"
#include "slam6d/bkd.h"
#include "shapes/convexplane.h"
#include "util.h"

using namespace std;

// Small little helper:
static inline double dot(Point p1, Point p2)
{
    double *p1d = new double[3];
    p1d[0] = p1.x;
    p1d[1] = p1.y;
    p1d[2] = p1.z;
    double *p2d = new double[3];
    p2d[0] = p2.x;
    p2d[1] = p2.y;
    p2d[2] = p2.z;
    double res = Dot(p1d, p2d);
    delete[] p1d;
    delete[] p2d;
    return res;
}

class NormalPlane
{
public:
    NormalPlane() {}

    /**
     * @brief: Constructs a Plane from the Convex Hull.
     * @param convex_hull: A vector of double* containing the points of the convex hull.
     */
    NormalPlane(vector<double*> &convex_hull);

    /**
     * @brief: Constructs a Plane from normal and reference point.
     * Note: Convex hull will be missing. Can be added at runtime using member functions.
     * @param n: The normal vector of the plane.
     * @param x: The reference point of the plane.
     */
    NormalPlane(double* n, double* x);

    /**
     * @brief: Constructs a Plane from normal and reference point.
     * Then, the Convex Hull is added without sanity checking for n and x.
     * @param n: The normal vector of the plane.
     * @param x: The reference point of the plane.
     * @param convex_hull: The vector containing the points of the convex hull.
     */
    NormalPlane(double* n, double* x, vector<double*> &convex_hull);

    /**
     * @brief: Fits a Plane through a cluster of Points.
     * If this constructor is used, all member attributes are well defined.
     * @param point_cluster: A vector of Points containing all points of the cluster.
     */
    //NormalPlane(vector<Point> &point_cluster);

    ~NormalPlane();

    /*
     * Plane Attributes. Not all of them are defined at runtime!
     * TODO: Make availability checks in member functions.
     */

    double n[3]; // normal
    double x[3]; // reference point
    double eigen[3]; // eigenvals of the plane (not always available!)
    double rho;
    char direction;

    // if created by a cluster, the normalplane stores all pts from the cluster
    vector<double*> hull; // convexhull

    vector<double*> all_pts; // all points
    BkdTree *search_tree; // all points in a search tree

    double& operator()(int); // get normal xyz at index
    double& operator[](int); // get reference xyz at index

    // Adds a 3D convexhull
    void addConvexHull(vector<double*> &ps);

    /**
     * @brief: Calculates the 3D convexhull from a given cluster.
     * Note: The function does what it says:
     * Calling this function does NOT create the all_pts attribute.
     * @param points: Vector containing all the points.
     */
    void calculateConvexHullFromCluster(vector<Point> &points);

    /**
     * @brief Calculates the amount of overlap between two planes.
     * @param p: Pointer to the plane to check overlap with
     * @return: Value between 0 and 1, representing the amount of overlap.
     */
    double overlap(NormalPlane* other);

    /**
     * @brief Calculates the hesse distance between to planes.
     * @param p: Pointer to the plane to calculate distance to.
     * @return: The distance.
     */
    double hesseDist2Plane(NormalPlane* other);

    /**
     * @brief Calculates the projected distance between to planes.
     * @param p: Pointer to the plane to calculate distance to.
     * @return: The distance.
     */
    double projDist2Plane(NormalPlane* other);

    /**
     * @brief Merges a plane with the current one.
     * @param p: Plane that should be merged to this.
     */
    void mergeWith(NormalPlane* other);

    void writePlane(string);

    // The part below converts the convex hull into an ordered point set.
    // The following static functions work with the point array to find minimum distances
    // to the Hull.

    Point* hullAsPointArr();

    /**
     * @brief Calculates a simple cubic bounding box around the convex hull
     * of the plane.
     * @param mins: Minimum x, y, z values for the bounding box.
     * @param maxs: Maximum x, y, z values for the bounding box.
     */
    void getBBox(double *mins, double *maxs);

    /**
     * @brief Calculates the minimum and maximum hesse distance between to planes.
     * Note: this function is order sensitive! this->getMinMaxHesseTo(other) may return
     * differently than other->getMinMaxHesseTo(this) !
     * @param p: Pointer to the plane to calculate distances to.
     * @param min: Reference to store the minimum distance.
     * @param max: Reference to store the maximum distance.
     */
    void getMinMaxHesseTo(NormalPlane* other, double &mind, double &maxd);

    /**
     * @brief Calculates the minimum and maximum hesse distance between to planes.
     * Note: this function is order sensitive! this->getMinMaxProjDistTo(other) may return
     * differently than other->getMinMaxProjDistTo(this) !
     * @param p: Pointer to the plane to calculate distances to.
     * @param min: Reference to store the minimum distance.
     * @param max: Reference to store the maximum distance.
     */
    void getMinMaxProjDistTo(NormalPlane* other, double &mind, double &maxd);

    static inline Point* vec2arr(vector<Point> &points)
    {
        size_t n = points.size();
        Point *arr = new Point[n];
        for ( size_t i = 0; i < n; ++i)
        {
            arr[i].x = points[i].x;
            arr[i].y = points[i].y;
            arr[i].z = points[i].z;
        }
        return arr;
    }

    // BEGIN COPYRIGHT

    // Copyright 2001, 2012, 2021 Dan Sunday
    // for the following code snippet:

    // isLeft(): tests if a point is Left|On|Right of an infinite line.
    //    Input:  three points P0, P1, and P2
    //    Return: >0 for P2 left of the line through P0 and P1
    //            =0 for P2  on the line
    //            <0 for P2  right of the line
    static inline int
    isLeft( Point, Point, Point );

    //===================================================================

    // cn_PnPoly(): crossing number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  0 = outside, 1 = inside
    // This code is patterned after [Franklin, 2000]
    static int
    cn_PnPoly( Point P, Point* V, int n );
    //===================================================================

    // wn_PnPoly(): winding number test for a point in a polygon
    //      Input:   P = a point,
    //               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
    //      Return:  wn = the winding number (=0 only when P is outside)
    static int
    wn_PnPoly( Point P, Point* V, int n );

    //===================================================================

    // END COPYRIGHT

    /**
    * @brief Calculates the nearest line segment of point <p> to <polygon>.
    * @param p - input: The point <p> that is checked.
    * @param poly - input: The <polygon> to check with.
    * @param p1 - output: start point of nearest line segment.
    * @param p2 - output: end point of nearest line segment.
    * @returns The minimum distance to the line segment.
    */
    static double nearestLineSegment(const Point& p, const Point* poly, int n, Point& p1, Point& p2);

    /**
     * @brief Converts a 3D point that lies on a plane to 2D.
     * @param p - input: The point to be converted.
     * @param dir - input: main direction of the plane normal.
     * @param out - output: Resulting point.
     */
    static void convert3Dto2D(Point p, char dir, Point& out);

    /**
     * @brief Converts multiple 3D points that lie on a plane into 2D.
     * @param ps - input: The points to be converted.
     * @param n - input: nr. of points.
     * @param dir - input: main direction of the plane normal.
     * @param out - output: Resulting points.
     */
    static void convert3Dto2D(Point* ps, int n, char dir, Point* out);

private:

    Point* _hull_parr;
};

typedef std::vector<NormalPlane*> Planes;

#endif // __NORMALPLANE_H_