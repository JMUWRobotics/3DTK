/** @file
 *  @brief Representation of the optimized k-d tree.
 *  @author Remus Dumitru. Jacobs University Bremen, Germany
 *  @author Corneliu-Claudiu Prodescu. Jacobs University Bremen, Germany
 *  @author Andreas Nuechter. Jacobs University Bremen, Germany
 *  @author Kai Lingemann. Inst. of CS, University of Osnabrueck, Germany
 *  @author Momchil Ivanov Ivanov. Jacobs University Bremen, Germany
 *  @author Igor Pener. Jacobs University Bremen, Germany
 *  @author Fabian Arzberger. Julius-Maximilian University Wuerzburg, Germany
 */

#ifndef __KD_TREE_IMPL_H__
#define __KD_TREE_IMPL_H__

#include "slam6d/kdparams.h"
#include "globals.icc"

#include <stdio.h>

#ifdef _MSC_VER
#if !defined _OPENMP && defined OPENMP
#define _OPENMP
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

struct Void {};

class PointCompare {
public:
    bool operator() (const std::pair<Point, double>& left,
				 const std::pair<Point, double>& right)
  {
    return left.second > right.second;
  }
};

/**
 * @brief The optimized k-d tree.
 *
 * A kD tree for points, with limited
 * capabilities (find nearest point to
 * a given point, or to a ray).
 *
 * You can hand over a whole scan of points and limit the operations on them by
 * specifying the indices of only a few of them. These indices are not necessarily
 * numeric though; they could also be the points themselves. This varies between the
 * supplied implementations (e.g. KDtree, KDtreeIndexed) and is a result of the
 * definition of these template parameters:
 *
 * PointData    the type of the input point data
 * AccessorData the type of indices
 * AccessorFunc retrieves data of type double[3], given an index of type
 *              AccessorData and the data of type PointData
 * PointType    the type that is stored in kdparams
 * ParamFunc    retrieves data of type PointType, given an index of type
 *              AccessorData and the data of type PointData
 **/
template<class PointData, class AccessorData, class AccessorFunc, class PointType, class ParamFunc>
class KDTreeImpl {
public:
  inline KDTreeImpl() { }

  virtual inline ~KDTreeImpl() {
    if (!npts) {
#ifdef WITH_OPENMP_KD
      omp_set_num_threads(OPENMP_NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
#endif
      for (int i = 0; i < 2; i++) {
        if (i == 0 && node.child1) delete node.child1;
        if (i == 1 && node.child2) delete node.child2;
      }
    } else {
      if (leaf.p) delete [] leaf.p;
    }
  }

  virtual void create(PointData pts, AccessorData *indices, size_t n,
                      unsigned int bucketSize = 20) {
    AccessorFunc point;

    if (n == 0) {
        throw std::runtime_error("cannot create kdtree with zero points");
    }

    // Find bbox and centroid
    double mins[3], maxs[3];
    double centroid[3];

    for (int i = 0; i < 3; i++) {
      // Initialize with first point
      mins[i] = point(pts, indices[0])[i];
      maxs[i] = point(pts, indices[0])[i];
      centroid[i] = point(pts, indices[0])[i];
    }

    for(size_t i = 1; i < n; i++) {
      for (int j = 0; j < 3; j++) {
        mins[j] = std::min(mins[j], point(pts, indices[i])[j]);
        maxs[j] = std::max(maxs[j], point(pts, indices[i])[j]);
        centroid[j] += point(pts, indices[i])[j];
      }
    }

    for (int i = 0; i < 3; i++) {
      centroid[i] /= n;
    }

    // Leaf nodes
    if ((n > 0) && (n <= bucketSize)) {
      npts = n;
      isleaf = true;
      leaf.p = new AccessorData[n];
      // fill leaf index array with indices
      for(size_t i = 0; i < n; ++i) {
        leaf.p[i] = indices[i];
      }
      return;
    }

    // Else, interior nodes
    npts = 0;
    isleaf = false;
    for (int i = 0; i < 3; i++) {
      node.center[i] = 0.5 * (mins[i]+maxs[i]);
    }

    node.dx = 0.5 * (maxs[0]-mins[0]);
    node.dy = 0.5 * (maxs[1]-mins[1]);
    node.dz = 0.5 * (maxs[2]-mins[2]);
    node.r = sqrt(sqr(node.dx) + sqr(node.dy) + sqr(node.dz));

    // Find longest axis
    if (node.dx > node.dy) {
      if (node.dx > node.dz) {
        node.splitaxis = 0;
      } else {
        node.splitaxis = 2;
      }
    } else {
      if (node.dy > node.dz) {
        node.splitaxis = 1;
      } else {
        node.splitaxis = 2;
      }
    }

    // Put points that were measured very closely together in the same bucket
    if ( fabs(std::max(std::max(node.dx,node.dy),node.dz)) < 0.01 ) {
      npts = n;
      isleaf = true;
      leaf.p = new AccessorData[n];
      // fill leaf index array with indices
      for(size_t i = 0; i < n; ++i) {
        leaf.p[i] = indices[i];
      }
      return;
    }

    // Partition

    // Old method, splitting at the center of the bbox
    // double splitval = node.center[node.splitaxis];

    // Now we split at the centroid (average of all points)
    node.splitval = centroid[node.splitaxis];

    AccessorData* left = indices,
                * right = indices + n - 1;
    while(true) {
      while(point(pts, *left)[node.splitaxis] < node.splitval)
        left++;
      while(point(pts, *right)[node.splitaxis] >= node.splitval)
        right--;
      if(right < left)
        break;
      std::swap(*left, *right);
    }

    // Build subtrees
    int i;
#ifdef WITH_OPENMP_KD                   // does anybody know the reason why this is slower ?? --Andreas
    omp_set_num_threads(OPENMP_NUM_THREADS);
    // we probably need to use a different OMP construct here (like "task") --Johannes B.
#pragma omp parallel for schedule(dynamic)
#endif
    for (i = 0; i < 2; i++) {
      if (i == 0) {
        node.child1 = new KDTreeImpl();
        node.child1->create(pts, indices, left - indices, bucketSize);
      }
      if (i == 1) {
        node.child2 = new KDTreeImpl();
        node.child2->create(pts, left, n - (left - indices), bucketSize);
      }
    }
  }

protected:
  /**
   * storing the parameters of the k-d tree, i.e., the current closest point,
   * the distance to the current closest point and the point itself.
   * These global variable are needed in this search.
   *
   * Padded in the parallel case.
   */
#ifdef _OPENMP
#ifdef __INTEL_COMPILER
  __declspec (align(16)) static KDParams params[MAX_OPENMP_NUM_THREADS];
#else
  static KDParams<PointType> params[MAX_OPENMP_NUM_THREADS];
#endif //__INTEL_COMPILER
#else
  static KDParams<PointType> params[MAX_OPENMP_NUM_THREADS];
#endif

  /**
   * number of points. If this is 0: intermediate node. If nonzero: leaf.
   */
  int npts;

  /**
   * When removing points, leaf nodes could get npts == 0. To prevent
   * treating them as intermediate nodes, we use this flag since unions wont
   * allow us to check whether a node is leaf or non-leaf
   */
  bool isleaf;

  /**
   * Cue the standard rant about anon unions but not structs in C++
   */
  union {
    /**
     * in case of internal node...
     */
    struct {
      double center[3]; ///< storing the center of the voxel (R^3)
      double dx,  ///< defining the voxel itself
	        dy,  ///< defining the voxel itself
	        dz,  ///< defining the voxel itself
             r;  ///< defining the voxel itself
      int splitaxis;   ///< defining the kind of splitaxis
      double splitval; ///< position of the split
      KDTreeImpl *child1;  ///< pointers to the childs
      KDTreeImpl *child2;  ///< pointers to the childs
    } node;
    /**
     * in case of leaf node ...
     */
    struct {
      /**
       * store the value itself
       * Here we store just a pointer to the data
       */
      AccessorData* p;
    } leaf;
  };

  void _CollectPts(const PointData& pts, int threadNum) const {
      AccessorFunc point;
      ParamFunc pointparam;

      // TODO: Here we need an implementation that stores shared_ptr<double>
      // That way, we dont need to allocate new memory everytime we merge trees in a bkd-forest.

      // TODO: Another idea: Write something like CollectPtsAndDelete(...) that would
      // automatically address new memory and free the old points memory. Ideally, we would
      // like to keep the original pointers untouched, while only moving references to them.
      if (npts) {
          for (int i = 0; i < npts; ++i)
              params[threadNum].collected_pts.push_back(pointparam(pts, leaf.p[i]));
              //params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));
          return;
      }

      node.child1->_CollectPts(pts, threadNum);
      node.child2->_CollectPts(pts, threadNum);
  }

  /**
   * @brief Removes a point from the tree by swapping.
   * @return How many points have been removed. Will be 0 or 1.
   */
  int _Remove(const PointData& pts, int threadNum) {
        AccessorFunc point;
        ParamFunc pointparam;
        int index_remove = -1;
        double closestd2 = __DBL_MAX__;

        // If leaf node (if nonzero)
        if (isleaf) {
            // Search for the (closest) point to be deleted
            for (int i = 0; i < npts; i++) {
                double d2 = Dist2( params[threadNum].p , point(pts, leaf.p[i]) );
                if (d2 < closestd2) {
                    closestd2 = d2;
                    params[threadNum].closest = pointparam(pts, leaf.p[i]);
                    index_remove = i;
                }
            }

            // Remove the (closest) point (if it is close enough)
            if (closestd2 < 0.000000001) {
                if (npts > 1 && index_remove != -1) {
                    // Swap elem to be removed with last elem
                    AccessorData *ptr2rmv = leaf.p + index_remove;
                    AccessorData *end = leaf.p + npts - 1;
                    std::swap(*ptr2rmv, *end);
                    // Exclude last elem in the future by decrementing nr of pts.
                    npts = npts - 1;
                    return 1; // we removed one point.
                }
                // only one point left...
                else if (npts == 1 && index_remove != -1) {
                    npts = npts - 1; // no need to swap this time
                    return 1; // we removed the last point.
                }
            }
            return 0; // we removed zero points
        }

        int removed = 0;
        // Else, If not leaf node (interior), traverse the tree recursivley
        double myd = node.splitval - params[threadNum].p[node.splitaxis];
        if (myd > 0.0) // go right
            removed += node.child1->_Remove(pts, threadNum);
        else if (myd < 0.0) // go left
            removed += node.child2->_Remove(pts, threadNum);
        else { // unsure, search both paths
            removed += node.child1->_Remove(pts, threadNum);
            removed += node.child2->_Remove(pts, threadNum);
        }
        return removed;
  }

  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the recursive case every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _FindClosest(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc   pointparam;

    // Leaf nodes
    if (isleaf) {
      for (int i = 0; i < npts; i++) {
        double myd2 = Dist2(params[threadNum].p, point(pts, leaf.p[i]));
        if (myd2 < params[threadNum].closest_d2) {
          params[threadNum].closest_d2 = myd2;
          params[threadNum].closest = pointparam(pts, leaf.p[i]);
        }
      }
      return;
    }

    // Quick check of whether to abort
    double approx_dist_bbox =
	 std::max(std::max(fabs(params[threadNum].p[0]-node.center[0])-node.dx,
		    fabs(params[threadNum].p[1]-node.center[1])-node.dy),
		fabs(params[threadNum].p[2]-node.center[2])-node.dz);
    if (approx_dist_bbox >= 0 &&
	   sqr(approx_dist_bbox) >= params[threadNum].closest_d2)
      return;

    // Recursive case
    double myd = node.splitval - params[threadNum].p[node.splitaxis];
    if (myd >= 0.0) {
      node.child1->_FindClosest(pts, threadNum);
      if (sqr(myd) < params[threadNum].closest_d2) {
        node.child2->_FindClosest(pts, threadNum);
      }
    } else {
      node.child2->_FindClosest(pts, threadNum);
      if (sqr(myd) < params[threadNum].closest_d2) {
        node.child1->_FindClosest(pts, threadNum);
      }
    }
  }

  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the check whether to abort every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _FindClosestAlongDir(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
      for (int i=0; i < npts; i++) {
        double p2p[] =  { params[threadNum].p[0] - point(pts, leaf.p[i])[0],
                          params[threadNum].p[1] - point(pts, leaf.p[i])[1],
                          params[threadNum].p[2] - point(pts, leaf.p[i])[2] };
        double myd2 = Len2(p2p) - sqr(Dot(p2p, params[threadNum].dir));
        if ((myd2 < params[threadNum].closest_d2)) {
          params[threadNum].closest_d2 = myd2;
          params[threadNum].closest = pointparam(pts, leaf.p[i]);
        }
      }
      return;
    }

    // Quick check of whether to abort
    double p2c[] = { params[threadNum].p[0] - node.center[0],
                     params[threadNum].p[1] - node.center[1],
                     params[threadNum].p[2] - node.center[2] };
    double myd2center = Len2(p2c) - sqr(Dot(p2c, params[threadNum].dir));
    if (myd2center > sqr(node.r + sqrt(params[threadNum].closest_d2)))
      return;

    // Recursive case
    if (params[threadNum].p[node.splitaxis] < node.splitval) {
      node.child1->_FindClosestAlongDir(pts, threadNum);
      node.child2->_FindClosestAlongDir(pts, threadNum);
    } else {
      node.child2->_FindClosestAlongDir(pts, threadNum);
      node.child1->_FindClosestAlongDir(pts, threadNum);
    }
  }

  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the check whether to abort every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _fixedRangeSearchBetween2Points(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	    for (int i = 0; i < npts; i++) {
        double p2p[] =  { params[threadNum].p[0] - point(pts, leaf.p[i])[0],
                          params[threadNum].p[1] - point(pts, leaf.p[i])[1],
                          params[threadNum].p[2] - point(pts, leaf.p[i])[2] };
        double myd2 = Len2(p2p) - sqr(Dot(p2p, params[threadNum].dir));
        if (myd2 < params[threadNum].closest_d2) {
          params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));
	      }
	    }
	    return;
    }

    // Quick check of whether to abort
    double c2c[] = { params[threadNum].p[0] - node.center[0],
                     params[threadNum].p[1] - node.center[1],
                     params[threadNum].p[2] - node.center[2] };

    double my_dist_2 = Len2(c2c); // Distance^2 camera node center
    double myd2center = my_dist_2 - sqr(Dot(c2c, params[threadNum].dir));
    //if (myd2center > (node.r2 + params[threadNum].closest_d2 + 2.0f * max(node.r2, params[threadNum].closest_d2)))

    if (myd2center > sqr(node.r + sqrt(params[threadNum].closest_d2)))
      return;
    //if (myd2center > (node.r2 + params[threadNum].closest_d2 + 2.0f * sqrt(node.r2) * sqrt(params[threadNum].closest_d2))) return;

    // check if not between points

    double p2c[] = { params[threadNum].p0[0] - node.center[0],
                     params[threadNum].p0[1] - node.center[1],
                     params[threadNum].p0[2] - node.center[2] };

    double distXP2 = Len2(p2c);
    if(params[threadNum].dist > distXP2 + node.r) return;

    if(params[threadNum].dist > sqrt(my_dist_2) + node.r) return;

    // Recursive case
    if (params[threadNum].p[node.splitaxis] < node.splitval) {
      node.child1->_fixedRangeSearchAlongDir(pts, threadNum);
      node.child2->_fixedRangeSearchAlongDir(pts, threadNum);
    } else {
      node.child2->_fixedRangeSearchAlongDir(pts, threadNum);
      node.child1->_fixedRangeSearchAlongDir(pts, threadNum);
    }

  }


  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the check whether to abort every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _fixedRangeSearchAlongDir(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	    for (int i = 0; i < npts; i++) {
        /*
        double p2pb[] =  { point(pts, leaf.p[i])[0] - params[threadNum].p[0],
                          point(pts, leaf.p[i])[1] - params[threadNum].p[1],
                          point(pts, leaf.p[i])[2] - params[threadNum].p[2]};
        double blub[3];
        Cross(p2pb, params[threadNum].dir, blub);
        double myd2b = Len2(blub) / Len2(params[threadNum].dir);
	      */

        double p2p[] =  { params[threadNum].p[0] - point(pts, leaf.p[i])[0],
                          params[threadNum].p[1] - point(pts, leaf.p[i])[1],
                          params[threadNum].p[2] - point(pts, leaf.p[i])[2] };
        double myd2 = Len2(p2p) - sqr(Dot(p2p, params[threadNum].dir));
        if (myd2 < params[threadNum].closest_d2) {
          params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));
	      }
	    }
	    return;
    }

    // Quick check of whether to abort
    double p2c[] = { params[threadNum].p[0] - node.center[0],
                     params[threadNum].p[1] - node.center[1],
                     params[threadNum].p[2] - node.center[2] };
    double myd2center = Len2(p2c) - sqr(Dot(p2c, params[threadNum].dir));
    //if (myd2center > (node.r2 + params[threadNum].closest_d2 + 2.0f * max(node.r2, params[threadNum].closest_d2)))
    if (myd2center > sqr(node.r + sqrt(params[threadNum].closest_d2)))
      return;

    // Recursive case
    if (params[threadNum].p[node.splitaxis] < node.splitval) {
      node.child1->_fixedRangeSearchAlongDir(pts, threadNum);
      node.child2->_fixedRangeSearchAlongDir(pts, threadNum);
    } else {
      node.child2->_fixedRangeSearchAlongDir(pts, threadNum);
      node.child1->_fixedRangeSearchAlongDir(pts, threadNum);
    }

  }

  /*
   * search for points inside the axis aligned bounding box given by p and p0
   * where p[0] < p0[0] && p[1] < p0[1] && p[2] < p0[2]
   */
  void _AABBSearch(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	 for (int i = 0; i < npts; i++) {
         double* tp = point(pts, leaf.p[i]);
         if (tp[0] >= params[threadNum].p[0] && tp[0] <= params[threadNum].p0[0]
          && tp[1] >= params[threadNum].p[1] && tp[1] <= params[threadNum].p0[1]
          && tp[2] >= params[threadNum].p[2] && tp[2] <= params[threadNum].p0[2]) {
             params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));
	   }
	 }
	 return;
    }

    // Quick check of whether to abort
    if (node.center[0]+node.dx < params[threadNum].p[0]
     || node.center[1]+node.dy < params[threadNum].p[1]
     || node.center[2]+node.dz < params[threadNum].p[2]
     || node.center[0]-node.dx > params[threadNum].p0[0]
     || node.center[1]-node.dy > params[threadNum].p0[1]
     || node.center[2]-node.dz > params[threadNum].p0[2])
        return;

    // Recursive case
    if (node.splitval > params[threadNum].p[node.splitaxis]) {
        node.child1->_AABBSearch(pts, threadNum);
        if (node.splitval < params[threadNum].p0[node.splitaxis]) {
            node.child2->_AABBSearch(pts, threadNum);
        }
    } else {
        node.child2->_AABBSearch(pts, threadNum);
    }
  }

  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the check whether to abort and the recursive
   *     case every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _FixedRangeSearch(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	 for (int i = 0; i < npts; i++) {
	   double myd2 = Dist2(params[threadNum].p, point(pts, leaf.p[i]));
	   if (myd2 < params[threadNum].closest_d2) {

		params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));

	   }
	 }
	 return;
    }

    // Quick check of whether to abort
    double approx_dist_bbox =
	 std::max(std::max(fabs(params[threadNum].p[0]-node.center[0])-node.dx,
		    fabs(params[threadNum].p[1]-node.center[1])-node.dy),
		fabs(params[threadNum].p[2]-node.center[2])-node.dz);
    if (approx_dist_bbox >= 0 &&
	   sqr(approx_dist_bbox) >= params[threadNum].closest_d2)
	 return;

    // Recursive case
    double myd = node.splitval - params[threadNum].p[node.splitaxis];
    if (myd >= 0.0) {
	 node.child1->_FixedRangeSearch(pts, threadNum);
	 if (sqr(myd) < params[threadNum].closest_d2) {
	   node.child2->_FixedRangeSearch(pts, threadNum);
	 }
    } else {
	 node.child2->_FixedRangeSearch(pts, threadNum);
	 if (sqr(myd) < params[threadNum].closest_d2) {
	   node.child1->_FixedRangeSearch(pts, threadNum);
	 }
    }
  }


  void _KNNSearch(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	 for (int i = 0; i < npts; i++) {
	   double myd2 = Dist2(params[threadNum].p, point(pts, leaf.p[i]));

        for (int j = 0; j < params[threadNum].k; j++)
            if (params[threadNum].distances[j] < 0.0f) {
                params[threadNum].closest_neighbors[j] = pointparam(pts, leaf.p[i]);
                params[threadNum].distances[j] = myd2;
                break;
            } else if (params[threadNum].distances[j] > myd2) {
                // move all other values one place up
                for (int l = params[threadNum].k - 1; l > j; --l) {
                    params[threadNum].closest_neighbors[l] = params[threadNum].closest_neighbors[l-1];
                    params[threadNum].distances[l] = params[threadNum].distances[l-1];
                }
                params[threadNum].closest_neighbors[j] = pointparam(pts, leaf.p[i]);
                params[threadNum].distances[j] = myd2;
                break;
            }
      }
      return;
    }

    int kN = params[threadNum].k-1;
	// FIXME comparing with zero doesn't work for indexed kdtree, we should
	// check the closest vector for the special negative values instead
     /**
       * @author Fabian Arzberger
       * Heres the fix to the FIXME above. Pls recheck somebody!
       * Instead of checking if closest_neighbors[] are set, we look at their
       * corresponding distances.
       */
    if (params[threadNum].distances[kN] != -1) {
        // Quick check of whether to abort
        double approx_dist_bbox
		= std::max(std::max(fabs(params[threadNum].p[0]-node.center[0])-node.dx,
				fabs(params[threadNum].p[1]-node.center[1])-node.dy),
			 fabs(params[threadNum].p[2]-node.center[2])-node.dz);
        if (approx_dist_bbox >= 0 &&
		  sqr(approx_dist_bbox) >= params[threadNum].distances[kN])
		return;
    }
    // Recursive case
    if (params[threadNum].p[node.splitaxis] < node.splitval) {
      node.child1->_KNNSearch(pts, threadNum);
      node.child2->_KNNSearch(pts, threadNum);
    } else {
      node.child2->_KNNSearch(pts, threadNum);
      node.child1->_KNNSearch(pts, threadNum);
    }
  }

  void _KNNRangeSearch(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    // Leaf nodes
    if (isleaf) {
	 for (int i = 0; i < npts; i++) {
	   double myd2 = Dist2(params[threadNum].p, point(pts, leaf.p[i]));
	   if (myd2 >= params[threadNum].closest_d2) {
		   continue;
	   }

        for (int j = 0; j < params[threadNum].k; j++)
            if (params[threadNum].distances[j] < 0.0f) {
                params[threadNum].closest_neighbors[j] = pointparam(pts, leaf.p[i]);
                params[threadNum].distances[j] = myd2;
                break;
            } else if (params[threadNum].distances[j] > myd2) {
                // move all other values one place up
                for (int l = params[threadNum].k - 1; l > j; --l) {
                    params[threadNum].closest_neighbors[l] = params[threadNum].closest_neighbors[l-1];
                    params[threadNum].distances[l] = params[threadNum].distances[l-1];
                }
                params[threadNum].closest_neighbors[j] = pointparam(pts, leaf.p[i]);
                params[threadNum].distances[j] = myd2;
                break;
            }
      }
      return;
    }

    int kN = params[threadNum].k-1;
	// Quick check of whether to abort
	double approx_dist_bbox =
		std::max(std::max(fabs(params[threadNum].p[0]-node.center[0])-node.dx,
					fabs(params[threadNum].p[1]-node.center[1])-node.dy),
				fabs(params[threadNum].p[2]-node.center[2])-node.dz);
	// FIXME comparing with zero doesn't work for indexed kdtree, we should
	// check the closest vector for the special negative values instead
    if (params[threadNum].closest_neighbors[kN] == 0) {
		if (approx_dist_bbox >= 0 &&
				sqr(approx_dist_bbox) >= params[threadNum].closest_d2)
			return;
	} else {
        if (approx_dist_bbox >= 0 &&
		  sqr(approx_dist_bbox) >= params[threadNum].distances[kN])
		return;
    }
    // Recursive case
    double myd = node.splitval - params[threadNum].p[node.splitaxis];
    if (myd >= 0.0) {
	 node.child1->_KNNRangeSearch(pts, threadNum);
	 if (sqr(myd) < params[threadNum].closest_d2) {
	   node.child2->_KNNRangeSearch(pts, threadNum);
	 }
    } else {
	 node.child2->_KNNRangeSearch(pts, threadNum);
	 if (sqr(myd) < params[threadNum].closest_d2) {
	   node.child1->_KNNRangeSearch(pts, threadNum);
	 }
    }
  }

  void _segmentSearch_all(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    double p2p[3], proj[3];
    double t, *comp;
    // Leaf nodes
    if (isleaf) {
        for (int i = 0; i < npts; i++) {
            p2p[0] = point(pts, leaf.p[i])[0] - params[threadNum].p[0];
            p2p[1] = point(pts, leaf.p[i])[1] - params[threadNum].p[1];
            p2p[2] = point(pts, leaf.p[i])[2] - params[threadNum].p[2];
            t = Dot(p2p, params[threadNum].segment_dir);
            if (t < 0.0) {
                // point is beyond point1 of the segment
                comp = params[threadNum].p;
            } else if (t > params[threadNum].segment_len2) {
                // point is beyond point2 of the segment
                comp = params[threadNum].p0;
            } else {
                // point is within the segment
                // calculate its projection onto the line
                proj[0] = params[threadNum].p[0] + t*params[threadNum].segment_n[0];
                proj[1] = params[threadNum].p[1] + t*params[threadNum].segment_n[1];
                proj[2] = params[threadNum].p[2] + t*params[threadNum].segment_n[2];
                comp = proj;
            }
            if (Dist2(comp,point(pts, leaf.p[i])) < params[threadNum].maxdist_d2) {
                params[threadNum].range_neighbors.push_back(pointparam(pts, leaf.p[i]));
            }
        }
        return;
    }

    // Quick check of whether to abort
    double approx_dist_bbox =
        std::max(std::max(fabs(params[threadNum].segment_center[0]-node.center[0])-node.dx,
                    fabs(params[threadNum].segment_center[1]-node.center[1])-node.dy),
                fabs(params[threadNum].segment_center[2]-node.center[2])-node.dz);
    if (approx_dist_bbox >= 0 &&
            sqr(approx_dist_bbox) >= params[threadNum].segment_r2)
        return;
    // Slower check of whether to abort
    p2p[0] = node.center[0] - params[threadNum].p[0];
    p2p[1] = node.center[1] - params[threadNum].p[1];
    p2p[2] = node.center[2] - params[threadNum].p[2];
    t = Dot(p2p, params[threadNum].segment_dir);
    if (t < 0.0) {
        // point is beyond point1 of the segment
        comp = params[threadNum].p;
    } else if (t > params[threadNum].segment_len2) {
        // point is beyond point2 of the segment
        comp = params[threadNum].p0;
    } else {
        // point is within the segment
        // calculate projection
        proj[0] = params[threadNum].p[0] + t*params[threadNum].segment_n[0];
        proj[1] = params[threadNum].p[1] + t*params[threadNum].segment_n[1];
        proj[2] = params[threadNum].p[2] + t*params[threadNum].segment_n[2];
        comp = proj;
    }
    if (Dist2(comp,node.center) > sqr(node.r+params[threadNum].maxdist_d))
        return;

    // Recursive case
    if (params[threadNum].p[node.splitaxis] < node.splitval) {
      node.child1->_segmentSearch_all(pts, threadNum);
      node.child2->_segmentSearch_all(pts, threadNum);
    } else {
      node.child2->_segmentSearch_all(pts, threadNum);
      node.child1->_segmentSearch_all(pts, threadNum);
    }

  }

  /*
   * TODO: benchmark what is faster:
   *   - squaring the distance in the check whether to abort and the recursive
   *     case every time?
   *   - or taking the square root once closest_d2 is updated?
   */
  void _segmentSearch_1NearestPoint(const PointData& pts, int threadNum) const {
    AccessorFunc point;
    ParamFunc pointparam;

    double p2p[3], proj[3];
    double t, newdist2;
    // Leaf nodes
    if (isleaf) {
        for (int i = 0; i < npts; i++) {
            p2p[0] = point(pts, leaf.p[i])[0] - params[threadNum].p[0];
            p2p[1] = point(pts, leaf.p[i])[1] - params[threadNum].p[1];
            p2p[2] = point(pts, leaf.p[i])[2] - params[threadNum].p[2];
            t = Dot(p2p, params[threadNum].segment_dir);
            if (t < 0.0) {
                // point is beyond point1 of the segment
                if (Dist2(params[threadNum].p,point(pts, leaf.p[i])) >= params[threadNum].maxdist_d2)
                    continue;
            } else if (t > params[threadNum].segment_len2) {
                // point is beyond point2 of the segment
                if (Dist2(params[threadNum].p0,point(pts, leaf.p[i])) >= params[threadNum].maxdist_d2)
                    continue;
            } else {
                // point is within the segment
                // calculate its projection onto the line
                proj[0] = params[threadNum].p[0] + t*params[threadNum].segment_n[0];
                proj[1] = params[threadNum].p[1] + t*params[threadNum].segment_n[1];
                proj[2] = params[threadNum].p[2] + t*params[threadNum].segment_n[2];
                if (Dist2(proj,point(pts, leaf.p[i])) >= params[threadNum].maxdist_d2)
                    continue;
            }
            newdist2 = Dist2(params[threadNum].p,point(pts, leaf.p[i]));
            if (newdist2 < params[threadNum].closest_d2) {
                params[threadNum].closest_d2 = newdist2;
                params[threadNum].closest = pointparam(pts, leaf.p[i]);
            }
        }
        return;
    }

    // Quick check of whether to abort (weeds out all nodes that are too far
    // away from the first point)
    double approx_dist_bbox =
        std::max(std::max(fabs(params[threadNum].p[0]-node.center[0])-node.dx,
                    fabs(params[threadNum].p[1]-node.center[1])-node.dy),
                fabs(params[threadNum].p[2]-node.center[2])-node.dz);
    if (approx_dist_bbox >= 0 &&
            sqr(approx_dist_bbox) >= params[threadNum].closest_d2)
        return;
    // Slower check of whether to abort (weeds out all nodes that are not in
    // the area to search)
    p2p[0] = node.center[0] - params[threadNum].p[0];
    p2p[1] = node.center[1] - params[threadNum].p[1];
    p2p[2] = node.center[2] - params[threadNum].p[2];
    t = Dot(p2p, params[threadNum].segment_dir);
    if (t < 0.0) {
        // point is beyond point1 of the segment
        if (Dist2(params[threadNum].p,node.center) > sqr(node.r+params[threadNum].maxdist_d))
            return;
    } else if (t > params[threadNum].segment_len2) {
        // point is beyond point2 of the segment
        if (Dist2(params[threadNum].p0,node.center) > sqr(node.r+params[threadNum].maxdist_d))
            return;
    } else {
        // point is within the segment
        // calculate projection
        proj[0] = params[threadNum].p[0] + t*params[threadNum].segment_n[0];
        proj[1] = params[threadNum].p[1] + t*params[threadNum].segment_n[1];
        proj[2] = params[threadNum].p[2] + t*params[threadNum].segment_n[2];
        if (Dist2(proj,node.center) > sqr(node.r+params[threadNum].maxdist_d))
            return;
    }

    // Recursive case
    double myd = node.splitval - params[threadNum].p[node.splitaxis];
    if (myd >= 0.0) {
      node.child1->_segmentSearch_1NearestPoint(pts, threadNum);
      if (sqr(myd) < params[threadNum].closest_d2) {
        node.child2->_segmentSearch_1NearestPoint(pts, threadNum);
      }
    } else {
      node.child2->_segmentSearch_1NearestPoint(pts, threadNum);
      if (sqr(myd) < params[threadNum].closest_d2) {
        node.child1->_segmentSearch_1NearestPoint(pts, threadNum);
      }
    }
  }
};


#endif
