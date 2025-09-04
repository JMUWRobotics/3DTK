#ifndef SHOWBOCTREE_H
#define SHOWBOCTREE_H


#include "slam6d/Boctree.h"
#include "linescan.h"
#include "lnnparams.h"

extern char amap[8][8];
extern char sequence2ci[8][256][8];

template <class T> class LBOctTree : public BOctTree<T> {

public:

  LBOctTree(vector<LineScan*> &linescans, T voxelSize ) {
    this->voxelSize = voxelSize;

    this->POINTDIM = 4; // 3 plus 1 for the linescan index

    BOctTree<T>::mins = new T[BOctTree<T>::POINTDIM];
    BOctTree<T>::maxs = new T[BOctTree<T>::POINTDIM];


    unsigned int n = 0;
    for (unsigned int j = 0; j < linescans.size(); j++) {
      n += linescans[j]->nrpts;
    }

    T **pts = new T*[n];

    unsigned int cn = 0;
    for (unsigned int j = 0; j < linescans.size(); j++) {
      for (int k = 0; k < linescans[j]->nrpts; k++) {
        pts[cn] = new T[BOctTree<T>::POINTDIM];
        pts[cn][0] = linescans[j]->points[k][0];
        pts[cn][1] = linescans[j]->points[k][1];
        pts[cn][2] = linescans[j]->points[k][2];
        pts[cn][3] = j;
        ::transform(pts[cn], linescans[j]->transMat);
        cn++;
      }
    }

    // initialising
    for (unsigned int i = 0; i < BOctTree<T>::POINTDIM; i++) {
      BOctTree<T>::mins[i] = pts[0][i];
      BOctTree<T>::maxs[i] = pts[0][i];
    }

    for (unsigned int i = 0; i < BOctTree<T>::POINTDIM; i++) {
      for (unsigned int j = 1; j < n; j++) {
        BOctTree<T>::mins[i] = min(BOctTree<T>::mins[i], (T)pts[j][i]);
        BOctTree<T>::maxs[i] = max(BOctTree<T>::maxs[i], (T)pts[j][i]);
      }
    }

    BOctTree<T>::center[0] = 0.5 * (BOctTree<T>::mins[0] + BOctTree<T>::maxs[0]);
    BOctTree<T>::center[1] = 0.5 * (BOctTree<T>::mins[1] + BOctTree<T>::maxs[1]);
    BOctTree<T>::center[2] = 0.5 * (BOctTree<T>::mins[2] + BOctTree<T>::maxs[2]);
    BOctTree<T>::size = max(max(0.5 * (BOctTree<T>::maxs[0] - BOctTree<T>::mins[0]), 0.5 * (BOctTree<T>::maxs[1] - BOctTree<T>::mins[1])), 0.5 * (BOctTree<T>::maxs[2] - BOctTree<T>::mins[2]));
    BOctTree<T>::size += 1.0; // for numerical reasons we increase size

    // calculate new buckets
    T newcenter[8][3];
    T sizeNew = BOctTree<T>::size / 2.0;

    for (unsigned char i = 0; i < 8; i++) {
      BOctTree<T>::childcenter(BOctTree<T>::center, newcenter[i], BOctTree<T>::size, i);
    }
    // set up values
    BOctTree<T>::uroot = new bitunion<T>();
    BOctTree<T>::root = &BOctTree<T>::uroot->node;

    BOctTree<T>::countPointsAndQueueFast(pts, n, newcenter, sizeNew, *BOctTree<T>::root, BOctTree<T>::center);
    BOctTree<T>::init();

    for (unsigned int i = 0; i < n; i++) {
      delete[] pts[i];
    }
    delete[] pts;

  }


  inline void findClosestInLeaf(bitunion<T> *node, int threadNum) {
    if (lparams[threadNum].count >= lparams[threadNum].max_count) return;
    lparams[threadNum].count++;
    T* points = node->getPoints();
    unsigned int length = node->getLength();
    for(unsigned int iterator = 0; iterator < length; iterator++ ) {

      if ( ( points[3] >= lparams[threadNum].begin1 && points[3] <= lparams[threadNum].end1 ) ||
           ( points[3] >= lparams[threadNum].begin2 && points[3] <= lparams[threadNum].end2 ) )

      {
        double myd2 = Dist2(lparams[threadNum].p, points);
        if (myd2 < lparams[threadNum].closest_d2) {
          lparams[threadNum].closest_d2 = myd2;
          lparams[threadNum].closest = points;
          if (myd2 <= 0.0001) {
            lparams[threadNum].closest_v = 0; // the search radius in units of voxelSize
          } else {
            lparams[threadNum].closest_v = sqrt(myd2) * BOctTree<T>::mult+1; // the search radius in units of voxelSize
          }
        }
      }
      points+=BOctTree<T>::POINTDIM;
    }
  }

  int FindClosestPoints(LineScan *l, double maxdist2, int mindist, int maxdist, int threadNum) {
    int count = 0;

      lparams[threadNum].begin1 = l->getIndex() - maxdist;
      lparams[threadNum].end1 = l->getIndex() - mindist;
      lparams[threadNum].begin2 = l->getIndex() + mindist;
      lparams[threadNum].end2 = l->getIndex() + maxdist;
      lparams[threadNum].closest_v = sqrt(maxdist2) * BOctTree<T>::mult+1; // the search radius in units of voxelSize

    for (unsigned int i = 0; i < l->nrQuery(); i++) {
      // TODO implement without copies of the points
      double *query = new double[3];

      query[0] = l->getQuery(i)[0];
      query[1] = l->getQuery(i)[1];
      query[2] = l->getQuery(i)[2];



      double *p = FindClosest(query, maxdist2, threadNum);
      if (p) {
        unsigned int index = p[3]; // get index of the linescan of this point
        //lppair pair = lppair(p, l->points[i]);
        double *closest = new double[3];
        closest[0] = p[0];
        closest[1] = p[1];
        closest[2] = p[2];
        //closest[3] = index;

        lppair pair = lppair(closest, query);
        l->addPtPair(pair, index);
        count++;
      }

    }
    return count;
  }

  /**
 * This function finds the closest point in the octree given a specified
 * radius. This implementation is quit complex, although it is already
 * simplified. The simplification incurs a significant loss in speed, as
 * several calculations have to be performed repeatedly and a high number of
 * unnecessary jumps are executed.
*/
  double *FindClosest(double *point, double maxdist2, int threadNum)
  {
    lparams[threadNum].closest = 0; // no point found currently
    lparams[threadNum].closest_d2 = maxdist2;
    lparams[threadNum].p = point;
    lparams[threadNum].x = (point[0] + BOctTree<T>::add[0]) * BOctTree<T>::mult;
    lparams[threadNum].y = (point[1] + BOctTree<T>::add[1]) * BOctTree<T>::mult;
    lparams[threadNum].z = (point[2] + BOctTree<T>::add[2]) * BOctTree<T>::mult;
    lparams[threadNum].count = 0;
    lparams[threadNum].max_count = 10000; // stop looking after this many buckets


    // box within bounds in voxel coordinates
    int xmin, ymin, zmin, xmax, ymax, zmax;
    xmin = max(lparams[threadNum].x-lparams[threadNum].closest_v, 0);
    ymin = max(lparams[threadNum].y-lparams[threadNum].closest_v, 0);
    zmin = max(lparams[threadNum].z-lparams[threadNum].closest_v, 0);

    xmax = min(lparams[threadNum].x+lparams[threadNum].closest_v, BOctTree<T>::largest_index);
    ymax = min(lparams[threadNum].y+lparams[threadNum].closest_v, BOctTree<T>::largest_index);
    zmax = min(lparams[threadNum].z+lparams[threadNum].closest_v, BOctTree<T>::largest_index);

    unsigned char depth = 0;
    unsigned int child_bit;
    unsigned int child_index_min;
    unsigned int child_index_max;

    bitunion<T> *node = &(*BOctTree<T>::uroot);

    int cx, cy, cz;

    child_bit = BOctTree<T>::child_bit_depth[depth];
    cx = BOctTree<T>::child_bit_depth[depth];
    cy = BOctTree<T>::child_bit_depth[depth];
    cz = BOctTree<T>::child_bit_depth[depth];

    while (true) { // find the first node where branching is required
      child_index_min = ((xmin & child_bit )!=0)  | (((ymin & child_bit )!=0 )<< 1) | (((zmin & child_bit )!=0) << 2);
      child_index_max = ((xmax & child_bit )!=0)  | (((ymax & child_bit )!=0 )<< 1) | (((zmax & child_bit )!=0) << 2);

      // if these are the same, go there
      // TODO: optimization: also traverse if only single child...
      if (child_index_min == child_index_max) {
        if (node->childIsLeaf(child_index_min) ) {  // luckily, no branching is required
          findClosestInLeaf(node->getChild(child_index_min), threadNum);
          return static_cast<double*>(lparams[threadNum].closest);
        } else {
          if (node->isValid(child_index_min) ) { // only descend when there is a child
            BOctTree<T>::childcenter(cx,cy,cz, cx,cy,cz, child_index_min, child_bit/2 );
            node = node->getChild(child_index_min);
            child_bit /= 2;
          } else {  // there is no child containing the bounding box => no point is close enough
            return 0;
          }
        }
      } else {
        // if min and max are not in the same child we must branch
        break;
      }
    }

    // node contains all box-within-bounds cells, now begin best bin first search
    _FindClosest(threadNum, node->node, child_bit/2, cx, cy, cz);
    return static_cast<double*>(lparams[threadNum].closest);
  }


  /**
   * This is the heavy duty search function doing most of the (theoretically unneccesary) work. The tree is recursively searched.
   * Depending on which of the 8 child-voxels is closer to the query point, the children are examined in a special order.
   * This order is defined in map, imap is its inverse and sequence2ci is a speedup structure for faster access to the child indices.
*/
  void _FindClosest(int threadNum, bitoct &node, int size, int x, int y, int z)
  {
    // Recursive case

    // compute which child is closest to the query point
    unsigned char child_index =  ((lparams[threadNum].x - x) >= 0) |
                                (((lparams[threadNum].y - y) >= 0) << 1) |
                                (((lparams[threadNum].z - z) >= 0) << 2);

    char *seq2ci = sequence2ci[child_index][node.valid];  // maps preference to index in children array
    char *mmap = amap[child_index];  // maps preference to area index

    bitunion<T> *children;
    bitoct::getChildren(node, children);
    int cx, cy, cz;
    cx = cy = cz = 0; // just to shut up the compiler warnings
    for (unsigned char i = 0; i < 8; i++) { // in order of preference
      child_index = mmap[i]; // the area index of the node
      if (  ( 1 << child_index ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(x,y,z, cx,cy,cz, child_index, size);
        if ( lparams[threadNum].closest_v == 0 ||  max(max(abs( cx - lparams[threadNum].x ),
                 abs( cy - lparams[threadNum].y )),
                 abs( cz - lparams[threadNum].z )) - size
        > lparams[threadNum].closest_v ) {
          continue;
        }
        // find the closest point in leaf seq2ci[i]
        if (  ( 1 << child_index ) & node.leaf ) {   // if ith node is leaf
          findClosestInLeaf( &children[seq2ci[i]], threadNum);
        } else { // recurse
          _FindClosest(threadNum, children[seq2ci[i]].node, size/2, cx, cy, cz);
        }
      }
    }
  }


  LNNParams lparams[100];
};


#endif
