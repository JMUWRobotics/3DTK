/*
 * kd implementation
 *
 * Copyright (C) by the 3DTK contributors
 *
 * Released under the GPL version 3.
 *
 */

/** @file
 *  @brief An optimized k-d tree implementation
 *  @author Remus Dumitru. Jacobs University Bremen, Germany
 *  @author Corneliu-Claudiu Prodescu. Jacobs University Bremen, Germany
 *  @author Andreas Nuechter. Jacobs University Bremen, Germany.
 *  @author Kai Lingemann. Inst. of CS, University of Osnabrueck, Germany.
 *  @author Thomas Escher Inst. of CS, University of Osnabrueck, Germany.
 *  @author Fabian Arzberger. Julius-Maximilian University Wuerzburg, Germany.
 */

#ifdef _MSC_VER
#define  _USE_MATH_DEFINES
#endif

#include "slam6d/kd.h"
#include "slam6d/globals.icc"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

// KDtree class static variables
template<class PointData, class AccessorData, class AccessorFunc, class PointType, class ParamFunc>
KDParams<PointType> KDTreeImpl<PointData, AccessorData, AccessorFunc, PointType, ParamFunc>::params[MAX_OPENMP_NUM_THREADS];

/**
 * Constructor
 *
 * Create a KD tree from the points pointed to by the array pts
 *
 * @param pts 3D array of points
 * @param n number of points
 */
KDtree::KDtree(double **pts, int n, int bucketSize)
{
    create(Void(), pts, n, bucketSize);
}

KDtree::~KDtree()
{
}

std::vector<double*> KDtree::CollectPts(int threadNum) const
{
    params[threadNum].collected_pts = std::vector<double*>(0);
    _CollectPts(Void(), threadNum);
    return params[threadNum].collected_pts;
}

int KDtree::Remove(double *_p, int threadNum)
{
    params[threadNum].closest = 0;
    params[threadNum].closest_d2 = __DBL_MAX__;
    params[threadNum].p = _p;
    return _Remove(Void(), threadNum);
}

/**
 * Finds the closest point within the tree,
 * wrt. the point given as first parameter.
 * @param _p point
 * @param maxdist2 maximal search distance.
 * @param threadNum Thread number, for parallelization
 * @return Pointer to the closest point
 */
double *KDtree::FindClosest(double *_p,
                            double maxdist2,
                            int threadNum) const
{
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  _FindClosest(Void(), threadNum);
  return params[threadNum].closest;
}

double *KDtree::FindClosestAlongDir(double *_p,
                                    double *_dir,
                                    double maxdist2,
                                    int threadNum) const
{
  params[threadNum].closest = NULL;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  params[threadNum].dir = _dir;
  _FindClosestAlongDir(Void(), threadNum);
  return params[threadNum].closest;
}

std::vector<Point> KDtree::kNearestNeighbors(double *_p,
                                        int _k,
                                        int threadNum) const
{
  std::vector<Point> result;
  params[threadNum].closest = 0;
  params[threadNum].p = _p;
  params[threadNum].k = _k;
  // todo fix this C/C++ mixture
  params[threadNum].closest_neighbors = (double **)calloc(_k,
                                                          sizeof(double *) );
  params[threadNum].distances = (double *)calloc(_k,
                                                 sizeof(double));
  // initialize distances to an invalid value to indicate unset neighbors
  for (int i = 0; i < _k; i++) {
      params[threadNum].distances[i] = -1.0;
  }

  _KNNSearch(Void(), threadNum);

  for (int i = 0; i < _k; i++) {
      // only push valid points
    if (params[threadNum].distances[i] >= 0.0f) {
    result.push_back(Point(params[threadNum].closest_neighbors[i][0],
                           params[threadNum].closest_neighbors[i][1],
                           params[threadNum].closest_neighbors[i][2]));
    }
  }

  free (params[threadNum].closest_neighbors);
  free (params[threadNum].distances);

  return result;
}

std::vector<Point> KDtree::kNearestRangeSearch(double *_p,
                                       int _k,
                                       double sqRad2,
                                       int threadNum) const
{
  std::vector<Point> result;
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = sqRad2;
  params[threadNum].p = _p;
  params[threadNum].k = _k;
  // todo fix this C/C++ mixture
  params[threadNum].closest_neighbors = (double **)calloc(_k,
                                                          sizeof(double *) );
  params[threadNum].distances = (double *)calloc(_k,
                                                 sizeof(double));
  // initialize distances to an invalid value to indicate unset neighbors
  for (int i = 0; i < _k; i++) {
      params[threadNum].distances[i] = -1.0;
  }
  _KNNRangeSearch(Void(), threadNum);

  for (int i = 0; i < _k; i++) {
      // only push valid points
    if (params[threadNum].distances[i] >= 0.0f) {
    result.push_back(Point(params[threadNum].closest_neighbors[i][0],
                           params[threadNum].closest_neighbors[i][1],
                           params[threadNum].closest_neighbors[i][2]));
    }
  }

  free (params[threadNum].closest_neighbors);
  free (params[threadNum].distances);

  return result;
}


std::vector<Point> KDtree::fixedRangeSearchBetween2Points(double *_p,
                      double *_p0,
                      double maxdist2,
                      int threadNum) const {
  std::vector<Point> result;
  params[threadNum].closest = _p0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  params[threadNum].dist = sqrt(Dist2(_p, _p0));

  double * _dir = new double[3];
  for(int i = 0; i < 3; i++) {
    _dir[i] = _p0[i] - _p[i];
  }

  Normalize3(_dir);

  params[threadNum].dir = _dir;
  params[threadNum].range_neighbors.clear();

  _fixedRangeSearchBetween2Points(Void(), threadNum);

  for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(Point(params[threadNum].range_neighbors[i][0],
                           params[threadNum].range_neighbors[i][1],
                           params[threadNum].range_neighbors[i][2]));
  }

  delete[] _dir;
  return result;
}


std::vector<Point> KDtree::fixedRangeSearchAlongDir(double *_p,
                      double *_dir,
                      double maxdist2,
                      int threadNum) const {
  std::vector<Point> result;
  params[threadNum].closest = NULL;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  params[threadNum].dir = _dir;
  params[threadNum].range_neighbors.clear();

  _fixedRangeSearchAlongDir(Void(), threadNum);

  for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(Point(params[threadNum].range_neighbors[i][0],
                           params[threadNum].range_neighbors[i][1],
                           params[threadNum].range_neighbors[i][2]));
  }

  return result;
}

std::vector<Point> KDtree::fixedRangeSearch(double *_p,
                                       double sqRad2,
                                       int threadNum) const
{
  std::vector<Point> result;
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = sqRad2;
  params[threadNum].p = _p;
  params[threadNum].range_neighbors.clear();
  _FixedRangeSearch(Void(), threadNum);

  for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(Point(params[threadNum].range_neighbors[i][0],
                           params[threadNum].range_neighbors[i][1],
                           params[threadNum].range_neighbors[i][2]));
  }

  return result;
}

std::vector<Point> KDtree::AABBSearch(double *_p,
                                 double* _p0,
                                 int threadNum) const
{
    if (_p[0] > _p0[0] || _p[1] > _p0[1] || _p[2] > _p0[2])
        throw std::logic_error("invalid bbox");
    std::vector<Point> result;
    params[threadNum].p = _p;
    params[threadNum].p0 = _p0;
    params[threadNum].range_neighbors.clear();
    _AABBSearch(Void(), threadNum);

    for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(Point(params[threadNum].range_neighbors[i][0],
                           params[threadNum].range_neighbors[i][1],
                           params[threadNum].range_neighbors[i][2]));
    }

    return result;
}

double *KDtree::segmentSearch_1NearestPoint(double *_p,
          double* _p0, double maxdist2, int threadNum) const
{
  params[threadNum].closest = 0;
  // the furthest a point can be away is the distance between the points
  // making the line segment plus maxdist
  params[threadNum].closest_d2 = sqr(sqrt(Dist2(_p,_p0))+sqrt(maxdist2));
  params[threadNum].maxdist_d2 = maxdist2;
  params[threadNum].maxdist_d = sqrt(maxdist2);
  params[threadNum].p = _p;
  params[threadNum].p0 = _p0;
  double *dir = new double[3]{_p0[0] - _p[0], _p0[1] - _p[1], _p0[2] - _p[2] };
  double len2 = Len2(dir);
  double *n = new double[3]{dir[0]/len2,dir[1]/len2,dir[2]/len2};
  params[threadNum].segment_dir = dir;
  params[threadNum].segment_len2 = len2;
  params[threadNum].segment_n = n;
  _segmentSearch_1NearestPoint(Void(), threadNum);
  delete[] dir;
  delete[] n;
  return params[threadNum].closest;
}
