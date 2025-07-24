/*
 * BruteForce Search implementation
 *
 * Copyright (C) Simon Engel, Michael Jechow
 *
 * Released under the GPL version 3.
 *
 */

#include "slam6d/BruteForceNotATree.h"
#include "slam6d/globals.icc"
#include <limits.h>

/**
 * Constructor
 *
 * copy given points to new array
 *
 * @param pts 3D array of points
 * @param n number of points
 */
BruteForceNotATree::BruteForceNotATree(double **_pts, size_t _n) {
  pts = (double **)new double[_n][3]; // stores all 3d points
  n   = _n;  // number of points

  // copy _pts into pts
  for (int i = 0; i < n; i++) {
    double *tmp = new double[3];  // create new 3d point
    tmp[0]      = _pts[i][0];
    tmp[1]      = _pts[i][1];
    tmp[2]      = _pts[i][2];
    pts[i]      = tmp;  // copy 3d point to pts
  }
}

/**
 * Destructor
 */
BruteForceNotATree::~BruteForceNotATree() {
  for (int i = 0; i < n; i++) {
    delete[] pts[i];
  }
  delete[] pts;
}

/**
 * Finds the closest point with brute force
 * wrt. the point given as first parameter.
 * @param _p point
 * @param maxdist2 maximal search distance.
 * @param threadNum Thread number, for parallelization
 * @return Pointer to the closest point
 */
double *BruteForceNotATree::FindClosest(double *_p, double maxdist2, int threadNum) const {
  size_t idx      = 0;
  double min_dist = __DBL_MAX__;

  for (int i = 0; i < n; i++) {
    // Calculate distance (without taking square root)
    double distance = (pts[i][0] - _p[0]) * (pts[i][0] - _p[0]) + (pts[i][1] - _p[1]) * (pts[i][1] - _p[1]) + (pts[i][2] - _p[2]) * (pts[i][2] - _p[2]);
    // if a shorter distance is found, update the index and distance of the closest point
    if (distance < min_dist) {
      min_dist = distance;
      idx      = i;
    }
  }
  
  if (Dist2(_p, pts[idx]) > maxdist2) return 0;
  return pts[idx];
}
