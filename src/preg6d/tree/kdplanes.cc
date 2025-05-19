#include "tree/kdplanes.h"

// KDtree class static variables
template<class PointData, class AccessorData, class AccessorFunc, class PointType, class ParamFunc>
KDParams<PointType> KDTreeImpl<PointData, AccessorData, AccessorFunc, PointType, ParamFunc>::params[MAX_OPENMP_NUM_THREADS];


KDtreePlanes::KDtreePlanes(PointPlane** pts, size_t n, int bucketSize)
{
    create(Void(), pts, n, bucketSize);
}

KDtreePlanes::KDtreePlanes(Planes planes, int bucketSize)
{
    int n = 0;
    for (const NormalPlane* plane : planes)
        n += plane->all_pts.size();

    PointPlane** pts = new PointPlane*[n];
    int j = 0;
    for (NormalPlane* plane : planes)
    {
        for(size_t i = 0; i < plane->all_pts.size(); ++i)
            pts[j++] = new PointPlane( plane->all_pts[i] , plane );
    }

    create(Void(), pts, n, bucketSize);
}

KDtreePlanes::~KDtreePlanes(){}

std::vector<PointPlane*> KDtreePlanes::CollectPts(int threadNum) const
{
    params[threadNum].collected_pts.clear();
    _CollectPts(Void(), threadNum);
    return params[threadNum].collected_pts;
}

int KDtreePlanes::Remove(double *_p, int threadNum)
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
PointPlane *KDtreePlanes::FindClosest(double *_p,
                            double maxdist2,
                            int threadNum) const
{
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  _FindClosest(Void(), threadNum);
  return params[threadNum].closest;
}

PointPlane *KDtreePlanes::FindClosestAlongDir(double *_p,
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

std::vector<PointPlane*> KDtreePlanes::kNearestNeighbors(double *_p,
                                        int _k,
                                        int threadNum) const
{
  std::vector<PointPlane*> result;
  params[threadNum].closest = 0;
  params[threadNum].p = _p;
  params[threadNum].k = _k;
  // todo fix this C/C++ mixture
  params[threadNum].closest_neighbors = (PointPlane **)calloc(_k,
                                                          sizeof(PointPlane *) );
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
        result.push_back(params[threadNum].closest_neighbors[i]);
    }
  }

  free (params[threadNum].closest_neighbors);
  free (params[threadNum].distances);

  return result;
}

std::vector<PointPlane*> KDtreePlanes::kNearestRangeSearch(double *_p,
                                       int _k,
                                       double sqRad2,
                                       int threadNum) const
{
  std::vector<PointPlane*> result;
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = sqRad2;
  params[threadNum].p = _p;
  params[threadNum].k = _k;
  // todo fix this C/C++ mixture
  params[threadNum].closest_neighbors = (PointPlane **)calloc(_k,
                                                          sizeof(PointPlane *) );
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
        result.push_back(params[threadNum].closest_neighbors[i]);
    }
  }

  free (params[threadNum].closest_neighbors);
  free (params[threadNum].distances);

  return result;
}

std::vector<PointPlane*> KDtreePlanes::fixedRangeSearchBetween2Points(double *_p,
                      double *_p0,
                      double maxdist2,
                      int threadNum) const {
  std::vector<PointPlane*> result;
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  params[threadNum].p0 = _p0;
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
    result.push_back(params[threadNum].range_neighbors[i]);
  }

  delete[] _dir;
  return result;
}

std::vector<PointPlane*> KDtreePlanes::fixedRangeSearchAlongDir(double *_p,
                      double *_dir,
                      double maxdist2,
                      int threadNum) const {
  std::vector<PointPlane*> result;
  params[threadNum].closest = NULL;
  params[threadNum].closest_d2 = maxdist2;
  params[threadNum].p = _p;
  params[threadNum].dir = _dir;
  params[threadNum].range_neighbors.clear();

  _fixedRangeSearchAlongDir(Void(), threadNum);

  for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(params[threadNum].range_neighbors[i]);
  }

  return result;
}

std::vector<PointPlane*> KDtreePlanes::fixedRangeSearch(double *_p,
                                       double sqRad2,
                                       int threadNum) const
{
  std::vector<PointPlane*> result;
  params[threadNum].closest = 0;
  params[threadNum].closest_d2 = sqRad2;
  params[threadNum].p = _p;
  params[threadNum].range_neighbors.clear();
  _FixedRangeSearch(Void(), threadNum);

  for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
    result.push_back(params[threadNum].range_neighbors[i]);
  }

  return result;
}

std::vector<PointPlane*> KDtreePlanes::AABBSearch(double *_p,
                                 double* _p0,
                                 int threadNum) const
{
    if (_p[0] > _p0[0] || _p[1] > _p0[1] || _p[2] > _p0[2])
        throw std::logic_error("invalid bbox");
    std::vector<PointPlane*> result;
    params[threadNum].p = _p;
    params[threadNum].p0 = _p0;
    params[threadNum].range_neighbors.clear();
    _AABBSearch(Void(), threadNum);

    for (size_t i = 0; i < params[threadNum].range_neighbors.size(); i++) {
        result.push_back(params[threadNum].range_neighbors[i]);
    }

    return result;
}

PointPlane *KDtreePlanes::segmentSearch_1NearestPoint(double *_p,
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