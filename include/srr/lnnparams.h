#ifndef __LNNPARAMS_H__
#define __LNNPARAMS_H__

struct LNNParams {
/**
   * pointer to the closest point.  size = 4 bytes of 32 bit machines
   */
  void *closest;

  /**
   * distance to the closest point. size = 8 bytes
   */
  double closest_d2;

  // distance to the closest point in voxels
  int closest_v;

  // location of the query point in voxel coordinates
  int x;
  int y;
  int z;

  /**
   * pointer to the point, size = 4 bytes of 32 bit machines
   */
  double *p;

  int begin1;
  int end1;
  int begin2;
  int end2;

  int count;
  int max_count;

};

#endif
