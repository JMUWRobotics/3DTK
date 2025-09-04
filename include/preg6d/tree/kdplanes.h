#ifndef __KDPLANES_H_
#define __KDPLANES_H_

#include "slam6d/kdparams.h"
#include "slam6d/searchTree.h"
#include "slam6d/kdTreeImpl.h"

#include "model/normalplane.h"

/*
 * Accociation between a point and corresponding plane pointer.
 * Usually, struct instances of this have 16 bytes on a 64bit operating system.
 */
struct PointPlane {
public:
    double* point;
    NormalPlane* plane;
    PointPlane(double *pt, NormalPlane *pl) : point(pt), plane(pl) {}
};

struct PointAccessor {
    inline double *operator() (Void, PointPlane* index) {
        return index->point;
    }
};

struct PTypeAccessor {
    inline PointPlane *operator() (Void, PointPlane* index) {
        return index;
    }
};

class KDtreePlanes : protected KDTreeImpl<Void, PointPlane*, PointAccessor, PointPlane*, PTypeAccessor>
{
protected:
    KDtreePlanes() {};

public:
    KDtreePlanes(PointPlane** pts, size_t n, int bucketSize = 20);

    KDtreePlanes(Planes planes, int bucketSize = 20);

    virtual ~KDtreePlanes();

    virtual int Remove(double *_p,
                     int threadNum = 0); // no const

    virtual std::vector<PointPlane*> CollectPts(int threadNum = 0) const;

    virtual PointPlane* FindClosest(double *_p,
                            double maxdist2,
                            int threadNum = 0) const;

    virtual PointPlane* FindClosestAlongDir(double *_p,
                                double *_dir,
                                double maxdist2,
                                int threadNum = 0) const;

    virtual std::vector<PointPlane*> fixedRangeSearchAlongDir(double *_p,
                    double *_dir,
                    double maxdist2,
                    int threadNum = 0) const;

    virtual std::vector<PointPlane*> fixedRangeSearchBetween2Points(double *_p,
                    double *_dir,
                    double maxdist2,
                    int threadNum = 0) const;

    virtual std::vector<PointPlane*> kNearestNeighbors(double *_p,
                                    int k,
                                    int threadNum = 0) const;

    virtual std::vector<PointPlane*> kNearestRangeSearch(double *_p,
								 int k,
								 double sqRad2,
								 int threadNum = 0) const;

    virtual std::vector<PointPlane*> fixedRangeSearch(double *_p,
                                    double sqRad2,
                                    int threadNum = 0) const;

    virtual std::vector<PointPlane*> AABBSearch(double *_p,
                                    double* _p0,
                                    int threadNum = 0) const;

    virtual PointPlane* segmentSearch_1NearestPoint(double *_p,
            double* _p0, double maxdist2, int threadNum) const;

};

#endif // __KDPLANES_H_

