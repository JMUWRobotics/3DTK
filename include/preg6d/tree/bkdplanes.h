#ifndef __BKDPLANES_H_
#define __BKDPLANES_H_

#include "tree/kdplanes.h"
#include <vector>

// This data structure stores trees.
// The trees themselfes do not have any pts stored.
// Instead, they operate on the "ForestElem" structure by reorganisation.
struct PForestElem
{
    // Create an empty Forest element, i.e. an empty tree (called 'sprout')
    PForestElem() : empty(true), nrpts(0) {} // empty constructor, we use c++ init list
    // Create a non-empty Forest element, i.e. a KD-tree (called 'tree')
    PForestElem(KDtreePlanes* t, size_t n) : empty(false), tree(t), nrpts(n) {}

    // Struct data
    bool empty;
    KDtreePlanes* tree;
    size_t nrpts;
};

class BkdTreePlanes : public KDtreePlanes
{
public:
    BkdTreePlanes(PointPlane** pts, size_t n, int bucketSize = 20);

    BkdTreePlanes(Planes planes, int bucketSize = 20);

    BkdTreePlanes(int bucketSize = 20);

    virtual ~BkdTreePlanes();

    // BKD specifics

    void insert(PointPlane*, int threadNum = 0);

    void insert(vector<PointPlane*> , int threadNum = 0);

    size_t size() const;

    std::string info() const;

    std::string _debug_info() const;

    friend std::ostream& operator<<(std::ostream&, const BkdTreePlanes&);

    // Inherited from KDtreePlanes

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

private:
        // Saves all the sub-trees (thus is called 'forest')
    std::vector<PForestElem> forest;
    std::vector<PointPlane*> buffer;
    int bucketSize;

    void mergeTreesLogarithmic(int, int, int threadNum = 0);

};

#endif // __BKDPLANES_H_