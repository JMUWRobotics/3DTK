#include "tree/bkdplanes.h"
#include <list>

BkdTreePlanes::BkdTreePlanes(int bucketsize)
{
    // Init buffer mem
    this->bucketSize = bucketsize;
    buffer.reserve(bucketsize); // keeps size at 0, reserves memory
}

BkdTreePlanes::BkdTreePlanes(PointPlane** pts, size_t n, int bucketsize)
{
    // nr of trees needed for n points with given bucketsize
    int nrtrees = std::max(1.0, ceil(log2(n/bucketsize)));

    // Init forrest
    /*
     * this way, all points would fit into the forest.
     * however, we add another tree and put all the points into
     * that large tree, thus all trees before are empty.
     * This way, point insertion starts optimaly
     */
    for(int i = 0; i < nrtrees; ++i) {
        PForestElem sprout;
        forest.push_back( sprout );
    }
    PForestElem tree(new KDtreePlanes(pts, n, bucketsize), n);
    forest.push_back( tree );
}

BkdTreePlanes::BkdTreePlanes(Planes planes, int bucketsize)
{
    int n = 0;
    for (const NormalPlane* plane : planes)
        n += plane->all_pts.size();

    // nr of trees needed for n points with given bucketsize
    int nrtrees = ceil(log2(n/bucketsize));
    for(int i = 0; i < nrtrees; ++i) {
        PForestElem sprout;
        forest.push_back( sprout );
    }
    PForestElem tree(new KDtreePlanes(planes, bucketsize), n);
    forest.push_back( tree );
}

BkdTreePlanes::~BkdTreePlanes()
{
    for (uint i = 0; i < forest.size(); ++i)
        forest.at(i).tree->~KDtreePlanes();
    forest.clear();
    buffer.clear();
}

std::vector<PointPlane*> BkdTreePlanes::CollectPts(int threadNum) const
{
    std::vector<PointPlane*> res;
    for (size_t i = 0; i < forest.size(); ++i) {
        if (forest.at(i).empty) continue;
        std::vector<PointPlane*> tree_res = forest.at(i).tree->CollectPts(threadNum);
        res.insert( std::end(res), std::begin(tree_res), std::end(tree_res) );
    }
    if (buffer.size() > 0)
        res.insert( std::end(res), std::begin(buffer), std::end(buffer) );
    return res;
}

/**
 * Neat for debugging.
 * Shows size of the bkd-tree, i.e. all sub-trees and how many
 * points they include (the buffer also).
 */
std::string BkdTreePlanes::info() const
{
    std::string info_str;
    std::stringstream info_ss;
    info_ss << "============= BKD INFO ===============" << std::endl;
    info_ss << "Buffer [" << std::to_string(buffer.size())
        << "/" << std::to_string(bucketSize) << "]" << std::endl;
    for (size_t i = 0; i < forest.size(); ++i) {
        info_ss << "Tree " << std::to_string(i) << ": "
            << std::to_string(forest.at(i).nrpts) << " pts." << std::endl;
    }
    info_ss << "======================================" << std::endl;
    info_str = info_ss.str();
    return info_str;
}

/**
 * Even more neat for debugging.
 * Shows size of the bkd-tree, i.e. all sub-trees and how many
 * points they include (the buffer also).
 * Furthermore, show all points inside the tree (use with care!)
 */
std::string BkdTreePlanes::_debug_info() const
{
    std::string info_str;
    std::stringstream info_ss;

    info_ss << "============== BKD DEBUG ==============" << std::endl;
    info_ss << "Buffer [" << std::to_string(buffer.size())
        << "/" << std::to_string(bucketSize) << "]" << std::endl;
    for (size_t i = 0; i < buffer.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            info_ss << buffer[i]->point[j] << " ";
        info_ss << std::endl;
    }

    for (size_t i = 0; i < forest.size(); ++i) {
        info_ss << "Tree " << std::to_string(i) << ": "
            << std::to_string(forest.at(i).nrpts) << " pts." << std::endl;
        if (forest.at(i).empty) continue;
        KDtreePlanes *t = forest.at(i).tree;
        std::vector<PointPlane*> data = t->CollectPts();

        for (size_t j = 0; j < data.size(); ++j) {
            for (int k = 0; k < 3; ++k)
                info_ss << data.at(j)->point[k] << " ";
            info_ss << std::endl;
        }
    }
    info_ss << "========================================" << std::endl;
    info_str = info_ss.str();
    return info_str;
}

size_t BkdTreePlanes::size() const
{
    size_t res = 0;
    res += buffer.size();
    for (size_t i = 0; i < forest.size(); ++i)
        if(!forest.at(i).empty)
            res += forest.at(i).nrpts;
    return res;
}

std::ostream& operator<<(std::ostream &out, const BkdTreePlanes &t)
{
    out << t.info().c_str();
    return out;
}

void BkdTreePlanes::insert(vector<PointPlane*> pts, int threadNum)
{
    for (PointPlane* p : pts)
        insert(p, threadNum);
}

/**
 *  @brief Insert point <pt> into the Bkd tree.
 *  Keeps the tree balanced and fully space utilized.
 */
void BkdTreePlanes::insert(PointPlane *pt, int threadNum)
{
    // Push point to the buffer
    #pragma omp critical
    buffer.push_back(pt);

    int b1 = buffer.size();
    int b2 = bucketSize;
    if (buffer.size() != bucketSize) return;
    // If the buffer is full
    if ( b1 == b2 ) {

        int nrnewpts = buffer.size();
        bool merged = false;

        // Search first free tree:
        // Iterate the forest
        for (int i = 0; i < forest.size(); ++i) {

            // Until you find an empty tree
            if ( forest.at(i).empty ) {
                // Merge all trees into the empty tree:
                mergeTreesLogarithmic(i, nrnewpts, threadNum);
                merged = true;
                break;
            }
            // If current tree is not empty, iterate further and store nrpts
            else
            {
                int npt = forest.at(i).nrpts;
                nrnewpts += npt;
            }
        }

        // All trees in the forest are already full
        if (!merged)
        {
            // add new elements
            int lastIndex = forest.size();
            PForestElem sprout; // just an empty dummy
            forest.push_back(sprout);
            // use logarithmic merging to get all prev trees into the new one
            mergeTreesLogarithmic( lastIndex, nrnewpts, threadNum );
        }

        // Clear buffer when point is inserted
        buffer.clear();
        //std::cout << "Cleared Buffer and finish." << std::endl;
        return;
    }
}

int BkdTreePlanes::Remove(double* pt, int threadNum)
{
    // Search if the point is inside the buffer
    for (size_t i = 0; i < buffer.size(); ++i) {
        double d2 = Dist2(pt, buffer[i]->point);
        if (d2 < 0.000000001) {
            buffer.erase( buffer.begin() + i );
            return 1;
        }
    }
    // Traverse the trees if the point is not inside the buffer
    for (size_t i = 0; i < forest.size(); ++i) {
        if (!forest.at(i).empty) {
            if (forest.at(i).tree->Remove(pt, threadNum)) {
                forest.at(i).nrpts--;
                if(forest.at(i).nrpts == 0)
                    forest.at(i).empty = true;
                return 1;
            }
        }
    }
    return 0;
}

PointPlane * BkdTreePlanes::FindClosest(double *_p,
                                double maxdist2,
                                int threadNum) const
{
    PointPlane *res = 0;
    PointPlane *tmp;
    double minDist2 = __DBL_MAX__;
    double dist2;
    // Search all the trees, pick best closest
    for (size_t i = 0; i < forest.size(); ++i) {

        // Skip the empty trees.
        if (forest.at(i).empty) continue;
        KDtreePlanes *tree = forest.at(i).tree;
        tmp = tree->FindClosest(_p, maxdist2, threadNum);
        dist2 = Dist2( tmp->point , _p );
        if ( dist2 < minDist2 ) {
            minDist2 = dist2;
            res = tmp;
        }
    }
    // Also search the buffer if there is any nearer point
    for (size_t i = 0; i < buffer.size(); ++i) {
        dist2 = Dist2( buffer[i]->point , _p);
        if (dist2 < minDist2) {
            minDist2 = dist2;
            res = buffer[i];
        }
    }
    return res;
}

PointPlane* BkdTreePlanes::FindClosestAlongDir(double *_p,
                                double *_dir,
                                double maxdist2,
                                int threadNum) const
{
    std::vector<PointPlane*> result;
    for ( const PForestElem& t : forest ) {
        PointPlane* tmp = t.tree->FindClosestAlongDir(_p, _dir, maxdist2, threadNum);
        result.push_back(tmp);
    }

    double min_dist = __DBL_MAX__, dist;
    PointPlane* min_p = 0;
    for(PointPlane* p : result)
    {
        dist = sqr(p->point[0] -_p[0])+sqr(p->point[1] -_p[1])+sqr(p->point[2] -_p[2]);
        if ( dist < min_dist )
        {
            min_dist = dist;
            min_p = p;
        }
    }
    return min_p;
}

std::vector<PointPlane*> BkdTreePlanes::fixedRangeSearchAlongDir(double *_p,
                                double *_dir,
                                double maxdist2,
                                int threadNum) const
{
    std::vector<PointPlane*> result;
    for ( const PForestElem& t : forest ) {
        std::vector<PointPlane*> tmp = t.tree->fixedRangeSearchAlongDir(_p, _dir, maxdist2, threadNum);
        result.insert( end(result), begin(tmp), end(tmp) );
    }
    return result;
}

std::vector<PointPlane*> BkdTreePlanes::fixedRangeSearchBetween2Points(double *_p,
                                double *_dir,
                                double maxdist2,
                                int threadNum) const
{
    std::vector<PointPlane*> result;
    for ( const PForestElem& t : forest ) {
        std::vector<PointPlane*> tmp = t.tree->fixedRangeSearchBetween2Points(_p, _dir, maxdist2, threadNum);
        result.insert( end(result), begin(tmp), end(tmp) );
    }
    return result;
}

std::vector<PointPlane*> BkdTreePlanes::kNearestNeighbors(double *_p,
								int k,
								int threadNum) const
{
    std::list<PointPlane*> candidates;
    // Include k nearest from all the trees
    for(size_t i = 0; i < forest.size(); ++i) {
        if (forest.at(i).empty) continue;
        KDtreePlanes *t = forest.at(i).tree;
        std::vector<PointPlane*> kNNs = t->kNearestNeighbors(_p, k, threadNum);
        candidates.insert( end(candidates), begin(kNNs), end(kNNs) );
    }
    // Include the buffer
    for(size_t i = 0; i < buffer.size(); ++i)
        candidates.push_back( buffer[i] );
    // Sort candidates by their distance to the query point _p
    candidates.sort(
        [&](PointPlane *p1, PointPlane *p2) -> bool {
            return sqr(p1->point[0] -_p[0])+sqr(p1->point[1] -_p[1])+sqr(p1->point[2] -_p[2])
                 < sqr(p2->point[0] -_p[0])+sqr(p2->point[1] -_p[1])+sqr(p2->point[2] -_p[2]);
        }
    );
    // Pick first k points
    std::vector<PointPlane*> result;
    std::list<PointPlane*>::iterator beg = candidates.begin();
    k = std::min(k, (int)candidates.size());
    for(int i = 0; i < k; ++i) {
        result.push_back( *beg );
        ++beg;
    }
    return result;
}

std::vector<PointPlane*> BkdTreePlanes::kNearestRangeSearch(double *_p,
								int k,
								double sqRad2,
								int threadNum) const
{
    //TODO: implement this.
    cout << "WARN: BkdTreePlanes::kNearestRangeSearch not implemented yet!" << endl;
    return std::vector<PointPlane*>(0);
}

std::vector<PointPlane*> BkdTreePlanes::fixedRangeSearch(double *_p,
								double sqRad2,
								int threadNum) const
{
    std::vector<PointPlane*> allPts;
    for (size_t i = 0; i < forest.size(); ++i) {
        if (forest.at(i).empty) continue;
        KDtreePlanes *t = forest.at(i).tree;
        std::vector<PointPlane*> pts = t->fixedRangeSearch(_p, sqRad2, threadNum);
        allPts.insert( end(allPts), begin(pts), end(pts) );
    }
    for (size_t i = 0; i < buffer.size(); ++i) {
        if (Dist2( buffer[i]->point , _p) <= sqRad2)
            allPts.push_back(buffer[i]);
    }
    return allPts;
}

/**
 * @brief Axis-aligned bounding box search.
 * @param _p: Minimum x, y, z values for the bounding box.
 * @param _p0: Maximum x, y, z values for the bouding box.
 * @return List of points that are inside the boudning box.
 */
std::vector<PointPlane*> BkdTreePlanes::AABBSearch(double *_p,
                                double* _p0,
                                int threadNum) const
{
    //TODO: implement this.
    cout << "WARN: BkdTreePlanes::AABBSearch not implemented yet!" << endl;
    return std::vector<PointPlane*>(0);
}

PointPlane * BkdTreePlanes::segmentSearch_1NearestPoint(double *_p,
                                double* _p0,
                                double maxdist2,
                                int threadNum) const
{
    //TODO: implement this.
    cout << "WARN: BkdTreePlanes::segmentSearch_1NearestPoint not implemented yet!" << endl;
    return 0;
}

/**
 *  @brief Merges the trees from index 0 up to index <index>.
 *  Take all the points that are in those trees, delete them, and setup
 *  a new tree with all the deleted points.
 */
void BkdTreePlanes::mergeTreesLogarithmic(int index, int nrpts, int threadNum)
{
    // setup new points memory
    PointPlane **newpts = new PointPlane*[nrpts];

    // copy the points
    int ptindex = 0;
    for (int k = 0; k < index; ++k) {

        if (forest.at(k).empty) continue;
        // else collect and copy pts.
        std::vector<PointPlane*> datak = forest.at(k).tree->CollectPts(threadNum);
        for (size_t q = 0; q < forest.at(k).nrpts; ++q) {
            newpts[ptindex] = datak[q];
            ptindex++;
        }

        // Free previous data fields
        forest.at(k).nrpts = 0;
        forest.at(k).empty = true;
    }

    // add the points in the buffer
    for (int i = 0; i < buffer.size(); ++i) {
        newpts[ptindex] = buffer.at(i);
        ptindex++;
    }

    // Allocate new points
    forest.at(index).tree = new KDtreePlanes(newpts, nrpts, bucketSize);
    forest.at(index).nrpts = nrpts;
    forest.at(index).empty = nrpts == 0;
    return;
}
