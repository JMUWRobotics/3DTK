/**
 * @file
 * @brief Representation of an octree for show
 * @author Jan Elseberg. Jacobs University Bremen gGmbH, Germany
 * @author Andreas Nuechter. Jacobs University Bremen gGmbH, Germany
 */

#ifndef SHOWBOCTREE_H
#define SHOWBOCTREE_H

#include "slam6d/Boctree.h"
#include "show/colormanager.h"
#include "show/scancolormanager.h"
#include "show/viewcull.h"
#include "show/colordisplay.h"
#include "slam6d/scan.h"

using namespace show;

/**
 * @brief Octree for show
 *
 * A cubic bounding box is calculated
 * from the given 3D points. Then it
 * is recusivly subdivided into smaller
 * subboxes
 *
 * It contains software culling functionalities
 */
template <class T>
class Show_BOctTree : public colordisplay
{
protected:
  BOctTree<T>* m_tree;

  unsigned long maxtargetpoints;

  unsigned int current_lod_mode;

  //! A copy of pointdim of the tree initialized from
  unsigned int POINTDIM;

  DataOcttree* m_cache_access;
  Scan* m_scan;
  ScanColorManager* scm;

  void init(ScanColorManager* _scm) {
    scm = _scm;
    setColorManager(0);
    if (scm) {
      scm->registerTree(this);
      scm->updateRanges(m_tree->getMins());
      scm->updateRanges(m_tree->getMaxs());
    }
    POINTDIM = m_tree->getPointdim();
    maxtargetpoints = maxTargetPoints(m_tree->getRoot());
    current_lod_mode = 0;
    m_cache_access = 0;
    m_scan = 0;
  }

public:
  //return the scm
  ScanColorManager* getScanColorManager()
  {
    return scm;
  }
  //! Create with tree held in cache, lock indefinitely
  Show_BOctTree(Scan* scan, DataOcttree* cache_access, ScanColorManager* scm = 0)
  {
    m_tree = &cache_access->get();
    init(scm);
    // save cache access and hold it until an unlock
    m_cache_access = cache_access;
    m_scan = scan;
  }

  //For Cylinder fitting
  std::vector<float*> RangeSearchAlongDir(double *point, double* dir, double maxDistR, double minDistR, double maxDistDir, int threadNum){
    return m_tree->RangeSearchAlongDir(point, dir, maxDistR, minDistR, maxDistDir, threadNum);
  }

  //! Retrieve a cache access object and update the tree pointer
  void lockCachedTree() {
    if(m_scan != 0 && m_cache_access == 0) {
      m_cache_access = new DataOcttree(m_scan->get("octtree"));
      m_tree = &(m_cache_access->get());
    }
  }

  //! Remove the data access
  void unlockCachedTree() {
    if(m_scan != 0 && m_cache_access != 0) {
      delete m_cache_access; m_cache_access = 0;
      //m_tree = 0;
    }
  }

  //! Create with already constructed tree and take ownership
  Show_BOctTree(BOctTree<T>* tree, ScanColorManager* scm = 0) :
    m_tree(tree)
  {
    init(scm);
  }

  //! Create tree by points
  template <class P>
  Show_BOctTree(P * const* pts, int n, T voxelSize, PointType pointtype = PointType(), ScanColorManager *scm = 0)
  {
    m_tree = new BOctTree<T>(pts, n, voxelSize, pointtype, true);
    init(scm);
  }

  //! Create tree by deserializing from file
  Show_BOctTree(std::string filename, ScanColorManager *scm = 0)
  {
    m_tree = new BOctTree<T>(filename);
    init(scm);
  }

  virtual ~Show_BOctTree() {
    // only delete cache access if created via this method
    if(m_cache_access) {
      delete m_cache_access;
    }
    // if not constructed by a cache-located tree, delete the owned tree
    if(!m_scan)
      delete m_tree;
  }

  int pow(int base, int exponent) {
    if (exponent == 0)
        return 1;
    else if (exponent % 2 == 0) {
        int half_pow = pow(base, exponent / 2);
        return half_pow * half_pow;
	} else
        return base * pow(base, exponent - 1);
  }
  BOctTree<T>* getTree() const { return m_tree; }

  void serialize(const std::string& filename) const { m_tree->serialize(filename); }

  unsigned int getMemorySize() const { return m_tree->getMemorySize(); }

  // virtual functions from colordisplay

  void selectRayBrushSize(std::set<T *> &points, int brushsize) {
    selectRayBS(points, m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize(), brushsize);
  }

  void selectRay(std::set<T *> &points, int depth = INT_MAX) {
    selectRay(points, m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize(), depth);
  }

  void selectRay(T * &point) {
    selectRay(point, m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize(), FLT_MAX);
  }

  void cycleLOD() {
    current_lod_mode = (current_lod_mode+1)%3;
  }

  void drawLOD(float ratio) {
    switch (current_lod_mode) {
      case 0:
        glBegin(GL_POINTS);
        displayOctTreeCulledLOD(maxtargetpoints * ratio, m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize());
        glEnd();
        break;
      case 1:
        glBegin(GL_POINTS);
        displayOctTreeCulledLOD2(ratio, m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize());
        glEnd();
        break;
      case 2:
#ifdef WITH_GLEE
        if (GLEE_ARB_point_parameters) {
          glPointParameterfARB(GL_POINT_SIZE_MIN_ARB, 1.0);
          glPointParameterfARB(GL_POINT_SIZE_MAX_ARB, 100000.0);
          GLfloat p[3] = {0.0, 0.0000, 0.0000005};
          glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB, p);
          displayOctTreeCPAllCulled(m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize(),
							 m_tree->getSize() / pow(2, std::min( (int)(ratio * m_tree->getMaxDepth()), (int)(m_tree->getMaxDepth() - 3))));
          p[0] = 1.0;
          p[2] = 0.0;
          glPointParameterfvARB(GL_POINT_DISTANCE_ATTENUATION_ARB, p);
        }
#endif
        break;
      default:
        break;
    }
  }

  void draw() {
    glBegin(GL_POINTS);
    displayOctTreeAllCulled(m_tree->getRoot(), m_tree->getCenter(), m_tree->getSize());
    glEnd();
  }

  // these are needed on windows so that show can be used as a library
  // we need to call the ExtractFrustum in a same process as the draw function on windows
  void extractFrustumAndDrawLOD(float ratio, short detail)
  {
	  ExtractFrustum(detail);
	  drawLOD(ratio);
  }
  void extractFrustumAndDraw(short detail)
  {
	  ExtractFrustum(detail);
	  draw();
  }

  // reroute center call (for recast from colordisplay to show_bocttree)
  void getCenter(double _center[3]) const {
    m_tree->getCenter(_center);
  }

protected:

  //! ?
  unsigned long maxTargetPoints(const bitoct &node) {
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    unsigned long max = 0;

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned long length = points[0].length;
          if (length > max) max = length;
        } else { // recurse
          unsigned long tp = maxTargetPoints(children->node);
          if (tp > max) max = tp;
        }
        ++children; // next child
      }
    }

    return max*POPCOUNT(node.valid);
  }

  void displayOctTreeAll(const bitoct &node) {
//    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point
          for(unsigned int iterator = 0; iterator < length; iterator++ ) {
            if(cm) cm->setColor(point);
            glVertex3f( point[0], point[1], point[2]);
            point+=POINTDIM;
          }
        } else { // recurse
          displayOctTreeAll( children->node);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeAllCulled(const bitoct &node, const T* center, T size ) {
    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeAll(node);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
       //   if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
              point+=POINTDIM;
            }
          //}
        } else { // recurse
          displayOctTreeAllCulled( children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeLOD2(float ratio, const bitoct &node, const T* center, T size ) {

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point

          int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point
          l = std::max((int)(l*l*ratio), 0);
          if (l > 1) {
            if ((int)length > l ) {
              T each = (T)POINTDIM * (T)((T)length/(T)l);
              T *p;
              int index;
              for(int iterator = 0; iterator < l; iterator++ ) {
                index = (T)iterator * each;
                p = point + index - index%POINTDIM;
                if(cm) cm->setColor(p);
                glVertex3f( p[0], p[1], p[2]);
              }
            } else if ((int)length <= l) {
              for(unsigned int iterator = 0; iterator < length; iterator++ ) {
                if(cm) cm->setColor(point);
                glVertex3f( point[0], point[1], point[2]);
                point+=POINTDIM;
              }
            } /* else if (l == 1) {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
            }*/
          } else {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
          }
        } else { // recurse
            int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point
            l = std::max((int)(l*l*ratio), 0);
            if (l > 0) {
              displayOctTreeCulledLOD2(ratio, children->node, ccenter, size/2.0);
            }
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCulledLOD2(float ratio, const bitoct &node, const T* center, T size ) {

    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeLOD2(ratio, node, center, size);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point

            int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point
            l = std::max((int)(l*l*ratio), 0);
            if (l != 0) {
              if ((int)length > l ) {
                T each = (T)POINTDIM * (T)((T)length/(T)l);
                T *p;
                int index;
                for(int iterator = 0; iterator < l; iterator++ ) {
                  index = (T)iterator * each;
                  p = point + index - index%POINTDIM;
                  if(cm) cm->setColor(p);
                  glVertex3f( p[0], p[1], p[2]);
                }
              } else if ((int)length <= l) {
                for(unsigned int iterator = 0; iterator < length; iterator++ ) {
                  if(cm) cm->setColor(point);
                  glVertex3f( point[0], point[1], point[2]);
                  point+=POINTDIM;
                }
              } else if (l == 1) {
                if(cm) cm->setColor(point);
                glVertex3f( point[0], point[1], point[2]);
              }
            }
          }
        } else { // recurse
          //int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point
          //l = max((int)(l*l*ratio), 0);
          //if (l > 0) {
            displayOctTreeCulledLOD2(ratio, children->node, ccenter, size/2.0);
          //}
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeLOD3(long targetpts, const bitoct &node, const T* center, T size ) {
    if (targetpts <= 0) return; // no need to display anything

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point
          int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point

          if ( l <= targetpts) {  // only a single pixel on screen only paint one point
            if(cm) cm->setColor(point);
            glVertex3f( point[0], point[1], point[2]);
          } else {
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
              point+=POINTDIM;
            }
          }
        } else { // recurse
          displayOctTreeLOD3(targetpts, children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCulledLOD3(long targetpts, const bitoct &node, const T* center, T size ) {
    if (targetpts <= 0) return; // no need to display anything

    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeLOD3(targetpts, node, center, size);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);


    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        //l = std::min( max((int)(l*l*ratio), 1), targetpts);

        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point

            int l = LOD2(ccenter[0], ccenter[1], ccenter[2], size/2.0);  // only a single pixel on screen only paint one point
          if (targetpts <= 0 ) std::cout << l << " " << targetpts << std::endl;
            if ( l <= targetpts) {  // only a single pixel on screen only paint one point
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
            } else {
              for(unsigned int iterator = 0; iterator < length; iterator++ ) {
                if(cm) cm->setColor(point);
                glVertex3f( point[0], point[1], point[2]);
                point+=POINTDIM;
              }
            }
          }

        } else { // recurse
          displayOctTreeCulledLOD3(targetpts, children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCulledLOD(long targetpts, const bitoct &node, const T* center, T size ) {
    if (targetpts <= 0) return; // no need to display anything

    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeLOD(targetpts, node, center, size);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    unsigned short nc = POPCOUNT(node.valid);
    long newtargetpts = targetpts;
    if (nc > 0) {
      newtargetpts = newtargetpts/nc;
      if (newtargetpts <= 0 ) return;
    }

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point

            if (length > 10 && !LOD(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {  // only a single pixel on screen only paint one point
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
            } else if (length <= newtargetpts) {        // more points requested than possible, plot all
              for(unsigned int iterator = 0; iterator < length; iterator++ ) {
                if(cm) cm->setColor(point);
                glVertex3f( point[0], point[1], point[2]);
                point+=POINTDIM;
              }
            } else {                         // select points to show
              // TODO smarter subselection of points here
              T each = (T)POINTDIM * (T)((T)length/(T)newtargetpts);
              T *p;
              int index;
              for(unsigned int iterator = 0; iterator < newtargetpts; iterator++ ) {
                index = (T)iterator * each;
                p = point + index - index%POINTDIM;
                if(cm) cm->setColor(p);
                glVertex3f( p[0], p[1], p[2]);
                //point += each;
              }
            }
          }

        } else { // recurse
          displayOctTreeCulledLOD(newtargetpts, children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeLOD(long targetpts, const bitoct &node, const T* center, T size ) {
    if (targetpts <= 0) return; // no need to display anything

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    unsigned short nc = POPCOUNT(node.valid);
    long newtargetpts = targetpts;
    if (nc > 0) {
      newtargetpts = newtargetpts/nc;
      if (newtargetpts <= 0 ) return;
    }


    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point
          if (length > 10 && !LOD(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {  // only a single pixel on screen only paint one point
            if(cm) cm->setColor(point);
            glVertex3f( point[0], point[1], point[2]);
          } else if (length <= newtargetpts) {        // more points requested than possible, plot all
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
              point+=POINTDIM;
            }
          } else {                         // select points to show
            // TODO smarter subselection of points here
            T each = (T)POINTDIM * (T)((T)length/(T)newtargetpts);
            T *p;
            int index;
            for(unsigned int iterator = 0; iterator < newtargetpts; iterator++ ) {
              index = (T)iterator * each;
              p = point + index - index%POINTDIM;
              if(cm) cm->setColor(p);
              glVertex3f( p[0], p[1], p[2]);
              //point += each;
            }
          }
        } else { // recurse
          displayOctTreeLOD(newtargetpts, children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
  }

  // this is because otherwise the following line is overriding selectRay from
  // the base class
  using colordisplay::selectRay;
  void selectRay(std::set<T *> &selpoints, const bitoct &node, const T* center, T size, int max_depth, int depth = 0) {
    if (depth < max_depth &&  !HitBoundingBox(center, size ))return;

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( !(depth+1  < max_depth) || HitBoundingBox(ccenter, size) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
                selpoints.insert(point);
              point+=POINTDIM;
            }
          }
        } else { // recurse
          selectRay( selpoints, children->node, ccenter, size/2.0, max_depth, depth+1);
        }
        ++children; // next child
      }
    }
  }

  void selectRayBS(std::set<T *> &selpoints, const bitoct &node, const T* center, T size, int brushsize) {
    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( HitBoundingBox(ccenter, size) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              //check if point is behind the nearplane and infront of the farplane
              float z_dis = RadialDistCamera(point);
              if (z_dis >= 0 && z_dis <= 1 && ScreenDist(point) < brushsize){
                selpoints.insert(point);
              }
            point+=POINTDIM;
            }
          }
        } else { // recurse
          selectRayBS( selpoints, children->node, ccenter, size/2.0, brushsize);
        }
        ++children; // next child
      }
    }
  }

  void selectRay(T * &selpoint, const bitoct &node, const T* center, T size, float min) {
    if (!HitBoundingBox(center, size ))return;

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( HitBoundingBox(ccenter, size) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              //check if point is behind the nearplane and infront of the farplane
              float z_dis = RadialDistCamera(point);
              T dist = RayDist(point);
              if (z_dis >= 0 && z_dis <= 1 && min > dist && ScreenDist(point) < 5) {
                selpoint = point;
                min = dist;
              }
              point+=POINTDIM;
            }
          }
        } else { // recurse
          selectRay( selpoint, children->node, ccenter, size/2.0, min);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCAllCulled(const bitoct &node, const T* center, T size, T minsize ) {
    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeCAll(node, center, size, minsize);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf || minsize > size ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            showCube(ccenter, size/2.0);
          }
        } else { // recurse
          displayOctTreeCAllCulled( children->node, ccenter, size/2.0, minsize);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCAll(const bitoct &node, const T* center, T size, T minsize ) {
    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf || minsize > size ) {   // if ith node is leaf get center
          showCube(ccenter, size/2.0);
        } else { // recurse
          displayOctTreeCAll( children->node, ccenter, size/2.0, minsize);
        }
        ++children; // next child
      }
    }
  }


  void displayOctTreeCPAllCulled(const bitoct &node, const T* center, T size, T minsize ) {
    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return;  // culled do not continue with this branch of the tree

    if (res == 2) { // if entirely within frustrum discontinue culling
      displayOctTreeCPAll(node, center, size, minsize);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if ( minsize > size ) {
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            if(cm) {
              if (( 1 << i ) & node.leaf )
                cm->setColor( &(children->getPoints()[1]) );
              else
                cm->setColor( m_tree->pickPoint(children->node) );
            }

            glPointSize(size/2.0);
            glBegin(GL_POINTS);
            glVertex3f( ccenter[0], ccenter[1], ccenter[2] );
            glEnd();
          }
        }else if ( ( 1 << i ) & node.leaf ) {
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            T *point = &(points[1].v);  // first point
            glPointSize(1.0);
            glBegin(GL_POINTS);
            for(unsigned int iterator = 0; iterator < length; iterator++ ) {
              if(cm) cm->setColor(point);
              glVertex3f( point[0], point[1], point[2]);
              point+=POINTDIM;
            }
            glEnd();
          }
        } else { // recurse
          displayOctTreeCPAllCulled( children->node, ccenter, size/2.0, minsize);
        }
        ++children; // next child
      }
    }
  }

  void displayOctTreeCPAll(const bitoct &node, const T* center, T size, T minsize ) {
    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if ( minsize > size ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            if(cm) {
              if (( 1 << i ) & node.leaf )
                cm->setColor( &(children->getPoints()[1]) );
              else
                cm->setColor( m_tree->pickPoint(children->node) );
            }

            glPointSize(size/2.0);
            glBegin(GL_POINTS);
            glVertex3f( ccenter[0], ccenter[1], ccenter[2] );
            glEnd();
          }
        }else if ( ( 1 << i ) & node.leaf ) {
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point
          glPointSize(1.0);
          glBegin(GL_POINTS);
          for(unsigned int iterator = 0; iterator < length; iterator++ ) {
            if(cm) cm->setColor(point);
            glVertex3f( point[0], point[1], point[2]);
            point+=POINTDIM;
          }
          glEnd();
        } else { // recurse
          displayOctTreeCPAll( children->node, ccenter, size/2.0, minsize);
        }
        ++children; // next child
      }
    }
  }

  void showCube(const T* center, T size) {
    glLineWidth(1.0);
    glBegin(GL_QUADS);      // draw a cube with 6 quads
    glColor3f(0.0f,1.0f,0.0f);      // Set The Color To Green
    glVertex3f(center[0] + size, center[1] + size, center[2] - size);
    glVertex3f(center[0] - size, center[1] + size, center[2] - size);
    glVertex3f(center[0] - size, center[1] + size, center[2] + size);
    glVertex3f(center[0] + size, center[1] + size, center[2] + size);
    glColor3f(1.0f,0.5f,0.0f);      // Set The Color To Orange

    glVertex3f(center[0] + size, center[1] - size, center[2] + size);
    glVertex3f(center[0] - size, center[1] - size, center[2] + size);
    glVertex3f(center[0] - size, center[1] - size, center[2] - size);
    glVertex3f(center[0] + size, center[1] - size, center[2] - size);

    glColor3f(1.0f,0.0f,0.0f);      // Set The Color To Red
    glVertex3f(center[0] + size, center[1] + size, center[2] + size);
    glVertex3f(center[0] - size, center[1] + size, center[2] + size);
    glVertex3f(center[0] - size, center[1] - size, center[2] + size);
    glVertex3f(center[0] + size, center[1] - size, center[2] + size);

    glColor3f(1.0f,1.0f,0.0f);      // Set The Color To Yellow
    glVertex3f(center[0] + size, center[1] - size, center[2] - size);
    glVertex3f(center[0] - size, center[1] - size, center[2] - size);
    glVertex3f(center[0] - size, center[1] + size, center[2] - size);
    glVertex3f(center[0] + size, center[1] + size, center[2] - size);

    glColor3f(0.0f,0.0f,1.0f);      // Set The Color To Blue
    glVertex3f(center[0] - size, center[1] + size, center[2] + size);
    glVertex3f(center[0] - size, center[1] + size, center[2] - size);
    glVertex3f(center[0] - size, center[1] - size, center[2] - size);
    glVertex3f(center[0] - size, center[1] - size, center[2] + size);

    glColor3f(1.0f,0.0f,1.0f);      // Set The Color To Violet
    glVertex3f(center[0] + size, center[1] + size, center[2] - size);
    glVertex3f(center[0] + size, center[1] + size, center[2] + size);
    glVertex3f(center[0] + size, center[1] - size, center[2] + size);
    glVertex3f(center[0] + size, center[1] - size, center[2] - size);

    glEnd();
  }


  void displayOctTreeAllCulled(const bitoct &node, const T* center, T size, float *frustum[6], unsigned char frustsize ) {
    float *new_frustum[6]; unsigned char counter = 0;
    for (unsigned char p = 0; p < frustsize; p++ ) {
      char res = PlaneAABB(center[0], center[1], center[2], size, frustum[p]);
      if (res == 0) { // cube is on the wrong side of the plane (not visible)
        return;
      } else if ( res == 1 ) {  // plane intersects this volume continue culling with this plane
        new_frustum[counter++] = frustum[p];
      } // other case is simply not to continue culling with the respective plane
    }
    if (counter == 0) { // if entirely within frustrum discontinue culling
      displayOctTreeAll(node);
      return;
    }

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          pointrep *points = children->getPointreps();
          unsigned int length = points[0].length;
          T *point = &(points[1].v);  // first point
          for(unsigned int iterator = 0; iterator < length; iterator++ ) {
            if(cm) cm->setColor(point);
            glVertex3f( point[0], point[1], point[2]);
            point+=POINTDIM;
          }
        } else { // recurse
          displayOctTreeAllCulled( children->node, ccenter, size/2.0, new_frustum, counter);
        }
        ++children; // next child
      }
    }
  }

  unsigned long int countVisiblePoints(const bitoct &node, const T* center, T size ) {
    unsigned long int result = 0;
    int res = CubeInFrustum2(center[0], center[1], center[2], size);
    if (res==0) return 0;  // culled do not continue with this branch of the tree

    T ccenter[3];
    bitunion<T> *children;
    bitoct::getChildren(node, children);

    for (short i = 0; i < 8; i++) {
      if (  ( 1 << i ) & node.valid ) {   // if ith node exists
        BOctTree<T>::childcenter(center, ccenter, size, i);  // childrens center
        if (  ( 1 << i ) & node.leaf ) {   // if ith node is leaf get center
          // check if leaf is visible
          if ( CubeInFrustum(ccenter[0], ccenter[1], ccenter[2], size/2.0) ) {
            pointrep *points = children->getPointreps();
            unsigned int length = points[0].length;
            result += length;
          }
        } else { // recurse
          result += countVisiblePoints( children->node, ccenter, size/2.0);
        }
        ++children; // next child
      }
    }
    return result;
  }
};

#endif
