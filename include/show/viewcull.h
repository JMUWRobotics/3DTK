#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef WITH_OPENGL
#ifdef _WIN32
#include <windows.h>
#include <GL/glu.h>
#include <GL/glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/glu.h>
#include <stdbool.h>
#else
#include <GL/glut.h>
#include <GL/glu.h>
#include <stdbool.h>
#endif
#else
#include "show/dummygl.h"
#endif

#ifndef __VIEWCULL_H__
#define __VIEWCULL_H__
namespace show{

/////////////////////////  Variable declarations...
/** The 6 planes of the viewing frustum */
extern float frustum[6][4];

/** the modelview * projection matrix to map a model point to an onscreen coordinate */
extern float matrix[16];
/** some useful variables for faster calculation of the pixel coordinates */
extern short VP[4];
/** a unit vector pointing to the right of the screen */
extern float right[3];
/** how much detail is shown, 0 means everything is plotted */
extern short DETAIL;

extern double SX, SY, SZ, EX, EY, EZ;
extern float origin[3], dir[3];   /*ray */
extern float dist;
extern int rayX,rayY;
extern float  rayVP[4];

#define NUMDIM  3
#define RIGHT 0
#define LEFT  1
#define MIDDLE  2

extern float minB[NUMDIM], maxB[NUMDIM];    /*box */
extern float coord[NUMDIM];       /* hit point */


template <class T> bool HitBoundingBox(const T center[3], T size )
{
  minB[0] = center[0] - size;
  minB[1] = center[1] - size;
  minB[2] = center[2] - size;

  maxB[0] = center[0] + size;
  maxB[1] = center[1] + size;
  maxB[2] = center[2] + size;

  bool inside = true;
  char quadrant[NUMDIM];
  int i;
  int whichPlane;
  float maxT[NUMDIM];
  float candidatePlane[NUMDIM];

  // Find candidate planes; this loop can be avoided if
  //  rays cast all from the eye(assume perspective view)
  for (i=0; i<NUMDIM; i++)
    if(origin[i] < minB[i]) {
      quadrant[i] = LEFT;
      candidatePlane[i] = minB[i];
      inside = false;
    }else if (origin[i] > maxB[i]) {
      quadrant[i] = RIGHT;
      candidatePlane[i] = maxB[i];
      inside = false;
    }else {
      candidatePlane[i] = 0.0;
      quadrant[i] = MIDDLE;
    }

  // Ray origin inside bounding box
  if(inside)  {
    return (true);
  }


  // Calculate T distances to candidate planes
  for (i = 0; i < NUMDIM; i++)
    if (quadrant[i] != MIDDLE && dir[i] !=0.)
      maxT[i] = (candidatePlane[i]-origin[i]) / dir[i];
    else
      maxT[i] = -1.;

  // Get largest of the maxT's for final choice of intersection
  whichPlane = 0;
  for (i = 1; i < NUMDIM; i++)
    if (maxT[whichPlane] < maxT[i])
      whichPlane = i;

  // Check final candidate actually inside box
  if (maxT[whichPlane] < 0.) return (false);
  for (i = 0; i < NUMDIM; i++)
    if (whichPlane != i) {
      coord[i] = origin[i] + maxT[whichPlane] *dir[i];
      if (coord[i] < minB[i] || coord[i] > maxB[i])
        return (false);
    } else {
      coord[i] = candidatePlane[i];
    }
  return (true);        // ray hits box
}

void calcRay(int x, int y, double znear, double zfar);

#include <float.h>
//calculates the z value of a point on the screen (0 = nearplane, 1 = farplane):
template <class T>
float RadialDistCamera(T *point){
  GLdouble modelMatrix[16];
  GLdouble projMatrix[16];
  int viewport[4];
  glGetDoublev(GL_MODELVIEW_MATRIX,modelMatrix);
  glGetDoublev(GL_PROJECTION_MATRIX,projMatrix);
  glGetIntegerv(GL_VIEWPORT,viewport);
  double winX, winY, winZ;
  gluProject(point[0], point[1], point[2], modelMatrix, projMatrix, viewport, &winX, &winY, &winZ);
  winY = viewport[3] - winY;
  return winZ;
}

template <class T>
float RayDist(T *point)
{
  return point[0] * dir[0] + point[1] * dir[1] + point[2] * dir[2] - dist;
}

template <class T>
short ScreenDist(T *point) {
  float pn[3];
  // x coordinate on screen, not normalized
  pn[0] = point[0] * matrix[0] + point[1] * matrix[4] + point[2] * matrix[8]  + matrix[12];
  pn[1] = point[0] * matrix[1] + point[1] * matrix[5] + point[2] * matrix[9]  + matrix[13];
  // normalization
  pn[2] = point[0] * matrix[3] + point[1] * matrix[7] + point[2] * matrix[11] + matrix[15];

  // normalized x coordinate on screen
  pn[0] /= pn[2];
  pn[1] /= pn[2];

  // true x coordinate in viewport coordinate system
  //Xi = pn[0]*VP[0] + VP[1];
  //fTempo[4]*0.5+0.5)*viewport[2]+viewport[0];

  float XX = ( (pn[0])*rayVP[0] + rayVP[1]);
  float YY = ( (pn[1])*rayVP[2] + rayVP[3]);

  short dx, dy;
  if (XX > rayX) dx = XX-rayX;
  else dx = rayX-XX;

  if (YY > rayY) dy = YY-rayY;
  else dy = rayY-YY;

  // for the benefit of visual studio's compiler cast to float
  return sqrt((float)(dx*dx + dy*dy));
}

void ExtractFrustum(short detail);
void ExtractFrustum(float *frust[6]);


bool CubeInFrustum( float x, float y, float z, float size );
int  CubeInFrustum2( float x, float y, float z, float size );
char PlaneAABB( float x, float y, float z, float size, float *plane );

void remViewport();
bool LOD(float x, float y, float z, float size);
int LOD2(float x, float y, float z, float size);

}
#endif
