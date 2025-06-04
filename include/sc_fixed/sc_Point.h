/**
 *  @file
 *  @brief Representation of a 3D point
 *  @author Tom Fleischmann, Jonas Wiesner, Yannik Winzer - University of Wuerzburg, Germany
 */

#ifndef __SC_POINT_H__
#define __SC_POINT_H__

#include <cmath>
#include "sc_fixed/sc_fixed_math.h"
#include <iostream>

#include <stdexcept>

/**
 * @brief Representation of a point in 3D space
 */
class sc_Point {

public:
  /**
   *	Default constructor
   */
  inline sc_Point();

  /**
   *	Copy constructor
   */
  inline sc_Point(const sc_Point& p);

  /**
   *	Constructor with an array, i.e., vecctor of coordinates
   */
  inline sc_Point(const double *p);

  /**
   *	Constructor with an array, i.e., vecctor of coordinates and color
   */
  inline sc_Point(const double *p, const char *c);

  /**
   *	Constructor with three double values
   */
  inline sc_Point(const double _x, const double _y, const double _z);

  /**
   *	Constructor with three double values and three color values
   */
  inline sc_Point(const double _x, const double _y, const double _z,
			const char _r, const char _g, const char _b);

  /**
   *	Constructor with three double values and three normal values
   */
  inline sc_Point(const double _x,  const double _y,  const double _z,
			const double _nx, const double _ny, const double _nz);

  inline sc_Point operator+(const sc_Point &p) const;
  inline sc_Point operator-(const sc_Point &p) const;
  inline sc_Point& operator-=(const sc_Point &p);
  inline sc_Point& operator+=(const sc_Point &p);
  inline bool operator!=(const sc_Point &p) const;
  inline bool operator==(const sc_Point &p) const;

  inline void transform(const double alignxf[16]);
  inline double distance(const sc_Point& p);
  inline friend std::ostream& operator<<(std::ostream& os, const sc_Point& p);
  inline friend std::istream& operator>>(std::istream& is, sc_Point& p);

  // also public; set/get functions not necessary here
  /// x coordinate in 3D space
  double x;
  /// y coordinate in 3D space
  double y;
  /// z coordinate in 3D space
  double z;
  /// normal x direction in 3D space
  double nx;
  /// normal x direction in 3D space
  double ny;
  /// normal x direction in 3D space
  double nz;
  /// additional information about the point, e.g., semantic
  ///  also used in veloscan for distiuguish moving or static
  int type;

  /////////////////////////for veloslam/////////////////////////////
  double rad;
  ///    tang in  cylindrical coordinates for veloscan
  double tan_theta;
  // point id in points for veloscan , you can use it find point.
  long point_id;
  /////////////////////////for veloslam/////////////////////////////

  // color information of the point between 0 and 255
  // rgb
  unsigned char rgb[3];

  float reflectance;
  float temperature;
  float amplitude;
  float deviation;

  static inline sc_Point cross(const sc_Point &X, const sc_Point &Y) {
    sc_Point res;
    res.x = X.y * Y.z - X.z * Y.y;
    res.y = X.z * Y.x - X.x * Y.z;
    res.z = X.x * Y.y - X.y * Y.x;
    return res;
  };

  static inline sc_Point norm(const sc_Point &p) {
    double l = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    sc_Point res(p.x/l, p.y/l, p.z/l);
    return res;
  };

};

inline sc_Point operator*(const double &v, const sc_Point &p) {
  sc_Point res;
  res.x = v * p.x;
  res.y = v * p.y;
  res.z = v * p.z;
  return res;
}



#include "sc_fixed/sc_Point.icc"

#endif
