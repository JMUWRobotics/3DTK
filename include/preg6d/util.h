#ifndef __UTIL_H_
#define __UTIL_H_

#include <sys/types.h>
#include <sys/stat.h>
#include "stdio.h"
#include <chrono>
#include <thread>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <math.h>
#include "slam6d/scan.h"
#include "slam6d/globals.icc"

/*
 * This code is from a stackoverflow thread.
 * https://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
 * Im sure many of you have also used this when working with rgb and hsv color spaces.
 */
typedef struct {
    double r;       // a fraction between 0 and 1
    double g;       // a fraction between 0 and 1
    double b;       // a fraction between 0 and 1
} rgb;

typedef struct {
    double h;       // angle in degrees
    double s;       // a fraction between 0 and 1
    double v;       // a fraction between 0 and 1
} hsv;

// hsv   rgb2hsv(rgb in);
// rgb   hsv2rgb(hsv in);

inline hsv rgb2hsv(rgb in)
{
    hsv         out;
    double      min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min  < in.b ? min  : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max  > in.b ? max  : in.b;

    out.v = max;                                // v
    delta = max - min;
    if (delta < 0.00001)
    {
        out.s = 0;
        out.h = 0; // undefined, maybe nan?
        return out;
    }
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        out.s = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0
        // s = 0, h is undefined
        out.s = 0.0;
        out.h = -1;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.h = ( in.g - in.b ) / delta;        // between yellow & magenta
    else
    if( in.g >= max )
        out.h = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
    else
        out.h = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.h *= 60.0;                              // degrees

    if( out.h < 0.0 )
        out.h += 360.0;

    return out;
}


inline rgb hsv2rgb(hsv in)
{
    double      hh, p, q, t, ff;
    long        i;
    rgb         out;

    if(in.s <= 0.0) {       // < is bogus, just shuts up warnings
        out.r = in.v;
        out.g = in.v;
        out.b = in.v;
        return out;
    }
    hh = in.h;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.v * (1.0 - in.s);
    q = in.v * (1.0 - (in.s * ff));
    t = in.v * (1.0 - (in.s * (1.0 - ff)));

    switch(i) {
    case 0:
        out.r = in.v;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.v;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.v;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = in.v;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.v;
        break;
    case 5:
    default:
        out.r = in.v;
        out.g = p;
        out.b = q;
        break;
    }
    return out;
}
// End copyright

// Check if a directory exists
inline int existsDir(const char* path)
{
    struct stat info;
    if (stat( path, &info ) != 0) return 0;
    else if ( info.st_mode & S_IFDIR ) return 1;
    else return 0;
}

inline bool rgbEquals(rgb c1, rgb c2)
{
    // yanky but works
    return
        fabs(c1.r - c2.r) < 0.01 &&
        fabs(c1.g - c2.g) < 0.01 &&
        fabs(c1.b - c2.b) < 0.01;
}

inline bool eigenValueOK( double *eigen, double thresh)
{
    // i.e. the ratio between extent of normal (z) direction
    // and extent in all (xyz) directions -> Measure of planarity.
    double den = (eigen[0] + eigen[1] + eigen[2]);
    if (den == 0.0) return false;
    return eigen[0]/den < thresh;
}

inline bool intervalOverlap(int x1, int x2, int y1, int y2)
{
    return x1 <= y2 && y1 <= x2;
}

/**
 * Small, simple timer class for very simple benchmarking
 */
template <class TimeT = std::chrono::milliseconds,
          class ClockT = std::chrono::steady_clock>
class Timer
{
    using timep_t = decltype(ClockT::now());
    timep_t _start, _end;
  public:
    void tick() {
        _end = timep_t{};
        _start = ClockT::now();
    }

    void tock() {
        _end = ClockT::now();
    }

    template <class TT = TimeT>
    TT duration() const {
        // Use gsl_Expects if your project supports it.
        assert(_end != timep_t{} && "Timer must toc before reading the time");
        return std::chrono::duration_cast<TT>(_end - _start);
    }
};

inline void debugCCCount(int **cccount, int nplanes, int nclusters)
{
    std::cout << "\t";
    for (int n = 0; n < nplanes; ++n) std::cout << "P" << n << "\t";
    std::cout << std::endl;
    for (int i = 0; i < nclusters; ++i) {
        std::cout << "C" << i << "\t";
        for (int j = 0; j < nplanes; ++j) {
            std::cout << cccount[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

inline void sortClustersBySize(std::vector<std::vector<size_t> > &clusters)
{
    // Reorganize cluster (sort according to size) in a list
    std::list< std::vector<size_t> > clusterlist(clusters.begin(), clusters.end());
    clusterlist.sort(
        [=] (const std::vector<size_t> &c1, const std::vector<size_t> &c2) -> bool {
            return c1.size() > c2.size();
        }
    );
    // Clear original cluster vector
    clusters.clear();
    // Convert ordered list back to vector
    clusters = std::vector<std::vector<size_t> >(clusterlist.begin(), clusterlist.end());
}

inline std::vector<size_t> clusterHist(std::vector<std::vector<size_t> > const &clusters)
{
    // Reorganize cluster (sort according to size)
    std::list< std::vector<size_t> > clusterlist(clusters.begin(), clusters.end());
    clusterlist.sort(
        [=] (const std::vector<size_t> &c1, const std::vector<size_t> &c2) -> bool {
            return c1.size() > c2.size();
        }
    );
    std::vector<size_t> histogram;
    for (std::vector<size_t> const &c : clusterlist) {
        histogram.push_back( c.size() );
    }
    return histogram;
}

inline void writeHist(std::vector<size_t> const &histogram, std::string path)
{
    std::ofstream file;
    file.open( path.c_str() );
    for (uint i = 0; i < histogram.size(); ++i)
        file << i+1 << " " << histogram[i] << "\n";
    file.flush();
    file.close();
}

enum InterpolType {
    LINEAR,
    QUADRATIC,
    CUBIC
};
/**
 * Handles growrate adaption.
 */
struct InterpolMap {

    double _imin; // Input, left border
    double _imax; // Input, right border
    double _omin; // Output, left border
    double _omax; // Output, right border
    InterpolType _t; // type of adaption

    InterpolMap(){} // dummy constructor. you should never use it

    // Use this constructor instead
    InterpolMap(double imin, double imax, double omin, double omax, InterpolType t)
        : _imin(imin), _imax(imax), _omin(omin), _omax(omax), _t(t) {}

    // Takes any input between imin and imax, spits out mapped output
    inline double operator()(double i) const {
        double tmp = (i - _imin) / (_imax - _imin);
        switch (_t)
        {
        case LINEAR:
            return _omin + (_omax - _omin)*tmp;
        case QUADRATIC:
            return _omin + (_omax - _omin)*sqr( tmp );
        case CUBIC:
            return _omin + (_omax - _omin)*tmp*tmp*tmp;
        default:
            cout << "ERR: this AdaptType enum has not been implemented." << endl;
            return -1;
        }
    }
};

/**
 * @brief Returns minimum distance from all points <pts> to local origin.
 * @param pts: Pointer to the points as double[nrpts][3] array.
 * @param nrpts: Number of points in the array.
 * @return Minimum distance.
 */
inline double minD2(double **pts, int nrpts)
{
    double mind2 = __DBL_MAX__;
    double origin[3] = {0, 0 ,0};
    for (int i = 0; i < nrpts; i++)
        mind2 = min( mind2, Dist2( pts[i], origin ) );
    return mind2;
}

/**
 * @brief Returns maximum distance from all points <pts> to local origin.
 * @param pts: Pointer to the points as double[nrpts][3] array.
 * @param nrpts: Number of points in the array.
 * @return Maximum distance.
 */
inline double maxD2(double **pts, int nrpts)
{
    double maxd2 = 0.0;
    double origin[3] = {0, 0 ,0};
    for (int i = 0; i < nrpts; i++)
        maxd2 = max( maxd2, Dist2( pts[i], origin ) );
    return maxd2;
}

/**
 * @brief Returns maximum distance from all points <pts> to origin <rPos>.
 * @param pts: Pointer to the points as double[nrpts][3] array.
 * @param nrpts: Number of points in the array.
 * @param rPos: Pointer to the origin as double[3] array.
 * @param min: Ref. to variable where min result gets stored.
 * @param max: Ref. to variable where max result gets stored.
 */
inline void minMaxD2(double **pts, int nrpts, double &mind2, double &maxd2)
{
    maxd2 = 0.0;
    mind2 = __DBL_MAX__;
    double d2;
    double origin[3] = {0, 0 ,0};
    for (int i = 0; i < nrpts; i++) {
        d2 = Dist2( pts[i], origin );
        maxd2 = max( maxd2, d2 );
        mind2 = min( mind2, d2 );
    }
}

/**
 * @brief
 */
inline void minMaxD2(DataXYZ xyz, double &mind2, double &maxd2) {
    size_t nrpts = xyz.size();
    maxd2 = 0.0;
    mind2 = __DBL_MAX__;
    double d2;
    double origin[3] = {0, 0 ,0};
    for (size_t i = 0; i < nrpts; i++) {
        d2 = Dist2( xyz[i], origin );
        maxd2 = max( maxd2, d2 );
        mind2 = min( mind2, d2 );
    }
}

inline double eigendiff(double *eig1, double *eig2) {
    return fabs(eig1[0] - eig2[0])
        + fabs(eig1[1] - eig2[1])
        + fabs(eig1[2] - eig2[2]);
}

inline void INFO_VEC3(std::string name, double* vec) {
    cout << name << ": " << vec[0] << " " << vec[1] << " " << vec[2];
}

inline void mean_centroid(std::vector<Point> pts, double* res)
{
    double sum[3] = {0, 0, 0};
    for (const Point &p : pts) {
        sum[0] += p.x;
        sum[1] += p.y;
        sum[2] += p.z;
    }
    for (int i = 0; i < 3; ++i) {
        sum[i] /= pts.size();
        res[i] = sum[i];
    }
};

#endif //__UTIL_H_