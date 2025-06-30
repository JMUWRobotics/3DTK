/** @file: ioplanes.h
 *
 * @brief: This class should never actually be instanced.
 * It has only static members, such as read() and write() which provide
 * IO functionality.
 *
 * When using read(), every plane will be stored in PlaneIO::allPlanes,
 * which is a vector of NormalPlane pointers (see normalplane.h)
 * The class can read multiple planes from entire directories, as specified
 * in the "normals.list" and "planes.list" files.
 *
 * @author Fabian Arzberger, JMU, Germany.
 *
 * Released under the GPL version 3.
 */

#ifndef __IO_CONV_PLANE_H_
#define __IO_CONV_PLANE_H_

#include "model/normalplane.h"
#include "slam6d/normals.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <sstream>

class PlaneIO {

public:
    static void read(std::string);
    static Planes allPlanes;

    static int write(Planes&, std::string&);

private:
    static int existsDir(const char* path)
    {
        struct stat info;
        if (stat( path, &info ) != 0) return 0;
        else if ( info.st_mode & S_IFDIR ) return 1;
        else return 0;
    }
    // wrapper around Dot product function
    static double dot(Point &p, double *p2)
    {
        double *p1 = new double[3];
        p1[0] = p.x;
        p1[1] = p.y;
        p1[2] = p.z;
        return Dot(p1, p2);
    }
};


#endif //__NORM_CONV_PLANE_H_