#include "io/ioplanes.h"

Planes PlaneIO::allPlanes;

void PlaneIO::read(std::string dir)
{
    string dir2 = dir;
    // check if the given path exists
    if ( !existsDir( dir.c_str() ) )
    {
        std::cout << "Using path: " << dir << std::endl;
        std::cout << "Invalid plane path. Doesnt exist." << std::endl;
        exit(0);
    }
    dir += "normals.list";
    dir2 += "planes.list";

    // Read the normals.list file to get all planeXXX.n files
    //std::vector<std::string> planefiles;
    std::vector<std::string> hullfiles;
    std::ifstream input, hinput;
    std::string line, hline;
    //input.open( dir.c_str() );
    hinput.open( dir2.c_str() );
    // while ( std::getline( input, line ) )
    // {
    //     cout << line << endl;
    //     cout.flush();
    //     line = line.substr( 7, line.size() ); // remove the "Normal " entry
    //     planefiles.push_back( line );
    // }
    while (std::getline(hinput, hline))
    {
        cout << hline << endl;
        cout.flush();
        hline = hline.substr( 6, hline.size() ); // Remove the "Plane " entry
        hullfiles.push_back( hline );
    }
    input.close();
    input.clear();
    hinput.close();
    hinput.clear();

    // open planeXXX.n files and normals and centerpoint
    for ( size_t i = 0; i < hullfiles.size(); ++i )
    {
        //const auto& path = planefiles[i];
        const auto& hpath = hullfiles[i];
        std::vector<double*> points = *new vector<double*>(0);
        double n[3] = {0 ,0 ,0};
        double x[3] = {0, 0, 0};
        //std::cout << "Reading " << path << " for normal information." << std::endl;
        std::cout << "Reading " << hpath << " for convex hull information." << std::endl;
        // extract points
        //input.open( path.c_str() );
        hinput.open ( hpath.c_str() );
        // Read normals and ref point.
        // while ( input.good() )
        // {
        //     try
        //     {
        //         input >> n[0] >> n[1] >> n[2];
        //         input >> x[0] >> x[1] >> x[2];
        //     }
        //     catch (...) { break; }
        // }
        // Read convex hull
        vector<Point> ppoints;
        while ( hinput.good() )
        {
            try
            {
                double *p = new double[3];
                hinput >> p[0] >> p[1] >> p[2];
                points.push_back(p);
                ppoints.push_back(Point(p));
            }
            catch (...) { break; }
        }
        input.close();
        input.clear();
        hinput.close();
        hinput.clear();

        cout << "Plane with " << points.size() << " pts" << endl;

        double eig[3] = {0, 0, 0};
        try {
            calculateNormal(ppoints, n, eig);
        } catch (NEWMAT::ConvergenceException &e) {
            cout << "Caught. n= " << n[0] << " " << n[1] << " " << n[2] << endl;
        }
        mean_centroid(ppoints, x);
        NormalPlane* new_plane = new NormalPlane( n, x, points );

        // Construct ConvexPlane objects
        allPlanes.push_back( new_plane );
    }

}

int PlaneIO::write(Planes &planes, std::string &d)
{
    std::string dir( d.c_str() );
    dir = dir + "hulls";
    #ifdef _WIN32
    int success = mkdir(dir.c_str());
    #else
    int success = mkdir(dir.c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
    #endif

    if(success != 0 && errno != EEXIST) {
        cerr << "Creating directory " << dir << " failed" << endl;
        return -1;
    }
    string listfile = dir + "/planes.list";
    ofstream out;
    out.open(listfile.c_str());
    for (uint i = 0; i < planes.size(); i++)
    {
        string planefile = dir + "/plane" + to_string(i, 3) + ".3d";
        planes.at(i)->writePlane(planefile);
        out << "Plane " << planefile << endl;
    }
    return 0;
}