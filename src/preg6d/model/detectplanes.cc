#include "model/detectplanes.h"

Detect::Detect(Scan* scan, std::string path, bool quiet, std::string config)
{
    this->scan = scan;
    if (path[ path.size() - 1] != '/')
        path += "/";
    this->outpath = path;
    detector = new Hough(scan, quiet, config);
}

void Detect::detect(plane_alg algo)
{
    if (!detector) {
        cout << "Attention! This detector has already been used for scan"
             << scan->getIdentifier() << endl;
        return;
    }
    if (algo >= RANSAC) {
        clusterRansac();
    } else switch (algo) {
        case RHT:
            detector->RHT();
            break;
        case SHT:
            detector->SHT();
            break;
        case PHT:
            detector->PHT();
            break;
        case PPHT:
            detector->PPHT();
            break;
        case APHT:
            detector->APHT();
            break;
        default:
            cout << "This should never happen." << endl;
            break;
    }
    detector->writeColoredPoints( outpath+"scan"+scan->getIdentifier()+".3d" );
    ofstream ofilepose(outpath+"scan"+scan->getIdentifier()+".pose", ios_base::out);
    ofilepose << scan->get_rPos()[0] << " "
              << scan->get_rPos()[1] << " "
              << scan->get_rPos()[2] << " ";
    ofilepose << deg(scan->get_rPosTheta()[0]) << " "
              << deg(scan->get_rPosTheta()[1]) << " "
              << deg(scan->get_rPosTheta()[2]) << endl;
}

void Detect::clusterRansac()
{
    DataXYZ xyz(Scan::allScans[0]->get("xyz"));
    RansacOctTree<double>* oct = new RansacOctTree<double>(PointerArray<double>(xyz).get(), xyz.size(), 50.0 );

    unsigned int stop = (unsigned int)(detector->allPoints->size()/100.0)*detector->myConfigFileHough.Get_MinSizeAllPoints();
    int counter = 0;
    while(detector->allPoints->size() > stop &&
            detector->planes.size() < (unsigned int)detector->myConfigFileHough.Get_MaxPlanes() &&
            counter++ < (int)detector->myConfigFileHough.Get_TrashMax())
    {

        vector<double *> points;
        CollisionPlane<double> * plane;
        plane = new CollisionPlane<double>(1.0); // 1.0 cm maxdist
        Ransac(*plane, oct, &points);

        cout << "nr points " << points.size() << endl;
        double nx,ny,nz,d;
        plane->getPlane(nx,ny,nz,d);
        if(!detector->quiet) cout << "DONE " << endl;

        if(!detector->quiet) cout << nx << " " << ny << " " << nz << " " << d << endl;

        double * normal = new double[4];
        normal[0] = nx;
        normal[1] = ny;
        normal[2] = nz;

        for(int i = points.size() - 1; i > -1; i--) {
            delete[] points[i];
        }

        detector->deletePoints(normal, -d);
        delete plane;
        delete normal;
    }
}

Detect::~Detect()
{
    delete detector;
}