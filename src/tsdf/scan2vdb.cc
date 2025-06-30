#include "slam6d/scan.h"
#include <cppunit/extensions/HelperMacros.h>
#include <openvdb/openvdb.h>
#include <openvdb/Exceptions.h>
#include <openvdb/Types.h>
#include <openvdb/tree/LeafNode.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/Types.h>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

class ParticleList
{
protected:
    struct Particle {
        openvdb::Vec3R p;
        openvdb::Real  r;
    };
    openvdb::Real _radiusScale;
    vector<Particle> _particleList;
public:
    typedef openvdb::Vec3R value_type;
    typedef openvdb::Vec3R PosType;
  
    ParticleList(openvdb::Real radiusScale=10) : _radiusScale(radiusScale) {}
  //    ParticleList(openvdb::Real radiusScale=1) : _radiusScale(radiusScale) {}  
    void add(const openvdb::Vec3R &p, const openvdb::Real &r) {
        Particle pa;
        pa.p = p;
        pa.r = r;
        _particleList.push_back(pa);
    }
    openvdb::Vec3R pos(int n) const { return _particleList[n].p; }
    openvdb::Real radius(int n) const { return _radiusScale * _particleList[n].r; }
    size_t size() const { return _particleList.size(); }
    void getPos(size_t n,  openvdb::Vec3R& pos) const { pos = _particleList[n].p; }
    void getPosRad(size_t n,  openvdb::Vec3R& pos, openvdb::Real& rad) const {
        pos = _particleList[n].p;
        rad = _radiusScale * _particleList[n].r;
    }
};

void validate(boost::any& v, const vector<string>& values, IOType*, int) {
    if (values.size() == 0)
        throw runtime_error("Invalid model specification");
    string arg = values.at(0);
    try {
        v = formatname_to_io_type(arg.c_str());
    } catch (...) { // runtime_error
        throw runtime_error("Format " + arg + " unknown.");
    }
}

int main(int argc, char **argv)
{
    string inDir;
    string outDir;
    IOType iotype;
    int startIndex;
    int endIndex;
    float voxelSize;
    float radiusSphere;

    po::options_description generic("Generic options");
    generic.add_options()
            ("help,h", "output this help message");

    po::options_description input("Input options");
    input.add_options()
            ("start,s", po::value<int>(&startIndex)->default_value(0),
             "start at scan <arg> (i.e., neglects the first <arg> scans) "
             "[ATTENTION: counting naturally starts with 0]")
            ("end,e", po::value<int>(&endIndex)->default_value(-1),
             "end after scan <arg>")
            ("radiusscale,r", po::value<float>(&radiusSphere)->default_value(1.0),
             "sets the radius of the spheres <arg>")
      
            ("format,f", po::value<IOType>(&iotype)->default_value(UOS),
             "using shared library <arg> for input. (chose F from {uos, uos_map, "
             "uos_rgb, uos_frames, uos_map_frames, old, rts, rts_map, ifp, "
             "riegl_txt, riegl_rgb, riegl_bin, zahn, ply})");

    po::options_description output("Output options");
    output.add_options()
            ("output,O", po::value<string>(&outDir)->default_value("."), "Output dir");

    po::options_description tsdf("TSDF options");
    tsdf.add_options()
            ("voxelsize,v", po::value<float>(&voxelSize)->default_value(0.1),
             "voxel size <arg>");

    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("input-dir", po::value<string>(&inDir), "input dir");

    po::options_description all;
    all.add(generic).add(input).add(output).add(tsdf).add(hidden);

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(input).add(output);

    po::positional_options_description pd;
    pd.add("input-dir", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(all).positional(pd).run(), vm);

    po::notify(vm);

    if (vm.count("help")) {
        cout << cmdline_options;
        return 0;
    }

    openvdb::initialize();

    const float halfWidth = 2.0f;
    openvdb::FloatGrid::Ptr ls = openvdb::createLevelSet<openvdb::FloatGrid>(voxelSize, halfWidth);
    openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> raster(*ls);

    Scan::openDirectory(false, inDir, iotype, startIndex, endIndex);
    if(Scan::allScans.size() == 0) {
        cerr << "No scans found. Did you use the correct format?" << endl;
        exit(-1);
    }

    for (ScanVector::iterator it = Scan::allScans.begin(); it != Scan::allScans.end(); ++it) {
        Scan* scan = *it;
        scan->setRangeFilter(-1, -1);
        scan->toGlobal();

        DataXYZ xyz(scan->get("xyz reduced"));

        ParticleList particleList(radiusSphere);

        for (uint i = 0; i < xyz.size(); i++) {
            particleList.add(openvdb::Vec3R(xyz[i][2] / 100.0, xyz[i][1] / 100.0, xyz[i][0] / 100.0), voxelSize * halfWidth);
        }

        cout << "Particles: " << particleList.size() << endl;

        raster.setGrainSize(4);
        raster.rasterizeSpheres(particleList);
        raster.finalize();

        ls->pruneGrid(0.001 * voxelSize);

        delete *it;
        *it = 0;
    }

    Scan::closeDirectory();

    cout << "Active voxel count: " << ls->activeVoxelCount() << endl;

    openvdb::io::File file(outDir + "/grid.vdb");
    openvdb::GridPtrVec grids;
    grids.push_back(ls);
    file.write(grids);
    file.close();

    return 0;
}
