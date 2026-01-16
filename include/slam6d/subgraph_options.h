#ifndef SGICP_OPTIONS
#define SGICP_OPTIONS

#include <string>
#include <boost/program_options.hpp>
#include "slam6d/io_types.h"

namespace po = boost::program_options;

// Boost needs to convert from string to IOType. Use this validation function.
void validate(boost::any& v, const std::vector<std::string>& values,
              IOType*, int) {
  if (values.size() == 0)
    throw std::runtime_error("Invalid model specification");
  std::string arg = values.at(0);
  try {
    v = formatname_to_io_type(arg.c_str());
  } catch (...) { // runtime_error
    throw std::runtime_error("Format " + arg + " unknown.");
  }
}

int parse_options(
    int argc,
    char** argv,
    std::string &scandir,
    int& start,
    int& end,
    IOType &type,
    bool &quiet,
    int &maxDist,
    int &minDist,
    int &octree,
    double &red,
    bool &scanserver,
    bool &use_frames,
    double &mdm,
    double &mdml,
    int& iter,
    int& iter_lum,
    double& epsICP,
    double& epsSLAM,
    size_t& subsize,
    int& clpairs,
    bool& meta,
    int& max_num_metascans,
    bool& icp_only
    )
{
    // Add program option descriptions and reference variables accordingly
    po::options_description generic("Generic options");
    po::options_description input("Program options");
    po::options_description reduction("Reduction options");
    po::options_description matching("Matching options");
    po::options_description hidden("Positional arguments");

    // Add program options in a curried way, using the overloaded ()-operator
    generic.add_options()
        ("help,h", "Display a very helpful message")
        ("quiet,q", po::bool_switch(&quiet)->default_value(false),
            "Supress any output (except segfault) of the program.");
    input.add_options()
        ("format,f", po::value<IOType>(&type)->default_value(UOS, "uos"),
            "chose input format from {uos, uosr, uos_rgb, ous_normal, old, xyz}")
        ("start,s", po::value<int>(&start)->default_value(0),
            "Skip the first <arg> scans.")
        ("end,e", po::value<int>(&end)->default_value(-1),
            "Stop at scan index <arg>")
        ("scanserver", po::value<bool>(&scanserver)->default_value(false),
            "Use the scanserver as an input method and handling of scan data")
        ("continue", po::bool_switch(&use_frames)->default_value(false),
            "Use pose specified in .frames file instead of .pose file.");
    reduction.add_options()
        ("max,m", po::value<int>(&maxDist)->default_value(-1),
            "neglegt all data points with a distance larger than <arg> 'units'")
        ("min,M", po::value<int>(&minDist)->default_value(-1),
            "neglegt all data points with a distance smaller than <arg> 'units'")
        ("reduce,r", po::value<double>(&red)->default_value(-1.0),
            "turns on octree based point reduction (voxel size= <arg>)")
        ("octree,O", po::value<int>(&octree)->default_value(0),
            "Use randomized octree based point reduction (pts per voxel=<arg>)");
    matching.add_options()
        ("dist,d", po::value<double>(&mdm)->default_value(25.0),
            "Minimum distance for considering valid correspondence for ICP")
        ("distSLAM,D", po::value<double>(&mdml)->default_value(5.0),
            "Minimum distance for considering valid correspondence for GraphSLAM")
        ("iter,i", po::value<int>(&iter)->default_value(100),
            "Maximum number of iterations.")
        ("iterSLAM,I", po::value<int>(&iter_lum)->default_value(50),
            "sets the maximal number of iterations for SLAM to <NR>"
            "(if not set, graphSLAM is not executed)")
        ("epsICP,E", po::value<double>(&epsICP)->default_value(0.00001),
            "Stop ICP matching if difference is smaller than <arg>")
        ("epsSLAM", po::value<double>(&epsSLAM)->default_value(0.05),
            "Stop Graph SLAM if difference is smaller than <arg>")
        ("clpairs,c", po::value<int>(&clpairs)->default_value(20),
            "Consider an edge in the graph if more than <arg> point pairs exist")
        ("size,S", po::value<size_t>(&subsize)->default_value(5),
            "Number of scans in the subgraphs.")
        ("metascan,2", po::bool_switch(&meta)->default_value(false),
            "Match current scan against a meta scan of all previous scans (default match against the last scan only)")
        ("maxmeta", po::value<int>(&max_num_metascans)->default_value(-1),
            "maximum nr of previous scans to combine to a metascan in scan matching. The total number will be <arg>*<size>.")
        ("icp", po::bool_switch(&icp_only)->default_value(false),
            "Subgraphs perform ICP matching, no graphSLAM");
    hidden.add_options()
        ("input-dir", po::value<string>(&scandir), "input dir");

    // All options together
    po::options_description alloptions;
    alloptions.add(generic).add(input).add(reduction).add(matching).add(hidden);

    // Only commandline, visible with --help
    po::options_description cmdoptions;
    cmdoptions.add(generic).add(input).add(reduction).add(matching);

    // positional argument for input directory
    po::positional_options_description pos;
    pos.add("input-dir", 1); // max 1 pos arg

    // Map and store option inputs to variables
    po::variables_map vars;
    po::store( po::command_line_parser(argc, argv).
                options(alloptions).
                positional(pos).
                run(),
                vars);

    // help display msg
    if ( vars.count("help") )
    {
        cout << cmdoptions;
        cout << endl;
        cout << "This is an experimental program which performs scan registration." << endl;
        cout << "The program puts S (--size) subsequent scans into a Graph." << endl;
        cout << "First, GraphSLAM is performed on these so called \"Sub-Graphs\" of size S." << endl;
        cout << "GraphSLAM or plain ICP is then performed on the subsequent Sub-Graphs." << endl;
        cout << endl << "Example usage:" << endl
            << "bin/sgicp dat/yourscans -f uosr -r 10 -O 1 -d 50 -S 5"
            << endl;
        exit(0);
    }
    po::notify(vars);

    // Add trailing directory slash if there is none. Works differently when compiling under Windows
#ifndef _MSC_VER
    if (scandir[ scandir.length()-1 ] != '/') scandir = scandir + "/";
#else
    if (scandir[ scandir.length()-1]  != '\\') scandir = scandir + "\\";
#endif

    return 0;
}

#endif // SGICP_OPTIONS
