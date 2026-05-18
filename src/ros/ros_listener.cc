#include <iostream>
#include <fstream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <omp.h>
#include <chrono>

#include <ros/ros.h>
#include <ros/spinner.h>
#include <std_msgs/Header.h>
#include <sensor_msgs/PointCloud2.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_listener.h>

#include <pcl_conversions/pcl_conversions.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/point_types.h>

#include <sys/types.h>
#include <sys/stat.h>

// checks if the given path exists.
int existsDir(const char* path)
{
    struct stat info;
    if (stat( path, &info ) != 0) return 0;
    else if ( info.st_mode & S_IFDIR ) return 1;
    else return 0;
}

using namespace std;
using namespace chrono;
namespace po = boost::program_options;

const char* PATH; // Current directory
const int PATH_CHAR_LEN = 5000; // what a shit define... maximum buffer size for output path
const double conv_rad2deg = 57.29577951308232;

bool skip_points = false;
bool export_intensity = false;
double min_dist2;
double max_dist2;
string map_frame = "";

uint seq; // this counts the scanfiles, starting from 000, 001, 002, ... , 999, 1000, 1001,...
bool firstLidarCallback = false;
bool verbose;

string lidar_frame = "";
tf::TransformListener* tf_listener = nullptr;

int parse_options (int argc, char** argv,
  string& outdir,
  string& lidar_topic,
  string& map_frame,
  bool& skip_points,
  bool& export_intensity,
  bool& verbose,
  double& minDist,
  double& maxDist
)
{
  po::options_description generic("Generic options");
  po::options_description input("Program options");
  po::options_description hidden("Hidden options");
  generic.add_options()
    ("help,h", "Display a very helpful message");
  input.add_options()
    ("lidar_topic,L", po::value<string>(&lidar_topic)->default_value("/lidar"),
     "Provide the name of the topic where LiDAR data gets published.\n"
     "Currently, only sensor_msgs::PointCloud2 format is supported.")
    ("map_frame,F", po::value<string>(&map_frame)->default_value("map"),
     "Name of the fixed (world) frame in the TF tree. The LiDAR frame will be\n"
     "looked up relative to this frame to obtain the scanner pose.")
    ("minDist,m", po::value<double>(&minDist)->default_value(0),
      "Ignore points closer to <arg> cm")
    ("maxDist,M", po::value<double>(&maxDist)->default_value(std::numeric_limits<double>::max()),
      "Ignore points further than <arg> cm")
    ("skip_points,s", po::bool_switch(&skip_points)->default_value(false),
    "Exports only pose data. Use if re-run on a bagfile where points have already been exported.")
    ("intensity,I", po::bool_switch(&export_intensity)->default_value(false),
    "Export intensity as a 4th column. If not set, only XYZ are exported.")
    ("verbose,v", po::bool_switch(&verbose)->default_value(false),
     "Makes this program talk more. Use to print debug information.");
  hidden.add_options()
    ("output-dir", po::value<string>(&outdir), "output-dir");
  // All options together
    po::options_description alloptions;
    alloptions.add(generic).add(input).add(hidden);

    // Only commandline, visible with --help
    po::options_description cmdoptions;
    cmdoptions.add(generic).add(input);

    // positional argument for input directory
    po::positional_options_description pos;
    pos.add("output-dir", 1); // max 1 pos arg

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
        cout << endl << "Example usage:" << endl
           << "\t bin/ros_listener dat/your/out/dir --lidar_topic=/livox/lidar --map_frame=map" << endl;
        exit(0);
    }
    po::notify(vars);

    // Add trailing directory slash if there is none.
    // Works differently when compiling under Windows
#ifndef _MSC_VER
    if (outdir[ outdir.length()-1 ] != '/') outdir = outdir + "/";
#else
    if (outdir[ outdir.length()-1]  != '\\') outdir = outdir + "\\";
#endif

    // return with success exit code
    return 0;
}

/**
 * Converts a right-hand-side matrix into a 3DTK matrix
 * @param *inMatrix pointer to matrix (double[16])
 * @param *outMatrix pointer to matrix (double[16])
 * @param scale used for unit conversion, default 100.0 for Riegl
 */
inline void to3DTKMat(const double *inMatrix,
				  double *outMatrix, float scale = 100.0)
{
    outMatrix[0] = inMatrix[5];
    outMatrix[1] = -inMatrix[9];
    outMatrix[2] = -inMatrix[1];
    outMatrix[3] = -inMatrix[13];
    outMatrix[4] = -inMatrix[6];
    outMatrix[5] = inMatrix[10];
    outMatrix[6] = inMatrix[2];
    outMatrix[7] = inMatrix[14];
    outMatrix[8] = -inMatrix[4];
    outMatrix[9] = inMatrix[8];
    outMatrix[10] = inMatrix[0];
    outMatrix[11] = inMatrix[12];
    outMatrix[12] = -scale*inMatrix[7];
    outMatrix[13] = scale*inMatrix[11];
    outMatrix[14] = scale*inMatrix[3];
    outMatrix[15] = inMatrix[15];
}

static inline void Matrix4ToEuler(const double *alignxf,
                                  double *rPosTheta,
                                  double *rPos = 0)
{
  double _trX, _trY;

  // Calculate Y-axis angle
  if(alignxf[0] > 0.0) {
    rPosTheta[1] = asin(alignxf[8]);
  } else {
    rPosTheta[1] = M_PI - asin(alignxf[8]);
  }

  double  C    =  cos( rPosTheta[1] );
  if ( fabs( C ) > 0.005 )  {                 // Gimbal lock?
    _trX      =  alignxf[10] / C;             // No, so get X-axis angle
    _trY      =  -alignxf[9] / C;
    rPosTheta[0]  = atan2( _trY, _trX );
    _trX      =  alignxf[0] / C;              // Get Z-axis angle
    _trY      = -alignxf[4] / C;
    rPosTheta[2]  = atan2( _trY, _trX );
  } else {                                    // Gimbal lock has occurred
    rPosTheta[0] = 0.0;                       // Set X-axis angle to zero
    _trX      =  alignxf[5];  //1                // And calculate Z-axis angle
    _trY      =  alignxf[1];  //2
    rPosTheta[2]  = atan2( _trY, _trX );
  }

  rPosTheta[0] = rPosTheta[0];
  rPosTheta[1] = rPosTheta[1];
  rPosTheta[2] = rPosTheta[2];

  if (rPos != 0) {
    rPos[0] = alignxf[12];
    rPos[1] = alignxf[13];
    rPos[2] = alignxf[14];
  }
}

void lidarMsgCallback(const sensor_msgs::PointCloud2::ConstPtr& msg)
{
  if (verbose) ROS_INFO("Callback LIDAR");

  // First Callback: record the LiDAR frame_id (done only once)
  if (!firstLidarCallback) {
    std::string info_string = std::string("Writing files to ")
                                    + std::string(PATH);
    ROS_INFO("%s", info_string.c_str());
    lidar_frame = msg->header.frame_id;
    ROS_INFO("LiDAR frame: \"%s\", map frame: \"%s\"", lidar_frame.c_str(), map_frame.c_str());
    firstLidarCallback = true;
  }

  // Opening 3d file to write into
  FILE* file_3d;
  char* file_name = new char[PATH_CHAR_LEN]();
  std::sprintf(file_name, "%sscan%03d.3d", PATH, seq);

  auto start_clock_pts = high_resolution_clock::now();
  if (!skip_points) {
    pcl::PCLPointCloud2 pcl_pc2;
    pcl_conversions::toPCL(*msg,pcl_pc2);
    bool has_intensity = false;
    for (size_t i = 0; i < msg->fields.size(); ++i) {
      if (msg->fields[i].name == "intensity") {
        has_intensity = true;
        break;
      }
    }

    // Writing the lidar data to the text file
    file_3d = fopen(file_name, "wb");;
    // This is a left handed coordinate system and we convert the values to cm.
    if (export_intensity && has_intensity) {
      pcl::PointCloud<pcl::PointXYZI>::Ptr temp_cloud(new pcl::PointCloud<pcl::PointXYZI>);
      pcl::fromPCLPointCloud2(pcl_pc2, *temp_cloud);
      for (size_t i = 0; i < temp_cloud->points.size(); ++i)
      {
          pcl::PointXYZI p = temp_cloud->points[i];
          double dist2 = 100*100*(p.x*p.x+p.y*p.y+p.z*p.z);
          if (std::isnan(dist2) || dist2 < min_dist2 || dist2 > max_dist2 )
            continue;
          // Convert to left handed and export XYZ + intensity
          fprintf(file_3d, "%lf %lf %lf %lf\n", 100*-p.y, 100*p.z, 100*p.x, p.intensity);
      }
    } else {
      if (export_intensity && !has_intensity) {
        ROS_WARN_THROTTLE(5.0, "Input cloud has no intensity field. Exporting XYZ only.");
      }
      pcl::PointCloud<pcl::PointXYZ>::Ptr temp_cloud(new pcl::PointCloud<pcl::PointXYZ>);
      pcl::fromPCLPointCloud2(pcl_pc2, *temp_cloud);
      for (size_t i = 0; i < temp_cloud->points.size(); ++i)
      {
          pcl::PointXYZ p = temp_cloud->points[i];
          double dist2 = 100*100*(p.x*p.x+p.y*p.y+p.z*p.z);
          if (std::isnan(dist2) || dist2 < min_dist2 || dist2 > max_dist2 )
            continue;
          // Convert to left handed and export XYZ only
          fprintf(file_3d, "%lf %lf %lf\n", 100*-p.y, 100*p.z, 100*p.x);
      }
    }
    fclose(file_3d);
  }
  auto stop_clock_pts = high_resolution_clock::now();

  // Look up the transform from the LiDAR frame into the fixed map frame via TF
  tf::StampedTransform lidar_to_map;
  try {
    tf_listener->waitForTransform(map_frame, lidar_frame, msg->header.stamp, ros::Duration(1.0));
    tf_listener->lookupTransform(map_frame, lidar_frame, msg->header.stamp, lidar_to_map);
  } catch (tf::TransformException& ex) {
    ROS_WARN("TF lookup failed, skipping scan %d: %s", seq, ex.what());
    delete[] file_name;
    return;
  }

  // Build homogeneous transformation matrix from TF result
  tf::Matrix3x3 m(lidar_to_map.getRotation());
  const double in_matrix[16] = {
    m[0][0], m[0][1], m[0][2], lidar_to_map.getOrigin().getX(),
    m[1][0], m[1][1], m[1][2], lidar_to_map.getOrigin().getY(),
    m[2][0], m[2][1], m[2][2], lidar_to_map.getOrigin().getZ(),
    0,       0,       0,       1
  };

  // Converting to left handed matrix for 3DTK (OpenGL style)
  double out_matrix[16], rPos[3], rPosTheta[16];
  to3DTKMat(in_matrix, out_matrix, 1);
  Matrix4ToEuler(out_matrix, rPosTheta, rPos);

  // Extracting Position and Orientaion
  double x = 100.0 * rPos[0]; // left-handed format in 3DTK uses cm instead of m
  double y = 100.0 * rPos[1];
  double z = 100.0 * rPos[2];
  double roll  = 1.0 * rPosTheta[0];
  double pitch = 1.0 * rPosTheta[1];
  double yaw   = 1.0 * rPosTheta[2];

  // Writing to the pose file
  // This is also a left hand coordinate system
  // - Thumb (x) to the right
  // - Pointy Finger (y) to  the top
  // - Middle Finger (z) into the room
  auto duration_pts = duration_cast<milliseconds>(stop_clock_pts - start_clock_pts);
  ROS_INFO("Writing %s ... [%d ms]", file_name, duration_pts);
  std::sprintf(file_name, "%sscan%03d.pose", PATH, seq);
  FILE* file_pose = fopen(file_name, "wb");
  fprintf(file_pose, "%lf %lf %lf %lf %lf %lf", x, y, z, roll*conv_rad2deg, pitch*conv_rad2deg, yaw*conv_rad2deg);
  fclose(file_pose);
  seq++;
  delete[] file_name;
}

int main(int argc, char **argv)
{
    // Program information
    ros::init(argc, argv, "file_writer");
    std::cout << "ros_listener - A program that listens to ROS topics" << std::endl;
    std::cout << "               and writes 3DTK files (UOSR-format) on disk." << std::endl;
    std::cout << "               Use -h for more information." << std::endl;
    std::cout << "Note: Remember that you need a running roscore" << std::endl;
    std::cout << "Tip of the day: Start this program first, it will wait for the topics." << std::endl;

    // Declaration of program parameters
    string outdir;
    string lidar_topic;
    double min_dist;
    double max_dist;

    // Definition and allocation of parameters
    parse_options(argc, argv,
        outdir, lidar_topic, map_frame, skip_points,
      export_intensity, verbose, min_dist, max_dist);

    min_dist2 = min_dist*min_dist;
    max_dist2 = max_dist*max_dist;

    // Create path on disk (if it does not exist yet)
    PATH = outdir.c_str();
    if ( !existsDir( PATH ) ) {
      ROS_WARN("\"%s\" does not exist, I will create it now.", PATH);
      ROS_INFO("Creating \"%s\".", PATH);
      if (!boost::filesystem::create_directory( PATH )) {
        ROS_ERROR("Could not create \"%s\". Aborting now.", PATH);
        return -1; // Unusual program end
      }
    } else ROS_INFO("\"%s\" exists, exporting there.", PATH);

    // Init iteration variable, count number of LiDAR callbacks
    seq = 0; // counts nrscans

    // Create the TF listener (must live for the duration of the program)
    tf::TransformListener listener;
    tf_listener = &listener;

    // Lidar Node Handle
    ros::NodeHandle nh_lidar;
    ros::Subscriber lidar_sub = nh_lidar.subscribe(lidar_topic, 100000, lidarMsgCallback);

    // Spin
    ros::spin();

    return 0; // normal program end
}
