/*
 * frames2tum implementation
 * Converts frames to tum format as in
 * https://github.com/MichaelGrupp/evo/wiki/Formats
 *
 *  timestamp x y z qx qy qz qw
 *
 * Copyright (C) Fabian Arzberger
 *
 * Released under the GPL version 3.
 *
 */
#include <fstream>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ios;
#include "slam6d/globals.icc"
#include <string.h>

#ifndef _MSC_VER
#include <unistd.h>
#else
#include "XGetopt.h"
#endif

#if WIN32
#define snprintf sprintf_s
#endif

int parseArgs(int argc,char **argv, char dir[255], int& start, int& end,int& sequence,bool& zero, bool& trustPose){
  start   = 0;
  end     = -1; // -1 indicates no limitation
  sequence=0;
  int  c;
  // from unistd.h
  extern char *optarg;
  extern int optind;

  cout << endl;
  while ((c = getopt (argc, argv, "s:e:q:zp")) != -1)
    switch (c)
   {
   case 's':
     start = atoi(optarg);
     if (start < 0) { cerr << "Error: Cannot start at a negative scan number.\n"; exit(1); }
     break;
   case 'e':
     end = atoi(optarg);
     if (end < 0)     { cerr << "Error: Cannot end at a negative scan number.\n"; exit(1); }
     if (end < start) { cerr << "Error: <end> cannot be smaller than <start>.\n"; exit(1); }
     break;

   case 'q':
     sequence = atoi(optarg);
     if (sequence < 0)     { cerr << "Error: Cannot end at a negative scan number.\n"; exit(1); }
     break;

   case 'z':
     zero = true;
     break;
   case 'p':
     trustPose = true;
     break;
   }

  if (optind != argc-1) {
    cerr << "\n*** Directory missing ***\n" << endl;
    cout << endl
	  << "Usage: " << argv[0] << "  [-s NR] [-e NR] directory" << endl << endl;

    cout << "  -s NR   start at scan NR (i.e., neglects the first NR scans)" << endl
       << "          [ATTENTION: counting starts with 0]" << endl
	  << "  -e NR   end after scan NR" << "" << endl
    << "  -z      Start position at zero (0,0,0)" << endl
	  << endl;
    cout << "Reads frame files from directory/scan???.frames and converts them to a single file trajectory.txt in the TUM pose file format." << endl;
    abort();
  }
  strncpy(dir,argv[optind],255);

#ifndef _MSC_VER
  if (dir[strlen(dir)-1] != '/') strcat(dir,"/");
#else
  if (dir[strlen(dir)-1] != '\\') strcat(dir,"\\");
#endif
  return 0;
}


int main(int argc, char **argv)
{
  int start = 0, end =-1;
  int sequence=0;
  char dir[255];
  bool zero = false;
  bool trustPose = false;
  parseArgs(argc, argv, dir, start,end,sequence,zero,trustPose);

  if (trustPose) cout << "Ignoring .frames !" << endl;

  int  fileCounter = start;
  char poseFileName[255];
  char frameFileName[255];
  ifstream frame_in;
  ifstream pose_in;
  FILE* pose_out;
  snprintf(poseFileName,255,"%strajectory.txt",dir);
  pose_out = fopen(poseFileName, "a");
  double tMatrix[17];
  double pose[7];
  double offset[3];
  bool init = true;

  for (;;) {
    if (end > -1 && fileCounter > end) break; // 'nuf read
    if (!trustPose) {
      snprintf(frameFileName,255,"%sscan%.3d.frames",dir,fileCounter);
      frame_in.open(frameFileName);
    }
    snprintf(poseFileName,255,"%sscan%.3d.pose",dir,fileCounter++);
    pose_in.open(poseFileName);
    if (!pose_in.good() && !frame_in.good()) break; // no more files in the directory

    cout << ".";
    cout.flush();

    // Get transformation from .frames
    while(frame_in.good()) {
      for (unsigned int i = 0; i < 17; frame_in >> tMatrix[i++]);
    }

    // Get timestamp from pose
    while(pose_in.good()) {
      for (unsigned int i = 0; i < 7; pose_in >> pose[i++]);
    }

    // Cleaning buffer
    if (!trustPose) {
      frame_in.close();
      frame_in.clear();
    }
    pose_in.close();
    pose_in.clear();

    // Convert to quaternion
    double quat[4];
    double t[3] = {1, 2, 3};
    if (trustPose) {
      double rPos[3] = {pose[0], pose[1], pose[2]};
      double rTheta[3] = {rad(pose[3]), rad(pose[4]), rad(pose[5])};
      double matrix4x4[16];
      EulerToMatrix4(rPos, rTheta, matrix4x4);
      Matrix4ToQuat(matrix4x4, quat, t);
    }
    else
      Matrix4ToQuat(tMatrix, quat, t);

    // Rotate the whole trajectory
    double tdummy[3] = {0,0,0};
    double rot90[3] = {rad(0.0), rad(0.0), rad(0.0)};
    double rot90Matrix[16];
    double newMatrix[16];
    QuatToMatrix4(quat, t, tMatrix);
    EulerToMatrix4(tdummy, rot90, rot90Matrix);
    MMult(rot90Matrix, tMatrix, newMatrix);
    Matrix4ToQuat(newMatrix, quat, t);

    if (zero) {
      if (init) {
        cout << "Offset set at" << endl;
        init = false;
        for (int i = 0; i < 3; ++i)
          offset[i] = t[i];
        cout << offset[0] << " " << offset[1] << " " << offset[2] << endl;
      }
      for (int i = 0; i < 3; ++i)
        t[i] = t[i] - offset[i];
    }

    // Convert left-handed to right-handed
    fprintf(pose_out, "%lf %lf %lf %lf %lf %lf %lf %lf\n", pose[6], t[2]/100.0, -t[0]/100.0, t[1]/100.0, -quat[3], quat[1], -quat[2], quat[0]);
  }
  fclose(pose_out);
  cout << " done writing trajectory.txt" << endl;
}
