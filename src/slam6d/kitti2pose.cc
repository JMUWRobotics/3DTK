/*
 * kitti2pose implementation
 *
 * Copyright (C) Shitong Du
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

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif
double deg(const double rad)
{
  return ( (rad * 360.0) / (2.0 * M_PI) );
}
 void Matrix4ToEuler( double *alignxf,
                                  double *rPosTheta,
                                  double *rPos = 0)
{

  double _trX, _trY;

  // Calculate Y-axis angle
  if(alignxf[0] > 0.0) {
    rPosTheta[1] = asin(alignxf[8]);
  } else {
    rPosTheta[1] = M_PI  - asin(alignxf[8]);
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

  rPosTheta[0] = deg(rPosTheta[0]);
  rPosTheta[1] = deg(rPosTheta[1]);
  rPosTheta[2] = deg(rPosTheta[2]);
 if (rPos != 0) {
    rPos[0] = alignxf[12];
    rPos[1] = alignxf[13];
    rPos[2] = alignxf[14];

  }
}

int parseArgs(int argc,char **argv, char dir[255], int& start, int& end,int& poseindex){
  start   = 0;
  end     = -1; // -1 indicates no limitation
  poseindex=0;
  int  c;
  // from unistd.h
  extern char *optarg;
  extern int optind;

  cout << endl;
  while ((c = getopt (argc, argv, "s:e:p:")) != -1)
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
   case 'p':
     poseindex = atoi(optarg);
     if (poseindex < 0)     { cerr << "Error: Cannot end at a negative scan number.\n"; exit(1); }
     if (poseindex > 10) { cerr << "Error: There are only 10 pose sequences.\n"; exit(1); }
     break;
   }

  if (optind != argc-1) {
    cerr << "\n*** Directory missing ***\n" << endl;
    cout << endl
	  << "Usage: " << argv[0] << "  [-s NR] [-e NR] directory" << endl << endl;

    cout << "  -s NR   start at scan NR (i.e., neglects the first NR scans)" << endl
       << "          [ATTENTION: counting starts with 0]" << endl
	  << "  -e NR   end after scan NR" << "" << endl
            << "  -p NR   read ground truth sequence NR" << "" << endl
	  << endl;
    cout << "Reads ground truth files from directory/??.txt and converts them to directory/scan???.pose." << endl;
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
  int start = 0, end = -1,poseindex=0;
  char dir[255];
  parseArgs(argc, argv, dir, start, end,poseindex);
  double inMatrix[16], tMatrix[16];
  double rPos[3], rPosTheta[16];
  int  fileCounter = start;
  char poseFileName[255];
  char truthFileName[255];

  ifstream pose_in;
  ofstream pose_out;
  snprintf(truthFileName,255,"%s%.2d.txt",dir,poseindex);
  pose_in.open(truthFileName);
  if (!pose_in.good())  // no more files in the directory
  { cerr << "Error: no more files in the directory.\n"; exit(1); }

  for (;;) {
    if (end > -1 && fileCounter > end) break; // 'nuf read

    snprintf(poseFileName,255,"%sscan%.3d.pose",dir,fileCounter++);

    cout << "Reading frame " << truthFileName << "..." << endl;

    bool validPose = true;

    for (int j = 0; j < 12; ++j){
      if (!(pose_in >> inMatrix[j])){
        validPose = false;
        break;
      }
    }
    if (!validPose){
      break;
    }

    inMatrix[12] = 0.0;
    inMatrix[13] = 0.0;
    inMatrix[14] = 0.0;
    inMatrix[15] = 1.0;


     tMatrix[ 0] = inMatrix[0];
     tMatrix[ 4] =-inMatrix[1];
     tMatrix[ 8] = inMatrix[2];
     tMatrix[ 12] = inMatrix[3];
     tMatrix[ 1] = -(inMatrix[4]);
     tMatrix[ 5] =(inMatrix[5]);
     tMatrix[ 9] = -(inMatrix[6]);
     tMatrix[ 13] = -(inMatrix[7]);
     tMatrix[2] = inMatrix[8];
     tMatrix[ 6] =- inMatrix[9];
     tMatrix[10] = inMatrix[10];
     tMatrix[14] =inMatrix[11];
     tMatrix[3] =inMatrix[12];
     tMatrix[7] =-inMatrix[13];;
     tMatrix[11] =inMatrix[14];
     tMatrix[15] =inMatrix[15];

     Matrix4ToEuler(tMatrix, rPosTheta, rPos);

    pose_out.open(poseFileName);

    cout << "Writing pose file... " << poseFileName << endl;

    for(int i = 0; i < 3; i++) {
      pose_out << rPos[i]*100.0 << " ";
    }
    pose_out << endl;

    for(int i = 0; i < 3; i++) {
      pose_out <<rPosTheta[i]<< " ";
    }
    pose_out << endl;

    pose_out.close();
    pose_out.clear();


    cout << " done." << endl;
  }
     pose_in.close();
     pose_in.clear();
}
