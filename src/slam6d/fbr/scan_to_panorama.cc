/*
 * scan_to_panorama implementation
 *
 * Copyright (C) Hamidreza Houshiar
 *
 * Released under the GPL version 3.
 *
 */

#include <stdio.h>
#include <fstream>
#include "slam6d/fbr/fbr_global.h"
#include "slam6d/fbr/scan_cv.h"
#include "slam6d/fbr/panorama.h"
#include "slam6d/fbr/feature.h"
#include "slam6d/fbr/feature_matcher.h"
#include "slam6d/fbr/registration.h"
#include "slam6d/fbr/feature_drawer.h"
#include "show/scancolormanager.h"
#include "show/show_Boctree.h"
#include "slam6d/point_type.h"
#include "show/show.h"
#include "slam6d/data_types.h"
#include "slam6d/Boctree.h"
#include "slam6d/basicScan.h"

using namespace std;
using namespace fbr;

struct information{
  string inDir;
  int start, end;
  IOType scanFormat;
  int panoramaWidth, panoramaHeight;
  scanner_type scannerType;
  double minReflectance, maxReflectance;
  int minHorizAngle, maxHorizAngle;
  int minVertAngle, maxVertAngle;
  panorama_map_method mapMethod;
  projection_method projectionMethod;
  panorama_format panoramaFormat;
  double panoramaFormatParam;
  int numberOfImages;
  double projectionParam;
  bool panoramaSizeOptimization;
  string  outDir;
  bool reflectance, color, range;
  bool normalizeRange;
  bool threeChannelRange;
  bool threeGrayscaleRange;
  bool oneGrayscaleRange;
  bool saveOct, loadOct;
} info;

void usage(int argc, char** argv){
  printf("\n");
  printf("USAGE: %s inDir -s start -e end \n", argv[0]);
  printf("\n");
  printf("\n");
  printf("\tOptions:\n");
  printf("\t\t-f scanFormat\t\t\t input scan file format [RIEGL_TXT|RXP|ALL SLAM6D SCAN_IO]\n");
  printf("\t\t-W panoramaWidth\t\t panorama image width\n");
  printf("\t\t-H panoramaHeight\t\t panorama image height\n");
  printf("\t\t-t scannerType \t\t\t scanner type[NONE | RIEGL | FARO | MANUAL]\n");
  printf("\t\t-b minReflectance \t\t Min Reflectance for manual reflectance normalization\n");
  printf("\t\t-B maxReflectance \t\t Max Reflectance for manual reflectance normalization\n");
  printf("\t\t-m minHorizAngle \t\t Scanner horizontal view minAngle \n");
  printf("\t\t-w maxHorizAngle \t\t Scanner horizontal view maxAngle \n");
  printf("\t\t-n minVertAngle \t\t Scanner vertical view minAngle \n");
  printf("\t\t-x maxVertAngle \t\t Scanner vertical view maxAngle \n");
  printf("\n");
  printf("\n");
  printf("\t\t-M mapMethod\t\t\t panorama map method [FARTHEST|EXTENDED|FULL]\n");
  printf("\t\t-p projectionMethod\t\t projection method [AZIMUTHAL|CONIC|CYLINDRICAL|EQUALAREACYLINDRICAL|EQUIRECTANGULAR|MERCATOR|MILLER|PANNINI|RECTILINEAR|STEREOGRAPHIC|ZAXIS]\n");
  printf("\t\t-F panoramaFormat\t\t panorama format [PNG|JPEG|JPEG2000|TIFF]\n");
  printf("\t\t-S panoramaFormatParam\t\t panorama format param panorama format related param mostly compression param\n");
  printf("\t\t-N numberOfImage\t\t number of Horizontal images used for some projections\n");
  printf("\t\t-P projectionParam\t\t special projection parameter (d for Pannini and r for stereographic)\n");
  printf("\t\t-i panoramaSizeOptimization \t Optimize the panorama image size based on projection \n");
  printf("\n");
  printf("\n");
  printf("\t\t-O outDir \t\t\t output directory if not stated same as input\n");
  printf("\t\t-R reflectance \t\t\t uses reflactance and creates the Reflectance image (for rxp input, the scanner type (-t) has to be set to RIEGL)\n");
  printf("\t\t-C color \t\t\t uses color and creates the Color image\n");
  printf("\t\t-A range \t\t\t creates the Range image\n");
  printf("\t\t-a normalizeiRange \t\t normalize the Range image to have values between 0--255\n");
  printf("\t\t-c threeChannelRange \t\t puts the range(meter*10000) value in a 3*8 bit rgb image\n");
  printf("\t\t-G oneGrayscaleRange \t\t puts the range value(meter*10000) in one CVU16 grayscale image\n");
  printf("\t\t-g threeGrayscaleRange \t\t puts the range(meter*10000) value in three CVU8 grayscale images\n");
  printf("\t\t-l loadOct \t\t\t load the Octtree\n");
  printf("\t\t-o saveOct \t\t\t save the Octtree\n");
  printf("\n");
  printf("\tExamples:\n");
  printf("\t\t%s -s 4 -e 4 -f rxp -t RIEGL -R -F PNG -p EQUIRECTANGULAR wü_city", argv[0]);
  exit(1);
}

void parssArgs(int argc, char** argv, information& info){

  //default values
  info.scanFormat = RIEGL_TXT;
  info.panoramaWidth = 3600;
  info.panoramaHeight = 1000;
  info.scannerType = NONE;
  info.minReflectance = -100;
  info.maxReflectance = 100;
  info.minHorizAngle = 0;
  info.maxHorizAngle = 360;
  info.minVertAngle = -40;
  info.maxVertAngle = 60;
  info.mapMethod = fbr::FARTHEST;
  info.projectionMethod = EQUIRECTANGULAR;
  info.panoramaFormat = PNG;
  info.panoramaFormatParam = 3;
  info.numberOfImages = 1;
  //depend on the projection method
  info.projectionParam = 0;
  info.panoramaSizeOptimization = false;
  info.outDir = "";
  info.reflectance = false;
  info.color = false;
  info.range = false;
  info.normalizeRange = false;
  info.threeChannelRange = false;
  info.threeGrayscaleRange = false;
  info.oneGrayscaleRange = false;
  info.saveOct = false;
  info.loadOct = false;

  int c;
  opterr = 0;
  //reade the command line and get the options
  while ((c = getopt (argc, argv, "aAb:B:cCe:f:F:gGH:ilm:M:n:N:oO:p:P:Rs:S:t:w:W:x:")) != -1)
    switch (c)
      {
      case 'a':
	info.normalizeRange = true;
	info.range = true;
	break;
      case 'A':
	info.range = true;
	break;
      case 'b':
	info.minReflectance = atof(optarg);
	break;
      case 'B':
	info.maxReflectance = atof(optarg);
	break;
      case 'c':
	info.threeChannelRange = true;
	info.range = true;
	break;
      case 'C':
	info.color = true;
	break;
      case 'e':
	info.end = atoi(optarg);
	break;
      case 'f':
	info.scanFormat = stringToScanFormat(optarg);
	break;
      case 'F':
	info.panoramaFormat = stringToPanoramaFormat(optarg);
	break;
      case 'g':
	info.threeGrayscaleRange = true;
	info.range = true;
	break;
      case 'G':
	info.oneGrayscaleRange = true;
	info.range = true;
	break;
      case 'H':
        info.panoramaHeight = atoi(optarg);
	break;
      case 'i':
	info.panoramaSizeOptimization = true;
	break;
      case 'l':
	info.loadOct = true;
	break;
      case 'm':
	info.minHorizAngle = atoi(optarg);
	break;
      case 'M':
	info.mapMethod = stringToPanoramaMapMethod(optarg);
	break;
      case 'n':
	info.minVertAngle = atoi(optarg);
	break;
      case 'N':
	info.numberOfImages = atoi(optarg);
	break;
      case 'o':
	info.saveOct = true;
	break;
      case 'O':
	info.outDir = optarg;
	break;
      case 'p':
	info.projectionMethod = stringToProjectionMethod(optarg);
	break;
      case 'P':
	info.projectionParam = atof(optarg);
	break;
      case 'R':
	info.reflectance = true;
	break;
      case 's':
	info.start = atoi(optarg);
	break;
      case 'S':
	info.panoramaFormatParam = atof(optarg);
	break;
      case 't':
	info.scannerType = stringToScannerType(optarg);
	break;
      case 'w':
	info.maxHorizAngle = atoi(optarg);
	break;
      case 'W':
	info.panoramaWidth = atoi(optarg);
	break;
      case 'x':
	info.maxVertAngle = atoi(optarg);
	break;
      case '?':
	cout<<"Unknown option character "<<optopt<<endl;
	usage(argc, argv);
	break;
      default:
	usage(argc, argv);
      }

  if(info.projectionMethod == PANNINI && info.projectionParam == 0){
    info.projectionParam = 1;
  }
  if(info.projectionMethod == STEREOGRAPHIC && info.projectionParam == 0){
    info.projectionParam = 2;
  }

  if (optind > argc - 1)
    {
      cout<<"Too few input arguments. At least inDir and two scan numbers are required."<<endl;
      usage(argc, argv);
    }

  info.inDir = argv[optind];
  if(info.outDir.empty()) info.outDir = info.inDir;
  if(info.outDir.back() != '/') info.outDir += '/';
}

void printInfo(information info){
  cout<<"Input Dir= "<<info.inDir<<endl;
  cout<<"Start= "<<info.start<<endl;
  cout<<"End= "<<info.end<<endl;
  cout<<"Scan Format= "<<scanFormatToString(info.scanFormat)<<endl;
  cout<<"Width= "<<info.panoramaWidth<<endl;
  cout<<"Height= "<<info.panoramaHeight<<endl;
  cout<<"scannerType= "<<info.scannerType<<endl;
  cout<<"minReflectance= "<<info.minReflectance<<endl;
  cout<<"maxReflectance= "<<info.maxReflectance<<endl;
  cout<<"minHorizAngle= "<<info.minHorizAngle<<endl;
  cout<<"maxHorizAngle= "<<info.maxHorizAngle<<endl;
  cout<<"minVertAngle= "<<info.minVertAngle<<endl;
  cout<<"maxVertAngle= "<<info.maxVertAngle<<endl;
  cout<<"Map Method= "<<panoramaMapMethodToString(info.mapMethod)<<endl;
  cout<<"Projection Method= "<<projectionMethodToString(info.projectionMethod)<<endl;
  cout<<"Panorama Format= "<<panoramaFormatToString(info.panoramaFormat)<<endl;
  cout<<"Panorama Format Param= "<<info.panoramaFormatParam<<endl;
  cout<<"Number of Images= "<<info.numberOfImages<<endl;
  cout<<"Projection Param= "<<info.projectionParam<<endl;
  cout<<"Panorama Size Optimization= "<<info.panoramaSizeOptimization<<endl;
  cout<<"Out Dir= "<<info.outDir<<endl;
  cout<<"Reflectance= "<<info.reflectance<<endl;
  cout<<"Color= "<<info.color<<endl;
  cout<<"Range= "<<info.range<<endl;
  cout<<"Normalize Range= "<<info.normalizeRange<<endl;
  cout<<"One Grayscale Range= "<<info.oneGrayscaleRange<<endl;
  cout<<"Three Channel Range= "<<info.threeChannelRange<<endl;
  cout<<"Three Gray Range= "<<info.threeGrayscaleRange<<endl;
  cout<<"Save Oct= "<<info.saveOct<<endl;
  cout<<"Load Oct= "<<info.loadOct<<endl;
  cout<<"-------------------------------"<<endl<<endl;
}



int run(int argc, char** argv){

  parssArgs(argc, argv, info);

  printInfo(info);

  //opendirectory of the scans
  bool scanserver = false;
  for(int s = info.start; s <= info.end; s++)
    {
      scan_cv scan(info.inDir, s, info.scanFormat, scanserver, info.scannerType, info.loadOct, info.saveOct, info.reflectance, info.color, -1, -1, info.minReflectance, info.maxReflectance);

      scan.convertScanToMat();

      //init the panorama
      fbr::panorama pImage;
      pImage.init(info.panoramaWidth, info.panoramaHeight, info.projectionMethod, info.numberOfImages, info.projectionParam, info.mapMethod, scan.getZMin(), scan.getZMax(), info.minHorizAngle, info.maxHorizAngle, info.minVertAngle, info.maxVertAngle, info.panoramaSizeOptimization, info.reflectance, info.range, info.color);

      //create panorama
      pImage.createPanorama(scan.getMatScan(), scan.getMatScanColor());

      //get the new panorama image size incase of optimized panorama size
      info.panoramaWidth = pImage.getImageWidth();
      info.panoramaHeight = pImage.getImageHeight();

      //write panorama to file
      string out;
      vector<int> panoramaFormatParams;
      switch(info.panoramaFormat){
      case PNG:
	panoramaFormatParams.push_back(cv::IMWRITE_PNG_COMPRESSION);
	panoramaFormatParams.push_back(info.panoramaFormatParam);
	break;
      case JPEG:
	panoramaFormatParams.push_back(cv::IMWRITE_JPEG_QUALITY);
	panoramaFormatParams.push_back(info.panoramaFormatParam);
	//panoramaFormatParams.push_back(CV_IMWRITE_JPEG_LUMA_QUALITY);
	//panoramaFormatParams.push_back(100);
	//panoramaFormatParams.push_back(CV_IMWRITE_JPEG_CHROMA_QUALITY);
	//panoramaFormatParams.push_back(99);
	break;
      case JPEG2000:
	break;
      case TIFF:
	break;
      }

      if(info.range == true)
	{
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_Range."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, pImage.getRangeImage(), panoramaFormatParams);
	}

      if(info.normalizeRange == true)
	{
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_NormalizedRange."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, pImage.getNormalizedRangeImage(), panoramaFormatParams);
	}

      if(info.threeChannelRange == true)
	{

	  //cout<<"ThreeChannelRange"<<endl;
	  //cv::Mat iThreeChannel24BitRange = pImage.getThreeChannel24BitRangeImage();
	  //for(int h = 0; h < info.panoramaHeight; h++)
	  //{
	  //for(int w = 0; w < info.panoramaWidth; w++)
	  //{
	  //float range = 0.0;
	  //unsigned char byte;
	  //byte = iThreeChannel24BitRange.at<cv::Vec3b>(h,w)[0];
	  //range+= byte<<16;
	  //byte = iThreeChannel24BitRange.at<cv::Vec3b>(h,w)[1];
	  //range+= byte<<8;
	  //byte = iThreeChannel24BitRange.at<cv::Vec3b>(h,w)[2];
	  //range+= byte;
	  //range/= 10000;
	  //cout<<h<<" "<<w<<" "<<range<<endl;
	  //}
	  //cout<<endl;
	  //}
	  //cout<<"---------------------"<<endl;

	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_ThreeChannel24BitRange."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, pImage.getThreeChannel24BitRangeImage(), panoramaFormatParams);
	}

     if(info.oneGrayscaleRange == true)
	{
	  cv::Mat range1;
	  pImage.getOneGrayscaleRangeImage(range1);


	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_OneGrayscaleRange."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, range1, panoramaFormatParams);

	}
      if(info.threeGrayscaleRange == true)
	{
	  cv::Mat range1, range2, range3;
	  pImage.getThreeGrayscaleRangeImages(range1, range2, range3);

	  //cout<<"ThreeGtauScaleRange"<<endl;
	  //for(int h = 0; h < info.panoramaHeight; h++)
	  //{
	  //for(int w = 0; w < info.panoramaWidth; w++)
	  //{
	  //float range = 0.0;
	  //unsigned char byte;
	  //byte = range1.at<uchar>(h,w);
	  //range+= byte<<16;
	  //byte = range2.at<uchar>(h,w);
	  //range+= byte<<8;
	  //byte = range3.at<uchar>(h,w);
	  //range+= byte;
	  //range/= 10000;
	  //cout<<h<<" "<<w<<" "<<range<<endl;
	  //}
	  //cout<<endl;
	  //}
	  //cout<<"---------------------"<<endl;

	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_ThreeGrayscaleRange_1."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, range1, panoramaFormatParams);
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_ThreeGrayscaleRange_2."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, range2, panoramaFormatParams);
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_ThreeGrayscaleRange_3."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, range3, panoramaFormatParams);
	}

      if(info.reflectance == true)
	{
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_Reflectance."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, pImage.getReflectanceImage(), panoramaFormatParams);
	}

      if(info.color == true)
	{
	  out = info.outDir+"scan"+to_string(s, 3)+"_"+projectionMethodToString(info.projectionMethod)+"_"+to_string(info.panoramaWidth)+"x"+to_string(info.panoramaHeight)+"_Color."+panoramaFormatToFileFormatString(info.panoramaFormat);
	  imwrite(out, pImage.getColorImage(), panoramaFormatParams);
	}
    }
    return 0;
}


