#include "continuousreg.h"
#include "linescan.h"
#include "lboctree.h"
#include <algorithm>
#include <scanio/helper.h>

#ifdef OPENMP
#include <omp.h>
#endif

using std::min;
using std::max;


vector <LineScan *> LineScan::allLineScans;
bool LineScan::outputFrames = false;
unsigned int LineScan::counter = 0;
unsigned int LineScan::mindist = 180;
unsigned int LineScan::maxdist = 5080;
std::string LineScan::error;

LineScan::~LineScan() {
  for (int i = 0; i < nrpts; i++) {
    delete[] points[i];
  }
  delete [] points;

}

LineScan::LineScan(Scan *scan) {
	index = counter++;
	identifier = std::string(scan->getIdentifier());
	basedir = scan->getPath();
	timestamp = scan->getTimeStamp();

	const double *scanRPos = scan->get_rPos();
	const double *scanRPosTheta = scan->get_rPosTheta();
	for (int i = 0; i < 3; i++) {
		rPos[i] = scanRPos[i];
		rPosTheta[i] = scanRPosTheta[i];
	}

	DataXYZ xyz(scan->get("xyz"));
	nrpts = xyz.size();
	if (nrpts != 0) {
		points = new double*[nrpts];

		for (int i = 0; i < nrpts; i++) {
			points[i] = new double[3];
			points[i][0] = xyz[i][0];
			points[i][1] = xyz[i][1];
			points[i][2] = xyz[i][2];
		}

		EulerToMatrix4(rPos, rPosTheta, orig_transMat);
		EulerToMatrix4(rPos, rPosTheta, transMat);
		allLineScans.push_back(this);
	}
	else {
		LineScan::error = "All points in Scan nr. " + identifier + " erased!";
	}
}

bool LineScan::applyCustomFilter(const Point& point, int &filterMode, int &nrOfParams, double* custFiltParams){
  bool filterTest = false;
  try {
      switch (filterMode)
      {
      case 0:
        // Custom Filter 0: symetrical, axis-parallel cuboid
        // all values inside the cuboid will be filtered (filterTest = false)
        // parameters: xFilterRange yFilterRange zFilterRange
        if (abs(point.x) > custFiltParams[0] || abs(point.y) > custFiltParams[1] || abs(point.z) > custFiltParams[2])
          filterTest = true;
        break;
      case 1:
        // Custom Filter 1: asymetrical axis-parallel cuboid
        // all values inside the cuboid will be filtered
        // parameters: xFilterRangeLow xFilterRangeHigh yFilterRangeLow yFilterRangeHigh zFilterRangeLow zFilterRangeHigh
        if (point.x < custFiltParams[0] || point.x > custFiltParams[1]
          || point.y < custFiltParams[2] || point.y > custFiltParams[3]
          || point.z < custFiltParams[4] || point.z > custFiltParams[5])
            filterTest = true;
          break;
      case 2:
          // As Custom Filter 1: asymetrical axis-parallel cuboid, with additional max range limitation
          // parameters: xFilterRangeLow xFilterRangeHigh yFilterRangeLow yFilterRangeHigh zFilterRangeLow zFilterRangeHigh maxRange
        if (point.x < custFiltParams[0] || point.x > custFiltParams[1]
          || point.y < custFiltParams[2] || point.y > custFiltParams[3]
          || point.z < custFiltParams[4] || point.z > custFiltParams[5]){
          if ((point.x * point.x + point.y * point.y + point.z * point.z) < (custFiltParams[6] * custFiltParams[6])){
            filterTest = true;
          }
          break;
        }
      case 10:
          // Custom Filter 10: symetrical, axis-parallel cuboid
          // all values outside the cuboid will be filtered (filterTest = false)
          // parameters: xFilterRange yFilterRange zFilterRange
          if (abs(point.x) < custFiltParams[0] && abs(point.y) < custFiltParams[1] && abs(point.z) < custFiltParams[2])
              filterTest = true;
          break;
      case 11:
          // Custom Filter 11: asymetrical axis-parallel cuboid
          // all values outside the cuboid will be filtered
          // parameters: xFilterRangeLow xFilterRangeHigh yFilterRangeLow yFilterRangeHigh zFilterRangeLow zFilterRangeHigh
          if (point.x > custFiltParams[0] && point.x < custFiltParams[1]
              && point.y > custFiltParams[2] && point.y < custFiltParams[3]
              && point.z > custFiltParams[4] && point.z < custFiltParams[5]) {
              filterTest = true;
          }
          break;
      case 20:
          // Custom Filter 20: two asymetrical axis-parallel cuboids, one inside the other
          // all values outside the first (outer) cuboid, and inside the second (inner) cuboid will be filtered
          // since the cuboids are supposed to overlap, only data in between the cuboids will remain
          // parameters: xFilterRangeLow1 xFilterRangeHigh1 yFilterRangeLow1 yFilterRangeHigh1 zFilterRangeLow1 zFilterRangeHigh1 xFilterRangeLow2 xFilterRangeHigh2 yFilterRangeLow2 yFilterRangeHigh2 zFilterRangeLow2 zFilterRangeHigh2

          // inside first cuboid..?
          if (point.x > custFiltParams[0] && point.x < custFiltParams[1]
              && point.y > custFiltParams[2] && point.y < custFiltParams[3]
              && point.z > custFiltParams[4] && point.z < custFiltParams[5]) {
              // inside first

              // ... but outside second..?
              if (point.x < custFiltParams[6] || point.x > custFiltParams[7]
                  || point.y < custFiltParams[8] || point.y > custFiltParams[9]
                  || point.z < custFiltParams[10] || point.z > custFiltParams[11]) {
                  filterTest = true;
              }
          }
          break;
      default:
          filterTest = true;
          break;
      }
  }
  catch (...){
      // Error occured - deactivate filter
      filterTest = true;
  }
  return filterTest;
}


bool LineScan::writeFrames(std::string basedir, std::string suffix, IOType type, bool append) {
    std::map<std::string, std::string> framesdata;
    for (unsigned int i = 0; i < LineScan::allLineScans.size(); i++) {
      framesdata[LineScan::allLineScans[i]->getFramesPath()] = LineScan::allLineScans[i]->getFramesData();
    }

    std::ios_base::openmode writemode = append ? std::ios_base::app : std::ios_base::out;

    bool state = write_multiple(framesdata,writemode);

    for (unsigned int i = 0; i < LineScan::allLineScans.size(); i++) {
      LineScan::allLineScans[i]->resetFramesData();
    }

    return state;
}


void LineScan::resetTransform() {
  Matrix4ToEuler(orig_transMat, rPosTheta, rPos);
  EulerToMatrix4(rPos, rPosTheta, transMat);
  EulerToMatrix4(rPos, rPosTheta, orig_transMat);

}

void LineScan::transformToEuler(double rP[3], double rPT[3], const Scan::AlgoType type, int islum)
{
  double tinv[16];
  double alignxf[16];
  M4inv(transMat, tinv);
  transform(tinv, Scan::INVALID, -1);
  EulerToMatrix4(rP, rPT, alignxf);
  transform(alignxf, type, islum);
}


void LineScan::transform(const double alignxf[16], const Scan::AlgoType type, int islum) {
  // update transformation
  double tempxf[16];
  MMult(alignxf, transMat, tempxf);
  memcpy(transMat, tempxf, sizeof(transMat));
  Matrix4ToEuler(transMat, rPosTheta, rPos);


  // store transformation
  int  found = 0;
  int end_loop = 0;

    switch (islum) {
      case -1:
        // write no tranformation
        break;
      case 0:
        end_loop = (int)allLineScans.size();
        for (int iter = 0; iter < end_loop; iter++) {
          if (allLineScans[iter]->sout.good()) {
            allLineScans[iter]->sout << allLineScans[iter]->transMat;
            if (allLineScans[iter] == this ) {
              found = iter;
              allLineScans[iter]->sout << type << endl;
            } else {
              if (found == 0) {
                allLineScans[iter]->sout << Scan::ICPINACTIVE << endl;
              } else {
                allLineScans[iter]->sout << Scan::INVALID << endl;
              }
            }
          } else {
            cerr << "ERROR: Cannot store frames." << endl;
            exit(1);
          }
        }
        break;
      case 1:
        if (sout.good()) {
          sout << transMat << type << endl;
        } else {
          cerr << "ERROR: Cannot store frames." << endl;
          exit(1);
        }
        break;
      case 2:
        end_loop = (int)allLineScans.size();
        for (int iter = 0; iter < end_loop; iter++) {
          if (allLineScans[iter] == this) {
            found = iter;
            if (sout.good()) {
              sout << transMat << type << endl;
            } else {
              cerr << "ERROR: Cannot store frames." << endl;
              exit(1);
            }

            if (allLineScans[0]->sout.good()) {
              allLineScans[0]->sout << allLineScans[0]->transMat << type << endl;
            } else {
              cerr << "ERROR: Cannot store frames." << endl;
              exit(1);
            }
            continue;
          }
          if (found != 0) {
            allLineScans[iter]->sout << allLineScans[iter]->transMat << Scan::INVALID << endl;
          }
        }
        break;
      case 3:
        break;
      default:
        cerr << "invalid point transformation mode" << endl;
    }

}


double *LineScan::FindClosest(double *p, double maxdist2, int threadNum) {
  //p is the query point to be found

  // TODO implement the correct search for corresponding points
  return p;
}


LScan::LScan(int _begin, int _end, int rep) {
  representative = min( max(0,rep), (int)LineScan::allLineScans.size()-1) ;
  begin = min( max(0,_begin), (int)LineScan::allLineScans.size()-1) ;
  end = min( max(0,_end), (int)LineScan::allLineScans.size()-1) ;
  scans.insert(scans.begin(), LineScan::allLineScans.begin() + begin, LineScan::allLineScans.begin() + end);
  memcpy(transMat, LineScan::allLineScans[representative]->transMat, 16*sizeof(double));
}

LScan::~LScan() {

}

void LScan::transform(const double alignxf[16], const Scan::AlgoType type, int islum) {
  cerr << "ERROR in LScan::transform()  Function not implemented!!!" << endl;

  double tempxf[16];
  MMult(alignxf, transMat, tempxf);
  memcpy(transMat, tempxf, sizeof(transMat));
  //Matrix4ToEuler(transMat, rPosTheta, rPos);


}
void LScan::writeAllPtPairs(vector<LScan *> &scans, int rnd, double max_dist_match2, string filename) {
  double centroid_m[3];
  double centroid_d[3];
  vector<PtPair> pairs;
  ofstream fout(filename.c_str());

  int k = 0; {
    for (unsigned int l = 0; l < scans.size(); l++) {
      getPtPairs(&pairs, scans[k], scans[l], 0, rnd, max_dist_match2, centroid_m, centroid_d);
      for (unsigned int j = 0; j < pairs.size(); j++) {
        fout <<
          pairs[j].p1.x << " " << pairs[j].p1.y << " " << pairs[j].p1.z
          << " " <<
          pairs[j].p2.x << " " << pairs[j].p2.y << " " << pairs[j].p2.z
          << endl;
      }
    }
  }

  fout.close();
  fout.clear();

}
void LScan::writePtPairs(LScan *LSource, LScan* LTarget, int rnd, double max_dist_match2, string filename) {
  double centroid_m[3];
  double centroid_d[3];
  vector<PtPair> pairs;
  ofstream fout(filename.c_str());
  getPtPairs(&pairs, LSource, LTarget, 0, rnd, max_dist_match2, centroid_m, centroid_d);
  for (unsigned int j = 0; j < pairs.size(); j++) {
    fout <<
      pairs[j].p1.x << " " << pairs[j].p1.y << " " << pairs[j].p1.z
      << " " <<
      pairs[j].p2.x << " " << pairs[j].p2.y << " " << pairs[j].p2.z
      << endl;
  }

  fout.close();
  fout.clear();


}

//#define GROUND_TRUTH_PAIRS

void LScan::getPtPairs(vector <PtPair> *pairs, LScan *Source, LScan* Target,  int thread_num,
              int rnd, double max_dist_match2, double *centroid_m, double *centroid_d) {
#ifdef GROUND_TRUTH_PAIRS
  centroid_m[0] = 0.0;
  centroid_m[1] = 0.0;
  centroid_m[2] = 0.0;
  centroid_d[0] = 0.0;
  centroid_d[1] = 0.0;
  centroid_d[2] = 0.0;
  int np = 0;

  int beg1,end1, beg2,end2;
  beg1 = Source->getBegin();
  end1 = Source->getEnd();
  beg2 = Target->getBegin();
  end2 = Target->getEnd();

  int s,e;
  e = -1;

  if (beg1 <= beg2) {
    s = beg2;
    if ( beg2 <= end1) {
      e = end1;
    } else if ( end2 <= end1) {
      e = end2;
    }
  } else {
    s = beg1;
    if ( beg1 <= end2) {
      e = end2;
    } else if ( end1 <= end2) {
      e = end1;
    }
  }

  for (int j = beg2; j <= end2; j++) {  // use beg2 and end2 if no links are i > 0
    LineScan *scan = LineScan::allLineScans[j];
    double **points = scan->getPoints();
    for (unsigned int i = 0; i < scan->getNrPoints(); i++) {
      double s[3];
      s[0] = points[i][0];
      s[1] = points[i][1];
      s[2] = points[i][2];
      double t[3];
      t[0] = points[i][0];
      t[1] = points[i][1];
      t[2] = points[i][2];
      ::transform(t, scan->transMat);
      ::transform(s, scan->ground_truth_mat);

      PtPair myPair(s, t);
      pairs->push_back(myPair);

      centroid_d[0] += t[0];
      centroid_d[1] += t[1];
      centroid_d[2] += t[2];
      centroid_m[0] += s[0];
      centroid_m[1] += s[1];
      centroid_m[2] += s[2];
      np++;
    }
  }
  centroid_d[0] /= np;
  centroid_d[1] /= np;
  centroid_d[2] /= np;
  centroid_m[0] /= np;
  centroid_m[1] /= np;
  centroid_m[2] /= np;


#else
  centroid_m[0] = 0.0;
  centroid_m[1] = 0.0;
  centroid_m[2] = 0.0;
  centroid_d[0] = 0.0;
  centroid_d[1] = 0.0;
  centroid_d[2] = 0.0;


  for (vector<LineScan*>::iterator it = Target->scans.begin(); it != Target->scans.end(); it++) {
    for (vector<LineScan*>::iterator its = LineScan::allLineScans.begin(); its != LineScan::allLineScans.end(); its++) {
      LineScan::getPtPairs(pairs, *its, *it, thread_num, rnd, max_dist_match2, centroid_m, centroid_d);
    }
  }
  ///////////////////

  centroid_m[0] /= pairs->size();
  centroid_m[1] /= pairs->size();
  centroid_m[2] /= pairs->size();
  centroid_d[0] /= pairs->size();
  centroid_d[1] /= pairs->size();
  centroid_d[2] /= pairs->size();
  return;
  ////////////////

  // this part is for adding a "stay where you are" weight
  int beg1,end1, beg2,end2;
  beg1 = Source->getBegin();
  end1 = Source->getEnd();
  beg2 = Target->getBegin();
  end2 = Target->getEnd();

  int s,e;
  e = -1;

  if (!(end1 < beg2 || end2 < beg1)) { // no overlap
    vector<int> ss;
    ss.push_back(beg1);
    ss.push_back(beg2);
    ss.push_back(end1);
    ss.push_back(end2);
    sort(ss.begin(), ss.end());
    s = ss[1];
    e = ss[2];

    unsigned int counter = 0;
    unsigned int max_counter = 1000;  // number of point pairs that "anchor" the scan
    for (int j = s; j < e; j++) {
      LineScan *scan = LineScan::allLineScans[j];
      double **points = scan->getPoints();
      for (int i = 0; i < scan->getNrPoints(); i++) {
        double s[3];
        s[0] = points[i][0];
        s[1] = points[i][1];
        s[2] = points[i][2];
        ::transform(s, scan->transMat);

        PtPair myPair(s, s);
        pairs->push_back(myPair);

        centroid_d[0] += s[0];
        centroid_d[1] += s[1];
        centroid_d[2] += s[2];
        centroid_m[0] += s[0];
        centroid_m[1] += s[1];
        centroid_m[2] += s[2];

        counter++;
        if (counter > max_counter) break;
      }
      if (counter > max_counter) break;
    }

  }
  centroid_m[0] /= pairs->size();
  centroid_m[1] /= pairs->size();
  centroid_m[2] /= pairs->size();
  centroid_d[0] /= pairs->size();
  centroid_d[1] /= pairs->size();
  centroid_d[2] /= pairs->size();
#endif

  return;
}

void LScan::getOdomPairs(vector <PtPair> *pairs, LScan *Source, LScan* Target,
              double *centroid_m, double *centroid_d, bool ground_truth) {
  centroid_m[0] = 0.0;
  centroid_m[1] = 0.0;
  centroid_m[2] = 0.0;
  centroid_d[0] = 0.0;
  centroid_d[1] = 0.0;
  centroid_d[2] = 0.0;

  LineScan *target = LineScan::allLineScans[Target->getRepresentative()];
  LineScan *source = LineScan::allLineScans[Source->getRepresentative()];


  double *orig_mat_prev;
  double *orig_mat;

  if (ground_truth) {
    orig_mat_prev = source->ground_truth_mat;
    orig_mat = target->ground_truth_mat;
  } else {
    orig_mat_prev = source->orig_transMat;
    orig_mat = target->orig_transMat;
  }

  double *mat_prev = source->transMat;

  double ompi[16];
  double odo[16];
  double glob2odo[16];


  M4inv(orig_mat_prev, ompi);
  MMult( orig_mat, ompi, odo);


if (ground_truth) {
  cout << "OMPI " << ompi << endl;
  cout << "MP   " << mat_prev << endl;
  cout << "G2OD " << glob2odo << endl;
}

  unsigned int max_counter = 9;  // number of point pairs that "anchor" the scan


  for (unsigned int i = 0; i < 100; i++) {
    double s[3], t[3];
    t[0] = s[0] = (rand() / (double)RAND_MAX) * 100.0;
    t[1] = s[1] = (rand() / (double)RAND_MAX) * 100.0;
    t[2] = s[2] = (rand() / (double)RAND_MAX) * 100.0;

    ::transform(s, source->transMat);
    ::transform(s, odo);
    ::transform(t, target->transMat);
    PtPair myPair(s, t);
    pairs->push_back(myPair);
    centroid_m[0] += s[0];
    centroid_m[1] += s[1];
    centroid_m[2] += s[2];
    centroid_d[0] += t[0];
    centroid_d[1] += t[1];
    centroid_d[2] += t[2];
  }

  for (unsigned int i = 0; i < max_counter/3; i++) {
    double s[3], t[3];
    t[0] = s[0] = 1.0;
    t[1] = s[1] = 0.0;
    t[2] = s[2] = 0.0;

    ::transform(s, source->transMat);
    ::transform(s, odo);
    ::transform(t, target->transMat);
    PtPair myPair(s, t);
    pairs->push_back(myPair);
    centroid_m[0] += s[0];
    centroid_m[1] += s[1];
    centroid_m[2] += s[2];
    centroid_d[0] += t[0];
    centroid_d[1] += t[1];
    centroid_d[2] += t[2];
  }
  for (unsigned int i = 0; i < max_counter/3; i++) {
    double s[3], t[3];
    t[0] = s[0] = 0.0;
    t[1] = s[1] = 1.0;
    t[2] = s[2] = 0.0;

    ::transform(s, source->transMat);
    ::transform(s, odo);
    ::transform(t, target->transMat);
    PtPair myPair(s, t);
    pairs->push_back(myPair);
    centroid_m[0] += s[0];
    centroid_m[1] += s[1];
    centroid_m[2] += s[2];
    centroid_d[0] += t[0];
    centroid_d[1] += t[1];
    centroid_d[2] += t[2];
  }
  for (unsigned int i = 0; i < max_counter/3; i++) {
    double s[3], t[3];
    t[0] = s[0] = 0.0;
    t[1] = s[1] = 0.0;
    t[2] = s[2] = 1.0;

    ::transform(s, source->transMat);
    ::transform(s, odo);
    ::transform(t, target->transMat);
    PtPair myPair(s, t);
    pairs->push_back(myPair);
    centroid_m[0] += s[0];
    centroid_m[1] += s[1];
    centroid_m[2] += s[2];
    centroid_d[0] += t[0];
    centroid_d[1] += t[1];
    centroid_d[2] += t[2];
  }

  centroid_m[0] /= pairs->size();
  centroid_m[1] /= pairs->size();
  centroid_m[2] /= pairs->size();
  centroid_d[0] /= pairs->size();
  centroid_d[1] /= pairs->size();
  centroid_d[2] /= pairs->size();

  return;
}


void LScan::transformToEuler(double rP[3], double rPT[3], const Scan::AlgoType type, int islum ) {
  // this is called by lum6dEuler
  EulerToMatrix4(rP, rPT, transMat);
}


void LScan::resetTransform() {
  double rPosTheta[3], rPos[3];
  Matrix4ToEuler(repscan->orig_transMat, rPosTheta, rPos);
  EulerToMatrix4(rPos, rPosTheta, transMat);

  for (unsigned int i = 0; i < scans.size(); i++) {
    scans[i]->resetTransform();
  }

}


void LineScan::getPtPairs(vector <PtPair> *pairs, LineScan *Source, LineScan* Target, int thread_num, int rnd, double max_dist_match2, double *centroid_m, double *centroid_d) {

  // check if anything should be done
  if ( ((int)Source->index - (int)maxdist <= (int)Target->index && (int)Target->index <= (int)Source->index - (int)mindist) ||
       ((int)Source->index + (int)mindist <= (int)Target->index && (int)Target->index <= (int)Source->index + (int)maxdist) ) {

    Mcorr::iterator corrs = Target->correspondences.find(Source->index);
    if (corrs == Target->correspondences.end() ) return;

    vlppair &corpts = corrs->second;

    // for each corresponding point
    for (unsigned int i = 0; i < corpts.size(); i++) {
      double *s = corpts[i].first;
      double *t = corpts[i].second;

      PtPair myPair(s, t);
      pairs->push_back(myPair);

      centroid_d[0] += t[0];
      centroid_d[1] += t[1];
      centroid_d[2] += t[2];
      centroid_m[0] += s[0];
      centroid_m[1] += s[1];
      centroid_m[2] += s[2];

    }


  }
}

void LineScan::findPtPairs(int i, double max_dist_match2) {
  if (i < 0 || i >= (int)allLineScans.size())
    return;

  vector<lppair> vpairs;

  LineScan *target = allLineScans[i];
  for (int ii = 0; ii < target->nrpts; ii++) {
    double *p = target->points[ii];

    double *closest = 0;

    double mdm2 = max_dist_match2;
    for (int j = 0; j < this->nrpts; j++) {
      double *q = this->points[j];
      double dist = Dist2(p, q);

      if (dist < mdm2) {
        closest = q;
        mdm2 = dist;
      }
    }
    if (closest) {
      vpairs.push_back(lppair(closest, p) );
    }
  }
  if (!vpairs.empty()) {
    correspondences[i] = vpairs;
  }

}

void LineScan::findPtPairs(double max_dist_match2) {
  for (int i = (int)(this->index) - (int)maxdist; i <= (int)(this->index) - (int)(mindist); i++) {
    findPtPairs(i, max_dist_match2);
//    cout << "  _ " << i << endl;
  }
  for (int i = (int)this->index + (int)mindist; i <= (int)this->index + (int)maxdist; i++) {
    findPtPairs(i, max_dist_match2);
//    cout << "   " << i << endl;
  }
}

void LineScan::computePtPairs(double max_dist_match2) {
// visual studio 2013 only supports OpenMP 2.5 which restricts the index
// variable in an OpenMP for statement to have signed integral type
// gcc >= 4.4 supports OpenMP 3.0 which allows unsigned int
#if defined(_WIN32) && defined(_OPENMP)
    omp_set_num_threads(OPENMP_NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < allLineScans.size(); i++) {
        allLineScans[i]->findPtPairs(max_dist_match2);
        cout << i << endl;

    }
#else
#ifdef _OPENMP
omp_set_num_threads(OPENMP_NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
#endif
  for (unsigned int i = 0; i < allLineScans.size(); i++) {
    allLineScans[i]->findPtPairs(max_dist_match2);
    cout << i << endl;

  }
#endif
}
void LineScan::computePtPairsOct(double max_dist_match2) {
  LBOctTree<double> *tree = new LBOctTree<double>(allLineScans, 10.0);
  cout << "tree created" << endl;
  int count = 0;

  vector<double*> query_pts;
  cout << "get random" << endl;
#ifdef _WIN32
  srand(0);
#else
  srandom(0);
#endif
  tree->GetOctTreeRandom(query_pts);
  for (unsigned int i = 0; i < allLineScans.size(); i++) {
    allLineScans[i]->clearQuery();
  }
  cout << "add queries" << endl;
  for (unsigned int i = 0; i < query_pts.size(); i++) {
    int index = query_pts[i][3];
    allLineScans[index]->addQuery(query_pts[i]);

  }

  cout << "find closest points" << endl;
  for (unsigned int i = 0; i < allLineScans.size(); i++) {
    int r = tree->FindClosestPoints(allLineScans[i], max_dist_match2, LineScan::mindist, LineScan::maxdist, 0);
    count += r;
  }

  cout << "deleting  tree" << endl;
  delete tree;
}

void LineScan::addPtPair(lppair &pair, unsigned int index) {
  correspondences[index].push_back(pair);
}

void LineScan::printPtPairs() {
  for (unsigned int l = 0; l < allLineScans.size(); l++) {
    Mcorr &correspondences = allLineScans[l]->correspondences;
    for (Mcorr::iterator it = correspondences.begin(); it!=correspondences.end(); it++) {
      for (unsigned int i = 0; i < it->second.size(); i++) {
        double *p = it->second[i].first;
        double *q = it->second[i].second;
        cerr << p[0] << " " << p[1] << " " << p[2] << " " << q[0] << " " << q[1] << " " << q[2]
          << endl;
      }

    }
  }
}

void LineScan::clearPtPairs() {

  for (Mcorr::iterator it = correspondences.begin(); it!=correspondences.end(); it++) {
    for (unsigned int i = 0; i < it->second.size(); i++) {
      double *q = it->second[i].second;
      double *p = it->second[i].first;
      delete[] q;
      delete[] p;
    }
  }

  correspondences.clear();
}
