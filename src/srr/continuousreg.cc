#include "continuousreg.h"
#include "slam6d/graphSlam6D.h"
#include "slam6d/Boctree.h"
#include "lboctree.h"
#include "lum6Deuler.h"
#include "ghelix6DQ2.h"
#include "slam6d/basicScan.h"
#include "srr/linescan.h"


/**
 * Creates a slam6D scan from the linescans begin-end, with the pose of index as reference coordinate system
 *
 */
Scan *joinLines(vector<LineScan*> &linescans, int begin, int end, int &index, double voxelsize, int nrpts)
{
  if (begin < 0) begin = 0;
  if (end >= (int)linescans.size()) end = linescans.size() -1;
  if (index < 0) index = 0;
  if (index >= (int)linescans.size()) index = linescans.size() -1;

  cout << "Create Scan between " << begin << " and " << end << "  with " << index << " as rep." << endl;
  vector<double *> pts;

  // transform points into coordinate system of index and store copies in pts
  double *transMat;
  double tinv[16];
  double rPos[3], rPosTheta[3];

  transMat = linescans[index]->transMat;
  Matrix4ToEuler(transMat, rPosTheta, rPos);
  M4inv(transMat, tinv);

  for (int i = begin; i <= end; i++) {
    transMat = linescans[i]->transMat;
    for (int j = 0; j < linescans[i]->nrpts; j++) {
      // TODO optimize
      double *p = new double[3];
      p[0] = linescans[i]->points[j][0];
      p[1] = linescans[i]->points[j][1];
      p[2] = linescans[i]->points[j][2];

      transform(p, transMat);
      transform(p, tinv);
      pts.push_back(p);
    }
  }


  if (voxelsize > 0.0) {
    vector<double*> tmp = pts;
    BOctTree<double> *oct = new BOctTree<double>(&(tmp[0]),tmp.size(), voxelsize);
    vector<double*> reduced;
    reduced.clear();
    if (nrpts > 0) {
      if (nrpts == 1) {
        oct->GetOctTreeRandom(reduced);
      }else {
        oct->GetOctTreeRandom(reduced, nrpts, false);
      }
    } else {
      oct->GetOctTreeCenter(reduced);
    }

    for (unsigned int i = 0; i < pts.size(); i++) {
      delete[] pts[i];
    }
    pts.clear();

    vector<double*> reduced2;
    for (unsigned int i = 0; i < reduced.size(); i++) {
      double *p = new double[3];
      p[0] = reduced[i][0];
      p[1] = reduced[i][1];
      p[2] = reduced[i][2];
      reduced2.push_back(p);
    }
    reduced.clear();

    delete oct;

    Scan *s = new BasicScan(rPos, rPosTheta, reduced2);
    for(size_t i=0;i<reduced2.size();i++) {
        delete[] reduced2[i];
    }
    reduced2.clear();
    return s;
  }

  Scan *s = new BasicScan(rPos,rPosTheta,pts);
  for(size_t i=0;i<pts.size();i++) {
    delete[] pts[i];
  }
  pts.clear();
  return s;
}


/**
 * Does ICP matching with a 3D scan created around linescan "index" (with "width" linescans in either direction)
 *
 * @param earliest is the index of the earliest linescan to match against
 * @param latest is the index of the latest linescan to match against
 *
 * Essentially we create 2 scans, that are rigidly matched using icp
 * the difference in the pose estimates of the "index" linescan is then distributed along the linescans in between
 *
 */
void preRegistration(vector<LineScan*> &linescans,
    int learliest, int llatest,
    int fearliest, int flatest,
    icp6D *icp, double voxelsize, int nr_octpts, int lindex, int findex) {

  // representative linescan for first scan
  if(findex < 0) {
    if(fearliest == 0) { // metascan case
      findex = flatest - (llatest - learliest)/2;
    } else {
      findex = fearliest + (flatest-fearliest)/2;
    }
  }
  if(lindex < 0) {
    lindex = learliest + (llatest-learliest)/2;
  }

  Scan *first = joinLines(linescans, fearliest, flatest, findex, voxelsize, nr_octpts);
  Scan *last  = joinLines(linescans, learliest, llatest, lindex, voxelsize, nr_octpts);

  vector<Scan*> scans;
  scans.push_back(first);
  scans.push_back(last);


#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif


    for (int i = 0; i < (int)scans.size(); i++) {
        scans[i]->setSearchTreeParameter(BOCTree);
        scans[i]->createSearchTree();
    }

  cout << "Doing icp..." << endl;
  icp->doICP(scans); // align scans to minimize odometry related errors

  cout << "DISTRIBUTE " << findex << " to " << lindex << endl;

  linearDistributeError(linescans, findex, lindex, last->get_transMat());

  // TODO delete Scans
  delete first;
  delete last;
}

void preRegistration(vector<LineScan*> &linescans,
    LScan* first, LScan* last,
    icp6D *icp, double voxelsize, int nr_octpts) {

  preRegistration(linescans,
      last->getBegin(),
      last->getEnd(),
      first->getBegin(),
      first->getEnd(),
      icp, voxelsize, nr_octpts,
      last->getRepresentative(),
      first->getRepresentative());
}

/**
 * Does slam6d matching with 3D scans
 *
 * @param earliest is the index of the earliest linescan to match against
 * @param latest is the index of the latest linescan to match against
 *
 * Essentially we create sequential scans, that are rigidly matched using slam6d
 * the differences in the pose estimates of the scans is then distributed along the linescans in between
 *
 */
void SemiRigidRegistration(vector<LineScan*> &linescans, graphSlam6DL *slam, Graph *&gr, int iterations, int slamiters, double mdm, int scaninterval, int scansize, double voxelsize, int nr_octpts)
{
  //////////////////////////
  double id[16];
  M4identity(id);
  for (unsigned int i = 0; i < linescans.size(); i++)  {
    linescans[i]->transform(id, Scan::ICP, 1 );
    linescans[i]->transform(id, Scan::ICP, 1 );
  }
  ///////////////////////
  vector<LScan*> scans;

  int nr_scans = linescans.size() / scaninterval + 1;


  for (int i = 0; i < nr_scans; i++) {
    int representative = i * scaninterval;
    int begin =  representative - scansize;
    int end = representative + scansize;
    cout << "create scan " << begin << " to " << end << " with " << representative << endl;
    LScan *s = new LScan(begin, end, representative);
    scans.push_back(s);
  }



  cout << "done with scans" << endl;


  // create slam graph
  if (gr == 0) {
    gr = new Graph();
    for (unsigned int i = 0; i < scans.size(); i++) {
      for (unsigned int j = 0; j < scans.size(); j++) {
        if (i==j) {
          continue;
        }
        if (i > j ) continue;
        gr->addLink(i, j);
      }
    }
  }
  cout << "done with graph" << endl;


  for (int iters = 0; iters < iterations; iters++) {
////////////////
    string fn = "trajectoryI" + to_string(iters,1) + ".txt";
    ofstream fout(fn.c_str());
  for (unsigned int i = 0; i < linescans.size(); i++)  {
    fout << i << " " << linescans[i]->rPos[2] << endl;
  }
  fout.close();
////////////////



  cout << "Iteration " << iters << " of semi rigid matching" << endl;
  cout << "finding ptpairs" << endl;

  LineScan::clearAllPtPairs();
  LineScan::computePtPairsOct(mdm*mdm);
  cout << "done with finding ptpairs" << endl;

  string fnc = "testcorr" + to_string(iters ,3) + ".3d";
  LScan::writeAllPtPairs(scans, 0, mdm*mdm, fnc.c_str());

    for (int i = 0; i < slamiters; i++) {
      cout << i << endl;
      double ret = slam->doGraphSlam6D(*gr, scans, 1, LineScan::allLineScans);
      if (ret < 0.001) {
        break;
      }
    }
  }

  cout << "deleting scans" << endl;

  for (unsigned int i = 0; i < scans.size(); i++)
    delete scans[i];
  scans.clear();

  delete gr;
}

double SSRR(vector<LineScan*> &linescans, graphSlam6D *slam, graphSlam6DL *slaml, Graph *&gr, icp6D *icp, int slamiters, int clpairs, double mdm, int scaninterval, int scansize, double voxelsize, int nr_octpts, volatile bool *abort_flag)
{
  std::cout << "create segment from linescans" << std::endl;

  if(abort_flag != nullptr && *abort_flag)
    return 0;

  std::vector<LSegment*> lsegments;
  LSegment* lsegment = new LSegment(0,linescans.size()-1, scaninterval, scansize);
  lsegments.push_back(lsegment);

  double ret = SSRR(linescans, lsegments, slam, slaml, gr, icp, slamiters,
               clpairs, mdm, voxelsize, nr_octpts, abort_flag);

  for(size_t i=0; i<lsegments.size();i++) {
    delete lsegments[i];
  }
  lsegments.clear();

  return ret;
}

double SSRR(vector<LineScan*> &linescans, vector<LSegment*> lsegments,
    graphSlam6D *slam, graphSlam6DL *slaml, Graph *&gr, icp6D *icp,
    int slamiters, int clpairs, double mdm, double voxelsize, int nr_octpts,
    volatile bool *abort_flag)
{

  vector<LScan*> lscans;
  vector<Scan*> sscans;
  vector<int> scan_reps;

  std::cout << "Segments: " << lsegments.size() << std::endl;
  for (size_t i = 0; i < lsegments.size() && (abort_flag == nullptr || *abort_flag == false); i++) {
    LSegment* lsegment = lsegments[i];
    vector<LScan*> &lsegment_lscans = lsegment->lscans;
    for(LScan* lscan : lsegment_lscans) {
      lscans.push_back(lscan);
      int representative = lscan->getRepresentative();
      Scan *ss = joinLines(linescans, lscan->getBegin(), lscan->getEnd(), representative, voxelsize, nr_octpts);
      sscans.push_back(ss);
      scan_reps.push_back(representative);
    }
  }

#ifdef _OPENMPxFIXME
#pragma omp parallel for schedule(dynamic)
#endif

  for (int i = 0; i < (int)lscans.size(); i++) {
#ifdef _OPENMPxFIXME
    if (abort_flag != nullptr && *abort_flag) continue;
#else
    if (abort_flag != nullptr && *abort_flag) break;
#endif
    sscans[i]->setSearchTreeParameter(BOCTree);
    sscans[i]->createSearchTree();
  }


  if (icp && (abort_flag == nullptr || *abort_flag == false))
    icp->doICP(sscans); // align scans to minimize odometry related errors


  if (abort_flag == nullptr || *abort_flag == false) {
    if (gr == 0) {
      gr = slam->computeGraph6Dautomatic(sscans, clpairs);
    }

    if (gr->getNrScans() <= 0) {
      throw runtime_error("graph contains zero scans");
    }

    cout << "Scans in Graph: " << gr->getNrScans() << " Linescans: " << (int)lscans.size() << endl;
    cout << "Links in Graph: " << gr->getNrLinks() << endl;
    cout << "do slam iters..." << endl;
  }
  for (int i = 0; i < slamiters && (abort_flag == nullptr || *abort_flag == false); i++) {
    cout << i << ": ";
    try {
	    double ret = slam->doGraphSlam6D(*gr, sscans, 1);
	    if (ret < 0.001) {
		    break;
	    }
    } catch (...) {
	    throw runtime_error("doGraphSlam6D failed");
    }
  }

  cout << "distributing errors..." << endl;

  double ret = 0;
  if(abort_flag == nullptr || *abort_flag == false) ret = slaml->doGraphSlam6D(*gr, lscans, 1, LineScan::allLineScans, &sscans, lsegments);

  cout << "deleting scans" << endl;

  // delete Scans
  for (unsigned int i = 0; i < sscans.size(); i++)
    delete sscans[i];
  sscans.clear();
  scan_reps.clear();
  lscans.clear();

  return ret/(LineScan::allLineScans.size());
}
