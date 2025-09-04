/**
 * @file
 * @brief The implementation of globally consistent scan matching algorithm by using APX correction
 * @author Andreas Nuechter. Institute of Computer Science, University of Osnabrueck, Germany.
 * @author Jan Elseberg. Institute of Computer Science, University of Osnabrueck, Germany.
 */

#ifdef _MSC_VER
#ifdef OPENMP
#define _OPENMP
#endif
#endif

#include "linescan.h"
#include "gapx6D.h"
#include "slam6d/icp6Dapx.h"
#include "suitesparse/cs.h"

#include <cfloat>
#include <fstream>
using std::flush;
#include <cstring>

#include "slam6d/globals.icc"

using std::ofstream;
using namespace NEWMAT;
/**
 * Constructor
 *
 * @param my_icp6Dminimizer Pointer to ICP minimization functor
 * @param mdm Maximum PtoP distance to which point pairs are collected for ICP
 * @param max_dist_match Maximum PtoP distance to which point pairs are collected for LUM
 * @param max_num_iterations Maximal number of iterations for ICP
 * @param quiet Suspesses all output to std out
 * @param meta Indicates if metascan matching has to be used
 * @param rnd Indicates if randomization has to be used
 * @param eP Extrapolate odometry?
 * @param anim Animate which frames?
 * @param epsilonICP Termination criterion for ICP
 * @param use_cache Shall we used cached k-d tree search
 * @param epsilonLUM Termination criterion for LUM
 */
gapx6DL::gapx6DL(icp6Dminimizer *my_icp6Dminimizer,
		   double mdm, double max_dist_match,
		   int max_num_iterations, bool quiet, bool meta, int rnd,
		   bool eP, int anim, double epsilonICP, bool use_cache, double epsilonLUM)
  : graphSlam6DL(my_icp6Dminimizer,
			 mdm, max_dist_match,
			 max_num_iterations, quiet, meta, rnd,
			 eP, anim, epsilonICP, use_cache, epsilonLUM)
{ }


/**
 * Destructor
 */
gapx6DL::~gapx6DL()
{
  delete my_icp;
}



double gapx6DL::genBAtransForLinkedPair( int firstScanNum, int secondScanNum,
							    double *centroids_m, double *centroids_d,
							    SymmetricMatrix *B, ColumnVector *A, ColumnVector &X)
{

  Point cm(centroids_m);
  Point cd(centroids_d);

  double x[3], dx[3] = {0.0, 0.0, 0.0};
  unsigned int vectorOffset;
  double alignxf[16];

  vectorOffset = (firstScanNum-1) * 3;

  if (firstScanNum != 0) {
    x[0] =  X(vectorOffset + 1);
    x[1] =  X(vectorOffset + 2);
    x[2] =  X(vectorOffset + 3);
  } else {
    x[0] = x[1] = x[2] = 0.0;
  }

  icp6D_APX::computeRt(x, dx, alignxf);

  cm.transform(alignxf);

  vectorOffset = (secondScanNum-1) * 3;

  x[0] =  X(vectorOffset + 1);
  x[1] =  X(vectorOffset + 2);
  x[2] =  X(vectorOffset + 3);

  icp6D_APX::computeRt(x, dx, alignxf);

  cd.transform(alignxf);

  ColumnVector Ak1(3);

  Ak1(1) = cm.x - cd.x;
  Ak1(2) = cm.y - cd.y;
  Ak1(3) = cm.z - cd.z;

#ifdef _OPENMP
#pragma omp critical (enterB)
#endif
  {
    if(firstScanNum != 0) {
	 A->Rows((firstScanNum-1)*3+1, (firstScanNum-1)*3+3) -= Ak1;
	 B->element(firstScanNum -1, firstScanNum  - 1)      += 1;
	 B->element(firstScanNum -1, secondScanNum - 1)      -= 1;
    }
    A->Rows((secondScanNum-1)*3+1, (secondScanNum-1)*3+3) += Ak1;
    B->element(secondScanNum - 1, secondScanNum - 1)      += 1;
  }    // of pragma omp critical

  return 1.0;

}


/**
 * This function generates the matrices B and Bd that are used for solving B * c = Bd.
 * This function has to be called once for every linked scan-pair.
 *
 * @param firstScanNum The number of the first scan of the linked scan-pair
 * @param secondScanNum The number of the second scan of the linked scan-pair
 * @param ptpairs Vector that holds all point-pairs for the actual scan-pair
 * @param B Matrix with dimension (6*(number of scans-1)) x (6 * (number of scans-1))
 * @param Bd Vector with dimension (6*(number of scans-1))
 * @return returns the sum of square distance
 */
double gapx6DL::genBArotForLinkedPair( int firstScanNum, int secondScanNum, vPtPair *ptpairs,
							  double *centroids_m, double *centroids_d,
                                     Matrix *B, ColumnVector *A)
{
  Matrix Mk(3,3), Dk(3,3);
  Matrix MkMkt(3,3), DkDkt(3,3), DkMkt(3,3), MkDkt(3,3);
  ColumnVector Ak1(3), Ak2(3);

  double p1x, p1y, p1z, p2x, p2y, p2z;
  double p1x2, p1y2, p1z2, p2x2, p2y2, p2z2;

  double p1xp1y, p1xp1z, p1yp1z, p1y2p1z2, p1x2p1z2, p1x2p1y2;
  double p2xp2y, p2xp2z, p2yp2z, p2y2p2z2, p2x2p2z2, p2x2p2y2;

  double p1xp2y, p1xp2z, p1yp2z;
  double p2xp1y, p2xp1z, p2yp1z;
  double p2zp1x, p2zp1y, p2yp1x, p1yp2x, p1zp2x, p1zp2y;

  double p1xp2xp1yp2y, p1xp2xp1zp2z, p1yp2yp1zp2z;

  int n = (*ptpairs).size();

  MkMkt = 0;
  DkDkt = 0;
  DkMkt = 0;
  MkDkt = 0;
  Ak1 = 0;
  Ak2 = 0;

  for (int i = 0; i < n; i++) {
    p1x = (*ptpairs)[i].p1.x - centroids_m[0]; // << korrekt
    p1y = (*ptpairs)[i].p1.y - centroids_m[1]; // << korrekt
    p1z = (*ptpairs)[i].p1.z - centroids_m[2]; // << korrekt

    p2x = (*ptpairs)[i].p2.x - centroids_m[0]; // << korrekt
    p2y = (*ptpairs)[i].p2.y - centroids_m[1]; // << korrekt
    p2z = (*ptpairs)[i].p2.z - centroids_m[2]; // << korrekt

    p1x2 = sqr(p1x);
    p1y2 = sqr(p1y);
    p1z2 = sqr(p1z);
    p2x2 = sqr(p2x);
    p2y2 = sqr(p2y);
    p2z2 = sqr(p2z);

    p1x2p1z2 = p1x2 + p1z2;
    p1x2p1y2 = p1x2 + p1y2;
    p1y2p1z2 = p1y2 + p1z2;

    p2x2p2z2 = p2x2 + p2z2;
    p2x2p2y2 = p2x2 + p2y2;
    p2y2p2z2 = p2y2 + p2z2;

    p1xp2xp1yp2y = p1x * p2x + p1y + p2y;
    p1xp2xp1zp2z = p1x * p2x + p1z + p2z;
    p1yp2yp1zp2z = p1y * p2y + p1z + p2z;

    p1xp1y = p1x * p1y;
    p1xp1z = p1x * p1z;
    p1yp1z = p1y * p1z;

    p2xp1y = p2x * p1y;
    p2xp1z = p2x * p1z;
    p2yp1x = p2y * p1x;
    p2yp1z = p2y * p1z;
    p2zp1x = p2z * p1x;
    p2zp1y = p2z * p1y;

    p1xp2y = p1x * p2y;
    p1xp2z = p1x * p2z;
    p1yp2z = p1y * p2z;
    p1yp2x = p1y * p2x;
    p1zp2x = p1z * p2x;
    p1zp2y = p1z * p2y;

    p2xp2y = p2x * p2y;
    p2xp2z = p2x * p2z;
    p2yp2z = p2y * p2z;

    MkMkt(1,1) += p1y2p1z2;
    MkMkt(1,2) -= p1xp1y;
    MkMkt(1,3) -= p1xp1z;
    MkMkt(2,1) -= p1xp1y;
    MkMkt(2,2) += p1x2p1z2;
    MkMkt(2,3) -= p1yp1z;
    MkMkt(3,1) -= p1xp1z;
    MkMkt(3,2) -= p1yp1z;
    MkMkt(3,3) += p1x2p1y2;

    DkDkt(1,1) += p2y2p2z2;
    DkDkt(1,2) -= p2xp2y;
    DkDkt(1,3) -= p2xp2z;
    DkDkt(2,1) -= p2xp2y;
    DkDkt(2,2) += p2x2p2z2;
    DkDkt(2,3) -= p2yp2z;
    DkDkt(3,1) -= p2xp2z;
    DkDkt(3,2) -= p2yp2z;
    DkDkt(3,3) += p2x2p2y2;

    MkDkt(1,1) += p1yp2yp1zp2z;
    MkDkt(1,2) -= p1yp2x;
    MkDkt(1,3) -= p1zp2x;
    MkDkt(2,1) -= p1yp2x;
    MkDkt(2,2) += p1xp2xp1zp2z;
    MkDkt(2,3) -= p1zp2y;
    MkDkt(3,1) -= p1zp2x;
    MkDkt(3,2) -= p1zp2y;
    MkDkt(3,3) += p1xp2xp1yp2y;

    DkMkt(1,1) += p1yp2yp1zp2z;
    DkMkt(1,2) -= p2yp1x;
    DkMkt(1,3) -= p2zp1x;
    DkMkt(2,1) -= p2yp1x;
    DkMkt(2,2) += p1xp2xp1zp2z;
    DkMkt(2,3) -= p2zp1y;
    DkMkt(3,1) -= p2zp1x;
    DkMkt(3,2) -= p2zp1y;
    DkMkt(3,3) += p1xp2xp1yp2y;

    Ak1(1) -= (p1z - p2z) * p2y - (p1y - p2y) * p2z;
    Ak1(2) -= (p1x - p2x) * p2z - (p1z - p2z) * p2x;
    Ak1(3) -= (p1y - p2y) * p2x - (p1x - p2x) * p2y;

    Ak2(1) += (p1z - p2z) * p1y - (p1y - p2y) * p1z;
    Ak2(2) += (p1x - p2x) * p1z - (p1z - p2z) * p1x;
    Ak2(3) += (p1y - p2y) * p1x - (p1x - p2x) * p1y;
  }

#ifdef _OPENMP
#pragma omp critical (enterB)
#endif
    {
	 if(firstScanNum != 0) {
	   A->Rows((firstScanNum-1)*3+1, (firstScanNum-1)*3+3) += Ak1;
 	   B->SubMatrix((firstScanNum-1)*3+1,(firstScanNum-1)*3+3, (firstScanNum-1)*3+1,(firstScanNum-1)*3+3) += MkMkt;

 	   B->SubMatrix((firstScanNum-1)*3+1,(firstScanNum-1)*3+3, (secondScanNum-1)*3+1,(secondScanNum-1)*3+3) +=  DkMkt;

 	   B->SubMatrix((secondScanNum-1)*3+1,(secondScanNum-1)*3+3, (firstScanNum-1)*3+1,(firstScanNum-1)*3+3) +=  MkDkt;

	 }
	 A->Rows((secondScanNum-1)*3+1, (secondScanNum-1)*3+3) += Ak2;
	 B->SubMatrix((secondScanNum-1)*3+1,(secondScanNum-1)*3+3, (secondScanNum-1)*3+1,(secondScanNum-1)*3+3) += DkDkt;

    }    // of pragma omp critical

  return 1.0;
}



/**
 * This function is used to match a set of laser scans with any minimally
 * connected Graph, using the globally consistent LUM-algorithm in 3D.
 *
 * @param gr Some Graph with no real subgraphs except for itself
 * @param allScans Contains all laser scans
 * @param nrIt The number of iterations the LUM-algorithm will run
 * @return Euclidian distance of all pose shifts
 */
double gapx6DL::doGraphSlam6D(Graph gr, vector <LScan *> allScans, int nrIt)
{

#ifdef WRITE_GRAPH_NET
  // for debug only:
  static int d = 0;
  cout << "writing graph.dat ....................................." << endl;
  d++;
  string gfilename = "graph_" + to_string(d, 3) + ".net";
  ofstream out(gfilename.c_str());
  for (int i=0; i < gr.getNrLinks(); i++) {
    int from = gr.getLink(i,0);
    int to = gr.getLink(i,1);
    // shouldn't be necessary, just in case a (out of date) graph file is loaded:
    if (from < (int)allScans.size() && to < (int)allScans.size()) {
      out << allScans[from]->get_rPos()[0] << " "
        << allScans[from]->get_rPos()[1] << " "
        << allScans[from]->get_rPos()[2] << endl
        << allScans[to  ]->get_rPos()[0] << " "
        << allScans[to  ]->get_rPos()[1] << " "
        << allScans[to  ]->get_rPos()[2] << endl << endl;
    }
  }
  out.close();
  out.clear();
#endif

  // the IdentityMatrix to transform some Scans with
  double id[16];
  M4identity(id);

  vPtPair **ptpairs = 0;                           // Contains sets of point pairs for all links
  double **centroids_m = 0, **centroids_d = 0;     // Contains centroids for all links
  Matrix B ( 3 * (gr.getNrScans()-1), 3 * (gr.getNrScans()-1) );
  SymmetricMatrix Bt ( gr.getNrScans()-1 ); Bt = 0;
  ColumnVector X( 3*(gr.getNrScans()-1) ); X = 0;
  ColumnVector T( 3*(gr.getNrScans()-1) ); T = 0;
  ColumnVector A( 3*(gr.getNrScans()-1) ); A = 0;

  B = 0.0;
  A = 0.0;

  double sum_position_diff = 0;
  double ret = DBL_MAX;

  for(int iteration = 0;
      iteration < nrIt && ret > epsilonLUM;
      iteration++) {
    sum_position_diff = 0;

    if (nrIt > 1) cout << "Iteration match " << iteration << endl;

    if (ptpairs != 0) {
      delete [] ptpairs;
    }
    if (centroids_m != 0) {
      delete [] centroids_m;
      delete [] centroids_d;
    }

    centroids_m = new double *[gr.getNrLinks()];
    centroids_d = new double *[gr.getNrLinks()];
    for (int i = 0; i < gr.getNrLinks(); i++){
      centroids_m[i] = new double[3];
      centroids_m[i][0] = centroids_m[i][1] = centroids_m[i][2] = 0.0;
      centroids_d[i] = new double[3];
      centroids_d[i][0] = centroids_d[i][1] = centroids_d[i][2] = 0.0;
    }
    ptpairs = new vPtPair*[gr.getNrLinks()];
    for (int i = 0; i < gr.getNrLinks(); i++) {
      ptpairs[i] = new vPtPair;
    }

    // Get all point pairs after ICP
    int end_loop = gr.getNrLinks();
#ifdef _OPENMP
    omp_set_num_threads(OPENMP_NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
#endif

    for (int i = 0; i < end_loop; i++) {
      cout << "P" << i << flush;
      LScan * FirstScan  = allScans[gr.getLink(i,0)];
      LScan * SecondScan = allScans[gr.getLink(i,1)];
#ifdef _OPENMP
      int thread_num = omp_get_thread_num();
#else
      int thread_num = 0;
#endif
      LScan::getPtPairs(ptpairs[i], FirstScan, SecondScan, thread_num, (int)my_icp->get_rnd(), (int)max_dist_match2_LUM, centroids_m[i], centroids_d[i]);

      // faulty network
      if (ptpairs[i]->size() <= 1) {
        cout << "Error: Link (" << gr.getLink(i,0)
          << " - " << gr.getLink(i, 1) << " ) is empty with "
          << ptpairs[i]->size() << " Corr. points. iteration = "
          << iteration << endl;
        //	   exit(1);

      } else

        sum_position_diff += genBArotForLinkedPair( gr.getLink(i,0), gr.getLink(i,1), ptpairs[i],
            centroids_m[i], centroids_d[i], &B, &A);

    }
    cout << " building rotation matrices done! " << endl;

    X = solveSparseCholesky(B, A);

    // TODO transformation bestimmen
    Bt = 0.0;
    A = 0.0;
    for ( int i = 0; i < end_loop; i++) {
      genBAtransForLinkedPair( gr.getLink(i,0), gr.getLink(i,1),
          centroids_m[i], centroids_d[i], &Bt, &A, X);
    }
    cout << " building translation matrices done! "<<endl;

    Bt = Bt.i();

    for (int i = 0; i < gr.getNrScans()-1; i++) {
      for (int j = 0; j < gr.getNrScans()-1; j++) {
        double bt = Bt.element(i, j);
        T.Rows(i*3 +1, i*3+3 ) += A.Rows(j*3 +1, j*3 +3 ) * bt;
      }
    }

    // delete ptPairs
    for (int i = 0; i < gr.getNrLinks(); i++) {
      ptpairs[i]->clear();
      delete (ptpairs[i]);
    }
    // delete centroids
    for (int i = 0; i < gr.getNrLinks(); i++){
      delete [] centroids_m[i];
      delete [] centroids_d[i];
    }

    ColumnVector t0(3), t(3), tlast(3);
    int vectorOffset;
    int loop_end = gr.getNrScans();
    double alignxf[16];

    for(int i = 1; i < loop_end; i++)
    {
      vectorOffset = (i-1) * 3;

      // Interpret results
      double x[3]  = { X(vectorOffset + 1), X(vectorOffset + 2), X(vectorOffset + 3) };
      double dx[3] = { T(vectorOffset + 1), T(vectorOffset + 2), T(vectorOffset + 3) };

      icp6D_APX::computeRt(x, dx, alignxf);

      // Update the Pose
      cout << "Old pose estimate, Scan " << i << endl;
      cout <<  "x: " << allScans[i]->get_rPos()[0]
        << " y: " << allScans[i]->get_rPos()[1]
        << " z: " << allScans[i]->get_rPos()[2]
        << " tx: " << allScans[i]->get_rPosTheta()[0]
        << " ty: " << allScans[i]->get_rPosTheta()[1]
        << " tz: " << allScans[i]->get_rPosTheta()[2]
        << endl;

      if (i < loop_end - 1) {
        allScans[i]->transform(alignxf, Scan::LUM, 1);
      } else {
        allScans[i]->transform(alignxf, Scan::LUM, 2);
      }

      cout <<  "x: " << allScans[i]->get_rPos()[0]
        << " y: " << allScans[i]->get_rPos()[1]
        << " z: " << allScans[i]->get_rPos()[2]
        << " tx: " << allScans[i]->get_rPosTheta()[0]
        << " ty: " << allScans[i]->get_rPosTheta()[1]
        << " tz: " << allScans[i]->get_rPosTheta()[2] << endl << endl;

      sum_position_diff += Len(dx);

    }
    cout << "Sum of Position differences = " << sum_position_diff << endl << endl;
    ret = (sum_position_diff / (double)gr.getNrScans());
  }

  delete [] centroids_m;
  delete [] centroids_d;

  delete [] ptpairs;
  ptpairs = 0;

  return ret;
}


