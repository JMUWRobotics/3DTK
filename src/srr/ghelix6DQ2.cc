/**
 * @file
 * @brief The implementation of globally consistent scan matching algorithm by using helix correction
 * @author Peter Schneider. Institute of Computer Science, University of Koblenz, Germany.
 */

#ifdef _MSC_VER
#ifdef OPENMP
#define _OPENMP
#endif
#endif

#include "linescan.h"
#include "ghelix6DQ2.h"
#include "slam6d/icp6Dhelix.h"
#include "suitesparse/cs.h"

#include <cfloat>
#include <fstream>
using std::flush;
#include <cstring>

#include "slam6d/globals.icc"

using std::ofstream;

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
ghelix6DQ2L::ghelix6DQ2L(icp6Dminimizer *my_icp6Dminimizer,
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
ghelix6DQ2L::~ghelix6DQ2L()
{
  delete my_icp;
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
double ghelix6DQ2L::genBBdForLinkedPair( int firstScanNum, int secondScanNum, vPtPair *ptpairs,
								Matrix *B, ColumnVector *Bd )
{
  double Btemp1[6][3];
  double Btemp2[6][3];
  memset(&Btemp1[0][0], 0, 18 * sizeof(double));
  memset(&Btemp2[0][0], 0, 18 * sizeof(double));
  double bd1[6];
  double bd2[6];
  memset(&bd1[0], 0, 6 * sizeof(double));
  memset(&bd2[0], 0, 6 * sizeof(double));

  double p1x, p1y, p1z, p2x, p2y, p2z, pDistX, pDistY, pDistZ, px1Sq, py1Sq, pz1Sq, px2Sq, py2Sq, pz2Sq;
  int n = (*ptpairs).size();
  double sum = 0;

  for (int i = 0; i < n; i++) {
    p1x = (*ptpairs)[i].p1.x;
    p1y = (*ptpairs)[i].p1.y;
    p1z = (*ptpairs)[i].p1.z;

    p2x = (*ptpairs)[i].p2.x;
    p2y = (*ptpairs)[i].p2.y;
    p2z = (*ptpairs)[i].p2.z;

    px1Sq = p1x * p1x;
    py1Sq = p1y * p1y;
    pz1Sq = p1z * p1z;

    px2Sq = p2x * p2x;
    py2Sq = p2y * p2y;
    pz2Sq = p2z * p2z;

    Btemp1[4][0] += -p2z;
//  Btemp1[3][1] += p2z;
    Btemp1[5][0] += p2y;
//  Btemp1[3][2] += -p2y;
    Btemp1[4][2] += p2x;
//  Btemp1[5][1] += -p2x;
    Btemp1[0][0] += pz2Sq + py2Sq;
    Btemp1[1][0] += p2y*-p2x;
    Btemp1[2][0] += -p2z*p2x;
    Btemp1[1][1] += pz2Sq + px2Sq;
    Btemp1[2][1] += p2z*-p2y;
    Btemp1[2][2] += px2Sq + py2Sq;

    pDistX = p1x - p2x;
    pDistY = p1y - p2y;
    pDistZ = p1z - p2z;

    bd1[0] += -p1z*pDistY + p1y*pDistZ;
    bd1[1] += p1z*pDistX - p1x*pDistZ;
    bd1[2] += -p1y*pDistX + p1x*pDistY;
    bd1[3] += pDistX;
    bd1[4] += pDistY;
    bd1[5] += pDistZ;

    bd2[0] += -p2z*-pDistY + p2y*-pDistZ;
    bd2[1] += p2z*-pDistX - p2x*-pDistZ;
    bd2[2] += -p2y*-pDistX + p2x*-pDistY;
    bd2[3] += -pDistX;
    bd2[4] += -pDistY;
    bd2[5] += -pDistZ;

    sum += pDistX*pDistX + pDistY*pDistY + pDistZ*pDistZ;
  }

#ifdef _OPENMP
  #pragma omp critical (enterB)
#endif
  {
  int matPlace1 = (firstScanNum-1) * 6;

  if(firstScanNum != 0)
  {
    (*B)(matPlace1+4,matPlace1+4) += n;
    (*B)(matPlace1+5,matPlace1+5) += n;
    (*B)(matPlace1+6,matPlace1+6) += n;

    (*B)(matPlace1+1,matPlace1+5) += Btemp1[4][0];
    (*B)(matPlace1+5,matPlace1+1) += Btemp1[4][0];
    (*B)(matPlace1+2,matPlace1+4) += -Btemp1[4][0];
    (*B)(matPlace1+4,matPlace1+2) += -Btemp1[4][0];
    (*B)(matPlace1+1,matPlace1+6) += Btemp1[5][0];
    (*B)(matPlace1+6,matPlace1+1) += Btemp1[5][0];
    (*B)(matPlace1+3,matPlace1+4) += -Btemp1[5][0];
    (*B)(matPlace1+4,matPlace1+3) += -Btemp1[5][0];
    (*B)(matPlace1+3,matPlace1+5) += Btemp1[4][2];
    (*B)(matPlace1+5,matPlace1+3) += Btemp1[4][2];
    (*B)(matPlace1+2,matPlace1+6) += -Btemp1[4][2];
    (*B)(matPlace1+6,matPlace1+2) += -Btemp1[4][2];
    (*B)(matPlace1+1,matPlace1+2) += Btemp1[1][0];
    (*B)(matPlace1+2,matPlace1+1) += Btemp1[1][0];
    (*B)(matPlace1+1,matPlace1+3) += Btemp1[2][0];
    (*B)(matPlace1+3,matPlace1+1) += Btemp1[2][0];
    (*B)(matPlace1+2,matPlace1+3) += Btemp1[2][1];
    (*B)(matPlace1+3,matPlace1+2) += Btemp1[2][1];
    (*B)(matPlace1+1,matPlace1+1) += Btemp1[0][0];
    (*B)(matPlace1+2,matPlace1+2) += Btemp1[1][1];
    (*B)(matPlace1+3,matPlace1+3) += Btemp1[2][2];

    (*Bd)(matPlace1+1) += bd1[0];
    (*Bd)(matPlace1+2) += bd1[1];
    (*Bd)(matPlace1+3) += bd1[2];
    (*Bd)(matPlace1+4) += bd1[3];
    (*Bd)(matPlace1+5) += bd1[4];
    (*Bd)(matPlace1+6) += bd1[5];
  }
//-------------------------------------------
  unsigned int matPlace2 = (secondScanNum-1) * 6;

  (*B)(matPlace2+4,matPlace2+4) += n;
  (*B)(matPlace2+5,matPlace2+5) += n;
  (*B)(matPlace2+6,matPlace2+6) += n;

  (*B)(matPlace2+1,matPlace2+5) += Btemp1[4][0];
  (*B)(matPlace2+5,matPlace2+1) += Btemp1[4][0];
  (*B)(matPlace2+2,matPlace2+4) += -Btemp1[4][0];
  (*B)(matPlace2+4,matPlace2+2) += -Btemp1[4][0];
  (*B)(matPlace2+1,matPlace2+6) += Btemp1[5][0];
  (*B)(matPlace2+6,matPlace2+1) += Btemp1[5][0];
  (*B)(matPlace2+3,matPlace2+4) += -Btemp1[5][0];
  (*B)(matPlace2+4,matPlace2+3) += -Btemp1[5][0];
  (*B)(matPlace2+3,matPlace2+5) += Btemp1[4][2];
  (*B)(matPlace2+5,matPlace2+3) += Btemp1[4][2];
  (*B)(matPlace2+2,matPlace2+6) += -Btemp1[4][2];
  (*B)(matPlace2+6,matPlace2+2) += -Btemp1[4][2];
  (*B)(matPlace2+1,matPlace2+2) += Btemp1[1][0];
  (*B)(matPlace2+2,matPlace2+1) += Btemp1[1][0];
  (*B)(matPlace2+1,matPlace2+3) += Btemp1[2][0];
  (*B)(matPlace2+3,matPlace2+1) += Btemp1[2][0];
  (*B)(matPlace2+2,matPlace2+3) += Btemp1[2][1];
  (*B)(matPlace2+3,matPlace2+2) += Btemp1[2][1];
  (*B)(matPlace2+1,matPlace2+1) += Btemp1[0][0];
  (*B)(matPlace2+2,matPlace2+2) += Btemp1[1][1];
  (*B)(matPlace2+3,matPlace2+3) += Btemp1[2][2];


  (*Bd)(matPlace2+1) += bd2[0];
  (*Bd)(matPlace2+2) += bd2[1];
  (*Bd)(matPlace2+3) += bd2[2];
  (*Bd)(matPlace2+4) += bd2[3];
  (*Bd)(matPlace2+5) += bd2[4];
  (*Bd)(matPlace2+6) += bd2[5];

//-------------------------------------------
  if( firstScanNum != 0)
  {
    (*B)(matPlace1+4,matPlace2+4) -= n;
    (*B)(matPlace1+5,matPlace2+5) -= n;
    (*B)(matPlace1+6,matPlace2+6) -= n;

    (*B)(matPlace1+1,matPlace2+5) += -Btemp1[4][0];
    (*B)(matPlace1+5,matPlace2+1) += -Btemp1[4][0];
    (*B)(matPlace1+2,matPlace2+4) += Btemp1[4][0];
    (*B)(matPlace1+4,matPlace2+2) += Btemp1[4][0];
    (*B)(matPlace1+1,matPlace2+6) += -Btemp1[5][0];
    (*B)(matPlace1+6,matPlace2+1) += -Btemp1[5][0];
    (*B)(matPlace1+3,matPlace2+4) += Btemp1[5][0];
    (*B)(matPlace1+4,matPlace2+3) += Btemp1[5][0];
    (*B)(matPlace1+3,matPlace2+5) += -Btemp1[4][2];
    (*B)(matPlace1+5,matPlace2+3) += -Btemp1[4][2];
    (*B)(matPlace1+2,matPlace2+6) += Btemp1[4][2];
    (*B)(matPlace1+6,matPlace2+2) += Btemp1[4][2];
    (*B)(matPlace1+1,matPlace2+2) += -Btemp1[1][0];
    (*B)(matPlace1+2,matPlace2+1) += -Btemp1[1][0];
    (*B)(matPlace1+1,matPlace2+3) += -Btemp1[2][0];
    (*B)(matPlace1+3,matPlace2+1) += -Btemp1[2][0];
    (*B)(matPlace1+2,matPlace2+3) += -Btemp1[2][1];
    (*B)(matPlace1+3,matPlace2+2) += -Btemp1[2][1];
    (*B)(matPlace1+1,matPlace2+1) += -Btemp1[0][0];
    (*B)(matPlace1+2,matPlace2+2) += -Btemp1[1][1];
    (*B)(matPlace1+3,matPlace2+3) += -Btemp1[2][2];
//----------------------------------------
    (*B)(matPlace2+4,matPlace1+4) -= n;
    (*B)(matPlace2+5,matPlace1+5) -= n;
    (*B)(matPlace2+6,matPlace1+6) -= n;

    (*B)(matPlace2+1,matPlace1+5) += -Btemp1[4][0];
    (*B)(matPlace2+5,matPlace1+1) += -Btemp1[4][0];
    (*B)(matPlace2+2,matPlace1+4) += Btemp1[4][0];
    (*B)(matPlace2+4,matPlace1+2) += Btemp1[4][0];
    (*B)(matPlace2+1,matPlace1+6) += -Btemp1[5][0];
    (*B)(matPlace2+6,matPlace1+1) += -Btemp1[5][0];
    (*B)(matPlace2+3,matPlace1+4) += Btemp1[5][0];
    (*B)(matPlace2+4,matPlace1+3) += Btemp1[5][0];
    (*B)(matPlace2+3,matPlace1+5) += -Btemp1[4][2];
    (*B)(matPlace2+5,matPlace1+3) += -Btemp1[4][2];
    (*B)(matPlace2+2,matPlace1+6) += Btemp1[4][2];
    (*B)(matPlace2+6,matPlace1+2) += Btemp1[4][2];
    (*B)(matPlace2+1,matPlace1+2) += -Btemp1[1][0];
    (*B)(matPlace2+2,matPlace1+1) += -Btemp1[1][0];
    (*B)(matPlace2+1,matPlace1+3) += -Btemp1[2][0];
    (*B)(matPlace2+3,matPlace1+1) += -Btemp1[2][0];
    (*B)(matPlace2+2,matPlace1+3) += -Btemp1[2][1];
    (*B)(matPlace2+3,matPlace1+2) += -Btemp1[2][1];
    (*B)(matPlace2+1,matPlace1+1) += -Btemp1[0][0];
    (*B)(matPlace2+2,matPlace1+2) += -Btemp1[1][1];
    (*B)(matPlace2+3,matPlace1+3) += -Btemp1[2][2];
  }
  }    // of pragma omp critical

  return sqrt( sum / (double) n );
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
double ghelix6DQ2L::doGraphSlam6D(Graph gr, vector <LScan *> allScans, int nrIt)
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

  vPtPair **ptpairs = 0;                 // Contains sets of point pairs for all links
  Matrix B ( 6 * (gr.getNrScans()-1), 6 * (gr.getNrScans()-1) );
  ColumnVector ccs( 6*(gr.getNrScans()-1) ), bd( 6*(gr.getNrScans()-1) );

  B = 0.0;
  bd = 0.0;

  double sum_position_diff = 0;
  double ret = DBL_MAX;

  for(int iteration = 0;
	 iteration < nrIt && ret > epsilonLUM;
	 iteration++) {
    sum_position_diff = 0;

    if (nrIt > 1) cout << "Iteration match " << iteration << endl;

    if (ptpairs != 0) delete [] ptpairs;
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

    for ( int i = 0; i < end_loop; i++) {
      cout << "P" << i << flush;
      LScan * FirstScan  = allScans[gr.getLink(i,0)];
      LScan * SecondScan = allScans[gr.getLink(i,1)];
#ifdef _OPENMP
      int thread_num = omp_get_thread_num();
#else
      int thread_num = 0;
#endif
	 double dummy_centroid_m[3];
	 double dummy_centroid_d[3];
   LScan::getPtPairs(ptpairs[i], FirstScan, SecondScan, thread_num, (int)my_icp->get_rnd(), (int)max_dist_match2_LUM, dummy_centroid_m, dummy_centroid_d);

      // faulty network
      if (ptpairs[i]->size() <= 1) {
	   cout << "Error: Link (" << gr.getLink(i,0)
		   << " - " << gr.getLink(i, 1) << " ) is empty with "
		   << ptpairs[i]->size() << " Corr. points. iteration = "
		   << iteration << endl;
	   //	   exit(1);

      } else
       // build the matrix B and vector bd
        genBBdForLinkedPair( gr.getLink(i,0), gr.getLink(i,1), ptpairs[i], &B, &bd );
    }
    cout <<" building matrices done! "<<endl;

    ccs = solveSparseCholesky(B, bd);

    // delete ptPairs
    for (int i = 0; i < gr.getNrLinks(); i++) {
	 ptpairs[i]->clear();
	 delete (ptpairs[i]);
    }

    ColumnVector t0(3), t(3), tlast(3);
    int vectorOffset;

    int loop_end = gr.getNrScans();
    double alignxfLum[16];
    for(int i = 1; i < loop_end; i++)
    {
      vectorOffset = (i-1) * 6;
      icp6D_HELIX::computeRt( &ccs, vectorOffset, alignxfLum);

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
        allScans[i]->transform(alignxfLum, Scan::LUM, 1);
      } else {
        allScans[i]->transform(alignxfLum, Scan::LUM, 2);
      }

	 cout <<  "x: " << allScans[i]->get_rPos()[0]
		 << " y: " << allScans[i]->get_rPos()[1]
		 << " z: " << allScans[i]->get_rPos()[2]
		 << " tx: " << allScans[i]->get_rPosTheta()[0]
		 << " ty: " << allScans[i]->get_rPosTheta()[1]
		 << " tz: " << allScans[i]->get_rPosTheta()[2] << endl << endl;

	 sum_position_diff += sqrt( sqr(alignxfLum[12]) + sqr(alignxfLum[13]) +  sqr(alignxfLum[14]));

    }
    cout << "Sum of Position differences = " << sum_position_diff << endl << endl;
    ret = (sum_position_diff / (double)gr.getNrScans());
  }
  delete [] ptpairs;
  ptpairs = 0;

  return ret;
}

