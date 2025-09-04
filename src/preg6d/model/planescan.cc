/** @file
 *  @brief this file contains the PlaneScan Class that represents a Scan
 *  and holds the information about Point-to-Plane correspondences.
 *
 *  @author Fabian Arzberger, JMU, Germany.
 *
 *  Released under the GPL version 3.
 */

#include <queue>

#include "scanio/helper.h"
#include "model/planescan.h"
#include "util.h"

// declare the static members
PlaneScans PlaneScan::allPlaneScans = PlaneScans(0);
unsigned int PlaneScan::idx_global = 0;
double PlaneScan::_eps_dist;
double PlaneScan::_eps_ppd;
bool PlaneScan::reduce;
Planes PlaneScan::allPlanes;
bool PlaneScan::use_correspondence_min;
bool PlaneScan::use_clustering;
int PlaneScan::k_neighbours;
double PlaneScan::d_growth;
double PlaneScan::d_growth_max_adapt;
double PlaneScan::d_thickness;
int PlaneScan::n_min_clusterpoints;
int PlaneScan::n_min_clusterpoints_far;
int PlaneScan::n_max_clusters;
double PlaneScan::percentile_filter;
double PlaneScan::eps_similarity;
string PlaneScan::cluster_outpath;
bool PlaneScan::read_normals;
bool PlaneScan::read_clusters;
bool PlaneScan::use_normal_cor;
size_t PlaneScan::min_scansize;
bool PlaneScan::export_cluster_rgb;
double PlaneScan::eigratio;
InterpolMap PlaneScan::adaptGrowth;
InterpolMap PlaneScan::adaptNrPts;
double PlaneScan::score_w_alpha;
double PlaneScan::score_w_overlap;
double PlaneScan::score_w_hesse;
double PlaneScan::score_w_eigen;
double PlaneScan::score_w_ppd;
bool PlaneScan::use_bruteforce;
int PlaneScan::_match_type;
plane_alg PlaneScan::plane_algo;

// Looks at two normals and returns if they are similar.
// Considers all quadrants, e.g. an angle between 180 deg. means equal normals.
bool PlaneScan::similar(double* n1, double* n2)
{
    return deg(angleBetweenNormals(n1, n2)) < eps_similarity;
}

// Most important constructor.
// Initializes a PlaneScan object from a scan, as read from the scan directory
// Infers additional information to the scan such as point normals or corresponding planes.
PlaneScan::PlaneScan(Scan *scan)
{
    // Read disk data
	idx_global += 1;
    idx = idx_global;
    basedir = scan->getPath();
	timestamp = scan->getTimeStamp();
    identifier = scan->getIdentifier();

    // Read pose data
	const double *scanRPos = scan->get_rPos();
	const double *scanRPosTheta = scan->get_rPosTheta();
	for (int i = 0; i < 3; i++)
    {
		rPos[i] = scanRPos[i];
		rPosTheta[i] = scanRPosTheta[i];
	}
    EulerToMatrix4(rPos, rPosTheta, orig_transMat);
	EulerToMatrix4(rPos, rPosTheta, transMat);

    // Get xyz point data
    string redstring = reduce ? " reduced show" : "";
	DataXYZ xyz(scan->get("xyz" + redstring));
    nrpts = xyz.size();

    // Read normals from scan files
    if (read_normals) {
        DataNormal scan_normals = scan->get("normal" + redstring);
        size_t nrnormals = scan_normals.size();
        normals = new double*[nrnormals];
        for (size_t i = 0; i < nrnormals; ++i) {
            normals[i] = new double[3];
            normals[i][0] = scan_normals[i][0]; //nx;
            normals[i][1] = scan_normals[i][1]; //ny;
            normals[i][2] = scan_normals[i][2]; //nz;
        }

    // Read clusters
    } else if (read_clusters) {
        string color_red_string = reduce ? "color reduced" : "rgb";
        DataRGB scan_color = scan->get(color_red_string);
        if (scan_color.size() != 0) {
            readClusters(scan_color);
        } else {
            DataType scan_type = scan->get( (reduce?"type reduced" : "type") );
            readClusters(scan_type);
        }
    }

    // Read point data, calculate normals if necessary
	if (nrpts != 0) {
		points = new double*[nrpts];

        // Store points
		for (size_t i = 0; i < nrpts; i++) {
			points[i] = new double[3];
			points[i][0] = xyz[i][0];
			points[i][1] = xyz[i][1];
			points[i][2] = xyz[i][2];
        }

        // Simple Normal Corresponances should be used, but normals are not read
        if (use_normal_cor && !read_normals) {
            cout << "Calc normals for scan" << identifier << endl;
            normals = new double*[nrpts];
            calcNormals(xyz);
        }

        // Clustering shall be used, but are not read
        if (use_clustering && !read_clusters)
        {
            // Calculate normals (if you havent read them already)
            if (!read_normals) {
                normals = new double*[nrpts];
                cout << "Calc normals for scan" << identifier << endl;
                calcNormals(xyz);
            }

            // Use normals to cluster points together
            calcLocalClusters(scan);
        }

    // Scan with index idx has no points, skip the scan.
	} else {
		cout << "Scan nr. " << idx << " skipped. Contains no data." << endl;
		empty = true;
	}
    matcher = 0;
    // Put the current PlaneScan object into the static collection
    #pragma omp critical
    allPlaneScans.push_back(this);
}

void PlaneScan::calcLocalClusters(Scan *scan)
{
    if (plane_algo == RG) {
        detectRG(); // no scan ptr needed.
    } else {
        Detect plane_detector(scan, cluster_outpath, true);
        plane_detector.detect( plane_algo );
        // TODO: delete points and normals and reallocate memory.
    }
}

void PlaneScan::reset() {
    for (int i = 0; i < 16; ++i)
        transMat[i] = orig_transMat[i];
    Matrix4ToEuler(orig_transMat, rPosTheta, rPos);
}

void PlaneScan::transform(const double* alignxf)
{
    double tempxf[16];
    MMult(alignxf, transMat, tempxf);
    memcpy(transMat, tempxf, sizeof(transMat));
    Matrix4ToEuler(transMat, rPosTheta, rPos);
}

void PlaneScan::mergeCoordinatesWithRoboterPosition(PlaneScan* prevScan, bool* dims)
{
    double tempMat[16], deltaMat[16];
    M4inv(prevScan->orig_transMat, tempMat);
    MMult(prevScan->transMat, tempMat, deltaMat);

    // apply delta transformation of the previous scan
    transform(deltaMat);

    // dimension lock
    if (dims != 0) {
        double rPos_orig[3], rPosTheta_orig[3];
        double rPos_new[3], rPosTheta_new[3];
        Matrix4ToEuler(orig_transMat, rPosTheta_orig, rPos_orig);
        Matrix4ToEuler(transMat, rPosTheta_new, rPos_new);
        // if(!dims[0]) rPosTheta_new[0] = rPosTheta_orig[0];
        // if(!dims[1]) rPosTheta_new[1] = rPosTheta_orig[1];
        // if(!dims[2]) rPosTheta_new[2] = rPosTheta_orig[2];
        if(!dims[3]) rPos_new[0] = rPos_orig[0];
        if(!dims[4]) rPos_new[1] = rPos_orig[1];
        if(!dims[5]) rPos_new[2] = rPos_orig[2];
        EulerToMatrix4(rPos_new, rPosTheta_new, transMat);
        Matrix4ToEuler(transMat, rPosTheta, rPos);
    }
}

// Calculates the Hesse distance from a transformed point "xyz" to plane "p"
double PlaneScan::dist2Plane(const double *xyz, const double *trans, const NormalPlane *p)
{
    // Convert normal and plane ref. point to coordinate system of the current scan
    double xyz_transformed[3];
    transform3(trans, xyz, xyz_transformed);
    double hesse = dist2Plane(xyz_transformed, p);
    return hesse;
}

// Calculates the Hesse distance from a point "xyz" to plane "p"
double PlaneScan::dist2Plane(const double *xyz, const NormalPlane *p)
{
    double v[3];
    v[0] = xyz[0] - p->x[0];
    v[1] = xyz[1] - p->x[1];
    v[2] = xyz[2] - p->x[2];
	return Dot(v, p->n);
}

// Calculates the polygon projection distance of a point to a plane.
double PlaneScan::projDist2Plane(const double *xyz, NormalPlane *p)
{
    double D = dist2Plane(xyz, p); // Hesse distance
    Point n(p->n); // normal vector
    Point Tp(xyz);
    Point Tp2d, p1, p2;
    Point projection( Tp - D*n );
    int polysize = p->hull.size();
    Point *polygon = p->hullAsPointArr();
    char direction = p->direction;
    Point polygon2d[polysize + 1]; // +1 because last point must be first point
    // Project everything in 2d space to check winding number there.
    NormalPlane::convert3Dto2D(
        projection, // input: 3D projection onto the corresponding plane.
        direction, // direction of plane normal.
        Tp2d // output: 2D projection of the 3D projection
    );
    NormalPlane::convert3Dto2D(
        polygon, // input: the 3D convex hull of the plane.
        polysize, // the size of the convex hull.
        direction, // the direction of plane normal.
        polygon2d // output: 2D projection of the convex hull.
    );
    // insert first point as last point again (close the polygon)
    polygon2d[polysize] = polygon2d[0];

    // // In the 2D representation, use winding number algorithm to check
    // // if the point is inside or outside of the polygon.
    int wn = NormalPlane::wn_PnPoly(Tp2d, polygon2d, polysize) ;
    bool outside = wn == 0;
    if (!outside) return 0.0;
    else return NormalPlane::nearestLineSegment(projection, polygon, polysize, p1, p2);
}

// Calculates the polygon projection distance of a transformed point to a plane.
double PlaneScan::projDist2Plane(const double *xyz, const double *trans, NormalPlane *p)
{
    // Transform point xyz
    double xyz_transformed[3];
    transform3(trans, xyz, xyz_transformed);
    double ppd = projDist2Plane(xyz_transformed, p);
    return ppd;
}

// Checks if a point "xyz" with transformation "trans" lies on the plane "p"
bool PlaneScan::isInPlane(double *xyz, double *trans, NormalPlane *p)
{
    double hesse = dist2Plane(xyz, trans, p);
    double ppd = projDist2Plane(xyz, trans, p);
    // Calculate distance to plane and apply threshold
    return fabs(hesse) < _eps_dist && fabs( ppd ) < _eps_ppd;
}

// Static setter that sets the distance threshold for plane correspondence.
// Use before constructing the object!
void PlaneScan::setEpsDist(double eps_dist)
{
	_eps_dist = eps_dist;
}

void PlaneScan::setEpsPPD(double eps_ppd)
{
    _eps_ppd = eps_ppd;
}

void PlaneScan::setUseCorrespondenceMin(bool val)
{
    use_correspondence_min = val;
}

void PlaneScan::setKNearest(int k)
{
    k_neighbours = k;
}

void PlaneScan::setUseClustering(bool val)
{
    use_clustering = val;
    if (use_clustering) {
        use_normal_cor = false;
    }
}

void PlaneScan::setPlaneAlgo(plane_alg val)
{
    plane_algo = val;
    if (!use_clustering)
        setUseClustering(true);
}

void PlaneScan::setGrowthThresh(double val)
{
    d_growth = val;
}

void PlaneScan::setMaxPlaneThickness(double val)
{
    d_thickness = val;
}

void PlaneScan::setEpsSimilarity(double val)
{
    eps_similarity = val;
}

void PlaneScan::setMinClusterPoints(int n)
{
    n_min_clusterpoints = n;
}

void PlaneScan::setMinClusterPointsFar(int n)
{
    n_min_clusterpoints_far = n;
}

void PlaneScan::setClusterOutputPath(string path)
{
    cluster_outpath = path;
}

void PlaneScan::setExportClusterOutputRGB(bool val)
{
    export_cluster_rgb = val;
}

void PlaneScan::setReadNormals(bool val)
{
    read_normals = val;
    if (read_normals && read_clusters) {
        read_clusters = false;
    }
}

void PlaneScan::setReadClusters(bool val)
{
    read_clusters = val;
    if (read_normals && read_clusters) {
        read_normals = false;
    }
}

void PlaneScan::setUseNormalCor(bool val) {
    use_normal_cor = val;
}

void PlaneScan::setMinScansize(size_t val) {
    min_scansize = val;
}

void PlaneScan::setMaxEigRatio(double val) {
    eigratio = val;
}

void PlaneScan::setGrowthMaxAdapt(double val) {
    d_growth_max_adapt = val;
}

void PlaneScan::setWeightAlpha(double val) {
    score_w_alpha = val;
}

void PlaneScan::setWeightOverlap(double val) {
    score_w_overlap = val;
}

void PlaneScan::setWeightHesse(double val) {
    score_w_hesse = val;
}

void PlaneScan::setWeightEigen(double val) {
    score_w_eigen = val;
}

void PlaneScan::setWeightPPD(double val) {
    score_w_ppd = val;
}

void PlaneScan::setUseRegionGrowingBruteForce(bool val) {
    use_bruteforce = val;
}

void PlaneScan::setPercentileFilter(double val) {
    percentile_filter = val;
}

void PlaneScan::setMaxClusters(int val) {
    n_max_clusters = val;
}

void PlaneScan::calcNormals(DataXYZ &xyz)
{
    vector<Point> ps; // points
    vector<Point> ns; // normals
    ps.reserve(xyz.size());
    ns.reserve(xyz.size());
    for(size_t j = 0; j < nrpts; j++)
    {
        double tmp[3] = { points[j][0], points[j][1], points[j][2] };
        // transform3(transMat, xyz[j], tmp);
        ps.push_back(Point( tmp ));
    }
    // Calc normals with k nearest neighbours
    double rPos_tmp[3] = {0, 0, 0};
    calculateNormalsIndexedKNN(ns, ps, k_neighbours, rPos_tmp);
    // Store normals in PlaneScan object
    for (size_t i = 0; i < nrpts; i++)
    {
        normals[i] = new double[3];
        normals[i][0] = ns.at(i).x;
        normals[i][1] = ns.at(i).y;
        normals[i][2] = ns.at(i).z;
    }
}

void PlaneScan::labelPoints(int match_type, bool quiet)
{
    if (!quiet) cout << "Labeling points on planes for scan" << identifier;
    correspondences.clear();
    global_matches.clear();
    global_mismatches.clear();
    point_pairs.clear();
    local_matches.clear();
    //delete matcher;
    if ( !matcher)
    {
        switch(match_type) {
            case 1:
            //single euklid approach (only where no ambiguities are present)
            matcher = new EuklidMatcher(this, false);
            break;
            case 2:
            //min euklid approach (take min distance at ambiguity)
            matcher = new EuklidMatcher(this, true);
            break;
            case 3:
            //min euklid approach with normal threshold
            matcher = new NormalMatcher(this, eps_similarity);
            break;
            case 4:
            //based on similarity score
            matcher = new ClusterMatcher(this, score_w_overlap,
                score_w_alpha, score_w_hesse, score_w_ppd, score_w_eigen, _eps_dist,
                _eps_ppd, eps_similarity, eigratio);
            break;
            case 5:
            //based on minimum energy
            matcher = new PlaneMatcher(this, _eps_dist, _eps_ppd, eps_similarity);
            break;
            default:
            cout << "This matcher has not been implemented." << endl;
            matcher = 0;
            break;
        }
    }
    matcher->match();

    if (!quiet)
    {
        int n = correspondences.size();
        cout << "... "
             << n << " correspondence" << (n == 1 ? "":"s")
             << " found." << endl;
    }
    _match_type = match_type;
}

vector<PointPlane*> PlaneScan::getGlobalPointPlanes()
{
    vector< PointPlane* > result;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        vector<double*> cluster;
        for (size_t j = 0; j < clusters[i].size(); j++)
        {
            double *p = new double[3];
            transform3(transMat, points[clusters[i][j]], p);
            cluster.push_back( p );
        }

        NormalPlane* plane = new NormalPlane(cluster);

        for (size_t j = 0; j < clusters[i].size(); j++)
            result.push_back( new PointPlane( points[clusters[i][j]], plane ) );

    }
    return result;
}

vector<PointPlane*> PlaneScan::getLocalPointPlanes()
{
    vector< PointPlane* > result;
    for (size_t i = 0; i < clusters.size(); i++)
    {
        vector<double*> cluster;
        for (size_t j = 0; j < clusters[i].size(); j++)
            cluster.push_back( points[clusters[i][j]] );

        NormalPlane* plane = new NormalPlane(cluster);

        for (size_t j = 0; j < clusters[i].size(); j++)
            result.push_back( new PointPlane( points[clusters[i][j]], plane ) );

    }
    return result;
}

PlaneScan::PlaneScan(Scan *scan, Planes &planes) : PlaneScan(scan)
{
    // Save plane data
    allPlanes = planes;
    // Label the points
	labelPoints(_match_type, false);
}

void PlaneScan::setPlanes(Planes &ps)
{
    allPlanes = ps;
}

// Return the point at index <n>
const double* PlaneScan::operator()(int n)
{
	return points[n];
}

bool PlaneScan::isEmpty()
{
    //
	return correspondences.size() <= min_scansize;
}

bool PlaneScan::isOkay()
{
	return !isEmpty();
}

size_t PlaneScan::index()
{
    return idx;
}

void PlaneScan::setReduce(bool val)
{
    reduce = val;
}

// instance-function, calculates squared point-2-plane error of one Scan
double PlaneScan::calcErr()
{
	double sum = 0.0;
	Correspondences::iterator it = correspondences.begin();
	for (; it != correspondences.end() ; ++it )
    {
        sum += sqr(dist2Plane(it->first, transMat, it->second));
    }

    // PointPairs::iterator pit = point_pairs.begin();
    // for (; pit != point_pairs.end(); ++pit) {
    //     double ptrans[3];
    //     transform3(transMat, pit->first, ptrans);
    //     sum += Dist2(ptrans , pit->second);
    // }

	return sqrt(sum / correspondences.size());
}

double PlaneScan::calcErrPtPairs()
{
    double sum = 0.0;
    PointPairs::iterator pit = point_pairs.begin();
    for (; pit != point_pairs.end(); ++pit) {
        double ptrans[3];
        transform3(transMat, pit->first, ptrans);
        sum += Dist2(ptrans , pit->second);
    }
    return sqrt(sum / correspondences.size());
}

// static function, calculates total point-2-plane distance of all Scans
double PlaneScan::totalError()
{
	double sum = 0.0;
	for ( const auto& pscan : allPlaneScans) sum += pscan->calcErr();
	return sum;
}

void PlaneScan::addFrame(Scan::AlgoType type)
{
    sout << transMat << type << endl;
}

string PlaneScan::getFrameData()
{
    return sout.str();
}

bool PlaneScan::writeFrames()
{
    map<std::string, std::string> framesdata;
    for (size_t i = 0; i < PlaneScan::allPlaneScans.size(); i++) {
        PlaneScan *pScan = PlaneScan::allPlaneScans[i];
        string framespath = boost::filesystem::path(pScan->basedir + "scan" + pScan->identifier + ".frames").string();
        framesdata[framespath] = pScan->getFrameData();
    }
    bool state = write_multiple(framesdata, ios_base::out);
    return state;
}

// Little helper to count lines in a string.
// Works for UNIX systems only.
unsigned int countLines(const string &s)
{
    int nlines = 0;
    for( char c : s ) if ( c == '\n' ) ++nlines;
    return nlines;
}

// Destructor
PlaneScan::~PlaneScan()
{
	allPlaneScans.clear();
	correspondences.clear();
}