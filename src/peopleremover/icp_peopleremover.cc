#include <peopleremover/common.h>

#include <boost/filesystem.hpp>

#include <slam6d/scan.h>
#include <slam6d/icp6D.h>
#include <slam6d/icp6Dsvd.h>

#ifdef WITH_MMAP_SCAN
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#endif

namespace po = boost::program_options;

// remove dynamic points
static void update_for_icp(Scan* scan, const std::set<size_t> &dynamic_indices) {
	DataXYZ xyz(scan->get("xyz"));
	std::vector<Point> static_pts;
	static_pts.reserve(xyz.size() - dynamic_indices.size());
	for (size_t i = 0; i < xyz.size(); ++i) {
		if (dynamic_indices.find(i) == dynamic_indices.end()) {
			const double *p = xyz[i];
			static_pts.emplace_back(p[0], p[1], p[2]);
		}
	}
    scan->clear("xyz");
    DataPointer dp = scan->create("xyz", sizeof(double)*3*static_pts.size());
    double *dptr = reinterpret_cast<double*>(dp.get_raw_pointer());
	for (size_t k = 0; k < static_pts.size(); ++k) {
		dptr[3*k]   = static_pts[k].x;
		dptr[3*k+1] = static_pts[k].y;
		dptr[3*k+2] = static_pts[k].z;
	}
}

int main(int argc, char *argv[])
{
	ssize_t start, end;
	IOType format;
	double fuzz;
	double voxel_size;
	size_t diff, normal_knearest;
	size_t cluster_size;
	normal_method_t normal_method;
	maxrange_method_t maxrange_method;
	std::string maskdir;
	std::string staticdir;
	std::string dir;
	bool no_subvoxel_accuracy;
	bool write_maxranges;
	int jobs;
	std::pair<size_t, double> reduce;

	size_t done = 0;

	int iterations = 3;
	int icp_iter	= 100;
	double icp_maxdist = 1.0;
#ifdef WITH_MMAP_SCAN
	std::string cachedir;
#endif

	parse_cmdline(argc, argv, start, end, format, fuzz, voxel_size, diff,
		      normal_knearest, cluster_size, normal_method,
		      maxrange_method, maskdir, staticdir, dir,
		      no_subvoxel_accuracy, write_maxranges, jobs, reduce//,
			  //iterations, icp_iter, icp_maxdist
#ifdef WITH_MMAP_SCAN
		      ,
		      cachedir
#endif
	);

	double voxel_diagonal = sqrt(3 * voxel_size * voxel_size);

	std::cerr << "dir: " << dir << std::endl;

#ifndef _MSC_VER
	double elapsed;
	struct timespec before, after;
	clock_gettime(CLOCK_MONOTONIC, &before);
#endif

#ifdef WITH_MMAP_SCAN
	Scan::openDirectory(false, dir, format, start, end, cachedir);
#else
	Scan::openDirectory(false, dir, format, start, end);
#endif
	if (Scan::allScans.size() == 0) {
		std::cerr << "No scans found. Did you use the correct format?"
			  << std::endl;
		exit(-1);
	}

	std::unordered_map<size_t, DataXYZ> points_by_slice;
	std::unordered_map<size_t, DataReflectance> reflectances_by_slice;
	std::unordered_map<size_t, DataRGB> rgb_by_slice;
	std::unordered_map<size_t, DataXYZ> orig_points_by_slice;
	std::unordered_map<
	    size_t, std::tuple<const double *, const double *, const double *>>
	    trajectory;
#ifdef WITH_MMAP_SCAN
	std::unordered_map<size_t, int> mmap_fds;
#endif
	std::cerr << "size: " << Scan::allScans.size() << std::endl;
	std::vector<size_t> scanorder;
	for (size_t id = 0; id < Scan::allScans.size(); ++id) {
		size_t i = id + start;
		scanorder.push_back(i);
		Scan *scan = Scan::allScans[id];
		/* The range filter must be set *before* transformAll() because
		 * otherwise, transformAll will move the point coordinates such
		 * that the rangeFilter doesn't filter the right points anymore.
		 * This in turn can lead to the situation that when retrieving
		 * the reflectance values, some reflectance values are not
		 * returned because when reading them in, those get filtered by
		 * the *original* point coordinates. This in turn can lead to
		 * the situation that the vector for xyz and reflectance are of
		 * different length.
		 */
		scan->setRangeFilter(-1, voxel_diagonal);
		DataXYZ xyz_orig(scan->get("xyz"));
		// copy points
		size_t raw_orig_data_size =
		    xyz_orig.size() * sizeof(double) * 3;
		unsigned char *xyz_orig_data;
#ifdef WITH_MMAP_SCAN
		if (!cachedir.empty()) {
			char filename[] = "ppl_XXXXXX";
			int fd = mkstemp(filename);
			if (fd == -1) {
				throw std::runtime_error(
				    "cannot create temporary file");
			}
			// by unlinking the file now, we make sure that there
			// are no leftover files even if the process is killed
			unlink(filename);
			int ret = fallocate(fd, 0, 0, raw_orig_data_size);
			if (ret == -1) {
				throw std::runtime_error("cannot fallocate");
			}
			xyz_orig_data = (unsigned char *)mmap(
			    NULL, raw_orig_data_size, PROT_READ | PROT_WRITE,
			    MAP_SHARED, fd, 0);
			if (xyz_orig_data == MAP_FAILED) {
				throw std::runtime_error(
				    std::string("cannot mmap: ") +
				    std::string(std::strerror(errno)) +
				    std::string(" ") + std::to_string(errno));
			}
			mmap_fds[i] = fd;
		} else {
#endif
			xyz_orig_data = new unsigned char[raw_orig_data_size];
#ifdef WITH_MMAP_SCAN
		}
#endif
		memcpy(xyz_orig_data, xyz_orig.get_raw_pointer(),
		       raw_orig_data_size);
		orig_points_by_slice[i] =
		    DataPointer(xyz_orig_data, raw_orig_data_size);
		/*
		// now that the original coordinates are saved, transform
		scan->transformAll(scan->get_transMatOrg());
		trajectory[i] =
		    std::make_tuple(scan->get_rPos(), scan->get_rPosTheta(),
				    scan->get_transMatOrg());
		DataXYZ xyz(scan->get("xyz"));
		*/
		DataReflectance refl(scan->get("reflectance"));
		if (refl.size() != 0) {
			if (xyz_orig.size() != refl.size()) {
				std::cerr << "refl mismatch: " << xyz_orig.size()
					  << " vs. " << refl.size()
					  << std::endl;
				exit(1);
			}
			reflectances_by_slice[i] = refl;
		}
		DataRGB rgb(scan->get("rgb"));
		if (rgb.size() != 0) {
			if (xyz_orig.size() != rgb.size()) {
				std::cerr << "rgb mismatch: " << xyz_orig.size()
					  << " vs. " << rgb.size() << std::endl;
				exit(1);
			}
			rgb_by_slice[i] = rgb;
		}
		//points_by_slice[i] = xyz;
		std::cerr << "number of points in scan " << i << ": "
			  << xyz_orig.size() << std::endl;
	}
#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &after);
	elapsed = (after.tv_sec - before.tv_sec);
	elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
	std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif
	std::set<struct voxel> last_free_voxels;
	std::unordered_map<struct voxel, std::set<size_t>> last_half_voxels;
	for (int iter = 0; iter < iterations; ++iter) {
		std::cerr << "\niteration " << iter+1 << " / " << iterations << "\n" << std::endl;
		// registration
		{
			icp6Dminimizer *icpMin = new icp6D_SVD(false);
			icp6D *my_icp = new icp6D(icpMin, icp_maxdist, icp_iter, false, false, -1, false, -1, 0.00001, simpleKD);
			std::cerr << "ICP..." << std::endl;
			my_icp->doICP(Scan::allScans, CLOSEST_POINT);
			delete my_icp; 
			delete icpMin;
		}
		// rebuild global points and trajectory
		points_by_slice.clear();
		trajectory.clear();
		for (size_t scan_id = 0; scan_id < Scan::allScans.size(); ++scan_id) {
			size_t ii = scan_id + start;
			DataXYZ orig = orig_points_by_slice[ii];
			size_t n = orig.size();
			double *glob = new double[3 * n];
			const double *transMat = Scan::allScans[scan_id]->get_transMatOrg();
			for (size_t k = 0; k < n; ++k) {
				double p[3] = {orig[k][0], orig[k][1], orig[k][2]};
				transform3(transMat, p);
				glob[3*k]   = p[0]; glob[3*k+1] = p[1]; glob[3*k+2] = p[2];
			}
			  points_by_slice[ii] = DataXYZ(DataPointer((unsigned char*)glob, n*sizeof(double)*3));
			trajectory[ii] = std::make_tuple(Scan::allScans[scan_id]->get_rPos(), Scan::allScans[scan_id]->get_rPosTheta(), transMat);
		}
	std::cerr << "calculate voxel occupation" << std::endl;

#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &before);
#endif
	std::unordered_map<struct voxel, std::set<size_t>>
	    voxel_occupied_by_slice;
	std::cerr << "0 %\r";
	std::cerr.flush();
	for (std::pair<size_t, DataXYZ> element : points_by_slice) {
		for (size_t i = 0; i < element.second.size(); ++i) {
			voxel_occupied_by_slice
			    [voxel_of_point(element.second[i], voxel_size)]
				.insert(element.first);
		}
		std::cerr << ((done + 1) * 100.0f / points_by_slice.size())
			  << " %\r";
		std::cerr.flush();
		done += 1;
	}
	std::cerr << std::endl;
#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &after);
	elapsed = (after.tv_sec - before.tv_sec);
	elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
	std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif

	if (voxel_occupied_by_slice.size() == 0) {
		std::cerr << "no voxel occupied" << std::endl;
		exit(1);
	}

	std::cerr << "occupied voxels: " << voxel_occupied_by_slice.size()
		  << std::endl;

	if (maxrange_method != NONE) {
		std::cerr << "compute maxranges and walk voxels" << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &before);
#endif
	} else {
		std::cerr << "walk voxels" << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &before);
#endif
	}
	std::cerr << "0 %\r";
	std::cerr.flush();

	std::set<struct voxel> free_voxels;

	/*
	 *  we need a separate variable to keep track of how many scans were
	 * done because they are not done in order when execution happens in
	 * parallel
	 */
	done = 0;
#ifdef _OPENMP
	omp_set_num_threads(jobs);
#pragma omp parallel for schedule(dynamic)
#endif
	for (
#if defined(_MSC_VER) and defined(_OPENMP)
	    // MSVC only supports OpenMP 2.5 where the counter must be signed
	    // There is also no ssize_t on non-POSIX platforms but sizeof(long)
	    // == sizeof(void*)
	    long
#else
	    size_t
#endif
		idx = 0;
	    idx < scanorder.size(); ++idx) {
		size_t i = scanorder[idx];
		// FIXME: this is just wasting memory
		std::vector<double> maxranges(
		    points_by_slice[i].size(),
		    std::numeric_limits<double>::infinity());
		std::vector<size_t> reduced;

		if (maxrange_method == NORMALS) {
			if (normal_method == KNEAREST_GLOBAL ||
			    normal_method == RANGE_GLOBAL) {
				exit(1);
			}
			std::cerr << "building spherical quad tree"
				  << std::endl;
			QuadTree qtree = QuadTree(orig_points_by_slice[i]);
			compute_maxranges(maxranges, qtree,
					  orig_points_by_slice[i],
					  normal_method, voxel_diagonal, fuzz);
			if (reduce.first != 0 && reduce.second != 0) {
				std::cerr << "reducing points..." << std::endl;
				reduced = qtree.reduce(
				    reduce.second, // angular diameter theta
				    reduce.first   // number of points
				);
			} else {
				for (size_t j = 0;
				     j < points_by_slice[i].size(); ++j) {
					reduced.push_back(j);
				}
			}
		} else if (maxrange_method == ONENEAREST) {
			exit(1);
		}
		// FIXME: fuzz also needs to be applied with maxrange_method ==
		// NONE

		if (write_maxranges) {
			std::cerr << "write maxranges" << std::endl;
			exit(1);
		}

		std::set<struct voxel> free;
		std::cerr << "shooting rays to " << reduced.size() << " points"
			  << std::endl;
		for (size_t j : reduced) {
			double p[3] = {points_by_slice[i][j][0],
				       points_by_slice[i][j][1],
				       points_by_slice[i][j][2]};
			if (maxranges[j] !=
			    std::numeric_limits<double>::infinity()) {
				double maxrange = maxranges[j];
				double r = Len(orig_points_by_slice[i][j]);
				double factor = maxrange / r;
				p[0] = orig_points_by_slice[i][j][0] * factor;
				p[1] = orig_points_by_slice[i][j][1] * factor;
				p[2] = orig_points_by_slice[i][j][2] * factor;
				transform3(std::get<2>(trajectory[i]), p);
			}
			struct visitor_args data = {};
			data.empty_voxels = &free;
			data.voxel_occupied_by_slice = &voxel_occupied_by_slice;
			data.current_slice = i;
			data.diff = diff;
			walk_voxels(std::get<0>(trajectory[i]), p, voxel_size,
				    visitor, &data);
		}
#ifdef _OPENMP
#pragma omp critical
#endif
		{
			for (struct voxel v : free) {
				free_voxels.insert(v);
			}
			std::cerr << ((done + 1) * 100.0f / scanorder.size())
				  << " %\r";
			std::cerr.flush();
			done += 1;
		}
	}
	std::cerr << std::endl;
#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &after);
	elapsed = (after.tv_sec - before.tv_sec);
	elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
	std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif

	std::cerr << "number of freed voxels: " << free_voxels.size() << " ("
		  << (free_voxels.size() * 100.0f /
		      voxel_occupied_by_slice.size())
		  << "% of occupied voxels)" << std::endl;

	if (cluster_size > 1) {
		std::cerr << "clustering voxels" << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &before);
#endif
		std::unordered_map<struct voxel, size_t> voxel_to_cluster;
		std::unordered_map<size_t, std::set<struct voxel>>
		    cluster_to_voxel;
		size_t i = 0;
		for (const struct voxel v : free_voxels) {
			std::set<size_t> neighbor_clusters;
			int vradius = 1;
			for (int i = 0 - vradius; i <= vradius; ++i) {
				for (int j = 0 - vradius; j <= vradius; ++j) {
					for (int k = 0 - vradius; k <= vradius;
					     ++k) {
						if (i == 0 && j == 0 &&
						    k == 0) {
							continue;
						}
						std::unordered_map<
						    struct voxel,
						    size_t>::iterator it2 =
						    voxel_to_cluster.find(
							voxel(v.x + i, v.y + j,
							      v.z + k));
						if (it2 !=
						    voxel_to_cluster.end()) {
							neighbor_clusters
							    .insert(
								it2->second);
						}
					}
				}
			}
			if (neighbor_clusters.size() == 0) {
				voxel_to_cluster[v] = i;
				cluster_to_voxel[i].insert(v);
				i += 1;
				continue;
			}
			if (neighbor_clusters.size() == 1) {
				size_t cluster = *(neighbor_clusters.begin());
				voxel_to_cluster[v] = cluster;
				cluster_to_voxel[cluster].insert(v);
				i += 1;
				continue;
			}
			size_t mincluster = *(neighbor_clusters.begin());
			for (size_t cluster : neighbor_clusters) {
				if (cluster == mincluster) {
					continue;
				}
				for (struct voxel v2 :
				     cluster_to_voxel[cluster]) {
					voxel_to_cluster[v2] = mincluster;
					cluster_to_voxel[mincluster].insert(v2);
				}
				cluster_to_voxel.erase(cluster);
			}
			voxel_to_cluster[v] = mincluster;
			cluster_to_voxel[mincluster].insert(v);
			i += 1;
		}
		for (std::pair<size_t, std::set<struct voxel>> p :
		     cluster_to_voxel) {
			if (p.second.size() >= cluster_size) {
				continue;
			}
			for (struct voxel voxel : p.second) {
				free_voxels.erase(voxel);
			}
		}
		std::cerr << "number of free voxels after clustering: "
			  << free_voxels.size() << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &after);
		elapsed = (after.tv_sec - before.tv_sec);
		elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
		std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif
	}

	std::unordered_map<struct voxel, std::set<size_t>> half_voxels;
	if (!no_subvoxel_accuracy) {
		std::cerr << "calculate half-free voxels" << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &before);
#endif
		for (struct voxel v : free_voxels) {
			std::set<struct voxel> neighbor_voxels;
			int vradius = 2;
			for (int i = 0 - vradius; i <= vradius; ++i) {
				for (int j = 0 - vradius; j <= vradius; ++j) {
					for (int k = 0 - vradius; k <= vradius;
					     ++k) {
						if (i == 0 && j == 0 &&
						    k == 0) {
							continue;
						}
						struct voxel neighbor = voxel(
						    v.x + i, v.y + j, v.z + k);
						if (free_voxels.find(
							neighbor) !=
						    free_voxels.end()) {
							continue;
						}
						if (voxel_occupied_by_slice
							.find(neighbor) ==
						    voxel_occupied_by_slice
							.end()) {
							continue;
						}
						for (size_t num :
						     voxel_occupied_by_slice
							 [v]) {
							half_voxels[neighbor]
							    .insert(num);
						}
					}
				}
			}
		}
		std::cerr << "number of half-free voxels: "
			  << half_voxels.size() << std::endl;
#ifndef _MSC_VER
		clock_gettime(CLOCK_MONOTONIC, &after);
		elapsed = (after.tv_sec - before.tv_sec);
		elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
		std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif
	}

#if 0
	std::cerr << "computing connectivity graph..." << std::endl;
	/*
	 * Output a graph where the nodes are the scan ids and edges connect nodes
	 * of scans that share at least one voxel
	 */
	{
		done = 0;
		std::cerr << "0 %\r";
		std::cerr.flush();
		/*
		 * go through all voxels and add all possible pairs of scan
		 * identifiers from a lower id to a higher id to a set
		 */
		std::map<std::pair<size_t, size_t>, size_t> result;
		for (std::pair<struct voxel, std::set<size_t>> el : voxel_occupied_by_slice) {
			/*
			 * the elements of std::set are sorted in ascending order, thus we
			 * do not need to sort them ourselves
			 */
			std::vector<size_t> slices(el.second.begin(), el.second.end());
			for (size_t i = 0; i < slices.size(); ++i) {
				for (size_t j = i+1; j < slices.size(); ++j) {
					result[std::make_pair(slices[i], slices[j])] += 1;
				}
			}
			std::cerr << ((done+1)*100.0f/voxel_occupied_by_slice.size()) << " %\r";
			std::cerr.flush();
			done += 1;
		}
		std::cerr << std::endl;
		size_t num_edges_complete_graph = scanorder.size()*(scanorder.size()-1)/2.0;
		std::cerr << "Scan pairs sharing the same voxel: " << result.size() << std::endl;
		std::cerr << "number edges complete graph: " << num_edges_complete_graph << std::endl;
		std::cerr << "graph G {" << std::endl;
		for (std::pair<std::pair<size_t, size_t>, size_t> el : result) {
			std::cout << "  " << el.first.first << " -- " << el.first.second << " [weight=" << el.second << "];" << std::endl;
		}
		std::cerr << "}" << std::endl;
	}
#endif
		last_free_voxels = free_voxels;
		last_half_voxels = half_voxels;
		std::unordered_map<size_t,std::set<size_t>> dyn_idx;
		for(size_t ii: scanorder) {
			for(size_t j=0;j<points_by_slice[ii].size();++j) {
				struct voxel v = voxel_of_point(points_by_slice[ii][j], voxel_size);
				if(free_voxels.find(v)!=free_voxels.end() ||
					(half_voxels.find(v)!=half_voxels.end() && half_voxels[v].count(ii)))
					dyn_idx[ii].insert(j);
			}
		}
		for(size_t scan_id=0; scan_id<Scan::allScans.size(); ++scan_id) {
			size_t ii = scan_id + start;
			update_for_icp(Scan::allScans[scan_id], dyn_idx[ii]);
		}
	}
	std::cerr << "write partitioning" << std::endl;
#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &before);
#endif

	/*
	 * we use a FILE object instead of an ofstream to write the data because
	 * it gives better performance (probably because it's buffered) and
	 * because of the availability of fprintf
	 */
	std::string outdir = dir + "/pplremover";
	boost::filesystem::create_directory(outdir);
	FILE *out_static = fopen((outdir + "/scan000.3d").c_str(), "wb");
	FILE *out_dynamic = fopen((outdir + "/scan001.3d").c_str(), "wb");
	done = 0;
	std::cerr << "0 %\r";
	std::cerr.flush();
	for (size_t i : scanorder) {
		std::unordered_map<size_t, DataReflectance>::const_iterator
		    refl_it = reflectances_by_slice.find(i);
		std::unordered_map<size_t, DataRGB>::const_iterator rgb_it =
		    rgb_by_slice.find(i);
		for (size_t j = 0; j < points_by_slice[i].size(); ++j) {
			FILE *out;
			struct voxel voxel =
			    voxel_of_point(points_by_slice[i][j], voxel_size);
			if (last_free_voxels.find(voxel) != last_free_voxels.end() ||
			    (last_half_voxels.find(voxel) != last_half_voxels.end() &&
			     last_half_voxels[voxel].find(i) !=
				 last_half_voxels[voxel].end())) {
				out = out_dynamic;
			} else {
				out = out_static;
			}
			// we print the mantissa with 13 hexadecimal digits
			// because the mantissa for double precision is 52 bits
			// long which is 6.5 bytes and thus 13 hexadecimal
			// digits
			//
			// According to
			// http://pubs.opengroup.org/onlinepubs/000095399/functions/printf.html
			// the %a conversion specifier automatically picks the
			// right precision to represent the value exactly.
			fprintf(out, "%a %a %a", points_by_slice[i][j][0],
				points_by_slice[i][j][1],
				points_by_slice[i][j][2]);
			if (refl_it != reflectances_by_slice.end()) {
				fprintf(out, " %a", refl_it->second[j]);
			}
			if (rgb_it != rgb_by_slice.end()) {
				// FIXME: rgb_it->second[j] is a RGB char array
				// not a single value
				fprintf(out, " %a", rgb_it->second[j]);
			}
			fprintf(out, "\n");
		}
		std::cerr << ((done + 1) * 100.0f / scanorder.size()) << " %\r";
		std::cerr.flush();
		done += 1;
	}
	std::cerr << std::endl;
	fclose(out_static);
	fclose(out_dynamic);

	FILE *pose_static = fopen((outdir + "/scan000.pose").c_str(), "wb");
	fprintf(pose_static, "0 0 0\n0 0 0\n");
	fclose(pose_static);
	FILE *pose_dynamic = fopen((outdir + "/scan001.pose").c_str(), "wb");
	fprintf(pose_dynamic, "0 0 0\n0 0 0\n");
	fclose(pose_dynamic);

	std::cerr << "write cleaned static scans" << std::endl;
	done = 0;
	std::cerr << "0 %\r";
	std::cerr.flush();
	if (staticdir == "") {
		staticdir = dir + "/pplremover/static";
	}
	boost::filesystem::create_directory(staticdir);
	for (std::pair<size_t, DataXYZ> element : points_by_slice) {
		size_t i = element.first;
		std::unordered_map<size_t, DataReflectance>::const_iterator
		    refl_it = reflectances_by_slice.find(i);
		std::ostringstream out;
		out << staticdir << "/scan" << std::setw(3) << std::setfill('0')
		    << i << ".3d";
		FILE *out_static = fopen(out.str().c_str(), "wb");
		if (out_static == NULL) {
			std::cerr << "cannot open" << out.str() << std::endl;
			exit(1);
		}
		for (size_t j = 0; j < element.second.size(); ++j) {
			int ret;
			double refl = 0;
			if (refl_it != reflectances_by_slice.end()) {
				refl = refl_it->second[j];
			}
			struct voxel voxel =
			    voxel_of_point(element.second[j], voxel_size);
			if (last_free_voxels.find(voxel) == last_free_voxels.end() &&
			    (last_half_voxels.find(voxel) == last_half_voxels.end() ||
			     last_half_voxels[voxel].find(i) ==
				 last_half_voxels[voxel].end())) {
				ret = fprintf(out_static, "%a %a %a %a\n",
					      orig_points_by_slice[i][j][0],
					      orig_points_by_slice[i][j][1],
					      orig_points_by_slice[i][j][2],
					      // 0.0f
					      refl);
			}
			if (ret < 0) {
				std::cerr << "failed to write to " << out.str()
					  << std::endl;
				exit(1);
			}
		}
		fclose(out_static);
		std::ostringstream out_pose;
		out_pose << staticdir << "/scan" << std::setw(3)
			 << std::setfill('0') << i << ".pose";
		FILE *pose = fopen(out_pose.str().c_str(), "wb");
		int ret =
		    fprintf(pose, "%.17f %.17f %.17f\n%.17f %.17f %.17f\n",
			    std::get<0>(trajectory[i])[0],
			    std::get<0>(trajectory[i])[1],
			    std::get<0>(trajectory[i])[2],
			    std::get<1>(trajectory[i])[0] * 180 / M_PI,
			    std::get<1>(trajectory[i])[1] * 180 / M_PI,
			    std::get<1>(trajectory[i])[2] * 180 / M_PI);
		if (ret < 0) {
			std::cerr << "failed to write to " << out_pose.str()
				  << std::endl;
			exit(1);
		}
		fclose(pose);
		std::cerr << ((done + 1) * 100.0f / points_by_slice.size())
			  << " %\r";
		std::cerr.flush();
		done += 1;
	}
	std::cerr << std::endl;

	std::cerr << "write masks" << std::endl;
	done = 0;
	std::cerr << "0 %\r";
	std::cerr.flush();

	if (maskdir == "") {
		maskdir = dir + "/pplremover/masks";
	}
	boost::filesystem::create_directory(maskdir);
	for (std::pair<size_t, DataXYZ> element : points_by_slice) {
		std::ostringstream out;
		out << maskdir << "/scan" << std::setw(3) << std::setfill('0')
		    << element.first << ".mask";
		FILE *out_mask = fopen(out.str().c_str(), "wb");
		if (out_mask == NULL) {
			std::cerr << "cannot open" << out.str() << std::endl;
			exit(1);
		}
		for (size_t i = 0; i < element.second.size(); ++i) {
			int ret;
			struct voxel voxel =
			    voxel_of_point(element.second[i], voxel_size);
			if (last_free_voxels.find(voxel) != last_free_voxels.end() ||
			    (last_half_voxels.find(voxel) != last_half_voxels.end() &&
			     last_half_voxels[voxel].find(element.first) !=
				 last_half_voxels[voxel].end())) {
				ret = fprintf(out_mask, "1\n");
			} else {
				ret = fprintf(out_mask, "0\n");
			}
			if (ret < 0) {
				std::cerr << "failed to write to " << out.str()
					  << std::endl;
				exit(1);
			}
		}
		fclose(out_mask);
		std::cerr << ((done + 1) * 100.0f / points_by_slice.size())
			  << " %\r";
		std::cerr.flush();
		done += 1;
	}
	std::cerr << std::endl;

#ifndef _MSC_VER
	clock_gettime(CLOCK_MONOTONIC, &after);
	elapsed = (after.tv_sec - before.tv_sec);
	elapsed += (after.tv_nsec - before.tv_nsec) / 1000000000.0;
	std::cerr << "took: " << elapsed << " seconds" << std::endl;
#endif

	std::cerr << "done" << std::endl;

	for (size_t idx = 0; idx < scanorder.size(); ++idx) {
		size_t i = scanorder[idx];
#ifdef WITH_MMAP_SCAN
		if (!cachedir.empty()) {
			int ret;
			ret = munmap(orig_points_by_slice[i].get_raw_pointer(),
				     orig_points_by_slice[i].size());
			if (ret != 0) {
				throw std::runtime_error("cannot munmap");
			}
			// since we called unlink() before, this also deletes
			// the file for good
			ret = close(mmap_fds[i]);
			if (ret != 0) {
				throw std::runtime_error("cannot close");
			}
		} else {
#endif
			delete orig_points_by_slice[i].get_raw_pointer();
#ifdef WITH_MMAP_SCAN
		}
#endif
	}

	Scan::closeDirectory();

	return 0;
}