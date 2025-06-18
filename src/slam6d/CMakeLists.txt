## e57 to scan tool
if(WITH_E57)
  message(STATUS "With e57 Support.")
  add_definitions(-DWITH_E57)
  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework CoreServices -framework CoreFoundation")
  endif()
  include_directories(${LIBE57_INCLUDE_DIRS})
  add_executable(e572scan e572scan.cpp)
  target_link_libraries(e572scan E57RefImpl ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
else()
  message(STATUS "Whitout e57 Support.")
endif()

if(WITH_METRICS)
  message(STATUS "With metrics in slam6d.")
  add_definitions(-DWITH_METRICS)
else()
  message(STATUS "Without metrics in slam6d.")
endif()

include_directories(${SUITESPARSE_INCLUDE_DIRS})
include_directories(${ANN_INCLUDE_DIRS})

### TOOLS

### SCAN_RED

if (WITH_OPENCV)
  add_executable(scan_red scan_red.cc ../scanio/writer.cc fbr/fbr_global.cc fbr/panorama.cc fbr/scan_cv.cc)
  target_include_directories(scan_red PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/rply-1.1.4)
  target_link_libraries(scan_red scan show ${ANN_LIBRARIES} fbr_cv_io fbr_panorama ${OpenCV_LIBS} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${SHOW_LIBS} rply spherical_quadtree)
endif()

### SCAN_DIFF
add_executable(scan_diff scan_diff.cc)
target_link_libraries(scan_diff scan ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
if (WIN32)
  target_link_libraries(scan_diff XGetopt)
endif()

add_executable(frame_to_graph frame_to_graph.cc)
add_executable(convergence convergence.cc)
add_executable(graph_balancer graph_balancer.cc)
add_executable(frames2riegl frames2riegl.cc)
add_executable(frames2tum frames2tum.cc)
add_executable(frames2pose frames2pose.cc)
add_executable(kitti2scan kitti2scan.cc)
add_executable(kitti2pose kitti2pose.cc)
add_executable(frames2kitti frames2kitti.cc)
add_executable(framesdiff2frames framesdiff2frames.cc)
add_executable(pose2frames pose2frames.cc)
add_executable(riegl2frames riegl2frames.cc)
add_executable(transformFrames transformFrames.cc)
add_executable(multFrames multFrames.cc)
add_executable(trajectoryLength trajectoryLength.cc)


#  add_executable(vigo23dtk vigo23dtk.cc)
#  add_executable(g23dtk g23dtk.cc)
add_executable(toGlobal toGlobal.cc)
add_executable(average6DoFposes average6DoFposes.cc)
add_executable(align sICP.cc)
add_executable(scan2scan_distance scan2scan_distance.cc)

if (WITH_OPENCV)
  add_executable(exportPoints exportPoints.cc ../scanio/writer.cc ../scanio/framesreader.cc)
  target_include_directories(exportPoints PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/rply-1.1.4)
  target_link_libraries(exportPoints scan ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${OpenCV_LIBS} rply)

  
  add_executable(condense condense.cc)
  target_include_directories(condense PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/rply-1.1.4)
  target_link_libraries(condense scan scanio ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${OpenCV_LIBS} rply)
  add_executable(atomize atomize.cc)
  target_include_directories(atomize PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/rply-1.1.4)
  target_link_libraries(atomize scan scanio ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${OpenCV_LIBS} rply)
endif()

if (WITH_SYSTEMC)
  add_executable(icpFixpoint icpFixpoint.cc ../scanio/writer.cc ../scanio/framesreader.cc ../sc_fixed/sc_ICP.cc ../sc_fixed/sc_ICPapx.cc)
  target_include_directories(icpFixpoint PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/rply-1.1.4)
  target_link_libraries(icpFixpoint scan ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${OpenCV_LIBS} ${SYSTEMC_LIBRARY} rply)
endif()

target_link_libraries(graph_balancer scan ${Boost_GRAPH_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} ${Boost_REGEX_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(align scan ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
target_link_libraries(transformFrames scan ${ANN_LIBRARIES} newmat ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(multFrames scan ${ANN_LIBRARIES} newmat ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(toGlobal scan)
target_link_libraries(convergence ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(pose2frames ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(frames2pose ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(frames2riegl ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(riegl2frames ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(scan2scan_distance scan ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(trajectoryLength ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY})

if (WIN32)
  target_link_libraries(frame_to_graph XGetopt)
  target_link_libraries(convergence XGetopt)
  target_link_libraries(graph_balancer XGetopt)
  target_link_libraries(kitti2scan XGetopt)
  target_link_libraries(kitti2pose XGetopt)
  target_link_libraries(frames2kitti XGetopt)
  target_link_libraries(framesdiff2frames XGetopt)
  target_link_libraries(toGlobal scan XGetopt)
endif()

target_link_libraries(average6DoFposes ${NEWMAT_LIBRARIES})

### SCANLIB

add_library(scan
        kd.cc             kdManaged.cc      kdMeta.cc         graphSlam6D.cc
        graph.cc          icp6D.cc          icp6Dapx.cc       icp6Dsvd.cc
        icp6Dortho.cc     icp6Dquat.cc      icp6Dhelix.cc     icp6Dlumeuler.cc
        icp6Dlumquat.cc   icp6Ddual.cc      lum6Deuler.cc     lum6Dquat.cc
        ghelix6DQ2.cc     gapx6D.cc         ann_kd.cc         elch6D.cc
        elch6Dquat.cc     elch6DunitQuat.cc elch6Dslerp.cc    elch6Deuler.cc
        point_type.cc     icp6Dquatscale.cc searchTree.cc     Boctree.cc
        scan.cc           basicScan.cc      managedScan.cc    metaScan.cc
        io_types.cc       io_utils.cc       pointfilter.cc    allocator.cc
        icp6Dnapx.cc      normals.cc        kdIndexed.cc      ../parsers/range_set_parser.cc
        bkd.cc            bkdIndexed.cc     BruteForceNotATree.cc
        )
set_property(TARGET scan PROPERTY POSITION_INDEPENDENT_CODE 1)
target_link_libraries(scan scanclient scanio ${ANN_LIBRARIES} ${NEWMAT_LIBRARIES} ${SUITESPARSE_LIBRARIES})

### SCAN2FEATURES

add_executable(scan2features scan2features.cc)

target_link_libraries(scan2features scan ${NEWMAT_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

### SLAM6D

add_executable(slam6D slam6D.cc)

target_link_libraries(slam6D scan ${NEWMAT_LIBRARIES} ${SUITESPARSE_LIBRARIES} ${ANN_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

#if(MSVC)
#	install(TARGETS slam6D RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/windows)
#endif()
