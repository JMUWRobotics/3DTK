if (WITH_OPENCV AND WITH_EIGEN3)
  include_directories(${APRILTAG_INCLUDE_DIRS})
  include_directories(${CCTag_INCLUDE_DIRS})
  include_directories(${PArUco_INCLUDE_DIRS})

  find_package(Eigen3 REQUIRED)
  include_directories(${EIGEN3_INCLUDE_DIRS})

  set(CALIBRATION_LIBRARY_SOURCES
      AprilTag.cc
      Chessboard.cc
      AprilTagDetector.cc
      ArucoDetector.cc
      ChessboardDetector.cc
      CirclesGridDetector.cc
      PinholeModel.cc
      DetectionFileHandler.cpp
      CCTagDetector.cc
      PArUcoDetector.cc
      )
  add_library(calibration ${CALIBRATION_LIBRARY_SOURCES})
  target_link_libraries(calibration apriltag CCTag paruco ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  add_executable(detect detect.cc )
  if (WIN32)
    target_link_libraries(detect winmm)
  endif()
  target_link_libraries(detect calibration)

  set(CAMCAL_SOURCES
      CalibrationPattern.cc
      PictureHandler.cc
      CalibrationToolbox.cc)
  add_executable(camcal camcal.cc ${CAMCAL_SOURCES})
  if (WIN32)
    target_link_libraries(camcal winmm)
  endif()
  target_link_libraries(camcal calibration apriltag ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  add_executable(undistortImage undistortImage.cc)
  target_link_libraries(undistortImage ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  add_executable(undistortPoints undistortPoints.cc)
  target_link_libraries(undistortPoints ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  add_executable(planeFromExtrinsics planeFromExtrinsics.cc)
  target_link_libraries(planeFromExtrinsics ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  add_executable(transformsFromExtrinsics transformsFromExtrinsics.cc)
  target_link_libraries(transformsFromExtrinsics ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  set(CALIBRATE_SOURCES
          Calibrator.cc)
  add_executable(calibrate calibrate.cc ${CALIBRATE_SOURCES})
  target_link_libraries(calibrate calibration ${OpenCV_LIBRARIES} ${Boost_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

  configure_file(createAprilGrid.py ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/createAprilGrid COPYONLY)
endif()
