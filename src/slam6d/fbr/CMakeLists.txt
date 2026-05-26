cmake_minimum_required(VERSION 3.14)

if (WITH_OPENCV)
  find_package(OpenCV REQUIRED)
  set(FBR_LIBS scan show ${OpenCV_LIBS})

  if(WITH_OPENCV_NONFREE)
    message(STATUS "Using non-free OpenCV functionality")
    add_definitions(-DWITH_OPENCV_NONFREE)
  else()
    message(STATUS "Not using non-free OpenCV functionality")
  endif()

  ### Export Static Libs
  set(FBR_IO_SRC fbr_global.cc scan_cv.cc)
  add_library(fbr_cv_io ${FBR_IO_SRC})
  target_link_libraries(fbr_cv_io ${FBR_LIBS})

  set(FBR_PANORAMA_SRC panorama.cc projection.cc)
  add_library(fbr_panorama ${FBR_PANORAMA_SRC} fbr_global.cc)
  target_link_libraries(fbr_panorama ${FBR_LIBS})

  if(OpenCV_VERSION VERSION_GREATER 2.1.0)

    set(FBR_FEATURE_SRC feature.cc)
    add_library(fbr_feature ${FBR_FEATURE_SRC})
    target_link_libraries(fbr_feature fbr_cv_io ${FBR_LIBS})

    set(FBR_FEATURE_MATCHER_SRC feature_matcher.cc feature.cc)
    add_library(fbr_feature_matcher ${FBR_FEATURE_MATCHER_SRC})
    target_link_libraries(fbr_feature_matcher fbr_cv_io ${FBR_LIBS})

    set(FBR_REGISTRATION_SRC registration.cc)
    add_library(fbr_registration ${FBR_REGISTRATION_SRC})
    target_link_libraries(fbr_registration fbr_cv_io ${FBR_LIBS})

    set(FBR_GEOREFRENCING_SRC geoRefrencing.cc)
    add_library(fbr_geoRefrencing ${FBR_GEOREFRENCING_SRC} fbr_global.cc)
    target_link_libraries(fbr_geoRefrencing ${FBR_LIBS})


    add_library(fbr scan_cv.cc panorama.cc projection.cc feature.cc feature_matcher.cc registration.cc fbr_global.cc feature_drawer.cc geoRefrencing.cc fbr_global.cc)
     

  add_library(scan_to_panorama_lib STATIC scan_to_panorama.cc) 
  target_link_libraries(scan_to_panorama_lib fbr ${FBR_LIBS} ) 
    target_link_libraries(fbr ${FBR_LIBS} scan ${ANN_LIBRARIES} ${OpenCV_LIBS} show)


    ###feature based registration
    add_executable(featurebasedregistration feature_based_registration.cc fbr_global.cc)
    #target_link_libraries(featurebasedregistration fbr_cv_io fbr_panorama fbr_feature fbr_feature_matcher fbr_registration ${FBR_LIBS})

    target_link_libraries(featurebasedregistration fbr ${FBR_LIBS} ${NEWMAT_LIBRARIES} show ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
    if (WIN32)
      target_link_libraries(featurebasedregistration XGetopt)
    endif()

    ###panorama matcher
    add_executable(panorama_matcher panorama_matcher.cc fbr_global.cc)

    target_link_libraries(panorama_matcher fbr ${Boost_SYSTEM_LIBRARY} ${SHOW_LIBS})

    if (WIN32)
      target_link_libraries(panorama_matcher XGetopt)
    endif()


    ###scan to panorama
    add_executable(scan_to_panorama scan_to_panorama.cc fbr_global.cc)

    target_link_libraries(scan_to_panorama fbr ${FBR_LIBS} ${NEWMAT_LIBRARIES} show ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

    if (WIN32)
	    target_link_libraries(scan_to_panorama XGetopt)
    endif()

    ###click tool
    add_executable(click_tool click_tool.cc fbr_global.cc)
    target_link_libraries(click_tool fbr ${NEWMAT_LIBRARIES} ${FBR_LIBS} show ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

    if (WIN32)
	    target_link_libraries(click_tool XGetopt)
    endif()

    ###panorama reader
    add_executable(panorama_to_scan panorama_to_scan.cc fbr_global.cc)

    target_link_libraries(panorama_to_scan fbr ${NEWMAT_LIBRARIES} show ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})

    if (WIN32)
	    target_link_libraries(panorama_to_scan XGetopt)
    endif()


  else()
    message("OpenCV Version > 2.2 required for FBR")
  endif()
endif()

  add_subdirectory(new_click_tool/MVP)


####TEST
#add_executable(fbr-icp fbr-icp.cc fbr_global.cc)

#if(UNIX)
#target_link_libraries(fbr-icp fbr newmat showstatic ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
#target_link_libraries(fbr-icp fbr newmat sparse showstatic ${SHOW_LIBS} ${OPENGL_LIBRARIES} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
#endif()

