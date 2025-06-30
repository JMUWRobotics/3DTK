if (WITH_OPENGL AND NOT APPLE)
  find_package(GLUT REQUIRED)
  find_package(OpenGL REQUIRED)
  # Compile glui library
  SET(GLUI_DIR ${PROJECT_SOURCE_DIR}/3rdparty/glui-2.36/)
  add_library(glui
   ${GLUI_DIR}algebra3.cpp
   ${GLUI_DIR}arcball.cpp
   ${GLUI_DIR}glui_add_controls.cpp
   ${GLUI_DIR}glui_bitmap_img_data.cpp
   ${GLUI_DIR}glui_bitmaps.cpp
   ${GLUI_DIR}glui_button.cpp
   ${GLUI_DIR}glui_checkbox.cpp
   ${GLUI_DIR}glui_column.cpp
   ${GLUI_DIR}glui_commandline.cpp
   ${GLUI_DIR}glui_control.cpp
   ${GLUI_DIR}glui.cpp
   ${GLUI_DIR}glui_edittext.cpp
   ${GLUI_DIR}glui_filebrowser.cpp
   ${GLUI_DIR}glui_listbox.cpp
   ${GLUI_DIR}glui_list.cpp
   ${GLUI_DIR}glui_mouse_iaction.cpp
   ${GLUI_DIR}glui_node.cpp
   ${GLUI_DIR}glui_panel.cpp
   ${GLUI_DIR}glui_radio.cpp
   ${GLUI_DIR}glui_rollout.cpp
   ${GLUI_DIR}glui_rotation.cpp
   ${GLUI_DIR}glui_scrollbar.cpp
   ${GLUI_DIR}glui_separator.cpp
   ${GLUI_DIR}glui_spinner.cpp
   ${GLUI_DIR}glui_statictext.cpp
   ${GLUI_DIR}glui_string.cpp
   ${GLUI_DIR}glui_textbox.cpp
   ${GLUI_DIR}glui_translation.cpp
   ${GLUI_DIR}glui_tree.cpp
   ${GLUI_DIR}glui_treepanel.cpp
   ${GLUI_DIR}glui_window.cpp
   ${GLUI_DIR}quaternion.cpp
   ${GLUI_DIR}viewmodel.cpp
  )
  target_include_directories(glui PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/glui-2.36/include/)
  target_include_directories(glui PUBLIC ${GLUT_INCLUDE_DIR})
  target_compile_definitions(glui PUBLIC "-DGLUI_NO_LIB_PRAGMA")

  target_link_libraries(glui ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES})
endif()

# Compile wxthings library
if (WITH_WXWIDGETS)
  find_package(wxWidgets REQUIRED core base gl)
  # wxWidgets include (this will do all the magic to configure everything)
  include( ${wxWidgets_USE_FILE})
  SET(WXTHINGS_DIR ${PROJECT_SOURCE_DIR}/3rdparty/wxthings/src/)
  add_library(wxthings STATIC 
    ${WXTHINGS_DIR}block.cpp
    ${WXTHINGS_DIR}bmpcombo.cpp
    ${WXTHINGS_DIR}dropdown.cpp
    # ${WXTHINGS_DIR}filebrws.cpp  # these two are not needed and would create unnecessary dependencies
    #${WXTHINGS_DIR}genergdi.cpp
    ${WXTHINGS_DIR}geometry.cpp
    ${WXTHINGS_DIR}menubtn.cpp
    ${WXTHINGS_DIR}optvalue.cpp
    ${WXTHINGS_DIR}range.cpp
    ${WXTHINGS_DIR}spinctld.cpp
    ${WXTHINGS_DIR}toggle.cpp
  )
  target_include_directories(wxthings PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/wxthings/include)
endif()

if (WITH_LASLIB)
# Compile las library
SET(LASLIB_DIR ${PROJECT_SOURCE_DIR}/3rdparty/lastools/laslib/src/)
FILE(GLOB LASLIB_SRCS "${LASLIB_DIR}/*.cpp" "${LASLIB_DIR}/*.hpp")
add_library(laslib ${LASLIB_SRCS})
target_compile_definitions(laslib PUBLIC "-DNDEBUG -DUNORDERED")
set_target_properties(laslib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(laslib PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/lastools/laslib/inc)

SET(LASSRC_DIR ${PROJECT_SOURCE_DIR}/3rdparty/lastools/src/)
foreach(prog lasinfo lasindex las2txt lasdiff)
  add_executable(${prog} "${LASSRC_DIR}/${prog}.cpp")
endforeach()
foreach(prog lasprecision txt2las las2las lasmerge laszip)
  add_executable(${prog} "${LASSRC_DIR}/${prog}.cpp" "${LASSRC_DIR}/geoprojectionconverter.cpp")
endforeach()
foreach(prog lasinfo lasindex las2txt lasdiff lasprecision txt2las las2las lasmerge laszip)
  target_compile_definitions(${prog} PUBLIC "-DNDEBUG -DUNORDERED")
  target_include_directories(${prog} PUBLIC ${PROJECT_SOURCE_DIR}/3rdparty/lastools/laslib/inc)
  target_link_libraries(${prog} laslib)
endforeach()

endif(WITH_LASLIB)

if (WITH_GLEE)
# Compile GLEE library 
IF (APPLE)
  FIND_LIBRARY(COREFOUNDATION_LIBRARY CoreFoundation )
ENDIF (APPLE)

SET(GLEE_DIR ${PROJECT_SOURCE_DIR}/3rdparty/glee/)

add_library(glee ${GLEE_DIR}GLee.c)

IF(APPLE)
  target_link_libraries(glee ${COREFOUNDATION_LIBRARY} ${OPENGL_LIBRARIES})
ENDIF(APPLE)

endif()

# Compile AprilTags library if required
if (WITH_SYSTEM_APRILTAG)
  find_package(apriltag REQUIRED)
  message(STATUS "Found system apriltag: ${apriltag_VERSION}")
  set(APRILTAG_INCLUDE_DIRS /usr/include/apriltag)
else()
  set(APRILTAG_DIR ${PROJECT_SOURCE_DIR}/3rdparty/apriltag-3.3.0/)
  file(GLOB TAG_FILES ${APRILTAG_DIR}/tag*.c)
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
      # otherwise clang takes about an hour to compile
      set_source_files_properties(SOURCE ${TAG_FILES} PROPERTIES COMPILE_FLAGS -O0)
  endif()
  set(APRILTAG_SRCS
    ${APRILTAG_DIR}/apriltag.c
    ${APRILTAG_DIR}/apriltag_pose.c
    ${APRILTAG_DIR}/apriltag_quad_thresh.c)
  aux_source_directory(${APRILTAG_DIR}/common COMMON_SRC)
  add_library(apriltag SHARED ${APRILTAG_SRCS} ${COMMON_SRC} ${TAG_FILES})
  target_include_directories(apriltag PRIVATE ${APRILTAG_DIR})
  if (MSVC)
      find_package(pthreads REQUIRED)
      target_include_directories(apriltag PUBLIC ${PThreads4W_INCLUDE_DIR})
      target_link_libraries(apriltag PThreads4W::PThreads4W winmm)
  else()
      find_package(Threads REQUIRED)
      target_link_libraries(apriltag Threads::Threads m)
  endif()
  set(APRILTAG_INCLUDE_DIRS ${APRILTAG_DIR} CACHE INTERNAL "")
endif()

add_subdirectory(rply-1.1.4)

if(WITH_CALIB)
  if(NOT EXISTS ${PROJECT_SOURCE_DIR}/3rdparty/cctag/CMakeLists.txt OR NOT EXISTS ${PROJECT_SOURCE_DIR}/3rdparty/paruco/CMakeLists.txt)
    message(FATAL_ERROR "CCTag submodule is missing. Please run `git submodule update --init --recursive` in the project root.")
  endif()

  set(CCTAG_SERIALIZE OFF)
  set(CCTAG_VISUAL_DEBUG OFF)
  set(CCTAG_NO_COUT ON)
  set(CCTAG_WITH_CUDA OFF)
  set(CCTAG_BUILD_APPS OFF)
  set(CCTAG_NVCC_WARNINGS OFF)
  set(CCTAG_EIGEN_MEMORY_ALIGNMENT OFF)
  set(CCTAG_USE_POSITION_INDEPENDENT_CODE ON)
  set(CCTAG_ENABLE_SIMD_AVX2 OFF)
  set(CCTAG_BUILD_TESTS OFF)
  set(CCTAG_BUILD_DOC OFF)
  set(BUILD_SHARED_LIBS OFF)
  add_subdirectory(cctag EXCLUDE_FROM_ALL)

  add_subdirectory(paruco EXCLUDE_FROM_ALL)
endif()

IF(WITH_E57)
  add_subdirectory(e57-3d-imgfmt)
ENDIF()
