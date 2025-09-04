cmake_minimum_required (VERSION 3.5)
project (3DTK)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/3rdparty/CMakeModules" ${CMAKE_MODULE_PATH})

# cmake no longer defines WIN32 on cygwin
set(CMAKE_LEGACY_CYGWIN_WIN32 0) # remove when cmake >= 2.8.4 is required

# On Windows, the symbols of a dynamic library have to be explicitly exported
# using __declspec(dllexport) or otherwise the library will not even be built.
# To avoid having to modify our headers, we just force the same behaviour as
# under Unix
if (MSVC)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

# Hide CMake variables
set (CMAKE_INSTALL_PREFIX "/usr/local" CACHE INTERNAL "" FORCE)
set (CMAKE_BUILD_TYPE "" CACHE INTERNAL "" FORCE)

# being able to set the output directory to a different one than the default
# (the source directory) is important for platforms like windows, where the
# produced binaries cannot be run if being placed in certain locations (like on
# a network drive)
set(OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}" CACHE PATH "The directory where the output will be placed into lib, obj and bin subdirectories (default: ${PROJECT_SOURCE_DIR})" )

# Set output directories for libraries and executables
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}/lib )
set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}/obj )
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${OUTPUT_DIRECTORY}/bin )
# Set output directories for multi-config builds (like with MSVC)
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIRECTORY}/lib )
  set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIRECTORY}/obj )
  set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${OUTPUT_DIRECTORY}/bin )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

include(CheckSymbolExists)
check_symbol_exists(mmap sys/mman.h HAVE_MMAP)
check_symbol_exists(mkstemp stdlib.h HAVE_MKSTEMP)
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
check_symbol_exists(fallocate fcntl.h HAVE_FALLOCATE)
set(CMAKE_REQUIRED_DEFINITIONS)
if(HAVE_MMAP AND HAVE_MKSTEMP AND HAVE_FALLOCATE)
  add_definitions(-DWITH_MMAP_SCAN)
endif()

#include_directories(OPENGL_INCLUDE_DIR)
if(WIN32)
  # Tells the config system not to automatically select  which libraries to
  # link against. Normally if a compiler supports #pragma lib, then the
  # correct library build variant will be automatically selected and linked
  # against, simply by the act of including one of that library's headers.
  # This macro turns that feature off.
  add_definitions(-DBOOST_ALL_NO_LIB)
  # Forces all libraries that have separate source, to be linked as dll's
  # rather than static libraries on Microsoft Windows (this macro is used to
  # turn on __declspec(dllimport) modifiers, so that the compiler knows which
  # symbols to look for in a dll rather than in a static library).
  add_definitions(-DBOOST_ALL_DYN_LINK)
endif()

if(MINGW)
  # https://github.com/msys2/MINGW-packages/issues/5233
  set(Boost_NO_BOOST_CMAKE ON)
endif()

set(Boost_ADDITIONAL_VERSIONS "1.42" "1.42.0" "1.44" "1.44.0" "1.45.0" "1.45" "1.46" "1.46.1" "1.47.0" "1.47" "1.48" "1.50" "1.52" "1.53" "1.55" "1.56")
if(WIN32)
  # for some unknown reason no one variant works on all windows platforms
  set(Boost_DEBUG 1)
  find_package( Boost COMPONENTS serialization graph regex filesystem system thread chrono date_time program_options system REQUIRED)
else()
  find_package( Boost COMPONENTS serialization graph regex filesystem system thread date_time program_options system REQUIRED)
endif()

link_directories(${BOOST_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})

find_package(SuiteSparse REQUIRED)

find_package(ANN)
if (NOT ANN_FOUND)
  add_subdirectory("3rdparty/ann")
  set(ANN_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/3rdparty/ann/ann_1.1.1_modified/include")
  set(ANN_LIBRARIES ann)
endif()

find_package(Newmat)
if(NOT NEWMAT_FOUND)
  add_subdirectory("3rdparty/newmat")
  set(NEWMAT_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/3rdparty/newmat/newmat-10")
  set(NEWMAT_LIBRARIES newmat)
endif()
include_directories(${NEWMAT_INCLUDE_DIRS})

if(APPLE)
  find_package(OpenCV REQUIRED)
else()
  if(WIN32)
    set(OpenCV_STATIC ON)
  endif()
endif()

if(EXISTS "${OpenCV_DIR}/OpenCVConfig.cmake")
  include("${OpenCV_DIR}/OpenCVConfig.cmake")
  set(ADDITIONAL_OPENCV_FLAGS
          "-DCV_MINOR_VERSION=${OpenCV_VERSION_MINOR} -DCV_MAJOR_VERSION=${OpenCV_VERSION_MAJOR}"
          CACHE STRING "OpenCV Version Defines)"
          )
  ## Include the standard CMake script
else()
  set(ADDITIONAL_OPENCV_FLAGS
          ""
          CACHE STRING "OpenCV Version Defines (BLUB)"
          )
endif()
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_OPENCV_FLAGS}")

#################################################
# Declare Options and modify build accordingly ##
#################################################

option(WITH_CGAL "Compile with CGAL support" OFF)
option(WITH_GMP "Compile with GMP support" OFF) # on Windows GMP is provided by MPIR
option(WITH_LIBZIP "Compile with libzip support" ON)
option(WITH_OPENGL "Compile with OpenGL support" ON)
option(WITH_OPENCV "Compile with OpenCV support" ON)
option(WITH_GLFW "Compile with GLFW support" ON)
option(WITH_FTGL "Compile with FTGL support" ON)
option(WITH_XMLRPC "Compile with XMLRPC support" ON)
option(WITH_EIGEN3 "Compile with eigen3 support" ON)
option(WITH_LIBCONFIG "Compile with libconfig support" ON)
option(WITH_ROS "Compile with ROS support" OFF)
option(WITH_PYTHON "Compile Python bindings" OFF)
option(WITH_QT "Compile tools relying on QT (qtshow)" ON)
option(WITH_WXWIDGETS "Compile with wxwidgets support" OFF)
option(WITH_IMGUI "Compile the imshow program" OFF)
option(WITH_OPENCV_NONFREE "Whether to use non-free (patent encumbered) OpenCV functionalities" OFF)
option(WITH_COMPACT_OCTREE "Whether to use the compact octree display ON/OFF" OFF)
option(WITH_GLEE "Whether to use OpenGL extensions, requires glee. ON/OFF" OFF)
option(WITH_LASLIB "Whether to build LASlib based scanio library" ON)
option(WITH_E57 "Build library for e57 file format support." OFF)
option(WITH_3DMOUSE "Support for 3DConnexion Space Mouse" OFF)
option(WITH_TSDF "Complie with TSDF (openVDB)" OFF)
option(WITH_CALIB "Complie with calibration" OFF)
option(WITH_SYSTEM_APRILTAG "Link against system-wide apriltag library" OFF)
option(WITH_SYSTEMC "Compile with SystemC library for fixedPoint arithmetic" ON)

# Adds imgui to show using opengl2 and glut
if(WITH_IMGUI)
  set(IMGUI_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/3rdparty/imgui/" "${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends" "${PROJECT_SOURCE_DIR}/3rdparty/imgui/imGuIZMO.quat")
  file(GLOB IMGUI_SRCS "${PROJECT_SOURCE_DIR}/3rdparty/imgui/*.cpp" "${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_opengl2.cpp" "${PROJECT_SOURCE_DIR}/3rdparty/imgui/backends/imgui_impl_glut.cpp" "${PROJECT_SOURCE_DIR}/3rdparty/imgui/imGuIZMO.quat/imGuIZMOquat.cpp")
endif()

if(WITH_3DMOUSE)
  #add_compile_definitions(SPACEMOUSE)
  if(APPLE)
    #find_library(3DCONNEXIONCLIENT_FRAMEWORK 3DconnexionClient)
    find_package(3DconnexionClient)
    include_directories(${3DCONNEXIONCLIENT_INCLUDE_DIR2})
    set(3DCONNEXION_LIB ${3DCONNEXIONCLIENT})
    set(SPNAV_LIB "")
    set(SPNAV_CONTROLLER "")
    #add_library(spnav STATIC IMPORTED)
    #set_property(TARGET spnav PROPERTY IMPORTED_LOCATION /usr/local/Cellar/libspnav/0.2.3/lib/libspnav.dylib)
    #include_directories(/usr/local/Cellar/libspnav/0.2.3/include)
  else()
    set(3DCONNEXION_LIB "")
    set(SPNAV_LIB "spnav")
    set(SPNAV_CONTROLLER "spnav_controller.cc")
    message(STATUS "With 3DMouse Support")
  endif()
else()
  set(3DCONNEXION_LIB "")
  set(SPNAV_LIB "")
  set(SPNAV_CONTROLLER "")
endif()

# cvblob is needed for thermo
if (WITH_OPENCV)
  find_package(OpenCV REQUIRED)
  if (OpenCV_VERSION VERSION_LESS 4)
  message(STATUS "building cvblob")
  include_directories(${OpenCV_INCLUDE_DIRS})
  add_subdirectory(3rdparty/cvblob)
  link_directories(${PROJECT_SOURCE_DIR}/3rdparty/cvblob)
  endif()
endif()

## CUDA accelerated collision detection
option(WITH_CUDA "Whether to build CUDA accelerated collision detection tools" OFF)
if(WITH_CUDA)
  find_package(CUDA REQUIRED)
  include_directories(${FOUND_CUDA_NVCC_INCLUDE})

  #owerrite users input - need to be fixed
  #select all NVIDIA GPU arch which support unified memory (CUDA toolkit >= 6.0) and arch>=30
  set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_30,code=sm_30;-gencode arch=compute_35,code=sm_35;-gencode arch=compute_50,code=sm_50;-gencode arch=compute_52,code=sm_52)

  message("CUDA_NVCC_FLAGS = ${CUDA_NVCC_FLAGS}")

  message(STATUS "With CUDA accelerated collision detection")
else()
  message(STATUS "Without CUDA accelerated collision detection")
endif()


## RivLib
option(WITH_RIVLIB "Whether the RIEGL rivlib is present ON/OFF" OFF)

if(WITH_RIVLIB)
  message(STATUS "Compiling a scan IO for RXP files")

  set(RIEGL_DIR ${PROJECT_SOURCE_DIR}/3rdparty/riegl/)

  if(UNIX)
    set(RiVLib_USE_STATIC_RUNTIME ON)
  endif()

  find_package(RiVLib QUIET COMPONENTS scanlib HINTS "${RIEGL_DIR}/cmake/")

  if(${RiVLib_FOUND})
    include_directories(${RiVLib_INCLUDE_DIRS})
  else()
    # TODO: Remove this if nobody is using the old RiVLib anymore.
    #       Change QUIET to REQUIRED in the find_package call for RiVLib above.

    message(STATUS "Cannot find current RiVLib. Trying to build scan IO for RXP files with old scanlib.")

    add_definitions(-DWITH_OLD_RIVLIB)

    find_package(LibXml2 REQUIRED)
    include_directories(${LIBXML2_INCLUDE_DIR})

    include_directories(${PROJECT_SOURCE_DIR}/3rdparty)
    if(WIN32)
      set(RiVLib_SCANLIB_LIBRARY ${RIEGL_DIR}libscanlib-mt.lib ${RIEGL_DIR}libctrllib-mt.lib ${RIEGL_DIR}libboost_system-mt-1_43_0-vns.lib ${LIBXML2_LIBRARIES})
    else()
      set(RiVLib_SCANLIB_LIBRARY ${RIEGL_DIR}libscanlib-mt-s.a ${RIEGL_DIR}libctrllib-mt-s.a ${RIEGL_DIR}libboost_system-mt-s-1_43_0-vns.a pthread ${LIBXML2_LIBRARIES})
    endif()
  endif()
else()
  message(STATUS "Do NOT compile a scan IO for RXP")
endif()

## Micro-Epsilon libraries
option(WITH_MICROEPSILONLIB "Whether the Micro-Epsilon libraries are present ON/OFF" OFF)

if(WITH_MICROEPSILONLIB)
  message(STATUS "Setting variables for Micro-Epsilon libraries")
  set(MICROEPSILON_PATH ${PROJECT_SOURCE_DIR}/3rdparty/microepsilon/)

  if (EXISTS "${MICROEPSILON_PATH}/libllt.so.0")
    if (EXISTS "${MICROEPSILON_PATH}/libmescan.so.0")
      set(MICROEPSILON_LIBRARY ${MICROEPSILON_PATH}libllt.so.0 ${MICROEPSILON_PATH}libmescan.so.0)
    else()
      message(STATUS "Cannot find Micro-Epsilon libraries. Provide libraries or disable WITH_MICROEPSILONLIB.")
    endif()
  else()
    message(STATUS "Cannot find Micro-Epsilon libraries. Provide libraries or disable WITH_MICROEPSILONLIB.")
  endif()
else()
  message(STATUS "Microepsilon libraries not present.")
endif()

option(WITH_OPENMP "Whether to use parallel processing capabilities of OPENMP. ON/OFF" ON)
if (WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  # For CMake < 3.9, we need to make the target ourselves
  if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX
            PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX
            PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)

  endif()
endif()

option(WITH_METRICS "Whether to use time metrics. ON/OFF" OFF)

#################################################
# OPERATING SYSTEM SPECIFIC BEHAVIOUR          ##
#################################################

## Special treatment for system specifics
if(APPLE)
  add_definitions(-Dfopen64=fopen)
endif()

## Multiple Cores
include(ProcessorCount)
ProcessorCount(PROCESSOR_COUNT)
if(NOT PROCESSOR_COUNT EQUAL 0)
  set(NUMBER_OF_CPUS "${PROCESSOR_COUNT}" CACHE STRING "The number of processors to use (default: ${PROCESSOR_COUNT})" )
else()
  set(NUMBER_OF_CPUS "1" CACHE STRING "The number of processors to use (default: 1)" )
endif()

# OPENMP_NUM_THREADS
if(WITH_OPENMP)
  message(STATUS "With OpenMP ")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMAX_OPENMP_NUM_THREADS=${NUMBER_OF_CPUS} -DOPENMP_NUM_THREADS=${NUMBER_OF_CPUS} ${OpenMP_CXX_FLAGS} -DOPENMP")
else()
  message(STATUS "Without OpenMP")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMAX_OPENMP_NUM_THREADS=1 -DOPENMP_NUM_THREADS=1")
endif()

# 3rdparty
if(WIN32)
  include_directories(${PROJECT_SOURCE_DIR}/3rdparty/windows/)
  link_directories(${PROJECT_SOURCE_DIR}/3rdparty/windows)
  link_directories(${BOOST_LIBRARYDIR})
  add_library(XGetopt ${PROJECT_SOURCE_DIR}/3rdparty/windows/XGetopt.cpp)
  set(CMAKE_STATIC_LIBRARY_SUFFIX "32.lib")
  # to be able to use numeric_limits<int>::max() and friends on windows
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNOMINMAX")
else()
  if (WITH_OPENGL)
    # Add include path for OpenGL without GL/-prefix
    # to avoid the include incompatibility between MACOS
    # and linux
    find_path(OPENGL_INC NAMES gl.h GL/gl.h PATHS /usr/include/GL)
    if (${OPENGL_INC})
      include_directories(${OPENGL_INC})
    else()
      message(STATUS "gl.h not found")
    endif()
  endif()
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
  # using Clang
  set(ADDITIONAL_CFLAGS "-O3 -Wall -Wno-write-strings -Wno-char-subscripts -Wno-unused-result" CACHE STRING "Additional flags given to the compiler (-O3 -Wall -Wno-write-strings)" )
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
  # using GCC
  set(ADDITIONAL_CFLAGS "-O3 -Wall -finline-functions -Wno-unused-but-set-variable -Wno-write-strings -Wno-char-subscripts -Wno-unused-result" CACHE STRING "Additional flags given to the compiler (-O3 -Wall -finline-functions -Wno-write-strings)" )
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL Intel)
  # using Intel C++
elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
  # using Visual Studio C++
  set(ADDITIONAL_CFLAGS "-O2" CACHE STRING "Additional flags given to the compiler ( -O2)" )
endif()

# Remove this conditional when CMake >= 3.1.0 is required.
# Older CMake is used in Debian Jessie (3.0.2) and Ubuntu Trusty (2.8.12.2).
#if(${CMAKE_VERSION} VERSION_LESS "3.1.0")
#	if (${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
#		# using Clang
#		set(ADDITIONAL_CFLAGS "${ADDITIONAL_CFLAGS} -std=c++0x")
#	elseif (${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
#		# using GCC
#		set(ADDITIONAL_CFLAGS "${ADDITIONAL_CFLAGS} -std=c++0x")
#	endif()
#else()
	# Specifying the C++ standard (-std=c++0x, -std=gnu++14, ...) by using
	# the CMAKE_CXX_FLAGS is not good cmake style (see Daniel Pfeifer in
	# https://youtu.be/bsXLMQ6WgIk).
	# By using CMAKE_CXX_STANDARD we also avoid having to implement
	# compiler specific command line arguments.
	set(CMAKE_CXX_STANDARD 17)
	# do not fall back to a previous standard
	set(CMAKE_CXX_STANDARD_REQUIRED ON)
	# enable compiler extensions like -std=gnu++11 on gcc
	# gnu extensions are needed for strdup, ftruncate, fseeko, ftello...
	set(CMAKE_CXX_EXTENSIONS ON)
#endif()

# Add OpenGL includes for MACOS if needed
# The OSX OpenGL frameworks natively supports freeglut extensions
if(APPLE)
  include_directories(/System/Library/Frameworks/GLUT.framework/Headers)
  include_directories(/System/Library/Frameworks/OpenGL.framework/Headers)
endif()

# hack to "circumvent" Debug and Release folders that are created under visual studio
# this is why the INSTALL target has to be used in visual studio
if(MSVC)
  install(DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release/ CONFIGURATIONS Release DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/3rdparty/windows/x64/ CONFIGURATIONS Release DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  else()
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/3rdparty/windows/ CONFIGURATIONS Release DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  endif()

  install(DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug/ CONFIGURATIONS Debug DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  if( CMAKE_SIZEOF_VOID_P EQUAL 8 )
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/3rdparty/windows/x64/ CONFIGURATIONS Debug DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  else()
    install(DIRECTORY ${PROJECT_SOURCE_DIR}/3rdparty/windows/ CONFIGURATIONS Debug DESTINATION ${PROJECT_SOURCE_DIR}/windows FILES_MATCHING PATTERN "*.dll" PATTERN "*.exe")
  endif()

endif()

#################################################
# Robot Operating System (ROS) Integration     ##
#################################################

# Export include directories if 3DTK is built within a catkin workspace.
# In that case, find_package(catkin) was already run and catkin_FOUND is true.
if( catkin_FOUND )
  catkin_package(
    INCLUDE_DIRS include addons/include
  )
endif()

#################################################
# GENERAL PROJECT SETTINGS                     ##
#################################################

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CFLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ADDITIONAL_CFLAGS}")

if(APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

# Set include and link dirs ...
include_directories(${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/3rdparty/wxthings/include/)
link_directories(${PROJECT_SOURCE_DIR}/obj)
link_directories(${PROJECT_SOURCE_DIR}/lib)

# src/show must be added first because it defines the cache variables
# SHOW_LIBS_* which are used by others
add_subdirectory(src/show)
add_subdirectory(src/slam6d)
add_subdirectory(src/scanio)
add_subdirectory(src/srr)
add_subdirectory(src/preg6d)
add_subdirectory(src/scanserver)
add_subdirectory(src/segmentation)
add_subdirectory(src/normals)
add_subdirectory(src/veloslam)
add_subdirectory(src/qtshow)
add_subdirectory(src/grid)
add_subdirectory(src/pmd)
add_subdirectory(src/shapes)
add_subdirectory(src/floorplan)
add_subdirectory(src/thermo)
add_subdirectory(src/slam6d/fbr)
add_subdirectory(src/scanner)
add_subdirectory(src/model)
add_subdirectory(src/collision)
add_subdirectory(src/peopleremover)
add_subdirectory(src/spherical_quadtree)
add_subdirectory(src/cuda)
add_subdirectory(src/ros)
add_subdirectory(src/tools)
add_subdirectory(src/gps)
add_subdirectory(src/curvefusion)
add_subdirectory(src/mesh)
add_subdirectory(src/detectCylinder)
add_subdirectory(src/sc_fixed)
# 3rdparty must come before src/calibration because it sets
# APRILTAG_INCLUDE_DIRS
add_subdirectory(3rdparty)
add_subdirectory(bindings)

if(WITH_TSDF)
  message(STATUS "With TSDF based on OpenVDB")
  add_subdirectory(src/tsdf)
endif()

if(WITH_CALIB)
  message(STATUS "With calibration.")
  add_subdirectory(src/calibration)
endif()

find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
enable_testing()
add_subdirectory(testing)


# SystemC integration

# Configurable constants for SystemC fixed-point value
set(FIXED_WORD_LENGTH 48 CACHE STRING "Word length for fixed-point type")
set(FIXED_INT_WORD_LENGTH 36 CACHE STRING "Integer word length for fixed-point type")
set(FIXED_HERON_ITERATIONS 3 CACHE STRING "Iterations ins heron algorithm for square root calculations for fixed point numbers")

if(WITH_SYSTEMC)
  message(STATUS "SystemC support enabled")
  find_package(SystemC REQUIRED)

  if(SystemC_FOUND)
    message(STATUS "SystemC installation found")
    include_directories(${SYSTEMC_INCLUDE_DIR})
    add_definitions(-DWITH_SYSTEMC)
    # Set configurable constants for SystemC fixed-point value as compiler definition
    target_compile_definitions(icpFixpoint PRIVATE
      FIXED_WORD_LENGTH=${FIXED_WORD_LENGTH}
      FIXED_INT_WORD_LENGTH=${FIXED_INT_WORD_LENGTH}
      FIXED_HERON_ITERATIONS=${FIXED_HERON_ITERATIONS}
    )
    set(SYSTEMC_LINK_LIBS ${SYSTEMC_LIBRARY})

  else()
    message(FATAL_ERROR "SystemC was enabled, but could not be found")
  endif()

  set(SYSTEMC_INCLUDE_DIR)

else()
  message(STATUS "SystemC support disabled")
endif()
# Dummy target with all header files
# This is a hint for some IDEs, such as Qt Creator, to show all headers in the project tree
file(GLOB_RECURSE 3DTK_HEADER_FILES "include/*.h")
add_custom_target(headers SOURCES ${3DTK_HEADER_FILES})

message(STATUS "Build environment is set up!")
