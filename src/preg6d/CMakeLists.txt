include_directories(${PROJECT_SOURCE_DIR}/include/preg6d)

include_directories(${SUITESPARSE_INCLUDE_DIRS})

include_directories(${NEWMAT_INCLUDE_DIRS})

include_directories(${ANN_INCLUDE_DIRS})

add_definitions(-DWITH_LIBZIP)

set(PREG6D_SRCS
planereg.cc
io/ioplanes.cc
model/planescan.cc
model/normalplane.cc
model/rg.cc
model/detectplanes.cc
match/matcher.cc
match/simplematchers.cc
match/clustermatcher.cc
match/planematcher.cc
opt/optimizer.cc
opt/adadelta6d.cc
opt/newtons6d.cc
opt/gaussnewton.cc
opt/svd.cc
tree/kdplanes.cc
tree/bkdplanes.cc
../slam6d/searchTree.cc
../scanio/framesreader.cc
)

set(POLYREG6D_SRCS
polyreg.cc
io/ioplanes.cc
model/planescan.cc
model/normalplane.cc
model/rg.cc
model/detectplanes.cc
match/matcher.cc
match/simplematchers.cc
match/clustermatcher.cc
match/planematcher.cc
opt/optimizer.cc
opt/adadelta6d.cc
opt/newtons6d.cc
opt/gaussnewton.cc
opt/svd.cc
tree/kdplanes.cc
tree/bkdplanes.cc
../slam6d/searchTree.cc
../scanio/framesreader.cc
)

set(CLUSTER_SRCS
cluster.cc
model/planescan.cc
model/normalplane.cc
model/rg.cc
model/detectplanes.cc
match/matcher.cc
match/simplematchers.cc
match/clustermatcher.cc
match/planematcher.cc
tree/kdplanes.cc
tree/bkdplanes.cc
../slam6d/searchTree.cc
)

add_executable(preg6D ${PREG6D_SRCS})
add_executable(polyreg ${POLYREG6D_SRCS})
add_executable(cluster ${CLUSTER_SRCS})

if(WIN32)
  target_link_libraries(preg6D shape scan scanio ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${ANN_LIBRARIES})
  target_link_libraries(polyreg shape scan scanio ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}${ANN_LIBRARIES})
  target_link_libraries(cluster shape scan scanio ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY}${ANN_LIBRARIES})
  else()
  target_link_libraries(preg6D shape scan scanio ${NEWMAT_LIBRARIES_SHARED} ${CXSPARSE_LIBRARIES} ${ANN_LIBRARIES_SHARED} dl ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${ANN_LIBRARIES})
  target_link_libraries(polyreg shape scan scanio ${NEWMAT_LIBRARIES_SHARED} ${CXSPARSE_LIBRARIES} ${ANN_LIBRARIES_SHARED} dl ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${ANN_LIBRARIES})
  target_link_libraries(cluster shape scan scanio ${NEWMAT_LIBRARIES_SHARED} ${CXSPARSE_LIBRARIES} ${ANN_LIBRARIES_SHARED} dl ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${ANN_LIBRARIES})
endif()
