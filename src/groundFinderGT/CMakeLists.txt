# CMakeLists.txt for groundFinderGT

add_executable(getNormalsGT getNormalsGT.cpp)

target_link_libraries(getNormalsGT 
    scan 
    scanclient
    scanio
    ${ANN_LIBRARIES} 
    ${NEWMAT_LIBRARIES}
    ${SUITESPARSE_LIBRARIES}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${Boost_SYSTEM_LIBRARY} 
    ${Boost_FILESYSTEM_LIBRARY}
)

if(WIN32)
    target_link_libraries(getNormalsGT XGetopt)
endif()

install(TARGETS getNormalsGT
    RUNTIME DESTINATION bin
)
