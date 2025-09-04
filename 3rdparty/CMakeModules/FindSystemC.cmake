# FindSystemC.cmake

find_path(SYSTEMC_INCLUDE_DIR
	NAMES systemc.h
	PATHS /usr/include
)

find_library(SYSTEMC_LIBRARY
	NAMES systemc
	PATHS /usr/lib /usr/lib/x86_64-linux-gnu
	
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SystemC DEFAULT_MSG
	SYSTEMC_LIBRARY	SYSTEMC_INCLUDE_DIR
)

mark_as_advanced(SYSTEMC_INCLUDE_DIR SYSTEMC_LIBRARY)

if(SystemC_FOUND)
  set(SystemC_INCLUDE_DIRS ${SYSTEMC_INCLUDE_DIR})
  set(SystemC_LIBRARIES ${SYSTEMC_LIBRARY})

  if(NOT TARGET SystemC::systemc)
    add_library(SystemC::systemc UNKNOWN IMPORTED)
    set_target_properties(SystemC::systemc PROPERTIES
      IMPORTED_LOCATION "${SYSTEMC_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${SYSTEMC_INCLUDE_DIR}"
    )
  endif()
endif()
