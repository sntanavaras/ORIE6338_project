cmake_minimum_required(VERSION 3.3)

project(EA_heuristic)

if(TARGET SCIP::SCIP)
  # find package by SCIP PATH
  find_package(SCIP CONFIG PATHS ${SCIP_BINARY_DIR} REQUIRED)
else()
  find_package(SCIP REQUIRED)
endif()

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(scip_ea_heur
  src/heur_EA.cpp
  src/cppmain.cpp)

target_link_libraries(scip_ea_heur ${SCIP_LIBRARIES})

if( TARGET examples )
      add_dependencies( examples scip_ea_heur )
endif()
