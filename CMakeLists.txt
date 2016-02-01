# CMAKE_COMPILER
# - Basic CMake setup
# Check version meets ou requirements
# Declare project, which will configure compiler for us
cmake_minimum_required(VERSION 2.8)
project(FalaisePlotModulePlugin)

# Use C++14
set(CMAKE_CXX_FLAGS "-W -Wall -std=c++14")

# Ensure our code can see the Falaise headers
#include_directories(${Falaise_INCLUDE_DIRS})
include_directories(${Falaise_BUILDPRODUCT_DIR}/include)
include_directories(${Falaise_BUILDPRODUCT_DIR}/include/falaise)
include_directories(${Bayeux_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../ParticleIdentification/source)

###########################################################################################
# - PlotModule modules:
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/source/falaise)

# - Headers:
list(APPEND FalaisePlotModulePlugin_HEADERS
  # source/falaise/snemo/analysis/control_plot_module.h
  source/falaise/snemo/analysis/vertices_plot_module.h
  source/falaise/snemo/analysis/halflife_limit_module.h
  # source/falaise/snemo/analysis/snemo_bfield_1e_module.h
  source/falaise/snemo/analysis/universal_plot_module.h
  )

# - Sources:
list(APPEND FalaisePlotModulePlugin_SOURCES
  # source/falaise/snemo/analysis/control_plot_module.cc
  source/falaise/snemo/analysis/vertices_plot_module.cc
  source/falaise/snemo/analysis/halflife_limit_module.cc
  source/falaise/snemo/analysis/halflife_limit_module.cc
  # source/falaise/snemo/analysis/snemo_bfield_1e_module.cc
  source/falaise/snemo/analysis/universal_plot_module.cc
  )

###########################################################################################

# Build a dynamic library from our sources
add_library(Falaise_PlotModule SHARED
  ${FalaisePlotModulePlugin_HEADERS}
  ${FalaisePlotModulePlugin_SOURCES})

target_link_libraries(Falaise_PlotModule Falaise Falaise_ParticleIdentification)

# Apple linker requires dynamic lookup of symbols, so we
# add link flags on this platform
if(APPLE)
  set_target_properties(Falaise_PlotModule
    PROPERTIES LINK_FLAGS "-undefined dynamic_lookup"
    )
endif()

# Install it:
install(TARGETS Falaise_PlotModule DESTINATION ${CMAKE_INSTALL_LIBDIR}/Falaise/modules)

# # Test support:
# option(FalaisePlotModulePlugin_ENABLE_TESTING "Build unit testing system for FalaisePlotModule" ON)
# if(FalaisePlotModulePlugin_ENABLE_TESTING)
#   enable_testing()
#   add_subdirectory(testing)
# endif()
