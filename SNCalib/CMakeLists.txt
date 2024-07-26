# Check cmake version meets our requirements
cmake_minimum_required(VERSION 3.3)

# Declare project, which will configure compiler for us
project(sncalib)

#Finds ROOT and Falaise packages and defines cmake variables
find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

include_directories(include ${Bayeux_INCLUDE_DIRS})

add_compile_definitions(Falaise_RESOURCE_DIR="${Falaise_RESOURCE_DIR}")
add_compile_definitions(FITTED_SPECTRA_DIR="${PROJECT_SOURCE_DIR}/build/Fits")

add_executable(sncalib src/sncalib.cc)

target_include_directories(sncalib PUBLIC "${PROJECT_SOURCE_DIR}/include")

target_sources(sncalib PRIVATE src/calibration_paremeter_finder.cc)

target_link_libraries(sncalib PUBLIC Falaise::Falaise ${ROOT_LIBRARIES})

configure_file("params.conf.in" "params.conf" @ONLY)