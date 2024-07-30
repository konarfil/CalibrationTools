project(sncalib)

find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

add_compile_definitions(Falaise_RESOURCE_DIR="${Falaise_RESOURCE_DIR}")
add_compile_definitions(FITTED_SPECTRA_DIR="${PROJECT_SOURCE_DIR}/build/Fits")

add_executable(sncalib src/sncalib.cc)

target_include_directories(sncalib PUBLIC "${PROJECT_SOURCE_DIR}/include")

target_sources(sncalib PRIVATE src/calibration_paremeter_finder.cc)

target_link_libraries(sncalib PUBLIC Falaise::Falaise ${ROOT_LIBRARIES})

configure_file("params.conf.in" "params.conf" @ONLY)