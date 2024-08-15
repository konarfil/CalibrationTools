project(CalibrationCutsModule)

find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

add_library(CalibrationCutsModule SHARED "src/calibration_cuts_module.cc" )
target_include_directories(CalibrationCutsModule PUBLIC "include")

target_link_libraries(CalibrationCutsModule PUBLIC Falaise::Falaise)  

install(TARGETS CalibrationCutsModule DESTINATION ${CMAKE_INSTALL_PLUGINDIR})

configure_file("calibration_cuts.conf.in" "calibration_cuts.conf" @ONLY)
configure_file("variant.profile.in" "variant.profile" @ONLY)
configure_file("source_positions.txt.in" "source_positions.txt" @ONLY)
