project(Charge2EnergyModule)

find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

add_library(Charge2EnergyModule SHARED "src/charge2energy_module.cc" )
target_include_directories(Charge2EnergyModule PUBLIC "include")

target_link_libraries(Charge2EnergyModule PUBLIC Falaise::Falaise ${ROOT_LIBRARIES} EnergyCorrectionCalculator)  

install(TARGETS Charge2EnergyModule DESTINATION ${CMAKE_INSTALL_PLUGINDIR})

configure_file("charge2energy.conf.in" "charge2energy.conf" @ONLY)
