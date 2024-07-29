# Check cmake version meets our requirements
cmake_minimum_required(VERSION 3.3)

# Declare project, which will configure compiler for us
project(Charge2EnergyModule)

#Finds ROOT and Falaise packages and defines cmake variables
find_package(Falaise REQUIRED)
find_package(ROOT    REQUIRED)

include_directories(include ${Bayeux_INCLUDE_DIRS})

add_compile_definitions(Falaise_RESOURCE_DIR="${Falaise_RESOURCE_DIR}")

file(GLOB SOURCES "src/*.cc" "src/dicts/*.cc")
add_library(Charge2EnergyModule SHARED ${SOURCES} )

target_link_libraries(Charge2EnergyModule PUBLIC Falaise::Falaise ${ROOT_LIBRARIES})  

install(TARGETS Charge2EnergyModule DESTINATION ${CMAKE_INSTALL_PLUGINDIR})

configure_file("charge2energy.conf.in" "charge2energy.conf" @ONLY)
