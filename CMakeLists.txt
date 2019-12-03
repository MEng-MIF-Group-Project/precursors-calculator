cmake_minimum_required(VERSION 3.12)

# set the project name and version
project(Samply VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# configure a header file to pass the version number only
configure_file("version.h.in" version.h)

# control where the static and shared libraries are built so that on windows
# we don't need to tinker with the path to run the executable
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

# find boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost COMPONENTS program_options filesystem regex) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    add_executable("${PROJECT_NAME}" samply.cxx) 
	target_compile_features("${PROJECT_NAME}" PRIVATE cxx_std_17)
    target_link_libraries("${PROJECT_NAME}" ${Boost_LIBRARIES})
endif()

# target_include_directories("${PROJECT_NAME}" PUBLIC "${PROJECT_BINARY_DIR}")
# target_include_directories("${PROJECT_NAME}" PUBLIC include)
