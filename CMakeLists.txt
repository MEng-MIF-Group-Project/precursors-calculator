cmake_minimum_required(VERSION 3.12)

# set the project name and version
project(samply VERSION 1.0)

set(CMAKE_TOOLCHAIN_FILE $ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake CACHE STRING "")
set(VCPKG_TARGET_TRIPLET "x64-windows" CACHE STRING "")

#set(Boost_NO_BOOST_CMAKE ON)

set(SOURCE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/elementdb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/input.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/reagent.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/reagentdb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/samply.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solutiondb.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/statedb.cpp
)

set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/defines.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/element.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/elementdb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/input.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/reagent.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/reagentdb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solution.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/solutiondb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/statedb.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/utils.h
)

set(CMAKE_BUILD_TYPE "Release")

set (Boost_DEBUG 1)
set (Boost_USE_STATIC_LIBS ON)
set (Boost_USE_MULTITHREADED ON)

set (BOOST_ROOT $ENV{VCPKG_ROOT}/installed/x64-windows/include)
set (BOOST_INCLUDEDIR $ENV{VCPKG_ROOT}/installed/x64-windows/include)
set (BOOST_LIBRARYDIR $ENV{VCPKG_ROOT}/installed/x64-windows/lib)

set (BOOST_LIBRARIES
    ${BOOST_LIBRARYDIR}/boost_regex-vc140-mt.lib
    ${BOOST_LIBRARYDIR}/boost_filesystem-vc140-mt.lib
    ${BOOST_LIBRARYDIR}/boost_program_options-vc140-mt.lib
    ${BOOST_LIBRARYDIR}/boost_system-vc140-mt.lib
)

set (EIGEN3_INCLUDEDIR $ENV{VCPKG_ROOT}/installed/x64-windows/include)

#find_package(Boost)
#find_package(Boost COMPONENTS program_options regex filestyem REQUIRED)



include_directories(${EIGEN3_INCLUDEDIR})
include_directories(${BOOST_INCLUDEDIR})

add_executable(${PROJECT_NAME} ${SOURCE} ${HEADERS}) 
#target_include_directories(${PROJECT_NAME} ${PROJECT_SOURCE_DIR})
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)
target_link_libraries(${PROJECT_NAME} ${BOOST_LIBRARIES})