cmake_minimum_required(VERSION 3.3.2 FATAL_ERROR)
cmake_policy(VERSION 3.3.2)

include(ExternalProject)

#TODO
# + add support for Apple systems (needs an apple part??)
# + make sure it correctly builds on the cobalt testbed
# + verify python build.

## project(lepton_injector VERSION 1.0.0 DESCRIPTION "lepton injector library")
project(lepton_injector VERSION 1.0.0 LANGUAGES CXX C)

SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_STANDARD 99)
## enable_language(cxx)
SET (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake
                       ${CMAKE_SOURCE_DIR}/cmake/Packages)

set(CMAKE_BUILD_TYPE Release CACHE STRING "")

INCLUDE (SuiteSparse)
INCLUDE (Python)

set(_FLAGS "-O2 -Wall -fPIC")
set(_FLAGS_DEBUG "-g -O0 -Wall -fPIC")
set(_FLAGS_RELEASE "-O2 -Wall -fPIC -s")

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(_FLAGS "${_FLAGS} -stdlib=libc++")
    set(_FLAGS_DEBUG "${_FLAGS_DEBUG} -stdlib=libc++")
    set(_FLAGS_RELEASE "${_FLAGS_RELEASE} -stdlib=libc++")
endif()

set(CMAKE_CXX_FLAGS ${_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${_FLAGS_DEBUG})
set(CMAKE_CXX_FLAGS_RELEASE ${_FLAGS_RELEASE})

if(CMAKE_VERSION VERSION_LESS 3.12.0)
    if(Python_ROOT_DIR)
        MESSAGE(WARNING "Python_ROOT_DIR is set but will be ignored by this version of CMake; set PYTHON_EXECUTABLE instead")
	endif()
	find_package(PythonInterp)
	find_package(PythonLibs)
	set(PYTHON_FOUND PYTHONLIBS_FOUND AND PYTHONINTERP_FOUND)
else()
    if(PYTHON_EXECUTABLE)
        MESSAGE(WARNING "PYTHON_EXECUTABLE is set but will be ignored by this version of CMake; set Python_ROOT_DIR instead")
	endif()
    find_package(Python COMPONENTS Interpreter Development)
    set(PYTHON_FOUND Python_Interpreter_FOUND AND Python_Development_FOUND)
    set(PYTHON_EXECUTABLE "${Python_EXECUTABLE}")
	set(PYTHON_INCLUDE_DIRS "${Python_INCLUDE_DIRS}")
	set(PYTHON_LIBRARIES "${Python_LIBRARIES}")
	set(PYTHON_VERSION_MAJOR "${Python_VERSION_MAJOR}")
	set(PYTHON_VERSION_MINOR "${Python_VERSION_MINOR}")
endif()

# find boost
find_package(Boost REQUIRED)
if(Boost_FOUND)
    MESSAGE("-- Found Boost at ${Boost_INCLUDE_DIRS}")
    include_directories(${Boost_INCLUDE_DIRS})
else()
    MESSAGE(FATAL_ERROR "Unable to find boost")
endif()

find_package(photospline REQUIRED)

if(photospline_FOUND)
    MESSAGE("-- photospline found")
endif()

include(GNUInstallDirs)

execute_process(
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
    COMMAND touch aclocal.m4 configure
    COMMAND touch Makefile.am Makefile.in aclocal.m4 m4/* rk.pc.in rk/Makefile.am rk/Makefile.in
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    ExternalProject_Add(
      rk
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
      CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env CXXFLAGS=-stdlib=libc++ ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk/configure --prefix=${CMAKE_BINARY_DIR}/extern/rk
      PREFIX ${CMAKE_BINARY_DIR}/extern/rk
      BUILD_COMMAND make
      BUILD_IN_SOURCE 0
    )
else()
    ExternalProject_Add(
      rk
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk
      CONFIGURE_COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/vendor/rk/configure --prefix=${CMAKE_BINARY_DIR}/extern/rk
      PREFIX ${CMAKE_BINARY_DIR}/extern/rk
      BUILD_COMMAND make
      BUILD_IN_SOURCE 0
    )
endif()

install(DIRECTORY "${CMAKE_BINARY_DIR}/extern/rk/include/rk"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    FILES_MATCHING PATTERN "*.hh")

install(DIRECTORY "${CMAKE_BINARY_DIR}/extern/rk/include/rk"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    FILES_MATCHING PATTERN "*.icc")

install(DIRECTORY "${CMAKE_BINARY_DIR}/extern/rk/lib"
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    FILES_MATCHING PATTERN "${CMAKE_SHARED_MODULE_PREFIX}rk${CMAKE_SHARED_LIBRARY_SUFFIX}*")

set(JUST_INSTALL_CEREAL ON CACHE INTERNAL "Cereal just install library.")
set(SKIP_PORTABILITY_TEST ON CACHE INTERNAL "Skip cereal portability tests.")
set(SKIP_PERFORMANCE_COMPARISON ON CACHE INTERNAL "Skip cereal performance comparison.")
add_subdirectory("${PROJECT_SOURCE_DIR}/vendor/cereal" "extern/cereal")
include_directories("${PROJECT_SOURCE_DIR}/vendor/cereal/include")

include(testing)

# add the pybind11 headers
add_subdirectory(vendor/pybind11)

add_subdirectory(vendor/delabella)

add_subdirectory(projects/utilities)
add_subdirectory(projects/serialization)
add_subdirectory(projects/math)
add_subdirectory(projects/dataclasses)
add_subdirectory(projects/geometry)
add_subdirectory(projects/detector)
add_subdirectory(projects/crosssections)
add_subdirectory(projects/distributions)
add_subdirectory(projects/injection)

add_library(LeptonInjector SHARED)
set_property(TARGET LeptonInjector PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(LeptonInjector
    LI_utilities
    LI_serialization
    LI_math
    LI_dataclasses
    LI_geometry
    LI_detector
    LI_crosssections
    LI_distributions
    LI_injection
)

install(TARGETS LeptonInjector
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

#add_subdirectory(tests)

configure_file("${PROJECT_SOURCE_DIR}/setup.py.in"
               "${PROJECT_BINARY_DIR}/setup.py")

MESSAGE("")
MESSAGE("Done!")
MESSAGE("Run `make` then `make install`")
