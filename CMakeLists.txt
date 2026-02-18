cmake_minimum_required (VERSION 3.10)
project (MLhydroAfterbuner LANGUAGES CXX)

set(libname "LibMLhydroAfterbuner")
set(exename "MLhydroAfterbuner")

# require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_HOME_DIRECTORY}")
set (CMAKE_CXX_FLAGS "-g -O3 -Wall")

INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

add_subdirectory (src)
