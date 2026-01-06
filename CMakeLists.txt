# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.30)
project(
  smithlab_cpp
  VERSION 1.2.0
  DESCRIPTION
  "smithlab_cpp"
  HOMEPAGE_URL https://github.com/smithlabcode/smithlab_cpp
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED off)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++17
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

include(GNUInstallDirs)

if(NOT PROJECT_IS_TOP_LEVEL)
  set(STATIC_ANALYSIS OFF)
endif()

if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
  -Werror
  -Wfatal-errors
)

file(GLOB cpp_files "*.cpp")
# exclude htslib_wrapper unless (below) USE_HTSLIB is set
list(FILTER cpp_files EXCLUDE REGEX "htslib_wrapper.cpp")

set(LIBRARY_OBJECTS "")
foreach(cpp_file ${cpp_files})
  get_filename_component(BASE_NAME ${cpp_file} NAME_WE)
  add_library(${BASE_NAME} OBJECT ${cpp_file})
  list(APPEND LIBRARY_OBJECTS ${BASE_NAME})
endforeach()

if(USE_HTSLIB)
  find_package(HTSLIB REQUIRED)
  find_package(Threads REQUIRED)
  add_library(htslib_wrapper OBJECT htslib_wrapper.cpp)
  list(APPEND LIBRARY_OBJECTS htslib_wrapper)
  target_link_libraries(htslib_wrapper PUBLIC
    HTSLIB::HTSLIB
    Threads::Threads
  )
endif()

add_library(smithlab_cpp)
target_include_directories(smithlab_cpp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(smithlab_cpp PUBLIC
  ${LIBRARY_OBJECTS}
)
