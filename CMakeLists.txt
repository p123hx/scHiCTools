cmake_minimum_required(VERSION 3.15)
project(scHiCtools)

set(CMAKE_CXX_STANDARD 14)
add_executable(scHiCtools main.cpp cmake-build-debug/load/straw.cpp cmake-build-debug/load/straw.h cmake-build-debug/load/ContactMaps.cpp cmake-build-debug/load/ContactMaps.h cmake-build-debug/load/load_hic_file.cpp cmake-build-debug/load/load_hic_file.h)