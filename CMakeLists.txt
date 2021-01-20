cmake_minimum_required(VERSION 3.15)
project(scHiCtools)
add_definitions(-DHAVE_CBLAS=1)
find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)
#[[load/straw.cpp load/straw.h]]
add_executable(scHiCtools embedding/reproducibility.cpp embedding/reproducibility.h load/ContactMaps.cpp load/ContactMaps.h load/load_hic_file.cpp load/load_hic_file.h embedding/embedding.cpp embedding/embedding.h test/test_reproducibility.cpp load/processing_utils.h load/processing_utils.cpp test/test_all.cpp)
if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif ()
target_link_libraries(scHiCtools xtensor xtensor::optimize -march=native xtensor-blas)


if (WIN32)
    find_package(OpenBLAS REQUIRED)
    set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif ()

message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})

target_link_libraries(scHiCtools ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_include_directories(scHiCtools PRIVATE ${xtensor_blas_INCLUDE_DIRS})
