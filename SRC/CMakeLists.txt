cmake_minimum_required(VERSION 3.16)

#set(MAKE_Fortran_COMPILER   /usr/bin/gfortran)
set(MAKE_Fortran_COMPILER  /opt/intel/oneapi/compiler/2023.2.0/linux/bin/intel64/ifort)

enable_language(Fortran)

Project(MSAT)

file(GLOB_RECURSE SRC_FILES ./*.f90)

add_executable(MSAT ${SRC_FILES})
find_package(LAPACK)

target_link_libraries(MSAT LAPACK::LAPACK)