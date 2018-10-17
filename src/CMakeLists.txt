cmake_minimum_required (VERSION 2.6)
project (SU3SimSuit)

enable_language(Fortran)

set(CMAKE_Fortran_COMPILER pgfortran)
set(CMAKE_LIBRARY_PATH_FLAG $PGI_HOME/lib)
set(CMAKE_LIBRARY_PATH_FLAG /opt/cuda/lib64)

set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
set(FERMION_EMT_DEPENDENCY ${MODULES_DIR}/fermion_prop.f90 ${MODULES_DIR}/lattice.f90 ${MODULES_DIR}/math.f90 ${MODULES_DIR}/IO.f90 ${MODULES_DIR}/types_params.f90 ${MODULES_DIR}/ziggurat.f90)

add_executable(fermion_tmunu.run fermion_tmunu.f90 ${FERMION_EMT_DEPENDENCY})

SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/fermion_prop.f90 PROPERTIES COMPILE_FLAGS "-Mcudalib=cublas")
#SET_TARGET_PROPERTIES(fermion_tmunu PROPERTIES LINK_FLAGS "-Mcudalib=cublas")

target_link_libraries(fermion_tmunu.run -lcudaforblas -lblas -lcublas -Mcudalib=cublas -lcudafor -Mcuda)
