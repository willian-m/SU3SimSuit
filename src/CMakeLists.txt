cmake_minimum_required (VERSION 2.6)
project (SU3SimSuit)

enable_language(Fortran)

#======================================================#
#Forces use of pgfortran compiler and set its lib paths#
#======================================================#
set(CMAKE_Fortran_COMPILER pgfortran)

#===================================================#
#Choose paths for libraries and location for modules#
#===================================================#
set(MODULES_DIR ${CMAKE_SOURCE_DIR}/modules)
set(FOX_PATH ${CMAKE_SOURCE_DIR}/lib/FoX)
set(MKLROOT /opt/intel/compilers_and_libraries_2019.2.187/linux/mkl)

find_library(FOXUTILS FoX_utils ${FOX_PATH}/objs/lib)
find_library(FOXCOMMON FoX_common ${FOX_PATH}/objs/lib)
find_library(FOXFSYS FoX_fsys ${FOX_PATH}/objs/lib)
find_library(FOXSAX FoX_sax ${FOX_PATH}/objs/lib)

link_directories(${MKLROOT}/lib/intel64)

MESSAGE(${FOX_PATH})
MESSAGE(${CMAKE_LIBRARY_PATH_FLAG})

#======================================================================#
#Create variables with the modules dependency chain of each application#
#======================================================================#
set(GEN_LAT_CONF_DEPENDENCY ${MODULES_DIR}/xml_parser.f90 ${MODULES_DIR}/heat_bath.f90  ${MODULES_DIR}/objects.f90 ${MODULES_DIR}/lattice.f90 ${MODULES_DIR}/IO.f90 ${MODULES_DIR}/math.f90 ${MODULES_DIR}/types_params.f90 ${MODULES_DIR}/ziggurat.f90) 
set(FERMION_EMT_DEPENDENCY ${MODULES_DIR}/fermion_prop.f90 ${MODULES_DIR}/lattice.f90 ${MODULES_DIR}/math.f90 ${MODULES_DIR}/IO.f90 ${MODULES_DIR}/types_params.f90 ${MODULES_DIR}/ziggurat.f90)
set(AVRG_PLAQUETTE_DEPENDENCY ${GEN_LAT_CONF_DEPENDENCY})
set(TMUNU_DEPENDENCY ${GEN_LAT_CONF_DEPENDENCY})
set(GREEN_FUNCTION_TMUNU_DEPENDENCY ${MODULES_DIR}/statistic.f90 ${MODULES_DIR}/xml_parser.f90 ${MODULES_DIR}/lattice.f90 ${MODULES_DIR}/IO.f90 ${MODULES_DIR}/math.f90 ${MODULES_DIR}/types_params.f90 ${MODULES_DIR}/ziggurat.f90)

#=====================================#
#Lists the applications to be compiled#
#=====================================#
add_executable(gen_lat_conf.run gen_lat_conf.f90 ${GEN_LAT_CONF_DEPENDENCY})
add_executable(avrg_plaquette.run avrg_plaquette.f90 ${AVRG_PLAQUETTE_DEPENDENCY})
add_executable(tmunu.run tmunu.f90 ${TMUNU_DEPENDENCY})
add_executable(green_function_tmunu.run green_function_tmunu.f90 ${GREEN_FUNCTION_TMUNU_DEPENDENCY})


#==========================#
#Set compile and link flags#
#==========================#
SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/xml_parser.f90 PROPERTIES COMPILE_FLAGS -I${FOX_PATH}/objs/finclude)
SET_SOURCE_FILES_PROPERTIES(gen_lat_conf.f90 PROPERTIES COMPILE_FLAGS -I${FOX_PATH}/objs/finclude)
SET_SOURCE_FILES_PROPERTIES(avrg_plaquette.f90 PROPERTIES COMPILE_FLAGS -I${FOX_PATH}/objs/finclude)
SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/statistic.f90  PROPERTIES COMPILE_FLAGS -i8)
SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/statistic.f90  APPEND_STRING PROPERTIES COMPILE_FLAGS -I${MKLROOT}/include)
SET_SOURCE_FILES_PROPERTIES(green_function_tmunu.f90 PROPERTIES COMPILE_FLAGS -I${MKLROOT}/include)
SET_SOURCE_FILES_PROPERTIES(green_function_tmunu.f90 APPEND_STRING PROPERTIES COMPILE_FLAGS -I${FOX_PATH}/objs/finclude)

target_link_libraries(gen_lat_conf.run ${FOXSAX} ${FOXFSYS} ${FOXCOMMON} ${FOXUTILS})
target_link_libraries(avrg_plaquette.run ${FOXSAX} ${FOXFSYS} ${FOXCOMMON} ${FOXUTILS})
target_link_libraries(tmunu.run ${FOXSAX} ${FOXFSYS} ${FOXCOMMON} ${FOXUTILS})
target_link_libraries(green_function_tmunu.run ${FOXSAX} ${FOXFSYS} ${FOXCOMMON} ${FOXUTILS} mkl_rt pthread m dl)

#======================================================#
#Unused flags commented out and kept here for reference#
#======================================================#
#SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/fermion_prop.f90 PROPERTIES COMPILE_FLAGS "-Mcudalib=cublas -Minfo=all,ccff -Mneginfo -ta=multicore")
#SET_SOURCE_FILES_PROPERTIES(${MODULES_DIR}/fermion_prop.f90 PROPERTIES COMPILE_FLAGS "-Mcudalib=cublas -Minfo=all,ccff -Mneginfo -ta=tesla,lineinfo")
#SET_PROPERTY(GLOBAL APPEND PROPERTY COMPILE_FLAGS "-Minfo=all,ccff -Mneginfo")
#SET_TARGET_PROPERTIES(fermion_tmunu.run PROPERTIES LINK_FLAGS "-acc")
#set(CMAKE_LIBRARY_PATH_FLAG ${CMAKE_LIBRARY_PATH_FLAG} "$PGI_HOME/lib" "opt/cuda/lib64" "${FOX_PATH}/objs/lib")
#set(CMAKE_LIBRARY_PATH_FLAG ${CMAKE_LIBRARY_PATH_FLAG} /opt/cuda/lib64)
#set(CMAKE_LIBRARY_PATH_FLAG ${CMAKE_LIBRARY_PATH_FLAG} ${FOX_PATH}/objs/lib)
#set(INCLUDE_DIRECTORIES ${INCLUDE_DIRECTORIES} ${FOX_PATH}/objs/finclude)
