#!/bin/bash
PWD=`pwd`
SRC=$PWD/src
FOX=$SRC/lib/FoX

FOXINCLUDE=`echo "-I$FOX/common -I$FOX/fsys -I$FOX/utils -I$FOX/sax"`
FOXOBJ=`echo "$FOX/sax/*.o $FOX/fsys/*.o $FOX/common/*.o $FOX/utils/*.o"`
#FOXLIB=-L$FOX/sax -lFoX_sax

pgfortran -Mcuda=cc20 -Mcudalib=cublas $FOXINCLUDE $FOXOBJ src/modules/types_params.f90 src/modules/ziggurat.f90 src/modules/math.f90 src/modules/IO.f90 src/modules/lattice.f90 src/modules/xml_parser.f90 src/modules/cuda_fermion_prop.f90 src/fermion_tmunu.f90 -o fermion_tmunu.run
