#!/bin/bash

#This script uses the following folder as root directory of the project
DIR=`pwd`
Ns=$1
Nt=$2
BETA=$3
THERM=$4
NSLICES=$5
NMC=$6
DELTA=$7
STEP=$8

#--------------------------------------------
mkdir -p $DIR/output
mkdir -p $DIR/output/lat_conf
cd $DIR/output/lat_conf
$DIR/bin/gen_lat_conf.run $Ns $Ns $Ns $Nt $BETA H $THERM $NSLICES $NMC $DELTA $STEP 
cd $DIR
