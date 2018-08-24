#!/bin/bash

#This script uses the following folder as root directory of the project
DIR=`pwd`

#--------------------------------------------
mkdir -p $DIR/output
mkdir -p $DIR/output/lat_conf
cd $DIR/output/lat_conf
$DIR/bin/gen_lat_conf.run $DIR/input.xml > $DIR/output/action.out
cd $DIR
