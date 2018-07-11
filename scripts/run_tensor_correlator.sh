#!/bin/bash

#This script uses the following folder as root directory of the project
DIR=`pwd`

THERM=$1
MEAS=$2
#--------------------------------------------
filepath=$DIR/output/tmunu

#Removes thermalization
cd $filepath
mkdir -p $filepath/pre_therm
mkdir -p $filepath/measurements_workspace
for file in  `find $filepath -name \*.dat | head -n $THERM`; do
   mv $file $filepath/pre_therm
done
for file in  `find $filepath -name \*.dat | head -n $MEAS`; do
   ln -s $file $filepath/measurements_workspace
done

#Gets lattice dimension
filename=`ls $DIR/output/lat_conf/ | head -n 1`
Ns=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}).*/\1/'`
Nt=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}){4}.*/\1/'`

#Create list of files and count them
find $filepath/measurements_workspace/* -name \*.dat > $DIR/output/list_files.in

NUM_FILES=`wc -l $DIR/output/list_files.in | sed -r 's/([0-9]*).*/\1/'`

cd $DIR/output
#Finally, we execute the program
$DIR/bin/tensor_correlator.run $Ns $Ns $Ns $Nt $DIR/output/list_files.in $NUM_FILES

#echo $Ns $Ns $Ns $Nt $DIR/output/list_files.in $NUM_FILES
#cat $DIR/output/list_files.in


#Reorganizes source dir
for file in `find $filepath/pre_therm/ -name \*.dat`; do
   mv $file $filepath
done
find $filepath/measurements_workspace -type l | xargs rm
cd $DIR
