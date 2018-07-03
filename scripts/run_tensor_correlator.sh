#!/bin/bash

#This script uses the following folder as root directory of the project
DIR=`pwd`

THERM=$1
#--------------------------------------------
filepath=$DIR/output/tmunu

#Removes thermalization
cd $filepath
rm -f `ls *.dat | head -n $THERM`
cd $DIR

#Gets lattice dimension
filename=`ls $DIR/output/lat_conf/ | head -n 1`
Ns=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}).*/\1/'`
Nt=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}){4}.*/\1/'`

#Create list of files and count them
ls $filepath > $DIR/output/list_files.in
NUM_FILES=`wc -l $DIR/output/list_files.in | sed -r 's/([0-9]*).*/\1/'`

#Finally, we execute the program
$DIR/bin/tensor_correlator.run $Ns $Ns $Ns $Nt $DIR/output/list_files.in $NUM_FILES
rm $filepath list_files.in