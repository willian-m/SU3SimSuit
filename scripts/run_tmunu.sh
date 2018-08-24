#!/bin/bash

#This script uses the following folder as root directory of the project
DIR=`pwd`

#--------------------------------------------
filename=`ls $DIR/output/lat_conf/ | head -n 1`

Ns=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}).*/\1/'`
Nt=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}){4}.*/\1/'`

mkdir $DIR/output/tmunu
cd $DIR/output/tmunu
COUNTER=1
for filename in $DIR/output/lat_conf/*; do
    echo 'Processing file'
    echo $filename
    $DIR/bin/tmunu.run $Ns $Ns $Ns $Nt $filename
    mv Tmunu.dat Tmunu`printf %08d $COUNTER`.dat
    let COUNTER=COUNTER+1  
done