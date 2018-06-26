#!/bin/bash
LC_NUMERIC=en_US.utf-8 #Needed to avoid user locale using a comma as decimal separator

#This script assumes it is being run from inside the project root directory
DIR=`pwd`
DIM1=$1
DIR1=$2
DIM2=$3
DIR2=$4

#--------------------------------------------
filePath=$DIR/output/lat_conf

#--------------------------------------------
cd $filePath
for filename in links*.dat; do
   Ns=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}).*/\1/'`
   Nt=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}){4}.*/\1/'`
   Nslices=`echo $filename | sed -r 's/.*nslices([0-9]{3}).*/\1/'`
   SWEEP=`echo $filename | sed -r 's/.*Sweep([0-9]{7}).*/\1/'`
   echo $Ns $Nt $Nslices $SWEEP $filename
   #$DIR/bin/wilson_loop_correlation.run $Ns $Ns $Ns $Nt $Nslices $DIM1 $DIR1 $DIM2 $SIR2 $filename
   #cp WilsonCorr_nt$Nt.dat `echo WilsonCorr_nt$(echo $Nt)_nslices$(echo $Nslices)_Sweep$(echo $SWEEP).dat`
done
cd $DIR
