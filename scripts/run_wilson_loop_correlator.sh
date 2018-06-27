#!/bin/bash
LC_NUMERIC=en_US.utf-8 #Needed to avoid user locale using a comma as decimal separator

#This script assumes it is being run from inside the project root directory
DIR=`pwd`
NSWEEPS=$1
NSUBSWEEPS=$2
STEPSIZE=$3
DIM1=$4
DIR1=$5
DIM2=$6
DIR2=$7

#--------------------------------------------
filePath=$DIR/output/lat_conf

#--------------------------------------------
cd $filePath
#Get first parameter as example
filename=`ls $DIR/output/lat_conf/ | head -n 1`

#Define parameters from filename
Ns=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}).*/\1/'`
Nt=`echo $filename | sed -r 's/[^0-9]*([0-9]{3}){4}.*/\1/'`
BETA=`echo $filename | sed -r 's/.*beta([0-9]\.[0-9]{2}).*/\1/'`
Nslices=`echo $filename | sed -r 's/.*nslices([0-9]{3}).*/\1/'`
filename="$filePath"'/links'"$Ns""$Ns""$Ns""$Nt"'beta'"$BETA"'nslices'"$Nslices"'Sweep'
$DIR/bin/wilson_loop_correlation.run $Ns $Ns $Ns $Nt $Nslices $NSWEEPS $NSUBSWEEPS $STEPSIZE $DIM1 $DIR1 $DIM2 $DIR2 $filename

   #cp WilsonCorr_nt$Nt.dat `echo WilsonCorr_nt$(echo $Nt)_nslices$(echo $Nslices)_Sweep$(echo $SWEEP).dat`
cd $DIR
