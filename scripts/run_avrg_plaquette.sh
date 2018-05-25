#!/bin/bash

#This script assumes it is being run from inside the project root directory
DIR=`pwd`
#--------------------------------------------
Ns=$1
Nt=$2
BETA=$3
COUNTER=$4
ENDCOUNTER=$5
STEP=$6


Nsf=`printf '%03d' $Ns`
Ntf=`printf '%03d' $Nt`
BETAf=`printf '%.2f' $BETA`
filePath=$DIR/output/lat_conf/links"$Nsf""$Nsf""$Nsf""$Ntf"beta"$BETAf"Sweep

#--------------------------------------------
cd $DIR/output
echo "#Avrg plaquette" > avrg_plaquette.out
while [  $COUNTER -le $ENDCOUNTER ]; do #For each configuration
   COUNTERf=`printf '%09d' $COUNTER`
   $DIR/bin/avrg_plaquette.run $Ns $Ns $Ns $Nt $filePath$COUNTERf.dat >> avrg_plaquette.out
   if [ $? -eq "1" ]; then
      echo  "Something went terribly wrong!!!"
      break
   fi
   echo "Computation for average plaquette in configuration "$COUNTER" finished."
   let COUNTER=COUNTER+STEP
done
cd $DIR
