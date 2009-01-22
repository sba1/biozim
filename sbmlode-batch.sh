#!/bin/bash


DRAWER=$1
SMBL2ODE=./sbml2ode 

if ["$DRAWER" -eq ""]
then
  DRAWER="."
fi

mkdir -p results

for filename in $DRAWER/*.xml
do
  echo "Processing $filename"
# echo $filename
  for i in `seq 1 10`
  do
#    echo `basename $filename .xml`.$i.txt
	${SMBL2ODE} --stochastic $filename --maxtime 100 >results/`basename $filename .xml`.$i.txt
  done
done
