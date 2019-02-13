#!/bin/bash
#
# Kanata Miyahana
# 2019/02/05
#

prog=ICDW_HMC
input=HMC_PARAM

mchi=2.23606797749979
gg=0.8
lam=0.1
nu=3.0
iseed=37812
ntherm=100
nskip=1
nsample=0

for nmd in 10 20 30 40 50 60 70 80 90
do

################################
# Create Input Parameter File
################################
#  jseed=`echo $iseed $nmd | awk '{printf("%d",$1+32*$2)}'`
  echo ${mchi}d0 ${gg}d0 ${lam}d0 ${nu}d0 > $input
  echo $ntherm $nskip $nsample >> $input
  echo $nmd $iseed >> $input

  output="nmd"${nmd}"_hamildiff.dat"

################################
# Run Program
################################
  echo "# Run program for nmd=${nmd}"
  ./ICDW_HMC > $output
  echo "# End program for nmd=${nmd}"

################################
# Statistical analysis on difference of Hamiltonian
################################

  echo "# Start statistical analysis for nmd=${nmd}"
  hdiff="nmd10to90_hamildiff.dat"

  ./jackbin $output >> $hdiff
  echo "# End statistical analysis for nmd=${nmd}"

done
