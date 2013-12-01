#! /bin/bash
# Runs a loop over multiple values of gamma for the newsvendor problem

for gamma in 0001 0002 0003 0005 0010 0015 0020 0025 0030 0033 0050 0067 0075 0100 0125 0133 0150 0167 0200 0225  0250 0300 0333 0375 0500 0667 0750 1000
do
   if [ $1 -ge $gamma ];
   then
      echo "./newsvendor -m $1 -g $gamma"
      ./newsvendor -m $1 -g $gamma
      echo "mv Newsvendor_MRP_soln.txt m$1_g$gamma.dat"
      mv Newsvendor_MRP_soln.txt m$1_g$gamma.dat
   fi
done
