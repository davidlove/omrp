#! /bin/bash
# Runs a loop over multiple values of gamma for the newsvendor problem

for gamma in 0001 0002 0003 0005 0010 0015 0020 0025 0030 0033 0050 0067 0075 0100 0125 0133 0150 0167 0200 0225  0250 0300 0333 0375 0500 0667 0750 1000 
do
   if [ $1 -ge $gamma ];
   then
      RUN="./sip -m $1 -g $gamma"
      echo ${RUN}
      eval ${RUN}
      MV_CMD="mv 10D_OMRP_soln.txt m$1_g$gamma.dat"
      echo ${MV_CMD}
      eval ${MV_CMD}
   fi
done

echo "Done running!"
./parser m${1}_g*.dat -o 10D_m${1}.dat
rm -f m${1}_g*.dat
