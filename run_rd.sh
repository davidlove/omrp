#! /bin/bash
# Runs a loop over multiple values of gamma for the newsvendor problem
# Arguments:
# $1: Problem Name
# $2: Value of m -- batch size
# $3: k -- number of nonoverlapping batches
# $4: Special designation of output file

for gamma in 0001 0005 0010 0020 0025 0033 0050 0067 0075 0100 0125 0133 0150 0167 0200 0225  0250 0300 0333 0350 0375 0400 0450 0500 0550 0600 0650 0667 0750 1000
do
   if [ $2 -ge $gamma ];
   then
      RUN="./decomp testproblems/${1} -s $2 -r 300 -k $3 -g $gamma | tee $1_m$2_$gamma.dat | grep Rep"
      echo ${RUN}
      eval ${RUN}
   fi
done

echo "Done!"
./parser/parser ${1}_m${2}_*.dat -o ${4}_${1}_m${2}.dat
