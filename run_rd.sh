#! /bin/bash
# Runs a loop over multiple values of gamma for the newsvendor problem
# Arguments:
# $1: Problem Name
# $2: Value of m -- batch size
# $3: k -- number of nonoverlapping batches
# $4: r -- number of replications
# $5: Special designation of output file

for gammabar in 5 10 20 33 50 75 100
do
   gamma=$[$gammabar*$2/100]
   RUN="./decomp testproblems/${1} -s $2 -r $4 -k $3 -g $gamma | tee $1_m$2_$gamma.dat | grep Rep"
   echo ${RUN}
   eval ${RUN}
done

echo "Done!"
./parser/parser ${1}_m${2}_*.dat -o ${5}_${1}_m${2}_r${4}.dat
