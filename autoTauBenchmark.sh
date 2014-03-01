#!/bin/bash

rm tauBenchmarkConvergence
rm tauBenchmarkConvergenceTable

for i in 1 2 3 4 5 6 7 8
do
  ./build/mc-mini paramFiles/tauBenchmark/tauBenchmark$(($i)) 1>> tauBenchmarkConvergence 2> /dev/null
done

echo     "U Velocity Convergence table" | tee -a tauBenchmarkConvergenceTable
echo -e  "Convergence rate\tResidual" | tee -a tauBenchmarkConvergenceTable
echo -ne "\t\t" | tee -a tauBenchmarkConvergenceTable
prevURes=0
for ures in `grep -v '^#' tauBenchmarkConvergence |
             awk '{print $1}' | 
             sed ':a;N;$!ba;s/\n/ /g'`
do
  ures=`echo ${ures} | sed -e 's/[eE]+*/\\*10\\^/' | bc -l`
  if [ $prevURes != 0 ]; then
    echo -n `echo "l($prevURes / $ures) / l(2)" | bc -l` | tee -a tauBenchmarkConvergenceTable
  fi
  echo -e "\t"$ures | tee -a tauBenchmarkConvergenceTable
  prevURes=$ures
done


echo -e  "\n\n" | tee -a tauBenchmarkConvergenceTable
echo     "V Velocity Convergence table" | tee -a tauBenchmarkConvergenceTable
echo -e  "Convergence rate\tResidual"  | tee -a tauBenchmarkConvergenceTable
echo -ne "\t\t" | tee -a tauBenchmarkConvergenceTable
prevVRes=0
for vres in `grep -v '^#' tauBenchmarkConvergence | 
             awk '{print $2}' | 
             sed ':a;N;$!ba;s/\n/ /g'`
do
  vres=`echo ${vres} | sed -e 's/[eE]+*/\\*10\\^/' | bc -l`
  if [ $prevVRes != 0 ]; then
    echo -n `echo "l($prevVRes / $vres) / l(2)" | bc -l` | tee -a tauBenchmarkConvergenceTable
  fi
  echo -e "\t"$vres | tee -a tauBenchmarkConvergenceTable
  prevVRes=$vres
done
