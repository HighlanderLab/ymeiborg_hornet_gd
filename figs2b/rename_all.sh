#!bin/bash
#for n in $(seq 1 10) ; do
  #cp slurm.job slurm${n}.job
#done

for n in $(seq 1 10) ; do
  cp fig3b.R fig3b_${n}.R
done
