#!/bin/bash
#SBATCH --job-name=TEST
#SBATCH --qos=cm_medium
#SBATCH --partition=catamount
#SBATCH --account=catamount
#SBATCH --nodes=4
#SBATCH --time=48:00:00

## Run command
module load mkl
module load python
module load numpy

ns=8
hdir=$(pwd)
for i in {0..5}
do
  for ((j=i; j<=5; j++))
  do
      for ((m=0; m<=ns; m++))
        do
          cwd='D'$i'/'$j'/'$m
          cd $cwd
          mpirun -np 64 -npernode 16 $EXEC  /global/home/users/iwinter/src/vasp.5.3/vasp > vasp.out
          cd $hdir
      done
  done
done
