#!/bin/bash
#SBATCH --job-name="N10MDFrho0315"
#SBATCH --output="N10.%j.%N.out"
#SBATCH --partition=short
#SBATCH --time=0-00:10:00
#SBATCH --error=joberr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --cpus-per-task=1
#SBATCH --export=ALL
#SBATCH --mail-user=aguerre2@uoregon.edu
#SBATCH --mail-type=end
#SBATCH --account=guenzagrp
BD=$PWD
cd $BD
module load slurm
module load racs-spack
module load python3
lmp < in.PEsoft > OUTsoft
cp FINALCONFIG RESTART
lmp < in.PEequil > OUTequil
cp FINALCONFIG RESTART
lmp read_restart FINALCONFIG
lmp write_data outdata
python3 backmap.py
cp lammps-AT-config backmap-round-0
exitlmp < in.PEsoft > OUTsoft
cp FINALCONFIG RESTART
lmp < in.PEequil > OUTequil
cp FINALCONFIG RESTART
lmp read_restart FINALCONFIG
lmp write_data outdata
python3 backmap.py
cp lammps-AT-config backmap-round-1
exit