#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --partition=medium
#SBATCH --mem=2000

export OMP_NUM_THREADS=1

scratch=/scratch/ivanchen/yusipov/os_lnd/$1
code_base=/home/ivanchen/yusipov/os_lnd/source/cpp/os_lnd/os_lnd
mkdir -p $scratch
mkdir -p $1
cd $scratch
cp $1/config.ini .

cat config.ini

srun $code_base/os_lnd.o

cp -r $scratch/* $1 # Better use a subdirectory of $HOME .
rm -r $scratch/*

