#!/bin/bash
#SBATCH --job-name=nanga-omp-130x130
#SBATCH --ntasks=4
#SBATCH --time=03:00:00
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu 
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --output=nanga-omp-130x130.%j

#SBATCH --mail-type=ALL
#SBATCH --mail-user=scott.carrion@tamu.edu

cd $HOME/NangaGPU_omp
ml GCC/6.4.0-2.28
ml CUDA/10.0.130
export PATH=$PATH:$SCRATCH/clang-offloading/install/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SCRATCH/clang-offloading/install/lib

/usr/bin/time build/skyview_omp 130x130
