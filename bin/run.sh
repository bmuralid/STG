#!/usr/bin/env bash
#

if [ ! -d "CPU" ]; then
  mkdir CPU
fi

time ../src/stg_cpu.x 1000 > res_cpu.out

mv uprime.dat vprime.dat wprime.dat sem.dat CPU

if [ ! -d "GPU" ]; then
  mkdir GPU
fi

time ../src/stg_gpu.x 1000 > res_gpu.out


mv uprime.dat vprime.dat wprime.dat sem.dat GPU


