
# Synthetic Turbulence Generator (STG)

This a C++ code for generating synthetic turbulence that can be used in scale resolving simulations such as Large Eddy Simulations (LES), Delayed Detached Eddy Simulation (DDES).

## Description

STG is a software tool for generating synthetic turbulence, compatible with scale-adaptive simulations like LES, DDES, or any zonal RANS/LES models. The generated turbulence data can be used to specify boundary conditions, such as inflow or at the RANS/LES interface, to accelerate the transition of the flow solution from RANS mode to LES mode.

## Getting Started

### Dependencies

* CUDA compatible C++ compiler (nvc++)

### Installing

* Makefiles are provided in `src` directory

### Executing program

* `bin` directory contains necessary inputs for a test case. Both the CPU and GPU version of STG can be executed from within the `bin` directory.

```
$ cd STG/bin 
$ ../src/stg_gpu.x 
```
On running, the code will generate the synthetic turbulence data in the `bin` directory. (`uprime.dat`, `vprime.dat`, `wprime.dat`).

For post-processing, the turbulence statistics can be compared using `post_proc.py` script in the `bin` directory.

## Authors

Balaji Muralidharan (balaji.murali@gmail.com)  

## Version History

* 0.1
    * Initial Release


