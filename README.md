# TSP with openMP and CUDA
This repository contains TSP implementation using plain C, openMP, CUDA and combined approach to test capabilities of parallelization.

# Content
This repository contains:
- tsp_plain.cu - source file of implementation in plain C
- tsp_openmp.cu - source file of implementation using openMP
- tsp_cuda.cu - source file of implementation using CUDA
- tsp_combined.cu - source file of implementation using openMP and CUDA
- testing.py - python script to test all versions against all test instances
- instances - directory with all instances to be tested
- bar graphs of performed tests

# Running
There is bash script "run.sh" that is used to compile and run chosen version.

Syntax:

./runsh [version] [instance_path]

Available versions:
- plain
- cuda
- openmp
- combined
