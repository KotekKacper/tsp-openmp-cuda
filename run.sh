#!/bin/bash

if [ "$1" == "plain" ]; then
    nvcc -Xcompiler -fopenmp ./tsp_plain.cu -o tsp_plain && ./tsp_plain "$2"
elif [ "$1" == "openmp" ]; then
    nvcc -Xcompiler -fopenmp ./tsp_openmp.cu -o tsp_openmp && ./tsp_openmp "$2"
elif [ "$1" == "cuda" ]; then
    nvcc -Xcompiler -fopenmp ./tsp_cuda.cu -o tsp_cuda && ./tsp_cuda "$2"
elif [ "$1" == "combined" ]; then
    nvcc -Xcompiler -fopenmp ./tsp_combined.cu -o tsp_combined && ./tsp_combined "$2"
else
    echo "Invalid argument. Please provide 'plain', 'openmp', 'cuda', or 'combined'."
fi
