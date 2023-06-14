#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include <cuda_runtime.h>

int N;
int **adj_matrix;
int *adj_matrix_flat;

typedef struct LINKED_LIST {
    int v;
    struct LINKED_LIST * parent;
} linked_list;

__device__ int global_adj_matrix_flat[500]; 

// Kernel function to save the global adjacency matrix flat on the GPU
__global__ void save_global_adj_matrix_flat(int* input_data, int data_size) {
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_id < data_size) {
        global_adj_matrix_flat[thread_id] = input_data[thread_id];
    }
}

// Kernel function to calculate the cost using CUDA
__global__ void calculate_cost_cuda(int* subset, int subset_size, int* cost_list, int k, int* combinations, int N) {
    int thread_id = threadIdx.x;
    if (thread_id < subset_size && thread_id != k) {
        int cost = combinations[subset[thread_id]] + global_adj_matrix_flat[subset[thread_id] * N + subset[k]];
        atomicMin(&cost_list[0], cost);
        if (cost == cost_list[0]) {
            cost_list[1] = subset[thread_id];
        }
    }
}
     
void generate_combinations(int level, int ***combinations, linked_list * x, int subset_size, int original_subset_size) {
    // Calculate the cost for combination when reached the end of it
    if (level==0) {
        // Create a subset array to store the vertices in the combination
        int *subset = (int *)malloc(original_subset_size * sizeof(int));
        int iter = 0;
        while (x) { 
            if(x->v > 0){
                subset[iter] = x->v;
                iter++;
            }
            x = x -> parent;
        }
        int bits = 0;
        for(int i = 0; i < original_subset_size; i++) {
            bits |= 1 << subset[i];
        }

        // Calculate the cost for each vertex in the combination and store the minimum cost and its parent
        #pragma omp parallel for
        for(int i=0; i<original_subset_size; i++) {
            int prev = bits & ~(1 << subset[i]);
            int cost_list[2] = {INT_MAX,0};

            // Paralization with CUDA
            // Alocating GPU memory
            int* d_subset;
            int* d_cost_list;
            int* d_combinations;
            cudaMalloc((void**)&d_subset, original_subset_size * sizeof(int));
            cudaMalloc((void**)&d_cost_list, 2 * sizeof(int));
            cudaMalloc((void**)&d_combinations, N * sizeof(int));
            // Copying data to GPU memory
            cudaMemcpy(d_subset, subset, original_subset_size * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_cost_list, cost_list, 2 * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_combinations, combinations[prev][0], N * sizeof(int), cudaMemcpyHostToDevice);
            // Running calculate_cost_cuda on GPU
            int threadsPerBlock = original_subset_size;
            int gridSize = 1;
            calculate_cost_cuda<<<gridSize, threadsPerBlock>>>(d_subset, original_subset_size, d_cost_list, i, d_combinations, N);
            // Loading data back from GPU
            cudaMemcpy(cost_list, d_cost_list, 2 * sizeof(int), cudaMemcpyDeviceToHost);
            // Releasing GPU memory
            cudaFree(d_subset);
            cudaFree(d_cost_list);
            cudaFree(d_combinations);

            // Store the minimum cost and its parent in the combinations array
            combinations[bits][0][subset[i]] = cost_list[0];
            combinations[bits][1][subset[i]] = cost_list[1];
        }
        free(subset);
    } else {
        linked_list l1, l2;
        // Generate combinations by considering two possibilities:
        // 1. Exclude the current level vertex
        // 2. Include the current level vertex

        #pragma omp task
        {
        if (level > subset_size) {
            // Exclude the current level vertex (set it to 0) and continue generating combinations
            l1.v = 0;
            l1.parent = x;
            generate_combinations(level - 1, combinations, &l1, subset_size, original_subset_size);
        }
        }
        #pragma omp task
        {
            if (subset_size > 0) {
                // Include the current level vertex and continue generating combinations with a reduced subset size
                l2.v = level;
                l2.parent = x;
                generate_combinations(level - 1, combinations, &l2, subset_size - 1, original_subset_size);
            }
        }
        #pragma omp taskwait  
    }
}

// TSP function to calculate the minimum cost path
int tsp(int *path) {
    int ***combinations;
    int combinations_d1 = pow(2, N) - 1;
    combinations = (int ***)malloc(combinations_d1 * N * 2 * sizeof(int));
    for(int i=0; i<combinations_d1;i++) {
        combinations[i] = (int **)malloc(N * 2 * sizeof(int));
        for(int j=0; j<2; j++) {
            combinations[i][j] = (int *)malloc(N * sizeof(int));
        }
    }

    // Initialize combinations for single vertices
    for(int i=1; i<N; i++) {
        combinations[1<<i][0][i] = adj_matrix_flat[0*N+i];
        combinations[1<<i][1][i] = 0;
    }

    // Generate combinations for subset sizes greater than 1
    for(int subset_size=2; subset_size<N; subset_size++){
        generate_combinations(N-1, combinations, NULL, subset_size, subset_size);
    }

    // Calculate optimal cost
    int bits = (pow(2, N) - 1) - 1;
    int optimum = INT_MAX;
    int parent;
    for (int i = 1; i < N; i++) {
        int cost = combinations[bits][0][i] + adj_matrix[i][0];
        if (cost < optimum) {
            optimum = cost;
            parent = i;
        }
    }

    // Retrieve the optimal path
    for (int i = N - 1; i > 0; i--) {
        path[i] = parent;
        int tmp = bits & ~(1 << parent);
        parent = combinations[bits][1][parent];
        bits = tmp;
    }

    path[0] = 0;

    // Freeing the memory
    free(combinations);

    return optimum;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("No instance filename given\n");
        return -1;
    }

    // Loading size and adjacency matrix from file
    FILE *instance_file;
    instance_file = fopen(argv[1], "r");
    if (instance_file == NULL) {
        printf("Error reading the file\n");
        return -1;
    }
    fscanf(instance_file, "%d", &N);
    adj_matrix = (int **)malloc(N * N * sizeof(int));
    for (int i = 0; i < N; i++) {
        adj_matrix[i] = (int *)malloc(N * sizeof(int));
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fscanf(instance_file, "%d", adj_matrix[i] + j);
        }
    }
    fclose(instance_file);

    adj_matrix_flat = (int *)malloc(N * N * sizeof(int));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            adj_matrix_flat[i * N + j] = adj_matrix[i][j];
        }
    }
    
    // Copying adj_matrix_flat to CUDA
    int data_size = N*N;
    int* d_input_data;
    cudaMalloc((void**)&d_input_data, data_size * sizeof(int));
    cudaMemcpy(d_input_data, adj_matrix_flat, data_size * sizeof(int), cudaMemcpyHostToDevice);
    int threadsPerBlock = 128;
    int gridSize = (threadsPerBlock + data_size - 1) / threadsPerBlock;
    save_global_adj_matrix_flat<<<gridSize, threadsPerBlock>>>(d_input_data, data_size);

    // Running tsp and measuring time
    int out_path[N];
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    int min_cost = tsp(out_path);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    // Printing output
    printf("Time: %f\n", cpu_time_used);
    printf("Path: ");
    for (int i = 0; i < N; i++)
        printf("%d ", out_path[i]);
    printf("\n");
    printf("Minimum cost: %d\n", min_cost);

    // Cleaning resources
    for (int i = 0; i < N; i++) {
        free(adj_matrix[i]);
    }
    free(adj_matrix);
    free(adj_matrix_flat);

    return 0;
}