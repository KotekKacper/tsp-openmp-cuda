#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <omp.h>

int N;                // Number of vertices
int **adj_matrix;     // Adjacency matrix representation of the graph

typedef struct LINKED_LIST {
    int v;
    struct LINKED_LIST *parent;
} linked_list;

void generate_combinations(int level, int ***combinations, linked_list *x, int subset_size, int original_subset_size) {
    // Calculate the cost for combination when reached the end of it
    if (level == 0) {
        // Create a subset array to store the vertices in the combination
        int *subset = (int *)malloc(original_subset_size * sizeof(int));
        int iter = 0;
        while (x) {
            if (x->v > 0) {
                subset[iter] = x->v;
                iter++;
            }
            x = x->parent;
        }
        int bits = 0;
        for (int i = 0; i < original_subset_size; i++) {
            bits |= 1 << subset[i];
        }

        // Calculate the cost for each vertex in the combination and store the minimum cost and its parent
        for (int i = 0; i < original_subset_size; i++) {
            int prev = bits & ~(1 << subset[i]);
            int cost_list[2] = {INT_MAX, INT_MAX};
            for (int j = 0; j < original_subset_size; j++) {
                int cost = combinations[prev][subset[j]][0] + adj_matrix[subset[j]][subset[i]];
                if (subset[j] != subset[i] && cost < cost_list[0]) {
                    cost_list[0] = cost;
                    cost_list[1] = subset[j];
                }
            }
            // Store the minimum cost and its parent in the combinations array
            memcpy(combinations[bits][subset[i]], cost_list, 2 * sizeof(int));
        }
        free(subset);
    } else {
        linked_list l1, l2;
        // Generate combinations by considering two possibilities:
        // 1. Exclude the current level vertex
        // 2. Include the current level vertex

        if (level > subset_size) {
            // Exclude the current level vertex (set it to 0) and continue generating combinations
            l1.v = 0;
            l1.parent = x;
            generate_combinations(level - 1, combinations, &l1, subset_size, original_subset_size);
        }
        if (subset_size > 0) {
            // Include the current level vertex and continue generating combinations with a reduced subset size
            l2.v = level;
            l2.parent = x;
            generate_combinations(level - 1, combinations, &l2, subset_size - 1, original_subset_size);
        }
    }
}


// TSP function to calculate the minimum cost path
int tsp(int *path) {
    int ***combinations;
    int combinations_d1 = pow(2, N) - 1;
    combinations = (int ***)malloc(combinations_d1 * N * 2 * sizeof(int));
    for (int i = 0; i < combinations_d1; i++) {
        combinations[i] = (int **)malloc(N * 2 * sizeof(int));
        for (int j = 0; j < N; j++) {
            combinations[i][j] = (int *)malloc(2 * sizeof(int));
        }
    }

    // Initialize combinations for single vertices
    for (int i = 1; i < N; i++) {
        combinations[1 << i][i][0] = adj_matrix[0][i];
        combinations[1 << i][i][1] = 0;
    }

    // Generate combinations for subset sizes greater than 1
    for (int subset_size = 2; subset_size < N; subset_size++) {
        generate_combinations(N - 1, combinations, NULL, subset_size, subset_size);
    }

    // Calculate optimal cost
    int bits = (pow(2, N) - 1) - 1;
    int optimum = INT_MAX;
    int parent;
    for (int i = 1; i < N; i++) {
        int cost = combinations[bits][i][0] + adj_matrix[i][0];
        if (cost < optimum) {
            optimum = cost;
            parent = i;
        }
    }

    // Retrieve the optimal path
    for (int i = N - 1; i > 0; i--) {
        path[i] = parent;
        int tmp = bits & ~(1 << parent);
        parent = combinations[bits][parent][1];
        bits = tmp;
    }

    path[0] = 0;

    // Freeing the memory
    for (int i = 0; i < combinations_d1; i++) {
        for (int j = 0; j < N; j++) {
            free(combinations[i][j]);
        }
        free(combinations[i]);
    }
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

    return 0;
}
