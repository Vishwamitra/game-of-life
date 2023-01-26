#include <stdio.h>  /*for printf()*/
#include <stdlib.h> /*for rand(),malloc(),free()*/
#include <string.h> /*for strcmp()*/
#include <unistd.h> /*for atoi()*/
#include <mpi.h>    /*for MPI functions*/
#include "glider.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define MESSAGE_TAG 0

// public variables to store common
// information for all the functions
uint8_t **next_generation;
int changed_cells;

uint8_t **local_matrix;
int local_cols;
int local_rows;
int thread_count;

local_rows = GLIDER_HEIGHT;
local_cols = GLIDER_WIDTH;
local_matrix = glider

               char *
               *allocate_memory(int rows, int columns)
{
    int i;
    char *data = malloc(rows * columns * sizeof(char));
    char **arr = malloc(rows * sizeof(char *));
    for (i = 0; i < rows; i++)
        arr[i] = &(data[i * columns]);

    return arr;
}
// to move to the next state
// by applying game rules
int get_next_state(int row, int col, int sum)
{
    int travelled_cells = 0;
    next_generation[i][j] = local_matrix[i][j]

        // If cell is alive = 1 and sum is less than 2
        // and greater than 3, set to dead
        if ((local_matrix[i][j] == '1') && (sum < 2 || sum > 3))
    {
        next_generation[i][j] = '0';
        travelled_cells++;
    }
    // If cell is dead and there are 3
    // consecutive active neighbours
    // let it remain alive = 1
    else if (local_matrix[i][j] == '0') && (sum == 3)
    {
        next_generation[i][j] = '1';
        travelled_cells++;
    }
    return travelled_cells;
}

void get_neighbour_sum(int x, int y, int *sum)
{

    *sum = 0;
    for (i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            sum += local_matrix[x + i][y + j];
        }
    }
    sum = sum - local_matrix[x][y];
}

void calculate_inner_matrix(void)
{

    // For all cells that require no communication at all
    int n;
    changed_cells = 0;
#pragma omp for schedule(static) reduction(+ \
                                           : changed_cells)
    for (n = 0; n < (local_cols - 2) * (local_rows - 2); ++n)
    {
        int i, j, sum;
        i = n / (local_rows - 2) + 2;
        j = n % (local_rows - 2) + 2;
        get_neighbour_sum(i, j, &sum);
        changed_cells += get_next_state(i, j, sum);
    }
}

void calculate_outer_matrix(void)
{
    int i, j;
    int sum;
    /*For all the border-cells*/
    for (i = 1; i <= local_cols; ++i)
        for (j = 1; j <= local_rows; ++j)
        {
            if (i == 1 || i == local_cols || j == 1 || j == local_rows)
            {
                get_neighbour_sum(i, j, &sum);
                get_next_state(i, j, sum);
            }
        }
}

void main()
{

    next_generation = allocate_memory(local_cols + 2, local_rows + 2);
    MPI_Request array_of_requests[16];
    MPI_Status array_of_statuses[16];
    MAX_GENS = 10;

    for (int gen = 0; gen < MAX_GENS; gen++)
    {
        // Start all requests [8 sends + 8 receives]
        MPI_Startall(16, array_of_requests);

        // Overlap communication [calculating inner matrix]
        calculate_inner_matrix();

        // Make sure all requests are completed
        MPI_Waitall(16, array_of_requests, array_of_statuses);

        // We are ready to calculate the outer matrix
        calculate_outer_matrix();

        temp = local_matrix;
        local_matrix = next_generation;
        next_generation = temp;
    }
}
