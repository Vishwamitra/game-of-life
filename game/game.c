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

        // Check if it has remained the same using a flag

        // if (termination_check(comm_2D, my_rank))
        // {
        //     if (!my_rank)
        //         printf("No change on %d generation\n", gen);
        //     break;
        // }

        // local matrix is set to next generation matrix
        temp = local_matrix;
        local_matrix = next_generation;
        next_generation = temp;
    }
}
// void find_neighbours(MPI_Comm comm_2D, int my_rank, int NPROWS, int NPCOLS, int *left, int *right, int *top, int *bottom, int *topleft, int *topright, int *bottomleft, int *bottomright)
// {

//     int source, dest, disp = 1;
//     int my_coords[2];
//     int corner_coords[2];
//     int corner_rank;

//     /*Finding top/bottom neighbours*/
//     MPI_Cart_shift(comm_2D, 0, disp, top, bottom);

//     /*Finding left/right neighbours*/
//     MPI_Cart_shift(comm_2D, 1, disp, left, right);

//     /*Finding top-right corner*/
//     MPI_Cart_coords(comm_2D, my_rank, 2, my_coords);
//     corner_coords[0] = my_coords[0] - 1;
//     corner_coords[1] = (my_coords[1] + 1) % NPCOLS;
//     if (corner_coords[0] < 0)
//         corner_coords[0] = NPROWS - 1;
//     MPI_Cart_rank(comm_2D, corner_coords, topright);

//     /*Finding top-left corner*/
//     MPI_Cart_coords(comm_2D, my_rank, 2, my_coords);
//     corner_coords[0] = my_coords[0] - 1;
//     corner_coords[1] = my_coords[1] - 1;
//     if (corner_coords[0] < 0)
//         corner_coords[0] = NPROWS - 1;
//     if (corner_coords[1] < 0)
//         corner_coords[1] = NPCOLS - 1;
//     MPI_Cart_rank(comm_2D, corner_coords, topleft);

//     /*Finding bottom-right corner*/
//     MPI_Cart_coords(comm_2D, my_rank, 2, my_coords);
//     corner_coords[0] = (my_coords[0] + 1) % NPROWS;
//     corner_coords[1] = (my_coords[1] + 1) % NPCOLS;
//     MPI_Cart_rank(comm_2D, corner_coords, bottomright);

//     /*Finding bottom-left corner*/
//     MPI_Cart_coords(comm_2D, my_rank, 2, my_coords);
//     corner_coords[0] = (my_coords[0] + 1) % NPROWS;
//     corner_coords[1] = my_coords[1] - 1;
//     if (corner_coords[1] < 0)
//         corner_coords[1] = NPCOLS - 1;
//     MPI_Cart_rank(comm_2D, corner_coords, bottomleft);
// }

// int termination_check(MPI_Comm comm_2D, int my_rank)
// {
//     int sum;
//     int changed;

//     if (changed_cells)
//         changed = 1;
//     else
//         changed = 0;

//     MPI_Allreduce(&changed, &sum, 1, MPI_INT, MPI_SUM, comm_2D);

//     if (sum == 0)
//         return 1;
//     else
//         return 0;
// }

// void print_local_matrix(void)
// {

//     int i, j;
//     for (i = 1; i <= local_cols; ++i)
//     {
//         for (j = 1; j <= local_rows; ++j)
//             printf("%c", local_matrix[i][j]);
//         printf("\n");
//     }
// }

// void print_neighbours(int my_rank, int left, int right, int top, int bottom, int topleft, int topright, int bottomleft, int bottomright)
// {

//     int rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     if (my_rank == rank)
//     {
//         printf("My_rank:%d\n============\n", my_rank);
//         printf("Left:%d\n", left);
//         printf("Right:%d\n", right);
//         printf("Top:%d\n", top);
//         printf("Bottom:%d\n", bottom);
//         printf("TopLeft:%d\n", topleft);
//         printf("TopRight:%d\n", topright);
//         printf("BottomLeft:%d\n", bottomleft);
//         printf("BottomRight:%d\n=============\n", bottomright);
//     }
// }

// void game(MPI_Comm comm_2D, int my_rank, int NPROWS, int NPCOLS, int MAX_GENS)
// {

//     int gen, i, j;
//     next_generation = allocate_memory(local_cols + 2, local_rows + 2);
//     char **temp;

//     /*16 requests , 16 statuses */
//     MPI_Request array_of_requests[16];
//     MPI_Status array_of_statuses[16];

//     /*Create 4 datatypes for sending*/
//     MPI_Datatype firstcolumn_send, firstrow_send, lastcolumn_send, lastrow_send;
//     create_datatype(&firstcolumn_send, 1, 1, local_cols, 1);
//     create_datatype(&firstrow_send, 1, 1, 1, local_rows);
//     create_datatype(&lastcolumn_send, 1, local_rows, local_cols, 1);
//     create_datatype(&lastrow_send, local_cols, 1, 1, local_rows);

//     /*Create 4 datatypes for receiving*/
//     MPI_Datatype firstcolumn_recv, firstrow_recv, lastcolumn_recv, lastrow_recv;
//     create_datatype(&firstcolumn_recv, 1, 0, local_cols, 1);
//     create_datatype(&firstrow_recv, 0, 1, 1, local_rows);
//     create_datatype(&lastcolumn_recv, 1, local_rows + 1, local_cols, 1);
//     create_datatype(&lastrow_recv, local_cols + 1, 1, 1, local_rows);

//     /*Find ranks of my 8 neighbours*/
//     int left, right, bottom, top, topleft, topright, bottomleft, bottomright;
//     find_neighbours(comm_2D, my_rank, NPROWS, NPCOLS, &left, &right, &top, &bottom, &topleft, &topright, &bottomleft, &bottomright);

//     MPI_Send_init(&(local_matrix[0][0]), 1, firstcolumn_send, left, MESSAGE_TAG, comm_2D, &array_of_requests[0]);
//     MPI_Send_init(&(local_matrix[0][0]), 1, firstrow_send, top, MESSAGE_TAG, comm_2D, &array_of_requests[1]);
//     MPI_Send_init(&(local_matrix[0][0]), 1, lastcolumn_send, right, MESSAGE_TAG, comm_2D, &array_of_requests[2]);
//     MPI_Send_init(&(local_matrix[0][0]), 1, lastrow_send, bottom, MESSAGE_TAG, comm_2D, &array_of_requests[3]);
//     MPI_Send_init(&(local_matrix[1][1]), 1, MPI_CHAR, topleft, MESSAGE_TAG, comm_2D, &array_of_requests[4]);
//     MPI_Send_init(&(local_matrix[1][local_rows]), 1, MPI_CHAR, topright, MESSAGE_TAG, comm_2D, &array_of_requests[5]);
//     MPI_Send_init(&(local_matrix[local_cols][local_rows]), 1, MPI_CHAR, bottomright, MESSAGE_TAG, comm_2D, &array_of_requests[6]);
//     MPI_Send_init(&(local_matrix[local_cols][1]), 1, MPI_CHAR, bottomleft, MESSAGE_TAG, comm_2D, &array_of_requests[7]);

//     MPI_Recv_init(&(local_matrix[0][0]), 1, firstcolumn_recv, left, MESSAGE_TAG, comm_2D, &array_of_requests[8]);
//     MPI_Recv_init(&(local_matrix[0][0]), 1, firstrow_recv, top, MESSAGE_TAG, comm_2D, &array_of_requests[9]);
//     MPI_Recv_init(&(local_matrix[0][0]), 1, lastcolumn_recv, right, MESSAGE_TAG, comm_2D, &array_of_requests[10]);
//     MPI_Recv_init(&(local_matrix[0][0]), 1, lastrow_recv, bottom, MESSAGE_TAG, comm_2D, &array_of_requests[11]);
//     MPI_Recv_init(&(local_matrix[0][0]), 1, MPI_CHAR, topleft, MESSAGE_TAG, comm_2D, &array_of_requests[12]);
//     MPI_Recv_init(&(local_matrix[0][local_rows + 1]), 1, MPI_CHAR, topright, MESSAGE_TAG, comm_2D, &array_of_requests[13]);
//     MPI_Recv_init(&(local_matrix[local_cols + 1][local_rows + 1]), 1, MPI_CHAR, bottomright, MESSAGE_TAG, comm_2D, &array_of_requests[14]);
//     MPI_Recv_init(&(local_matrix[local_cols + 1][0]), 1, MPI_CHAR, bottomleft, MESSAGE_TAG, comm_2D, &array_of_requests[15]);

// for (gen = 0; gen < MAX_GENS; gen++)
// {
//     // Start all requests [8 sends + 8 receives]
//     MPI_Startall(16, array_of_requests);

//     // Overlap communication [calculating inner matrix]
//     calculate_inner_matrix();

//     // Make sure all requests are completed
//     MPI_Waitall(16, array_of_requests, array_of_statuses);

//     // We are ready to calculate the outer matrix
//     calculate_outer_matrix();

//     // Check if it has remained the same using a flag

//     if (termination_check(comm_2D, my_rank))
//     {
//         if (!my_rank)
//             printf("No change on %d generation\n", gen);
//         break;
//     }

//     // local matrix is set to next generation matrix
//     temp = local_matrix;
//     local_matrix = next_generation;
//     next_generation = temp;
// }

//     // Free resources
//     MPI_Type_free(&firstcolumn_send);
//     MPI_Type_free(&firstrow_send);
//     MPI_Type_free(&lastcolumn_send);
//     MPI_Type_free(&lastrow_send);

//     MPI_Type_free(&firstcolumn_recv);
//     MPI_Type_free(&firstrow_recv);
//     MPI_Type_free(&lastcolumn_recv);
//     MPI_Type_free(&lastrow_recv);

//     free(next_generation[0]);
//     free(next_generation);
// }

// int main(int argc, char **argv)
// {
//     /************************************************************************************************************/
//     int size, rank, i, j, proc, COLS, ROWS, MAX_GENS;
//     int flag1 = 0, flag2 = 0, flag3 = 0, flag4 = 0, flag5 = 0;
//     double local_start, local_finish, local_elapsed, elapsed;
//     char *filename;

//     MPI_Init(&argc, &argv);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

//     for (i = 0; i < argc; ++i)
//     {
//         if (!strcmp("-n", argv[i]))
//         {
//             ROWS = atoi(argv[i + 1]);
//             flag1 = 1;
//         }
//         else if (!strcmp("-m", argv[i]))
//         {
//             COLS = atoi(argv[i + 1]);
//             flag2 = 1;
//         }
//         else if (!strcmp("-max", argv[i]))
//         {

//             MAX_GENS = atoi(argv[i + 1]);
//             flag3 = 1;
//         }

//         else if (!strcmp("-f", argv[i]))
//         {
//             filename = malloc((strlen(argv[i + 1] + 1)) * sizeof(char));
//             strcpy(filename, argv[i + 1]);
//             flag4 = 1;
//         }
//         else if (!strcmp("-t", argv[i]))
//         {
//             thread_count = atoi(argv[i + 1]);
//             flag5 = 1;
//         }
//     }

//     /*Sanity Check*/
//     if (!flag1 || !flag2 || !flag3 || !flag5)
//     {
//         if (rank == 0)
//             printf("Usage:mpiexec [-n <NoPROCESSES>] [-f <machine_file>] ./gol -n <ROWS> -m <COLUMNS> -max <MAX_GENS> -t <threads>[-f <inputfile>] \nExiting...\n\n");
//         MPI_Finalize();
//         exit(1);
//     }

//     /************************************************************************************************************/
//     /* Setup virtual 2D topology*/
//     int dims[2] = {0, 0};
//     MPI_Dims_create(size, 2, dims);
//     int periods[2] = {1, 1}; /*Periodicity in both dimensions*/
//     int my_coords[2];
//     MPI_Comm comm_2D;
//     MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_2D);
//     MPI_Cart_coords(comm_2D, rank, 2, my_coords);

//     const int NPROWS = dims[0]; /* Number of 'block' rows */
//     const int NPCOLS = dims[1]; /* Number of 'block' cols */
//     int *num_rows;              /* Number of rows for the i-th process [local_cols]*/
//     int *num_cols;              /* Number of columns for the i-th process [local_rows]*/
//     int *extent;                /* Extent for the i-th process [local_extent]*/
//     int *disps;                 /*Displacement for the i-th process [local_disp]*/

//     int local_disp;
//     int local_extent;

//     /************************************************************************************************************/
//     // char name[MPI_MAX_PROCESSOR_NAME];
//     // int namelen = MPI_MAX_PROCESSOR_NAME;
//     // MPI_Get_processor_name(name,&namelen);
//     // printf("Process:%d | Machine:%s\n",my_rank+1,name);

//     /*Start timing*/
//     MPI_Barrier(comm_2D);
//     local_start = MPI_Wtime();

//     /* Calculate rows,cols,displacement,extent for each process */
//     if (rank == 0)
//     {

//         num_rows = (int *)malloc(size * sizeof(int));
//         num_cols = (int *)malloc(size * sizeof(int));

//         calculate_rows_columns(num_rows, num_cols, size, ROWS, COLS, NPROWS, NPCOLS);

//         if (flag4)
//         {
//             disps = (int *)malloc(size * sizeof(int));
//             extent = (int *)malloc(size * sizeof(int));

//             calculate_disp(disps, num_rows, num_cols, COLS, NPROWS, NPCOLS);
//             calculate_extent(extent, num_cols, NPROWS, NPCOLS);
//         }
//     }

//     /************************************************************************************************************/

//     // Scatter dimensions,displacement,extent of each process
//     MPI_Scatter(num_rows, 1, MPI_INT, &local_cols, 1, MPI_INT, 0, comm_2D);
//     MPI_Scatter(num_cols, 1, MPI_INT, &local_rows, 1, MPI_INT, 0, comm_2D);

//     if (rank == 0)
//     {
//         free(num_rows);
//         free(num_cols);
//     }

//     if (flag4)
//     {
//         MPI_Scatter(disps, 1, MPI_INT, &local_disp, 1, MPI_INT, 0, comm_2D);
//         MPI_Scatter(extent, 1, MPI_INT, &local_extent, 1, MPI_INT, 0, comm_2D);

//         if (rank == 0)
//         {
//             free(disps);
//             free(extent);
//         }
//     }

//     /************************************************************************************************************/

//     local_matrix = allocate_memory(local_cols + 2, local_rows + 2);

//     /*If a file has been provided*/
//     if (flag4)
//         read_file(filename, local_disp, local_extent, rank);
//     /*If no file has been provided*/
//     /*Each process will fill its own submatrix randomly*/
//     else
//     {
//         int n;
//         // # pragma omp parallel for num_threads(thread_count) schedule(dynamic)
//         for (n = 0; n < local_cols * local_rows; ++n)
//         {
//             int i, j;
//             i = n / local_rows + 1;
//             j = n % local_rows + 1;
//             if (rand() % 2)
//                 local_matrix[i][j] = '1';
//             else
//                 local_matrix[i][j] = '0';
//         }
//     }

//     /*Each process will start the game*/
//     game(comm_2D, rank, NPROWS, NPCOLS, MAX_GENS);

//     /*Getting finish time*/
//     local_finish = MPI_Wtime();
//     local_elapsed = local_finish - local_start;
//     MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, comm_2D);

//     if (!rank)
//         printf("Elapsed time:%.3f seconds\n", elapsed);

//     free(local_matrix[0]);
//     free(local_matrix);

//     MPI_Finalize();
// }

// void calculate_rows_columns(int *num_rows, int *num_cols, int size, int ROWS, int COLS, int NPROWS, int NPCOLS)
// {
//     int i, j;
//     for (i = 0; i < size; i++)
//     {
//         num_rows[i] = ROWS / NPROWS;
//         num_cols[i] = COLS / NPCOLS;
//     }
//     for (i = 0; i < (ROWS % NPROWS); i++)
//     {
//         for (j = 0; j < NPCOLS; j++)
//         {
//             num_rows[i * NPCOLS + j]++;
//         }
//     }
//     for (i = 0; i < (COLS % NPCOLS); i++)
//     {
//         for (j = 0; j < NPROWS; j++)
//         {
//             num_cols[i + NPROWS * j]++;
//         }
//     }
// }

// void calculate_disp(int *disps, int *num_rows, int *num_cols, int COLS, int NPROWS, int NPCOLS)
// {

//     int i, j;

//     for (i = 0; i < NPROWS; i++)
//         for (j = 0; j < NPCOLS; j++)
//             if (j == 0)
//             {
//                 int row;
//                 disps[i * NPCOLS + j] = 0;
//                 /*For all rows above me*/
//                 for (row = 1; row <= i; row++)
//                 {
//                     disps[i * NPCOLS + j] += (COLS + 1) * num_rows[i * NPCOLS + j - row * NPCOLS];
//                 }
//             }
//             else
//                 /*Just add num_cols of the left process to its displacement*/
//                 disps[i * NPCOLS + j] = disps[i * NPCOLS + j - 1] + num_cols[i * NPCOLS + j - 1];
// }

// void calculate_extent(int *extent, int *num_cols, int NPROWS, int NPCOLS)
// {
//     int i, j;

//     for (i = 0; i < NPROWS; i++)
//         for (j = 0; j < NPCOLS; j++)
//         {
//             int current = i * NPCOLS + j;
//             int block;

//             /*Add newline to the extent*/
//             /*Add num_cols of this process*/
//             extent[current] = 1 + num_cols[current];

//             /*For all blocks on the left of me*/
//             for (block = i * NPCOLS; block < current; block++)
//                 extent[current] += num_cols[block];
//             /*For all blocks on the right */
//             for (block = current + 1; block < i * NPCOLS + NPCOLS; block++)
//                 extent[current] += num_cols[block];
//         }
// }

// void read_file(char *filename, int local_disp, int local_extent, int rank)
// {

//     MPI_File fh;
//     int i, j, error;

//     /*Open the file*/
//     MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
//     error = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

//     if (error != MPI_SUCCESS)
//     {
//         /*Only Process 0 will print the error message */
//         if (rank == 0)
//         {
//             char error_string[BUFSIZ];
//             int length_of_error_string;

//             MPI_Error_string(error, error_string, &length_of_error_string);
//             fprintf(stderr, "%s\nExiting...\n\n", error_string);
//         }

//         /*Every process is about to exit*/
//         MPI_Finalize();
//         exit(1);
//     }

//     int lb = 0;
//     int count;
//     int buffsize = local_cols * local_rows;
//     char buff[buffsize];
//     MPI_Datatype etype, filetype, contig;
//     etype = MPI_CHAR;

//     /*Create contiguous datatype of local_rows chars*/
//     MPI_Type_contiguous(local_rows, MPI_CHAR, &contig);

//     /*Extend this datatype*/
//     MPI_Type_create_resized(contig, lb, (MPI_Aint)(local_extent) * sizeof(char), &filetype);

//     /*Commit this datatype*/
//     MPI_Type_commit(&filetype);

//     /*Each process will now have its own file view*/
//     MPI_File_set_view(fh, (MPI_Offset)(local_disp) * sizeof(char), etype, filetype, "native", MPI_INFO_NULL);

//     /*At first we read all local_cols*local_rows characters in a temp buffer*/
//     MPI_Status status;

//     MPI_File_read_all(fh, buff, buffsize, MPI_CHAR, &status);
//     MPI_Get_count(&status, MPI_CHAR, &count);

//     if (count != buffsize)
//     {
//         fprintf(stderr, "Read:%d instead of %d\n", count, buffsize);
//         fflush(stderr);
//         MPI_Finalize();
//         exit(1);
//     }
//     else
//     {
//         /*We initialize our local matrix with the characters in the buffer*/
//         int k = 0;
//         for (i = 1; i <= local_cols; ++i)
//             for (j = 1; j <= local_rows; ++j)
//             {
//                 local_matrix[i][j] = buff[k];
//                 k++;
//             }
//     }

//     MPI_File_close(&fh);
// }