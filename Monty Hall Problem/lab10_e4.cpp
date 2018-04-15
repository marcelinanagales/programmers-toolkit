# include <mpi.h>               // OpenMPI Example 4
# include <iostream>            // cout 
# include <stdlib.h>            // rand


// function for exercise 4
double openMPI_monteHall(int num, int cores){

//  srand(time(NULL));              // for the random chosen and picked

    double winPercent = 0;
    int chosen;  
    int choice;
    int empty;
    int new_choice = 0;

    #pragma omp parallel \
    shared( num) \
    private( i ) 

    #pragma omp for reduction(+: winPercent) \
    num_threads(threadnum)
    for(int i = 0; i < num; i++){
        chosen = rand() % 3;  
        choice = rand() % 3;
        empty = chosen;

        // reveal loss door //
        while (empty == chosen || empty == choice){
            empty = rand() % 3;
        }

        if (choice == chosen){
            winPercent++;
        }

    }

    return winPercent/num;

}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
           " out of %d processors\n",
           processor_name, world_rank, world_size);

    // Finalize the MPI environment.
    MPI_Finalize();
}