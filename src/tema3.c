#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

# define NR_CLUSTERS 4
# define CLUSTERS_TAG 0

void print_topology(int rank, int *nrOfProc, int *topology) {
    char print_top[200];
    sprintf(print_top,"%d -> ",rank);
    int index = 0;
    for(int k = 0; k < NR_CLUSTERS; k++) {
        char coord[5];
        sprintf(coord,"%d:",k);
        strcat(print_top,coord);
        for (int i = index; i < index + nrOfProc[k]; i++) {
            char el[5];
            sprintf(el,"%d", topology[i]);
            strcat(print_top,el);
            if (i + 1 == index + nrOfProc[k]){
                strcat(print_top," ");
            } else {
                strcat(print_top,",");
            }
        }
        index += nrOfProc[k];
    }
    printf("%s\n",print_top);
}

int *send_vector(int rank, int index, int iterations, int *nrOfProc, int *v, int *workers ,int N) {
    MPI_Status status;
    for (int i = 0; i < nrOfProc[rank]; i++) {
        MPI_Send(&N,1, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        MPI_Send(v,N, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        MPI_Send(&iterations,1, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        MPI_Send(&index,1, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",rank,workers[i]);
        MPI_Recv(&v[index], iterations, MPI_INT, workers[i], 0, MPI_COMM_WORLD, &status);
        index += iterations;
    }
    return v;
}

/* Sends the worker 1 element to multiply if the rest is not 0*/
int *send_rest(int rank, int *v, int *workers,int *nrOfProc, int N, int *count, int *rest){
    MPI_Status status;
    for (int i = 0; i < nrOfProc[rank]; i++) {
        MPI_Send(v,N, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        MPI_Send(rest,1, MPI_INT, workers[i], 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",rank,workers[i]);
        if (*rest > 0) {
            MPI_Recv(v, N, MPI_INT, workers[i], 0, MPI_COMM_WORLD, &status);
            (*count)++;
        }
        (*rest)--;
    }
    return v;
}

int main (int argc, char *argv[])
{
    int numtasks, rank;
    int nrOfProc[NR_CLUSTERS];
    char filename[50];
    int workers[50];  // stores the workers (only for roots)
    int *topology; // stores the topology
    int root;  // stores the parent of a worker (for roots it will be themselves, until they transmit it to their workers)
     int N;
     int iterations;
     int rest;
     int index;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);


    for (int i = 0; i < NR_CLUSTERS; i++) {
        
        if (rank == i) {
            // Sets the root
            root = i;

            //converts i to string;
            char index[10];
            sprintf(index,"%d",i);

            // Form the file for each cluster (ex: cluster0.txt for RANK 0)
            strcpy(filename,"cluster");
            strcat(filename,index);
            strcat(filename,".txt");

            FILE *in = fopen(filename,"rt");
            fscanf(in,"%d",&nrOfProc[i]);

            for (int j = 0; j < nrOfProc[rank]; j++) {
                fscanf(in,"%d",&workers[j]);
            }
            fclose(in);
        }
    }

    // Transmits each other the number of workers they have
    MPI_Status status;
    if (rank == 1) {

        MPI_Send(&nrOfProc[1], 1, MPI_INT, 2, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",1,2);
        MPI_Recv(&nrOfProc[0], 4, MPI_INT, 2, 1, MPI_COMM_WORLD, &status);

    } else if (rank == 2) {

        MPI_Recv(&nrOfProc[1], 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&nrOfProc[1], 2, MPI_INT, 3, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",2,3);
        MPI_Recv(&nrOfProc[0], 4, MPI_INT, 3, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&nrOfProc[0], 4, MPI_INT, 1, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",2,1);

    } else if (rank == 3) {

        MPI_Recv(&nrOfProc[1], 2, MPI_INT, 2, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&nrOfProc[1], 3, MPI_INT, 0, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,0);
        MPI_Recv(&nrOfProc[0], 4, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&nrOfProc[0], 4, MPI_INT, 2, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,2);

    } else if (rank == 0) {

        MPI_Recv(&nrOfProc[1], 3, MPI_INT, 3, 1, MPI_COMM_WORLD, &status);
        MPI_Send(&nrOfProc[0], 4, MPI_INT, 3, 1, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",0,3);
    }

    // /* This will store the whole topology
    //    from 0 to nrofProc[0] will be the workers of 0, and so on*/
    //if (rank <= 4) {
     topology = calloc(100, sizeof(int));
    //}
    // // Transmit the topology to all coordonators and all the workers
    if (rank == 1) {
        // Copies his workers to the topology
        for (int i = nrOfProc[0], x = 0; i < nrOfProc[0] + nrOfProc[1]; i++, x++) {
            topology[i] = workers[x];
        }

        MPI_Send(&topology[nrOfProc[0]], nrOfProc[1], MPI_INT, 2, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",1,2);
        MPI_Recv(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

        print_topology(rank,nrOfProc,topology);

        // Send the topology to his workers
        for (int i = 0; i < nrOfProc[1]; i++) {
            MPI_Send(&root, 1, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            printf("M(%d,%d)\n",1,workers[i]);
            MPI_Send(&nrOfProc[0], 4, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            MPI_Send(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
        }

    } else if (rank == 2) {
        // Copies his workers to the topology
        for (int i = nrOfProc[0] + nrOfProc[1], x = 0; i < nrOfProc[0] + nrOfProc[1] + nrOfProc[2]; i++, x++) {
            topology[i] = workers[x];
        }

        MPI_Recv(&topology[nrOfProc[0]], nrOfProc[1], MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&topology[nrOfProc[0]], nrOfProc[1] + nrOfProc[2], MPI_INT, 3, 0, MPI_COMM_WORLD);
        MPI_Recv(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 3, 0, MPI_COMM_WORLD, &status);

        print_topology(rank,nrOfProc,topology);

        // Send the topology to his workers
        for (int i = 0; i < nrOfProc[2]; i++) {
            MPI_Send(&root, 1, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            printf("M(%d,%d)\n",2,workers[i]);
            MPI_Send(&nrOfProc[0], 4, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            MPI_Send(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
        }

        // Send the full topology to 1
        MPI_Send(&topology[0],nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",2,1);
        
    } else if (rank == 3) {
        // Copies his workers to the topology
        for (int i = nrOfProc[0] + nrOfProc[1] + nrOfProc[2], x = 0; i < nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]; i++, x++) {
            topology[i] = workers[x];
        }
        
        MPI_Recv(&topology[nrOfProc[0]], nrOfProc[1] + nrOfProc[2], MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
        MPI_Send(&topology[nrOfProc[0]], nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,0);
        MPI_Recv(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        print_topology(rank,nrOfProc,topology);

        // Send the topology to his workers
        for (int i = 0; i < nrOfProc[3]; i++) {
            MPI_Send(&root, 1, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            printf("M(%d,%d)\n",3,workers[i]);
            MPI_Send(&nrOfProc[0], 4, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            MPI_Send(&topology[0], nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
        }
        
        // Send the full topology to 2
        MPI_Send(&topology[0],nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 2, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,2);
        
    } else if (rank == 0) {
        // Copies his workers to the topology
        for (int i = 0, x = 0; i < nrOfProc[0]; i++, x++) {
            topology[i] = workers[x];
        }

        // Receives the full topology and prints it
        MPI_Recv(&topology[nrOfProc[0]], nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
        print_topology(rank,nrOfProc,topology);

        // Send the root and the topology to his workers
        for (int i = 0; i < nrOfProc[0]; i++) {
            MPI_Send(&root, 1, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            printf("M(%d,%d)\n",0,workers[i]);
            MPI_Send(&nrOfProc[0], 4, MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
            MPI_Send(topology, nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, workers[i], workers[i], MPI_COMM_WORLD);
        }

        // Send the topology to cluster 3
        MPI_Send(&topology[0],nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, 3, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",0,3);
    }

    // Transmiting the topology to the workers by tag
    if (rank >= 4) {
        MPI_Recv(&root, 1, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
         MPI_Recv(&nrOfProc[0], 4, MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
        // topology = calloc(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], sizeof(int));
         MPI_Recv(topology, nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3], MPI_INT, MPI_ANY_SOURCE, rank, MPI_COMM_WORLD, &status);
         print_topology(rank,nrOfProc,topology);
    }

    // Transform the vector (starting from rank 0)
    if (rank == 0) {
        root = -1;
        int *v = malloc (100 * sizeof(int));

        N = atoi(argv[1]);
        for (int k = 0; k < N; k++) {
            v[k] = N - k - 1;
        }
        iterations = N/(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);
        index = 0;

        MPI_Send(&N,1, MPI_INT, 3, 0, MPI_COMM_WORLD);
        MPI_Send(v,N, MPI_INT, 3, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",0,3);

        v = send_vector(rank,index,iterations,nrOfProc,v,workers,N);

        // Receives the vector from 2
        int count;
        MPI_Recv(&rest, 1, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&count, 1, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&v[nrOfProc[0] * iterations],iterations * (nrOfProc[1] + nrOfProc[2] + nrOfProc[3]) + count , MPI_INT, 3, 0, MPI_COMM_WORLD, &status);

        v = send_rest(rank,v,workers,nrOfProc,N,&count,&rest);

        // Prints the final vector:
        printf("Rezultat:");
        for (int k = 0; k < N; k++) {
            printf(" %d",v[k]);
        }
        printf("\n");

        free(v);
    } else if (rank == 3) {
        root = -1;
        int *v = malloc (100 * sizeof(int));

        MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(v, N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

        iterations = N/(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);
        rest = N%(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);
        index = nrOfProc[0] * iterations;

        MPI_Send(&N,1, MPI_INT, 2, 0, MPI_COMM_WORLD);
        MPI_Send(v,N, MPI_INT, 2, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,2);

        v = send_vector(rank,index,iterations,nrOfProc,v,workers,N);

        // Receives the vector from 2
        int count;
        MPI_Recv(&rest, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&count, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&v[(nrOfProc[0] + nrOfProc[3]) * iterations], iterations * (nrOfProc[1] + nrOfProc[2]) + count, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

        v = send_rest(rank,v,workers,nrOfProc,N,&count,&rest);

        // Sends the number of additional itteration ant the vector
        MPI_Send(&rest,1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&count,1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&v[index],iterations * (nrOfProc[1] + nrOfProc[2] + nrOfProc[3]) + count, MPI_INT, 0, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",3,0);

        free(v);
    } else if (rank == 2) {
        root = -1;
        int *v = malloc (100 * sizeof(int));

        MPI_Recv(&N, 1, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(v, N, MPI_INT, 3, 0, MPI_COMM_WORLD, &status);

        iterations = N/(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);
        index = (nrOfProc[0] + nrOfProc[3]) * iterations;

        MPI_Send(&N,1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        MPI_Send(v,N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",2,1);

        v = send_vector(rank,index,iterations,nrOfProc,v,workers,N);

        // Receives the vector from 1
        int count;
        MPI_Recv(&rest, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&count, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&v[(nrOfProc[0] + nrOfProc[2] + nrOfProc[3]) * iterations], iterations * nrOfProc[1] + count, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);

        v = send_rest(rank,v,workers,nrOfProc,N,&count,&rest);

        // Sends the number of additional itteration ant the vector
        MPI_Send(&rest,1, MPI_INT, 3, 0, MPI_COMM_WORLD);
        MPI_Send(&count,1, MPI_INT, 3, 0, MPI_COMM_WORLD);
        MPI_Send(&v[index],iterations * (nrOfProc[1] + nrOfProc[2]) + count, MPI_INT, 3, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",2,3);
        
        free(v);
    } else if (rank == 1) {
        root = -1;
        int *v = malloc (100 * sizeof(int));

        MPI_Recv(&N, 1, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(v, N, MPI_INT, 2, 0, MPI_COMM_WORLD, &status);

        iterations = N/(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);
        index = (nrOfProc[0] + nrOfProc[2] + nrOfProc[3]) * iterations;
        rest = N%(nrOfProc[0] + nrOfProc[1] + nrOfProc[2] + nrOfProc[3]);

        // Sends the vector to the workers
        v = send_vector(rank,index,iterations,nrOfProc,v,workers,N);
        
        // Split one operation to the workers until rest is 0
        int count = 0;
        v = send_rest(rank,v,workers,nrOfProc,N,&count,&rest);

        // Sends the number of additional itteration ant the vector
        MPI_Send(&rest,1, MPI_INT, 2, 0, MPI_COMM_WORLD);
        MPI_Send(&count,1, MPI_INT, 2, 0, MPI_COMM_WORLD);
        MPI_Send(&v[index],iterations * nrOfProc[1] + count, MPI_INT, 2, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",1,2);
        free(v);
    }

    // Here workers will multiply the vector by certain parts
    if (rank >= 4) {
        int *v = malloc (100 * sizeof(int));

        MPI_Recv(&N, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(v, N, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&iterations, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&index, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);

        for (int i = index; i < index + iterations; i++) {
            v[i] *= 5;
        }
        MPI_Send(&v[index],iterations, MPI_INT, root, 0, MPI_COMM_WORLD);
        printf("M(%d,%d)\n",rank,root);

        MPI_Recv(v, N, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&rest, 1, MPI_INT, root, 0, MPI_COMM_WORLD, &status);
        if (rest > 0){
            v[N - rest] *= 5;
            MPI_Send(v,N, MPI_INT, root, 0, MPI_COMM_WORLD);
            printf("M(%d,%d)\n",rank,root);
        }

        free(v);
    } 

    free(topology);
    MPI_Finalize();
}