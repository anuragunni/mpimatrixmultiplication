#define N 10000
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include "mpi.h"
#include <omp.h>

int main(int argc, char *argv[])
{
    int rank, size, dest, src, tt;
    double a11,a12,a21,a22,b11,b12,b21,b22;
    double cc[2][N],dum[2][N];
    double *a[N];
    double *b[N];
    double *c[N];
    double timestamp,timestamp2;
    double p1,p2,p3,p4,p5,p6,p7;
    MPI_Init(&argc, &argv);
    for(int k ; k<N ; k++){
    a[k] = (double*)malloc(N*sizeof(double));
    b[k] = (double*)malloc(N*sizeof(double));
    c[k] = (double*)malloc(N*sizeof(double));}
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
                a[i][j] = i%10;
                b[i][j] = j%10;
                c[i][j] = 0;
        }
    }
    for(int i = 0; i<2; i++){
        for(int j=0 ; j<N ; j++){
                dum[i][j]=0;
        }
    }
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    MPI_Request req[size];
        if (rank == 0){
                int l = 0, k = 1, t = 0, fl = 0, flag = 1;
                timestamp = MPI_Wtime();
                for (int tl = 1 ; tl < size; tl++){
                MPI_Isend(&(b[0][0]), N*N, MPI_DOUBLE, tl,99, MPI_COMM_WORLD, &req[tl]);}
                for (int tl = 1; tl < size; tl++){
                                         MPI_Wait(&req[tl], &status);}
                while(l<N){
                        for (int i=0; i<2; i++){
                                for (int j=0; j<N; j++){
                                        dum[i][j] = a[l][j];
                                }
                                l = l+1;
                        }
                        dest = (k % size);
                        MPI_Isend(&(dum[0][0]), 2*N, MPI_DOUBLE, dest, 99,  MPI_COMM_WORLD, &req[dest]);
                        k = k + 1;
                        if (dest == size-1 || l == N){
                                if (dest == size-1){ fl = size-1; for (int k = 1; k <=fl ; k++){MPI_Wait(&req[k], &status);}}
                                        else{fl = dest;  for (int k = 1; k <=fl ; k++){MPI_Wait(&req[k], &status);}}
                                for (int k = 1; k <= fl; k++){
                                src = (k % size);
                                MPI_Irecv(&(cc[0][0]), 2*N, MPI_DOUBLE, src, 99,  MPI_COMM_WORLD, &req[src]);
                                MPI_Wait(&req[src], &status);
                                for (int i=0; i<2; i++){
                                        for (int j=0; j<N; j++){
                                                c[t][j] = cc[i][j];
                                        }
                                        t = t + 1;
                                }
                                }
                        k = 1;}
                }
                timestamp2 = MPI_Wtime();
                printf("Total Time for Computation %f \n", timestamp2-timestamp);
}


        if (rank != 0){
                MPI_Irecv(&(b[0][0]), N*N, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &req[rank]);
                MPI_Wait(&req[rank], &status);
                int K = N;
                int k = 0, iter = 0, flag = 0;
                iter = K/(2*(size-1)) + 1;
                flag = (K/2) % (size-1);
                while(iter  > 1 || rank <= flag){
                        int i = 0, j = 0, l = 0, r = 0, col = 0;
                        for (int p = 0 ; p < 2; p++){
                                for (int q = 0; q < N; q++){
                                        cc[p][q] = 0;}}
                        MPI_Irecv(&(dum[0][0]), 2*N, MPI_DOUBLE, 0, 99,  MPI_COMM_WORLD, &req[rank]);
                        MPI_Wait(&req[rank], &status);

                        while(l < N){
                                #pragma omp parallel for private(j,k) shared(dum,cc) num_threads(20)
                                for(k = 0; k<N; k=k+2){
                                                a11=dum[i][j];
                                                a12=dum[i][j+1];
                                                a21=dum[i+1][j];
                                                a22=dum[i+1][j+1];
                                                b11=b[k][l];
                                                b12=b[k][l+1];
                                                b21=b[k+1][l];
                                                b22=b[k+1][l+1];
                                                p1 = (a11+a22)*(b11+b22);
                                                p2 = (a21+a22)*b11;
                                                p3 = (a11)*(b12 - b22);
                                                p4 = a22*(b21 - b11);
                                                p5 = (a11 + a12)*b22;
                                                p6 = (a21 - a11)*(b11 + b12);
                                                p7 = (a12 - a22)*(b21 + b22);
                                                cc[r][col] = cc[r][col] + p1 + p4 - p5 + p7;
                                                cc[r][col+1] = cc[r][col+1] + p3 + p5;
                                                cc[r+1][col] = cc[r+1][col] + p2 + p4;
                                                cc[r+1][col+1] = cc[r+1][col+1] + p1 + p3 - p2 + p6;
                                                j=j+2;
                                }
                                l = l + 2;
                                col = col + 2;
                                j = 0;
                        }
                        MPI_Isend(&(cc[0][0]), 2*N, MPI_DOUBLE, 0, 99,  MPI_COMM_WORLD, &req[rank]);
                        MPI_Wait(&req[rank], &status);
                        iter = iter - 1;
                        if(iter == 0)
                                 { flag = 0; }
                }
                }

   MPI_Finalize();
}

