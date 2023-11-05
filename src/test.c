#include <stdio.h>
#include "mpi.h"

#define m 8
#define n 12

 int main(int argc, char* argv[]) {

    int ProcNum, ProcRank, RecvRank;
    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    double data[m][n];
    double data_flattened[m*n];

    if ( ProcRank == 0 ){

        int i;

        for (i = 0;i<m;i++) {
            int j;
            for (j = 0;j<n;j++) {
                data[i][j] = 1;
            }
        }

        for (i = 0; i < m; i++) {
            int j;
            for (j = 0; j < n; j++) {
                data_flattened[i * n + j] = data[i][j];
            }
        }
        
    }

    MPI_Bcast(data_flattened, m * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int i;
    for (i = 0; i < m; i++) {
        int j;
        for (j = 0; j < n; j++) {
            data[i][j] = data_flattened[i * n + j];
        }
    }
    
    int left, right;

    if (ProcRank == 0) left = 0;
    else left = (n / ProcNum) * ProcRank - 1;

    if (ProcRank == ProcNum - 1) right = n - 1;
    else right = (n / ProcNum) * (ProcRank + 1);

    printf("%d : left = %d, right = %d\n", ProcRank, left, right);

    double ProcSum = 0.0;
    double TotalSum;

    for (i = 0; i < m; i++) {
        int j;
        for (j = left; j <= right; j++) {
            ProcSum += data[i][j];
        }
    }

    if (ProcRank == 0) {

        TotalSum = ProcSum; // часть суммы, подсчитанная потоком 0
        printf("%d: %lf\n", 0, ProcSum);

        for (i=1; i < ProcNum; i++ ) {
            MPI_Recv(&ProcSum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
            printf("%d: %lf\n", i, ProcSum);
            TotalSum = TotalSum + ProcSum;
        }

    }

    else
        // все потоки отсылают свои частичные суммы потоку 0
        MPI_Send(&ProcSum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    

    // вывод результата
    if ( ProcRank == 0 )
        printf("TotalSum: %lf\n", TotalSum);



    MPI_Finalize();
    return 0;
 } 