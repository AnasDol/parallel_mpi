#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define T 100 // постоянная температура в правом верхнем углу
#define T3 50 // начальная температура внутри области

#define EPSILON 0.01  // Допустимая погрешность
#define MAX_ITER 1000   // Максимальное количество итераций

#define min(a,b) (((a) < (b)) ? (a) : (b))

void initialize(double* temperature, int m, int n, double R);
double get_new_temp(double* emperature, int m, int n, int current_row, int current_col, double R);
double compute_temp_one_step(double* temperature, int m, int n, double R, int ProcRank, int ProcNum);
int left(int m, int n, int ProcRank, int ProcNum) { return (ProcRank == 0) ? 1 : (n / ProcNum) * ProcRank; }
int right(int m, int n, int ProcRank, int ProcNum) { return (ProcRank == ProcNum-1) ? n-1 : (n / ProcNum) * (ProcRank + 1) - 1; }

int main(int argc, char* argv[]) {

    int ProcNum, ProcRank;

    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int n, m;
    double R;

    FILE* logfile;

    double start = 0, finish = 0;

    if (ProcRank == 0) {
        if (argc < 2) {
            printf("Enter log filename as command prompt argument\n");
            return 0;
        }

        char filename[255];
        strncpy(filename, argv[1], 255);

        logfile = fopen(filename, "w+");
        if (!logfile) {
            printf("File not found\n");
            return 0;
        }

        printf("Row number (>=2): ");
        if (!(scanf("%d", &m)==1 && m>=2)) {
            printf("Input error\n");
            return 0;
        }

        printf("Column number (>=2): ");
        if (!(scanf("%d", &n)==1 && n>=2)) {
            printf("Input error\n");
            return 0;
        }

        printf("Cut radius (1<R<%.2lf): ", (double)min(m,n)/2);
        if (!(scanf("%lf", &R)==1 && R>=1 && R<=(double)min(m,n)/2)) {
            printf("Input error\n");
            return 0;
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&R, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* temperature = (double*)malloc(n * m * sizeof(double));
    double* flattened_temperature = (double*)malloc(n * m * sizeof(double));

    int iter = 0;

    if (ProcRank == 0) {

        initialize(temperature, m, n, R);

        if (m <= 30 && n <= 15) {
            printf("Initial state:");
            for (int i = 0; i< n * m; i++) {
                if (i%n == 0) printf("\n");
                printf("%7.3lf ", temperature[i]);
            }
        }

        if (logfile != NULL) {
            fprintf(logfile, "[%d] ", iter);
            for (int i = 0; i < n; i++) {
                fprintf(logfile, "%.3lf ", temperature[i]);
            }
            fprintf(logfile, "\n");
        }

        start = MPI_Wtime();
        
    }

    MPI_Barrier(MPI_COMM_WORLD);

    double max_diff = 1.0;
    double diff;

    while (max_diff > EPSILON) {

        MPI_Bcast(temperature, m*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        diff = compute_temp_one_step(temperature, m, n, R, ProcRank, ProcNum);

        MPI_Barrier(MPI_COMM_WORLD);

        if (ProcRank == 0) {

            iter++;

            max_diff = diff;

            for (int i = 1; i < ProcNum; i++) {

                MPI_Recv(flattened_temperature, m*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
                MPI_Recv(&diff, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
            
                max_diff += diff;

                for (int j = 0; j < m; j++) {
                    for (int k = left(m, n, i, ProcNum); k <= right(m, n, i, ProcNum); k++) {
                        temperature[j * n + k] = flattened_temperature[j * n + k];
                    }
                }

            }

            //printf("\n[%d, max_diff = %lf]: \n", iter, max_diff);

            // for (int i = 0; i < m; i++) {
            //     for (int j = 0; j < n; j++) {
            //         printf("%7.3lf ", temperature[i][j]);
            //     }
            //     printf("\n");
            // }

            if (logfile != NULL) {
                fprintf(logfile, "[%d, diff = %lf] ", iter, max_diff);
                for (int i = 0; i < n; i++) {
                    fprintf(logfile, "%.3lf ", temperature[i]);
                }
                fprintf(logfile, "\n");
            }

            finish = MPI_Wtime();
            printf("Time: %lf\n", finish - start);
            return 0;
            

            // Вывод верхней границы
            // for (int i = 0; i < m; i++) {
            //     printf("%-5.3lf ", flattenedTemperature[i]);
            // }

        } else {

            MPI_Send(temperature, m*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&diff, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&max_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }  

    if (ProcRank == 0) {
        finish = MPI_Wtime();
        printf("\nOutput:\n");
        for (int i = 0; i < n; i++) {
            printf("%.3lf\t", temperature[i]);
        }
        printf("\nIntermediate results is saved in the log file.\n");
        printf("Time: %lf\n", finish - start);
        printf("Count: %d\n", iter);
    } 
    

    free(temperature);
    free(flattened_temperature);

    // Завершаем работу с MPI
    MPI_Finalize();

    return 0;
}

void initialize(double* temperature, int m, int n, double R) {
    // Инициализация значений температур в начальный момент времени
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (j == n - 1) {
                temperature[i * n + j] = T - (double)T/(m-1)*i;  // равномерно убывающая температура на правой границе
            } 
            else if (j == 0 || i == m-1) {
                temperature[i * n + j] = 0;  // нулевая температура на левой и нижней границах
            }
            else if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) {
                temperature[i * n + j] = 0;
            } 
            else {
                temperature[i * n + j] = T3; // начальная температура внутри области
            }
        }
    }
}

double get_new_temp(double* temperature, int m, int n, int current_row, int current_col, double R) {

    if (current_col == 0 || current_col == n-1 || current_row == m-1) return temperature[current_row * n + current_col]; // постоянная температура
    
    else if (pow(current_row-round((double)(m/2-1)), 2) + pow(current_col-round((double)(n/2)), 2) <= R*R) return 0.0;

    double newTemp = 0.0;
    int count = 0;

    for (int i = current_row-1; i<=current_row+1;i++) {
        if (i < 0 || i > m-1) continue;
        for(int j = current_col-1;j<=current_col+1;j++) {
            if (j < 0 || j > n-1) continue;
            if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) > R*R) {
                newTemp += temperature[i * n + j];
                count++;
            }
        }
    }

    if (count != 0) return newTemp/count;
    else temperature[current_row * n + current_col];

}


double compute_temp_one_step(double* temperature, int m, int n, double R, int ProcRank, int ProcNum) {

    double diff = 0;

    double* temp = (double*)malloc(m * n * sizeof(double));

    for (int i = 0; i < m * n; i++) {
        temp[i] = temperature[i];
    }

    for (int i = 0; i < m; i++) {
        for (int j = left(m, n, ProcRank, ProcNum); j <= right(m, n, ProcRank, ProcNum); j++) { 
            double new_temp = get_new_temp(temperature, m, n, i, j, R); 
            diff += fabs(new_temp - temp[i * n + j]);
            temp[i * n + j] = new_temp;
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = left(m, n, ProcRank, ProcNum); j <= right(m, n, ProcRank, ProcNum); j++) {
            temperature[i * n + j] = temp[i * n + j];
        }
    }

    return diff;

}

