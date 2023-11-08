#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>

#define T 100 // постоянная температура в правом верхнем углу
#define T3 50 // начальная температура внутри области

#define EPSILON 0.001  // Допустимая погрешность
#define MAX_ITER 1000   // Максимальное количество итераций

#define min(a,b) (((a) < (b)) ? (a) : (b))

void initialize(double** temperature, int m, int n);
double get_new_temp(double** temperature, int rows, int cols, int current_row, int current_col, double R);
double compute_temp_one_step(double** temperature, int m, int n, double R, int ProcRank, int ProcNum);
int left(int m, int n, int ProcRank, int ProcNum) { return (ProcRank == 0) ? 1 : (n / ProcNum) * ProcRank; }
int right(int m, int n, int ProcRank, int ProcNum) { return (ProcRank == ProcNum-1) ? n-1 : (n / ProcNum) * (ProcRank + 1) - 1; }

int main(int argc, char* argv[]) {

    int ProcNum, ProcRank, RecvRank;

    MPI_Status Status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int n, m;
    double R;

    FILE* logfile;

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

    //printf("%d: m = %d, n = %d, R = %lf\n", ProcRank, m, n, R);

    double** temperature = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        temperature[i] = (double*)malloc(n * sizeof(double));
    }

    double* flattenedTemperature = (double*)malloc(n * m * sizeof(double));

    int iter = 0;

    if (ProcRank == 0) {

        initialize(temperature, m, n);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                flattenedTemperature[i * n + j] = temperature[i][j];
            }
        }

        printf("Initial state:\n");
        for (int i = 0; i< n * m; i++) {
            if (i%n == 0) printf("\n");
            printf("%.2lf ", flattenedTemperature[i]);
        }

        if (logfile != NULL) {
            fprintf(logfile, "[%d] ", iter);
            for (int i = 0; i < n; i++) {
                fprintf(logfile, "%.2lf ", temperature[0][i]);
            }
            fprintf(logfile, "\n");
        }
        
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(flattenedTemperature, m*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double max_diff = 1.0;
    double diff;

    while (max_diff > EPSILON) {

        MPI_Bcast(flattenedTemperature, m*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                temperature[i][j] = flattenedTemperature[i * n + j];
            }
        }

        diff = compute_temp_one_step(temperature, m, n, R, ProcRank, ProcNum);

        MPI_Barrier(MPI_COMM_WORLD);

        if (ProcRank == 0) {

            printf("\n\n[%d]:\n", iter);
            iter++;

            max_diff = diff;

            for (int i = 1; i < ProcNum; i++) {

                MPI_Recv(flattenedTemperature, m*n, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
                MPI_Recv(&diff, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
            
                if (diff > max_diff) max_diff = diff;

                for (int j = 0; j < m; j++) {
                    for (int k = left(m, n, i, ProcNum); k <= right(m, n, i, ProcNum); k++) {
                        temperature[j][k] = flattenedTemperature[j * n + k];
                    }
                }

            }

            if (logfile != NULL) {
                fprintf(logfile, "[%d] ", iter);
                for (int i = 0; i < n; i++) {
                    fprintf(logfile, "%.2lf ", temperature[0][i]);
                }
                fprintf(logfile, "\n");
            }

            for (int i = 0; i < m; i++) {
                int j;
                for (j = 0; j < n; j++) {
                    flattenedTemperature[i * n + j] = temperature[i][j];
                }
            }

            for (int i = 0; i < m*n; i++) {
                if (i%n==0) printf("\n");
                printf("%-5.2lf ", flattenedTemperature[i]);
            }

        } else {

            for (int i = 0; i < m; i++) {
                for (int j = 0; j < n; j++) {
                    flattenedTemperature[i * n + j] = temperature[i][j];
                }
            }

            MPI_Send(flattenedTemperature, m*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            MPI_Send(&diff, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&max_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }    
    
    for (int i = 0;i<m;i++) {
        free(temperature[i]);
    }
    free(temperature);

    free(flattenedTemperature);

    // Завершаем работу с MPI
    MPI_Finalize();

    return 0;
}

void initialize(double** temperature, int rows, int cols) {
    // Инициализация значений температур в начальный момент времени
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (j == cols - 1) {
                temperature[i][j] = T - (double)T/(rows-1)*i;  // равномерно убывающая температура на правой границе
            } 
            else if (j == 0 || i == rows-1) {
                temperature[i][j] = 0;  // нулевая температура на левой и нижней границах
            } else {
                temperature[i][j] = T3; // начальная температура внутри области
            }
        }
    }
}

double get_new_temp(double** temperature, int rows, int cols, int current_row, int current_col, double R) {

    if (current_col == 0 || current_col == cols-1 || current_row == rows-1) return temperature[current_row][current_col]; // постоянная температура
    
    else if (pow(current_row-round((double)(rows/2-1)), 2) + pow(current_col-round((double)(cols/2)), 2) <= R*R) return 0.0;

    double newTemp = 0.0;
    int count = 0;

    for (int i = current_row-1; i<=current_row+1;i++) {
        if (i < 0 || i > rows-1) continue;
        for(int j = current_col-1;j<=current_col+1;j++) {
            if (j < 0 || j > cols-1) continue;
            if (pow(i-round((double)(rows/2-1)), 2) + pow(j-round((double)(cols/2)), 2) > R*R) {
                newTemp += temperature[i][j];
                count++;
            }
        }
    }

    if (count != 0) return newTemp/count;
    else temperature[current_row][current_col];

}

double compute_temp_one_step(double** temperature, int m, int n, double R, int ProcRank, int ProcNum) {

    double diff = 0;

    double** new_temp = (double**)malloc(m * sizeof(double*));
    int i;
    for (i = 0; i < m; i++) {
        new_temp[i] = (double*)malloc(n * sizeof(double));
    }

    for (i = 0; i < m; i++) {
        int j;
        for (j = 0; j < n; j++) {
            new_temp[i][j] = temperature[i][j];
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = left(m, n, ProcRank, ProcNum); j <= right(m, n, ProcRank, ProcNum); j++) {   
            new_temp[i][j] = get_new_temp(temperature, m, n, i, j, R);
        }
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0;j < n; j++) {
            diff += fabs(new_temp[i][j] - temperature[i][j]);
            temperature[i][j] = new_temp[i][j];
        }
    }

    return diff;

}

