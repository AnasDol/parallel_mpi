#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>

#define T 100 // постоянная температура в правом верхнем углу
#define T3 50 // начальная температура внутри области

#define EPSILON 0.01  // Допустимая погрешность
#define MAX_ITER 1000   // Максимальное количество итераций

#define min(a,b) (((a) < (b)) ? (a) : (b))

void initialize(double** temperature, int m, int n, double R);
double get_new_temp(double** temperature, int m, int n, int current_row, int current_col, double R);
int compute_temperature(double** temperature, int m, int n, double R, FILE* logfile);

int main(int argc, char* argv[]) {

    if (argc < 2) {
        printf("Enter log filename as command prompt argument\n");
        return 0;
    }

    char filename[255];
    strncpy(filename, argv[1], 255);

    FILE* logfile = fopen(filename, "w+");
    if (!logfile) {
        printf("File not found\n");
        return 0;
    }

    printf("Row number (>=2): ");
    int m; // размер области по вертикали
    if (!(scanf("%d", &m)==1 && m>=2)) {
        printf("Input error\n");
        return 0;
    }

    printf("Column number (>=2): ");
    int n; // размер области по горизонтали
    if (!(scanf("%d", &n)==1 && n>=2)) {
        printf("Input error\n");
        return 0;
    }

    printf("Cut radius (1<R<%.2lf): ", (double)min(m,n)/2);
    double R; // радиус круглого выреза
    if (!(scanf("%lf", &R)==1 && R>=1 && R<=(double)min(m,n)/2)) {
        printf("Input error\n");
        return 0;
    }

    double** temperature = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        temperature[i] = (double*)malloc(n * sizeof(double));
    }

    initialize(temperature, m, n, R);

    if (m <= 30 && n <= 15) {
        printf("Initial state:\n");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                printf("%7.3lf ", temperature[i][j]);
            }
            printf("\n");
        }
    }

    //double start = omp_get_wtime();

    int count = compute_temperature(temperature, m, n, R, logfile);

    // double finish = omp_get_wtime();

    // printf("\nOutput:\n");
    // for (int i = 0; i < n; i++) {
    //     printf("%.3lf\t", temperature[0][i]);
    // }
    // printf("\nIntermediate results is saved in the log file.\n", count);
    // printf("Time: %lf\n", finish - start);
    // printf("Count: %d\n", count);

    // Вывод распределения температур
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) printf("   \t");
    //         else printf("%.2f\t", temperature[i][j]);
    //     }
    //     printf("\n");
    // }


    for (int i = 0;i<m;i++) {
        free(temperature[i]);
    }
    free(temperature);

    fclose(logfile);
    
    return 0;
}

void initialize(double** temperature, int m, int n, double R) {
    // Инициализация значений температур в начальный момент времени
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (j == n - 1) {
                temperature[i][j] = T - (double)T/(m-1)*i;  // равномерно убывающая температура на правой границе
            } 
            else if (j == 0 || i == m-1) {
                temperature[i][j] = 0;  // нулевая температура на левой и нижней границах
            } 
            else if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) <= R*R) {
                temperature[i][j] = 0;
            } 
            else {
                temperature[i][j] = T3; // начальная температура внутри области
            }
        }
    }
}

double get_new_temp(double** temperature, int m, int n, int current_row, int current_col, double R) {

    if (current_col == 0 || current_col == n-1 || current_row == m-1) return temperature[current_row][current_col]; // постоянная температура
    else if (pow(current_row-round((double)(m/2-1)), 2) + pow(current_col-round((double)(n/2)), 2) <= R*R) return 0.0;

    double newTemp = 0.0;
    int count = 0;

    for (int i = current_row-1; i<=current_row+1;i++) {
        if (i < 0 || i > m-1) continue;
        for(int j = current_col-1;j<=current_col+1;j++) {
            if (j < 0 || j > n-1) continue;
            if (pow(i-round((double)(m/2-1)), 2) + pow(j-round((double)(n/2)), 2) > R*R) {
                newTemp += temperature[i][j];
                count++;
            }
        }
    }

    if (count != 0) return newTemp/count;
    else return temperature[current_row][current_col];

}

int compute_temperature(double** temperature, int m, int n, double R, FILE* logfile) {

    // Метод последовательных приближений
    double diff = 1.0;  // Разница между текущей и предыдущей итерациями
    int iter = 0;       // Счетчик итераций

    double** new_temp = (double**)malloc(m * sizeof(double*));
    for (int i = 0; i < m; i++) {
        new_temp[i] = (double*)malloc(n * sizeof(double));
    }

    if (logfile != NULL) {
            fprintf(logfile, "[%d] ", iter);
            for (int i = 0; i < n; i++) {
                fprintf(logfile, "%.3lf ", temperature[0][i]);
            }
            fprintf(logfile, "\n");
        }

    double start = omp_get_wtime();

    while (diff > EPSILON) {

        iter++;

        diff = 0.0;

        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                
                new_temp[i][j] = get_new_temp(temperature, m, n, i, j, R);
            }
        }

        for (int i = 0; i < m; i++) {
            for (int j = 0;j < n; j++) {
                diff += fabs(new_temp[i][j] - temperature[i][j]);
                temperature[i][j] = new_temp[i][j];
                
            }
        }

        if (logfile != NULL) {
            fprintf(logfile, "[%d, diff = %lf] ", iter, diff);
            for (int i = 0; i < n; i++) {
                fprintf(logfile, "%.3lf ", temperature[0][i]);
            }
            fprintf(logfile, "\n");
        }

        double finish = omp_get_wtime();
        printf("Time: %lf\n", finish - start);
        return 0;


    }

    for (int i = 0;i<m;i++) {
        free(new_temp[i]);
    }
    free(new_temp);

    return iter;
    
}